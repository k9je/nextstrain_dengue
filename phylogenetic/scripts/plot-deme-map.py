#!/usr/bin/env python3
"""
Draw a publication-quality phylogeographic map from a Nextstrain/Auspice JSON.

Unlike Auspice, which scales deme circles continuously by tip count (so one
heavily-sequenced deme swamps the map), this draws:

  * one PIE per deme, sized from a 3-level BIN of its tip count, with slices
    coloured by Nextclade major lineage;
  * transmission lines whose widths come from a 3-level BIN of the number of
    inferred deme transitions, with deliberately exaggerated contrast.

Transmission lines are counted the same way Auspice draws them: a branch whose
parent and child carry different reconstructed `deme` values is one event.

Because dengue demes cluster tightly in Southeast Asia, overlapping pies are
pushed apart by a short force-directed relaxation and tied back to their true
centroid with a leader line and a dot. Labels are then placed by testing
candidate positions against the real rendered text extents.
"""
import argparse
import json
from collections import Counter, defaultdict
from itertools import combinations

import matplotlib

matplotlib.use("Agg")
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.patheffects as pe
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.patches import Circle, FancyArrowPatch, Wedge

# --- binning ---------------------------------------------------------------
# Pie size: tip-count bins. Radii (degrees) are deliberately NOT
# area-proportional -- that is the whole point of binning.
SIZE_BINS = [(1, 9, 1.7, "1-9"), (10, 99, 3.0, "10-99"), (100, 10**9, 4.6, "100+")]

# Line width: transition-count bins, with exaggerated contrast.
LINE_BINS = [
    (1, 2, 0.6, 0.26, False, "1-2"),
    (3, 9, 2.4, 0.58, True, "3-9"),
    (10, 10**9, 5.6, 0.88, True, "10+"),
]

LINE_COLOR = "#1f1f1f"
# Focus deme highlighting: imports red, exports blue.
IN_COLOR = "#d62728"
OUT_COLOR = "#1f77b4"
FOCUS_MIN_WIDTH = 1.5   # keep highlighted lines visible even in the 1-2 bin
FOCUS_MIN_ALPHA = 0.92

PIE_PAD = 2.6          # extra degrees kept between pies (room for labels)
MAX_SHIFT = 30.0       # max degrees a pie may be displaced from its centroid
PUSH_GAIN = 0.35       # how hard overlapping pies repel
SPRING_GAIN = 0.0030   # how hard a pie is pulled back to its true centroid
LEADER_MIN = 0.9       # draw a leader line past this displacement

# Map extent and layout (inches). The axes box is sized to the map's own
# aspect ratio so cartopy does not letterbox it inside a taller box. The
# latitude window is derived from the data so no pie is clipped at an edge.
LAT_PAD = 6.0
LAT_MIN_BOUNDS = (-80.0, -42.0)   # (hard floor, latest possible top of range)
LAT_MAX_BOUNDS = (52.0, 84.0)
FIG_WIDTH = 16.5
AX_WIDTH_FRAC = 0.984
LEGEND_BAND = 2.0
TITLE_BAND = 0.5


def bin_for(value, bins):
    for b in bins:
        if b[0] <= value <= b[1]:
            return b
    return bins[-1]


def walk_tree(root, deme_key, lineage_key):
    """Return (tips-per-deme-per-lineage, directed transition counts)."""
    per_deme = defaultdict(Counter)
    transitions = Counter()

    def recurse(node, parent_deme):
        attrs = node.get("node_attrs", {})
        deme = attrs.get(deme_key, {}).get("value")
        children = node.get("children")
        if parent_deme and deme and parent_deme != deme:
            transitions[(parent_deme, deme)] += 1
        if not children:
            if deme:
                lineage = attrs.get(lineage_key, {}).get("value") or "unassigned"
                per_deme[deme][lineage] += 1
        for child in children or []:
            recurse(child, deme or parent_deme)

    recurse(root, None)
    return per_deme, transitions


def lineage_colors(meta, lineage_key, lineages):
    """Reuse the dataset's own colour scale so the figure matches Auspice."""
    colors = {}
    for coloring in meta.get("colorings", []):
        if coloring["key"] == lineage_key:
            for entry in coloring.get("scale", []) or []:
                colors[entry[0]] = entry[1]
    fallback = plt.get_cmap("tab20")
    for i, lineage in enumerate(sorted(lineages)):
        colors.setdefault(lineage, matplotlib.colors.to_hex(fallback(i % 20)))
    colors["unassigned"] = "#bdbdbd"
    return colors


def deme_coordinates(meta, deme_key):
    for resolution in meta.get("geo_resolutions", []):
        if resolution["key"] == deme_key:
            return {
                name: (float(v["longitude"]), float(v["latitude"]))
                for name, v in resolution["demes"].items()
            }
    raise SystemExit(f"No geo_resolution {deme_key!r} in this dataset")


def relax_positions(true_pos, radii, pad=PIE_PAD, max_shift=MAX_SHIFT,
                    iterations=1500):
    """Push overlapping pies apart, with a weak spring back to the centroid."""
    pos = {k: np.array(v, dtype=float) for k, v in true_pos.items()}
    names = list(pos)
    for _ in range(iterations):
        moved = False
        for a, b in combinations(names, 2):
            delta = pos[b] - pos[a]
            dist = float(np.hypot(*delta))
            want = radii[a] + radii[b] + pad
            if dist < 1e-6:
                delta = np.array([0.5, 0.5])
                dist = float(np.hypot(*delta))
            if dist < want:
                push = (want - dist) * PUSH_GAIN
                unit = delta / dist
                pos[a] -= unit * push
                pos[b] += unit * push
                moved = True
        # weak spring home, and hard cap on displacement
        for name in names:
            home = np.array(true_pos[name], dtype=float)
            pos[name] += (home - pos[name]) * SPRING_GAIN
            offset = pos[name] - home
            dist = float(np.hypot(*offset))
            if dist > max_shift:
                pos[name] = home + offset / dist * max_shift
        if not moved:
            break
    return pos


def wrapped_offsets(x1, x2):
    """Longitude offsets needed to draw a line that may cross the map edge."""
    if x2 - x1 > 180:
        return [(0.0, -360.0), (360.0, 0.0)]
    if x2 - x1 < -180:
        return [(0.0, 360.0), (-360.0, 0.0)]
    return [(0.0, 0.0)]


def place_labels(ax, fig, entries, pie_display, fontsize):
    """Greedily place labels at the first candidate offset that does not collide.

    `entries` is a list of (name, text, data_xy, radius_deg) sorted so the
    most important labels are placed first. Collision testing uses the real
    rendered text extents in display coordinates.
    """
    renderer = fig.canvas.get_renderer()
    placed = []
    results = []

    # candidate directions, in order of preference
    directions = [
        (0, -1), (0, 1), (1, 0), (-1, 0),
        (0.85, -0.85), (-0.85, -0.85), (0.85, 0.85), (-0.85, 0.85),
    ]

    for name, text, (dx, dy), radius_deg in entries:
        cx, cy = ax.transData.transform((dx, dy))
        redge = abs(ax.transData.transform((dx + radius_deg, dy))[0] - cx)

        probe = ax.text(0, 0, text, fontsize=fontsize, ha="center", va="center")
        bbox = probe.get_window_extent(renderer)
        tw, th = bbox.width, bbox.height
        probe.remove()

        best, best_score = None, None
        for scale in (1.0, 1.5, 2.1, 2.9, 3.8):
            for ux, uy in directions:
                gap = redge + 8 + (th / 2 if uy else tw / 2) * 0.55
                px = cx + ux * gap * scale
                py = cy + uy * gap * scale
                rect = (px - tw / 2, py - th / 2, px + tw / 2, py + th / 2)
                score = 0.0
                for ox, oy, orad in pie_display:
                    nx = max(rect[0], min(ox, rect[2]))
                    ny = max(rect[1], min(oy, rect[3]))
                    if np.hypot(nx - ox, ny - oy) < orad:
                        score += orad
                for prect in placed:
                    ox1 = min(rect[2], prect[2]) - max(rect[0], prect[0])
                    oy1 = min(rect[3], prect[3]) - max(rect[1], prect[1])
                    if ox1 > 0 and oy1 > 0:
                        score += (ox1 * oy1) ** 0.5
                score += 0.35 * (scale - 1.0) * (tw + th) / 2
                if best_score is None or score < best_score:
                    best_score, best = score, (px, py, rect)
                if score == 0:
                    break
            if best_score == 0:
                break

        px, py, rect = best
        placed.append(rect)
        inv = ax.transData.inverted()
        results.append((name, text, inv.transform((px, py)), (dx, dy), redge))
    return results


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--auspice-json", required=True)
    parser.add_argument("--output", required=True, help="output path; .pdf also written")
    parser.add_argument("--title", default=None)
    parser.add_argument("--deme-key", default="deme")
    parser.add_argument("--lineage-key", default="major_lineage")
    parser.add_argument("--central-longitude", type=float, default=None,
                        help="default 150 globally, or the centre of --extent")
    parser.add_argument(
        "--extent", nargs=4, type=float, default=None,
        metavar=("LON_MIN", "LON_MAX", "LAT_MIN", "LAT_MAX"),
        help="crop to this box instead of drawing the whole world",
    )
    parser.add_argument(
        "--pie-scale", type=float, default=1.0,
        help="multiply all pie radii (and the spacing/displacement that "
             "depend on them). Use well below 1 when cropping, since radii "
             "are in degrees and a cropped map spans far fewer of them.",
    )
    parser.add_argument("--line-scale", type=float, default=1.0,
                        help="multiply all transmission-line widths")
    parser.add_argument("--fig-width", type=float, default=FIG_WIDTH)
    parser.add_argument("--coastlines", default=None, choices=["110m", "50m", "10m"],
                        help="default 110m globally, 50m when cropped")
    parser.add_argument("--label-fontsize", type=float, default=6.6)
    parser.add_argument("--min-line-count", type=int, default=1)
    parser.add_argument(
        "--focus-deme", default="Singapore",
        help="highlight transitions INTO this deme in red and OUT of it in "
             "blue; pass an empty string to disable",
    )
    parser.add_argument("--no-labels", action="store_true")
    args = parser.parse_args()

    with open(args.auspice_json, encoding="utf-8") as fh:
        data = json.load(fh)
    meta = data["meta"]
    per_deme, transitions = walk_tree(data["tree"], args.deme_key, args.lineage_key)
    coords = deme_coordinates(meta, args.deme_key)
    per_deme = {d: c for d, c in per_deme.items() if d in coords}

    all_lineages = sorted({lin for c in per_deme.values() for lin in c})
    colors = lineage_colors(meta, args.lineage_key, all_lineages)

    # The lineage legend drives how deep the legend band has to be: DENV-1 and
    # DENV-2 have 24 Nextclade major lineages where DENV-4 has 7, and a fixed
    # band would push the legend off the bottom of the figure.
    lineage_ncol = max(3, min(5, -(-len(all_lineages) // 6)))
    lineage_rows = -(-len(all_lineages) // lineage_ncol)
    legend_band = max(LEGEND_BAND, 0.55 + lineage_rows * 0.21 + 0.50)

    if args.central_longitude is not None:
        central_longitude = args.central_longitude
    elif args.extent:
        central_longitude = (args.extent[0] + args.extent[1]) / 2.0
    else:
        central_longitude = 150.0
    projection = ccrs.PlateCarree(central_longitude=central_longitude)

    # Pie radii are in degrees, so a cropped map needs much smaller ones.
    scale = args.pie_scale

    # Lay the pies out before building the canvas, so the latitude window can
    # be sized to whatever the relaxation produced.
    totals = {d: sum(c.values()) for d, c in per_deme.items()}
    radii = {d: bin_for(n, SIZE_BINS)[2] * scale for d, n in totals.items()}
    true_pos = {
        d: projection.transform_point(*coords[d], ccrs.PlateCarree())
        for d in per_deme
    }
    pos = relax_positions(true_pos, radii, pad=PIE_PAD * scale,
                          max_shift=MAX_SHIFT * scale)
    leader_min = LEADER_MIN * scale

    if args.extent:
        lon_min, lon_max, lat_min, lat_max = args.extent
        lon_span = lon_max - lon_min
    else:
        lon_min = lon_max = None
        lon_span = 360.0
        tops = [pos[d][1] + radii[d] for d in pos] + [t[1] for t in true_pos.values()]
        bottoms = [pos[d][1] - radii[d] for d in pos] + [t[1] for t in true_pos.values()]
        lat_max = min(LAT_MAX_BOUNDS[1], max(LAT_MAX_BOUNDS[0], max(tops) + LAT_PAD))
        lat_min = max(LAT_MIN_BOUNDS[0], min(LAT_MIN_BOUNDS[1], min(bottoms) - LAT_PAD))

    # Size the canvas to the map's aspect ratio, plus bands for title/legends.
    ax_width_in = args.fig_width * AX_WIDTH_FRAC
    ax_height_in = ax_width_in / (lon_span / (lat_max - lat_min))
    fig_height = ax_height_in + legend_band + TITLE_BAND
    legend_y = legend_band / fig_height
    fig = plt.figure(figsize=(args.fig_width, fig_height))
    ax = fig.add_axes(
        [(1 - AX_WIDTH_FRAC) / 2, legend_y, AX_WIDTH_FRAC,
         ax_height_in / fig_height],
        projection=projection,
    )
    resolution = args.coastlines or ("50m" if args.extent else "110m")
    ax.add_feature(cfeature.LAND.with_scale(resolution), facecolor="#f4f3ee",
                   edgecolor="none", zorder=0)
    ax.add_feature(cfeature.OCEAN.with_scale(resolution), facecolor="#ffffff",
                   edgecolor="none", zorder=0)
    ax.add_feature(cfeature.COASTLINE.with_scale(resolution), linewidth=0.35,
                   edgecolor="#9a9a9a", zorder=1)
    ax.add_feature(cfeature.BORDERS.with_scale(resolution), linewidth=0.22,
                   edgecolor="#c9c9c9", zorder=1)
    if args.extent:
        ax.set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())
    else:
        ax.set_global()
        ax.set_ylim(lat_min, lat_max)
    ax.spines["geo"].set_edgecolor("#999999")
    ax.spines["geo"].set_linewidth(0.6)

    # --- transmission lines (under the pies) -------------------------------
    focus = args.focus_deme or None
    if focus and focus not in pos:
        print(f"  NOTE: focus deme {focus!r} is not on this tree; no highlighting")
        focus = None

    def styled(source, target, width, alpha):
        """Colour, alpha, width and z-order for one transmission line."""
        if focus and target == focus:
            return IN_COLOR, max(alpha, FOCUS_MIN_ALPHA), max(width, FOCUS_MIN_WIDTH), 3.7
        if focus and source == focus:
            return OUT_COLOR, max(alpha, FOCUS_MIN_ALPHA), max(width, FOCUS_MIN_WIDTH), 3.6
        return LINE_COLOR, alpha, width, 2 + width / 10

    drawn = Counter()
    focus_counts = Counter()
    # Draw background lines first, then the highlighted ones on top.
    ordered = sorted(
        transitions.items(),
        key=lambda kv: (focus in kv[0] if focus else False, kv[1]),
    )
    for (source, target), count in ordered:
        if count < args.min_line_count or source not in pos or target not in pos:
            continue
        _lo, _hi, width, alpha, arrow, label = bin_for(count, LINE_BINS)
        drawn[label] += 1
        color, alpha, width, zorder = styled(source, target, width, alpha)
        width *= args.line_scale
        if focus and target == focus:
            focus_counts["in"] += count
        elif focus and source == focus:
            focus_counts["out"] += count
        (x1, y1), (x2, y2) = pos[source], pos[target]
        for off1, off2 in wrapped_offsets(x1, x2):
            ax.add_patch(FancyArrowPatch(
                (x1 + off1, y1), (x2 + off2, y2),
                connectionstyle="arc3,rad=0.16",
                arrowstyle="-|>" if arrow else "-",
                mutation_scale=8 + 2.0 * width,
                linewidth=width, color=color, alpha=alpha,
                zorder=zorder, capstyle="round", joinstyle="round",
                shrinkA=2.5, shrinkB=2.5,
            ))

    # --- leader lines back to true centroids --------------------------------
    for deme in per_deme:
        (tx, ty), (px, py) = true_pos[deme], pos[deme]
        if np.hypot(px - tx, py - ty) > leader_min:
            ax.plot([tx, px], [ty, py], color="#7a7a7a", linewidth=0.5,
                    alpha=0.75, zorder=3, solid_capstyle="round")
            ax.plot([tx], [ty], marker="o", markersize=1.7, color="#5a5a5a",
                    zorder=3.5, linestyle="")

    # --- pies ---------------------------------------------------------------
    for deme in sorted(per_deme, key=lambda d: -totals[d]):
        counts, total, radius = per_deme[deme], totals[deme], radii[deme]
        x, y = pos[deme]
        ax.add_patch(Circle((x, y), radius * 1.1, facecolor="white",
                            edgecolor="none", alpha=0.8, zorder=4))
        start = 90.0
        for lineage, n in sorted(counts.items(), key=lambda kv: (-kv[1], kv[0])):
            extent = 360.0 * n / total
            ax.add_patch(Wedge((x, y), radius, start - extent, start,
                               facecolor=colors[lineage], edgecolor="white",
                               linewidth=0.4, zorder=5))
            start -= extent
        ax.add_patch(Circle((x, y), radius, facecolor="none",
                            edgecolor="#4a4a4a", linewidth=0.5, zorder=6))

    # --- labels -------------------------------------------------------------
    # When cropped, only label demes actually inside the frame; the rest are
    # off-map and only contribute lines running off the edge.
    x_lo, x_hi = ax.get_xlim()
    y_lo, y_hi = ax.get_ylim()
    visible = [
        d for d in per_deme
        if x_lo <= pos[d][0] <= x_hi and y_lo <= pos[d][1] <= y_hi
    ]

    if not args.no_labels:
        fig.canvas.draw()
        pie_display = []
        for deme in visible:
            cx, cy = ax.transData.transform(pos[deme])
            redge = abs(ax.transData.transform(
                (pos[deme][0] + radii[deme], pos[deme][1]))[0] - cx)
            pie_display.append((cx, cy, redge))

        entries = [
            (d, f"{d} ({totals[d]})", pos[d], radii[d])
            for d in sorted(visible, key=lambda d: -totals[d])
        ]
        for name, text, (lx, ly), (ax_, ay), redge in place_labels(
                ax, fig, entries, pie_display, fontsize=args.label_fontsize):
            ax.text(lx, ly, text, ha="center", va="center", fontsize=args.label_fontsize,
                    color="#141414", zorder=8,
                    path_effects=[pe.withStroke(linewidth=2.1, foreground="white")])
            # short connector when the label sits well away from its pie
            dist = np.hypot(lx - ax_, ly - ay)
            if dist > radii[name] * 2.6:
                ux, uy = (lx - ax_) / dist, (ly - ay) / dist
                ax.plot([ax_ + ux * radii[name], lx - ux * 1.4],
                        [ay + uy * radii[name], ly - uy * 1.4],
                        color="#8a8a8a", linewidth=0.42, alpha=0.8, zorder=7)

    # --- legends, in a band below the map -----------------------------------
    lineage_handles = [
        Line2D([], [], marker="o", linestyle="", markersize=7.5,
               markerfacecolor=colors[l], markeredgecolor="white",
               markeredgewidth=0.5, label=l)
        for l in all_lineages
    ]
    leg1 = fig.legend(handles=lineage_handles, title="Nextclade major lineage",
                      loc="upper left", bbox_to_anchor=(0.035, legend_y - 0.015),
                      frameon=True, fontsize=8.2, title_fontsize=9.2, ncol=lineage_ncol,
                      borderpad=0.8, handletextpad=0.6, columnspacing=1.3)

    size_handles = [
        Line2D([], [], marker="o", linestyle="", markerfacecolor="#cfcfcf",
               markeredgecolor="#4a4a4a", markeredgewidth=0.55,
               markersize=3.2 + 1.5 * r, label=f"{lbl} tips")
        for (_lo, _hi, r, lbl) in SIZE_BINS
    ]
    leg2 = fig.legend(handles=size_handles, title="Tips per deme (binned)",
                      loc="upper left", bbox_to_anchor=(0.40, legend_y - 0.015),
                      frameon=True, fontsize=8.2, title_fontsize=9.2,
                      borderpad=0.9, labelspacing=1.5, handletextpad=1.2)

    line_handles = [
        Line2D([], [], color=LINE_COLOR, linewidth=w * args.line_scale, alpha=a,
               label=f"{lbl} events  (n={drawn.get(lbl, 0)})")
        for (_lo, _hi, w, a, _arrow, lbl) in LINE_BINS
    ]
    leg3 = fig.legend(handles=line_handles, title="Inferred deme transitions (binned)",
                      loc="upper left", bbox_to_anchor=(0.635, legend_y - 0.015),
                      frameon=True, fontsize=8.2, title_fontsize=9.2,
                      borderpad=0.9, labelspacing=1.4, handlelength=3.6)

    legends = [leg1, leg2, leg3]
    if focus:
        focus_handles = [
            Line2D([], [], color=IN_COLOR, linewidth=3.0, alpha=FOCUS_MIN_ALPHA,
                   label=f"into {focus}  ({focus_counts['in']} events)"),
            Line2D([], [], color=OUT_COLOR, linewidth=3.0, alpha=FOCUS_MIN_ALPHA,
                   label=f"out of {focus}  ({focus_counts['out']} events)"),
        ]
        leg4 = fig.legend(handles=focus_handles, title=f"{focus} focus",
                          loc="upper left", bbox_to_anchor=(0.845, legend_y - 0.015),
                          frameon=True, fontsize=8.2, title_fontsize=9.2,
                          borderpad=0.9, labelspacing=1.4, handlelength=3.6)
        legends.append(leg4)

    for leg in legends:
        leg.get_frame().set_edgecolor("#bbbbbb")
        leg.get_frame().set_linewidth(0.6)

    fig.text(0.035, 0.022,
             "Pie area is binned, not proportional. Grey leader lines join a "
             "displaced pie to its true centroid (dot).\nArrowheads show "
             "inferred direction for bins of 3+ events.",
             fontsize=7.4, color="#555555", va="bottom")

    if args.title:
        ax.set_title(args.title, fontsize=14, pad=11)

    fig.savefig(args.output, dpi=330, facecolor="white")
    pdf = args.output.rsplit(".", 1)[0] + ".pdf"
    fig.savefig(pdf, facecolor="white")

    print(f"wrote {args.output} and {pdf}")
    print(f"  demes {len(per_deme)} | tips {sum(totals.values())} | "
          f"lineages {len(all_lineages)}")
    print(f"  lines drawn by bin: {dict(drawn)}  (total {sum(drawn.values())})")
    print(f"  demes by size bin: "
          f"{dict(Counter(bin_for(n, SIZE_BINS)[3] for n in totals.values()))}")
    if focus:
        print(f"  {focus}: {focus_counts['in']} events in (red), "
              f"{focus_counts['out']} events out (blue)")
    shifted = sum(1 for d in per_deme
                  if np.hypot(*(np.array(pos[d]) - np.array(true_pos[d]))) > leader_min)
    print(f"  pies displaced for legibility: {shifted}/{len(per_deme)}")


if __name__ == "__main__":
    main()

# Phylogeographic demes and transmission lines

## Where things live

| Path | What |
|---|---|
| `analysis/figures/` | the rendered figures (vector **PDF**, tracked in git) |
| `analysis/README.md` | this document |
| `phylogenetic/` | the Snakemake workflow, including all the deme changes |

The workflow changes are **not** in this folder and cannot be: Snakemake loads
`rules/deme.smk`, `defaults/deme_overrides.tsv` and `scripts/*.py` by fixed
relative path, and several of the changes are edits to upstream files
(`Snakefile`, `rules/export.smk`, `rules/merge_inputs.smk`,
`rules/prepare_sequences.smk`, `defaults/config_dengue.yaml`,
`defaults/color_orderings.tsv`). Moving them would stop the workflow running.

**Every path below is relative to `phylogenetic/`.**

The high-resolution PNGs are deliberately not tracked (60 MB); regenerate the
whole set, PNG and PDF, with:

```bash
bash phylogenetic/scripts/make-all-figures.sh
```

This fork adds a configurable **deme** trait so Auspice draws **transmission
lines** on the map between geographic units of our choosing, rather than only
between Nextstrain's stock `country` / `region` columns.

## What a deme is here

| Where | Deme |
|---|---|
| Inside Asia | the **country** (Thailand, Vietnam, India, …) |
| Outside Asia | the **region** (Africa, Europe, Oceania, South America) |
| Nextstrain's "North America" | split into **North America**, **Central America**, **Caribbean** |

Nextstrain's region vocabulary has only six values and no "Central America":
its `North America` lumps the USA/Canada/Mexico together with all of Central
America and the whole Caribbean. `defaults/deme_overrides.tsv` splits those
three apart. Trinidad and Tobago and French Guiana are deliberately *not*
reassigned — Nextstrain places them in South America and they stay there.

Records with a missing region (or a missing country inside Asia) get `deme=?`
and are dropped by `augur filter`.

The same applies to Nextstrain's **"Oceans and Seas"** region, which holds
records with `country = "Pacific Ocean"` / `"Atlantic Ocean"` and a *WHO
region* as the division (9 records in DENV-1, 3 in DENV-2, 4 in DENV-3; none
in DENV-4). These have no usable location, so they are treated as missing
rather than given a meaningless centroid in open ocean. See
`NON_GEOGRAPHIC_REGIONS` in `scripts/add-deme.py`.

## Why transmission lines need this

Auspice draws a transmission line between two map locations when the
*internal* nodes of the tree carry a geographic attribute and that attribute
has coordinates. Three things must therefore line up, and all three are wired
up here:

1. `augur traits` must reconstruct `deme` on internal nodes
   — `traits.traits_columns` in `defaults/config_dengue.yaml`.
2. `deme` must be listed in `geo_resolutions`
   — `rules/export.smk`, `prepare_auspice_config`.
3. **Every** deme must have coordinates, supplied to `augur export v2`
   via `--lat-longs` — `rules/deme.smk`, `deme_lat_longs`.

Point 3 is the easy one to get wrong: a deme with no coordinates is silently
dropped from the map along with every transmission line touching it. For that
reason `scripts/make-deme-lat-longs.py` **fails the build** rather than skip a
deme it cannot place.

## Files

| File | Purpose |
|---|---|
| `scripts/add-deme.py` | adds the `deme` column to the metadata |
| `scripts/make-deme-lat-longs.py` | emits a `deme`-keyed lat/longs TSV; errors on any deme it cannot place |
| `defaults/deme_overrides.tsv` | country → deme overrides (the Americas split) |
| `defaults/lat_longs_deme_extra.tsv` | centroids Augur lacks (Caribbean, Central America, East Timor, Borneo, …) |
| `rules/deme.smk` | the two Snakemake rules |
| `build-configs/denv4-deme-phylogeo/config.yaml` | the DENV-4 pilot |

`add_deme` deliberately produces the canonical
`results/{serotype}/metadata.tsv` that every downstream rule already reads
(`merge_inputs.smk` now writes `metadata_raw.tsv` instead), so filter, refine,
traits, export and tip_frequencies all see the `deme` column unchanged.

## The two datasets

Each build writes **two** Auspice datasets from the same tree and the same
`augur traits` run:

| Dataset | Opens on |
|---|---|
| `auspice/dengue_{serotype}_{gene}.json` | the deme map (country in Asia, region elsewhere) |
| `auspice/dengue_{serotype}_{gene}_country.json` | the standard country-everywhere map |

Both contain `deme`, `country` and `region` as geo-resolutions, so either can
be switched to the other view from the Auspice sidebar; they differ only in
which resolution they open on.

## Publication figures

`scripts/plot-deme-map.py` renders a static map from any of the exported
Auspice JSONs, for use in a paper or slide. It exists because Auspice scales
deme circles continuously by tip count, which is unreadable when one deme has
880 tips and another has 1.

* **Pies**, one per deme, sliced by Nextclade major lineage, reusing the
  dataset's own colour scale so the figure matches Auspice.
* **Three binned pie sizes** — 1-9, 10-99, 100+ tips. Deliberately *not*
  area-proportional.
* **Three binned line widths** — 1-2, 3-9, 10+ inferred transitions, with
  exaggerated contrast; arrowheads on the 3+ bins only.
* **Focus highlighting** — `--focus-deme` (default `Singapore`) draws imports
  red and exports blue.
* Crowded pies are pushed apart by a short force-directed relaxation and tied
  back to their true centroid by a grey leader line and dot; labels are placed
  by testing candidate offsets against the real rendered text extents.

```bash
python scripts/plot-deme-map.py \
    --auspice-json auspice/dengue_denv4_E.json \
    --output figures/denv4_E_deme_map.png \
    --title "DENV-4 global phylogeography (E gene)"
```

Pass `--deme-key country` for the country-everywhere version, `--focus-deme ""`
to disable highlighting, and `--min-line-count N` to drop weakly supported
lines. Each run writes both a 330 dpi PNG and a vector PDF. The bins are fixed
constants at the top of the script so figures stay comparable across builds;
a bin can legitimately come out empty (the genome build has no 10+ line).

To render every figure for every serotype that has been built:

```bash
bash scripts/make-all-figures.sh [focus-deme]     # default focus: Singapore
```

That writes three figures per serotype x gene (global deme map, global country
map, South & Southeast Asia crop) — 24 in all for the four serotypes.

The plotting dependencies (cartopy/matplotlib) are **not** in the workflow
environment; they live in a separate `denvfig` conda env.

One caveat on the **country-everywhere** maps: DENV-1 and DENV-2 reach 96-99
country demes, and over half the pies have to be displaced to stop them
overlapping, so the drawn positions drift well away from the true centroids.
The deme maps (27-35 demes, 0-12 displaced) are the more trustworthy geography;
treat the country versions as supplementary.

## Running it

```bash
cd phylogenetic
snakemake --configfile build-configs/denv4-deme-phylogeo/config.yaml -c 32 \
    --set-threads align=32 tree=32
```

To run another serotype, copy the pilot config and change `serotypes`
(`denv1` / `denv2` / `denv3` / `all`), or add `E` to `genes`. The deme logic is
serotype-independent; `defaults/color_orderings.tsv` already lists demes for
all four serotypes, and `assign-colors.py` drops any value absent from the
metadata.

#!/usr/bin/env bash
# Render every phylogeographic figure for every serotype that has been built.
#
# For each serotype x gene we produce:
#   <serotype>_<gene>_deme_map        global, demes = country in Asia / region outside
#   <serotype>_<gene>_country_map     global, demes = country everywhere
#   <serotype>_<gene>_deme_map_sea    cropped to South & Southeast Asia
#
# Every map highlights transitions into the focus deme (default Singapore) in
# red and out of it in blue. Each run writes a 330 dpi PNG and a vector PDF.
#
# Usage:  bash scripts/make-all-figures.sh [focus-deme]
set -uo pipefail

PYTHON=${PYTHON:-/hpc/home/ansellim/micromamba/envs/denvfig/bin/python}
FOCUS=${1:-Singapore}
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
cd "$HERE"
# Figures land in the repo-level analysis/ folder, not under phylogenetic/.
OUT=${OUT:-"$(cd "$HERE/.." && pwd)/analysis/figures"}
mkdir -p "$OUT"

# South & Southeast Asia crop. Pie radii are in degrees, so a cropped map
# needs a much smaller scale than the global one.
# Pies and fonts are both 1.5x the original crop values: big enough to read,
# small enough that the legend band stays shallow and the map keeps most of
# the figure. The legend band measures and re-packs itself to suit.
SEA_EXTENT=(72 142 -12 38)
SEA_ARGS=(--pie-scale 0.24 --line-scale 1.15 --label-fontsize 8.5 --font-scale 1.5)

declare -A GENE_LABEL=([genome]="whole genome" [E]="E gene")

failures=0
for serotype in denv1 denv2 denv3 denv4; do
  for gene in genome E; do
    json="auspice/dengue_${serotype}_${gene}.json"
    [[ -f $json ]] || { echo "SKIP  $json (not built)"; continue; }
    upper=${serotype/denv/DENV-}
    glabel=${GENE_LABEL[$gene]}

    echo "=== $upper $gene ==="
    $PYTHON -I scripts/plot-deme-map.py --auspice-json "$json" \
      --focus-deme "$FOCUS" \
      --output "$OUT/${serotype}_${gene}_deme_map.png" \
      --title "$upper global phylogeography ($glabel): demes = country in Asia, region elsewhere" \
      2>&1 | grep -vi downloadwarning | sed 's/^/    /' || ((failures++))

    $PYTHON -I scripts/plot-deme-map.py --auspice-json "$json" \
      --focus-deme "$FOCUS" --deme-key country \
      --output "$OUT/${serotype}_${gene}_country_map.png" \
      --title "$upper global phylogeography ($glabel): demes = country everywhere" \
      2>&1 | grep -vi downloadwarning | sed 's/^/    /' || ((failures++))

    $PYTHON -I scripts/plot-deme-map.py --auspice-json "$json" \
      --focus-deme "$FOCUS" \
      --extent "${SEA_EXTENT[@]}" "${SEA_ARGS[@]}" \
      --output "$OUT/${serotype}_${gene}_deme_map_sea.png" \
      --title "$upper phylogeography ($glabel), South & Southeast Asia: demes = country in Asia, region elsewhere" \
      2>&1 | grep -vi downloadwarning | sed 's/^/    /' || ((failures++))
  done
done

echo
echo "=== wrote $(ls "$OUT"/*.png 2>/dev/null | wc -l) PNGs / $(ls "$OUT"/*.pdf 2>/dev/null | wc -l) PDFs into $OUT; $failures failure(s) ==="

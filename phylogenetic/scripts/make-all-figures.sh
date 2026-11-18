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
mkdir -p figures

# South & Southeast Asia crop. Pie radii are in degrees, so a cropped map
# needs a much smaller scale than the global one.
# Pies are 4x the original crop size and every font/legend is 3x, for use at
# slide/poster scale; the legend band measures and re-packs itself to suit.
SEA_EXTENT=(72 142 -12 38)
SEA_ARGS=(--pie-scale 0.64 --line-scale 1.15 --label-fontsize 8.5 --font-scale 3)

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
      --output "figures/${serotype}_${gene}_deme_map.png" \
      --title "$upper global phylogeography ($glabel): demes = country in Asia, region elsewhere" \
      2>&1 | grep -vi downloadwarning | sed 's/^/    /' || ((failures++))

    $PYTHON -I scripts/plot-deme-map.py --auspice-json "$json" \
      --focus-deme "$FOCUS" --deme-key country \
      --output "figures/${serotype}_${gene}_country_map.png" \
      --title "$upper global phylogeography ($glabel): demes = country everywhere" \
      2>&1 | grep -vi downloadwarning | sed 's/^/    /' || ((failures++))

    $PYTHON -I scripts/plot-deme-map.py --auspice-json "$json" \
      --focus-deme "$FOCUS" \
      --extent "${SEA_EXTENT[@]}" "${SEA_ARGS[@]}" \
      --output "figures/${serotype}_${gene}_deme_map_sea.png" \
      --title "$upper phylogeography ($glabel), South & Southeast Asia: demes = country in Asia, region elsewhere" \
      2>&1 | grep -vi downloadwarning | sed 's/^/    /' || ((failures++))
  done
done

echo
echo "=== wrote $(ls figures/*.png 2>/dev/null | wc -l) PNGs / $(ls figures/*.pdf 2>/dev/null | wc -l) PDFs; $failures failure(s) ==="

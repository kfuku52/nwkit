#!/usr/bin/env bash
set -euo pipefail
# Run from the NWKIT repository root. All data in this example are synthetic.
example=examples/radte/species-uncertainty
result_dir=${1:-output/pdf/species-uncertainty}
mkdir -p "$result_dir"
common=(--gene-tree "$example/gene.nwk" --species-tree "$example/species.nwk"
        --species-map-tsv "$example/species-map.tsv" --reconcile lca --max-age 30
        --rate-sd 0.3 --species-node-intervals-tsv "$example/species-intervals.tsv")
python -m nwkit.cli radte "${common[@]}" --uncertainty profile \
  --out-prefix "$result_dir/fixed" --figure-out "$result_dir/fixed.pdf"
python -m nwkit.cli radte "${common[@]}" --uncertainty profile \
  --species-node-bounds-tsv "$example/species-bounds.tsv" \
  --out-prefix "$result_dir/bounded" --figure-out "$result_dir/bounded.pdf"
python -m nwkit.cli radte "${common[@]}" --uncertainty input-ensemble \
  --species-tree-ensemble "$example/species-samples.nwk" \
  --ensemble-within-uncertainty profile \
  --out-prefix "$result_dir/ensemble" --figure-out "$result_dir/ensemble.pdf"
python -m nwkit.cli radte-compare --fixed-prefix "$result_dir/fixed" \
  --bounded-prefix "$result_dir/bounded" --ensemble-prefix "$result_dir/ensemble" \
  --species-tree "$example/species.nwk" --out-prefix "$result_dir/comparison"

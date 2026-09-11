#!/usr/bin/env bash
set -euo pipefail
repo=/Users/kf/repos/nwkit
review=/work/reviews/joint-engine-optimization-2026-09-11
runtime=(docker run --rm --entrypoint python
  -e OPENBLAS_NUM_THREADS=1 -e OMP_NUM_THREADS=1 -e MKL_NUM_THREADS=1
  -v "$repo:/work" -v "$repo/tools:/bench/tools:ro" -w /tmp)
for engine in baseline optimized; do
  library=source
  if [[ "$engine" == baseline ]]; then library=source-before; fi
  "${runtime[@]}" -e "PYTHONPATH=$review/$library" c8c69075b5d5 \
    /bench/tools/benchmark_joint_engines.py --engine "$engine" \
    --output "$review/kernel-$engine.json"
done
"${runtime[@]}" -e "PYTHONPATH=$review/source-before" c8c69075b5d5 \
  /bench/tools/benchmark_ten_shifts_missing.py \
  --job '{"truth":"shared","mode":"shared","replicate":0,"traits":2,"missing_rate":0.2,"timeout":1800,"part":"timing"}' \
  --output "$review/baseline-current-shared-0.json"
"${runtime[@]}" -e "PYTHONPATH=$review/source" c8c69075b5d5 \
  /bench/tools/run_ten_shift_benchmark.py --output "$review/timing-final-v2" \
  --part timing --replicates 1 --workers 1 --timeout 1800 --traits 2
cp "$repo/reviews/joint-engine-optimization-2026-09-11/timing-final-v2/shared-0-shared.json" \
   "$repo/reviews/joint-engine-optimization-2026-09-11/candidate-current-shared-0.json"
"${runtime[@]}" -e "PYTHONPATH=$review/source" c8c69075b5d5 \
  /bench/tools/run_ten_shift_benchmark.py --output "$review/accuracy-final" \
  --part accuracy --replicates 5 --workers 8 --timeout 1800 --traits 2

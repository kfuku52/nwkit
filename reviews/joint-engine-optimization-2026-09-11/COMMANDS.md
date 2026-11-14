# Reproduction

The baseline is the exact pre-implementation package in `source-before/nwkit/`; the final candidate is in `source/nwkit/`. Both run in GeneGalleon image `c8c69075b5d5` with the same dependencies, on Apple M2 Max with one BLAS/OpenMP thread. Native source hashes and dependency versions are recorded in the kernel and search manifests.

Mount benchmark helpers at `/bench/tools`, separately from either native package. A legacy helper prepends its parent to Python's search path; separating this mount ensures it cannot silently select a different checkout.

Baseline kernel, with a new output filename:

```bash
docker run --rm --entrypoint python \
  -e OPENBLAS_NUM_THREADS=1 -e OMP_NUM_THREADS=1 -e MKL_NUM_THREADS=1 \
  -e PYTHONPATH=/work/reviews/joint-engine-optimization-2026-09-11/source-before \
  -v /Users/kf/repos/nwkit:/work \
  -v /Users/kf/repos/nwkit/reviews/joint-engine-optimization-2026-09-11/source/tools:/bench/tools:ro \
  -w /tmp c8c69075b5d5 /bench/tools/benchmark_joint_engines.py \
  --engine baseline --output /work/reviews/joint-engine-optimization-2026-09-11/reproduction-kernel-baseline.json
```

After it finishes, repeat with `--engine optimized` and `PYTHONPATH` ending in `/source`, using another output filename. Each case has three warmups and seven batches of five evaluations. Outputs include likelihood, coefficients, coefficient covariance, all tip predictions, timings, and process peak RSS.

Candidate search replay:

```bash
docker run --rm --entrypoint python \
  -e OPENBLAS_NUM_THREADS=1 -e OMP_NUM_THREADS=1 -e MKL_NUM_THREADS=1 \
  -e PYTHONPATH=/work/reviews/joint-engine-optimization-2026-09-11/source \
  -v /Users/kf/repos/nwkit:/work \
  -v /Users/kf/repos/nwkit/reviews/joint-engine-optimization-2026-09-11/source/tools:/bench/tools:ro \
  -w /tmp c8c69075b5d5 /bench/tools/run_ten_shift_benchmark.py \
  --output /work/reviews/joint-engine-optimization-2026-09-11/reproduction-timing \
  --part timing --replicates 1 --workers 1 --timeout 1800 --traits 2
```

For detection outcomes, use `--part accuracy --replicates 5 --workers 8` with a separate directory. Parallel accuracy times are not used for speed ratios. The data generator, seeds, observation errors, missingness, and search budgets match the earlier report in `../ten-shifts-errors-missing-2026-09-11/`.

For the single pre-implementation search, use the baseline `PYTHONPATH` and call `benchmark_ten_shifts_missing.py` directly:

```text
--job '{"truth":"shared","mode":"shared","replicate":0,"traits":2,"missing_rate":0.2,"timeout":1800,"part":"timing"}' --output /work/reviews/joint-engine-optimization-2026-09-11/reproduction-baseline-shared-0.json
```

The executed command sequence is preserved in `run_isolated.sh` and its log in `isolated-final.log`. The `candidate-current-shared-0.json` record is a copy of the shared/shared replicate-0 result in `timing-final-v2/`, not another independent timing repetition. Final detection outcomes are in `accuracy-final/`.

Rebuild the audit from this checkout with:

```bash
python tools/validate_joint_engine_benchmark.py \
  reviews/joint-engine-optimization-2026-09-11 \
  reviews/ten-shifts-errors-missing-2026-09-11
```

`source-v1/`, `accuracy/`, `timing-final-v1/` and files/directories prefixed `preflight` are development evidence, excluded from final measurements. They document the squared-noise zero derivative, finite-difference truncation bias, and the earlier eager geometry preparation. Kernel provenance validation explicitly checks the pre-implementation and final source trees.

Validation used a temporary container venv with system site packages for quality/build tools. The scientific runtime dependencies were retained. `check-full.log`, `check-focused.log`, `baseline-failures.log`, `check-distribution.log` and `check-package-final.log` record what ran. `package_check.py` stages distribution checks in a separate container directory to preserve existing workspace build outputs.

# Native SHIFT: 1,000 tips and up to 100 shifts

This development target is separate from the production-adoption gates in
[NATIVE_SHIFT_VALIDATION.md](NATIVE_SHIFT_VALIDATION.md). A successful fixed-
covariance search does not establish calibrated selection or support-bootstrap
runtime, statistical error control, or replacement readiness for GeneGalleon.

## First implementation step

Candidate likelihood profiles now retain an orthogonal reduction of the joint
whitened response and candidate design. Subsequent layouts operate on at most
`candidate_pool + 1` rows instead of the full observed-tip count. Including the
response preserves its component outside the candidate design span. The method
does not use normal equations. Numerical rank checks retain the original observed-
tip count in their tolerance, including with missing observations.

The reduction changes the representation used to rank candidates; it does not
change candidate budgets, covariance models, final unpenalized fits, or the
bootstrap procedure. Full refits and structural checks still have their own costs. Structural checks
are shared between traits only when both alpha and the observed-tip mask match;
different masks or alpha values retain their separate checks.

## Reproducible development measurements

The search benchmark accepts explicit budgets, so its requested shift count need
not fit the small default candidate pool. For example:

```sh
OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONPATH=. \
python tools/benchmark_native_search.py \
  --tips 1000 --traits 2 --shifts 100 --candidate-pool 128 \
  --refit-budget 220 --screening-budget 20000 --beam-width 1 \
  --no-convergence --output /tmp/native-search-100.json
```

This is an explicit **fixed-covariance, distinct-regime, single-search** workload.
The output records `largest_fitted_shift_count`; a requested cap alone is not
proof that the search reached it. `--convergence` enables shared-regime search,
and `--fit-covariance` enables covariance estimation. Measure these workloads
separately; neither is implied by the command above. Existing benchmark defaults
remain unchanged. Output paths must be new.

For a repeated candidate-ranking measurement with 100 shifts, 128 candidate
branches, two traits and partial missingness:

```sh
OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONPATH=. \
python tools/benchmark_native_profiles.py --output /tmp/native-profiles.json
```

Use `--shape pectinate` for the unbalanced-tree case and `--traits 4` for four
traits. Candidate terminal branches are observed in every trait; missing data
remain on other tips, so every timed candidate must have a finite score. The benchmark records setup
separately, one warmup and three scoring repetitions, every likelihood for
comparison, and peak process RSS. Compare the same script, inputs, dependencies
and thread limits between source snapshots.

## Remaining scale work

- Structural rank checks and final covariance estimation remain outside the
  compressed candidate kernel. Profile full searches before selecting the next
  optimization.
- Screening constructs dense tip-by-branch matrices. Its memory budget is a
  screening guard, not a bound on peak RAM of the complete analysis.
- Candidate, scoring and refit budgets must permit the requested complexity.
  Shared-regime search has many more candidates than distinct-regime search.
- Calibration repeats the complete configured search; support bootstrap nests
  calibrated selections. Single-search timings cannot be presented as complete
  workflow timings.
- Validate 1/4-trait, balanced/pectinate and estimated-covariance workloads,
  convergence failures, recovery and false selection before production adoption.

Measured results and verification are recorded in
[the development report](reviews/native-scale-2026-09-10/REPORT.md).

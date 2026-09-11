# Ten shifts with observation error and missing coordinates

Frozen formal design, 2026-09-11. Resource preflight uses replicate 999 and is excluded from formal estimates. Formal design was frozen after observing preflight per-layout resource costs, before running formal searches.

- Frozen GeneGalleon Docker image: sha256:c8c69075b5d591499bf3195ec912ef5b5e1aa1cd163f03864639ed6a669e9429. Archived native source from the preceding comparison. Docker results do not establish SIF compatibility.
- 100-tip balanced rooted ultrametric tree, height 1. Sixteen disjoint clades of 6 or 7 tips are eligible. Ten are sampled without replacement for each dataset; the other six supply background observations.
- Two traits, full diffusion covariance with correlation 0.8. Marginal latent tip variance 1 for each trait. Generating alpha either [1,1] or [0.25,4]. Holding diffusion correlation does not hold tip correlation when alpha differs.
- Ten independent regime shifts, no convergence. Each clade receives a random unit direction with expected tip displacement of Euclidean length 3 relative to background. The same clades, directions, mask and random seed are used for both generating-alpha scenarios. These are strong shifts in small clades; no claim of general power across effect sizes or trees.
- Known independent sampling SE 0.1 and extra independent measurement variance 0.04 per trait. Fit estimates extra measurement variance as well as alpha, full process covariance and means.
- MCAR missingness: independent probability 0.2 per tip/trait coordinate. No imputation; mask passed to simulation and likelihood. Observed fraction is recorded per dataset.
- Both fitting alpha models receive identical observations and missingness, with recorded data hashes. Alpha is estimated, never fixed to truth. Exact zero and infinity modes and optimizer acceptance follow native defaults without relaxation.
- Search maximum 12 shifts, candidate pool 32, refit budget 26, beam 2, screening budget 10000, lasso iterations 200. Truth is not injected. Equal budgets permit different model-selected candidate pools. This is a budget-limited search, not exhaustive inference.
- AIC and BIC select among the search's retained maximum-likelihood layouts. These are information-criterion detection results, not bootstrap-calibrated significance or controlled 5% false-positive tests.
- Exact branch TP/FP/FN, precision, recall, F1, estimated count and exact-set recovery. All true sets contain ten shifts; null false-positive rate is not assessed. Failure and timeout rates stay in the denominator for operational reporting; paired accuracy explicitly identifies successfully completed pairs.
- Statistical jobs may execute in parallel and their elapsed times are not isolated speed benchmarks. Separate serial runs provide paired timing; profiled runs are excluded from timing ratios. Single-thread BLAS/OpenMP/MKL.

## Formal replication and resource limits

- Accuracy: five independently generated datasets per generating-alpha scenario (seeds with replicate 0 through 4), ten datasets and twenty searches total. This is an exploratory benchmark; five datasets per scenario do not support precise population-level power estimates. Ten true shifts per dataset are not ten independent simulation replicates.
- Accuracy controller: eight concurrent subprocesses, BLAS/OpenMP/MKL each one thread. Per-search wall-time limit 1800 seconds. Include every outcome, and do not rerun failed cases selectively.
- Isolated timing: replay replicate 0 from both generating-alpha scenarios, both fitting models, one search each, serial with no profilers or competing benchmark jobs. Four end-to-end times on two paired datasets; not a repeated-run median or population-level speed estimate. Preflight has warmed libraries and the filesystem; fixture generation excluded from search time.
- Kernel proposal microbenchmark: three warmups then seven batches of ten likelihood evaluations, both generating-alpha datasets, null and true ten-shift layouts, shared/unequal/zero/infinite alpha. Compare log likelihood, mean coefficients and their covariance to 1e-8. This validates tested likelihood evaluations, not optimizer/search equivalence or total search speedup.

## Interpretation details

NWKIT joint BIC uses the number of tip vectors with at least one observed coordinate as sample size. Its parameter count includes shift locations, regime means, covariance, alpha and estimated measurement variance. We use this existing implementation without redefining BIC for this benchmark.

The dense-kernel proposal computes the same tested log likelihood, mean coefficients, coefficient covariance and predictions, but is not a drop-in production implementation: full interface validation, all supported root models, and optimization/search equivalence remain to be implemented and verified. Its quadratic storage and cubic dense factorization are unsuitable as an unconditional replacement for tree pruning. The timing comparison is restricted to this small workload and excludes one-time tree shared-time preparation, whose time is separately recorded.

Profiling used separate replicate 999 and a 600-second cap. Both profiles are partial searches, with profiling overhead and overlapping concurrent work. Their cumulative times identify expensive call paths; overlapping categories must not be added, and their wall times are not unprofiled model-speed comparisons.

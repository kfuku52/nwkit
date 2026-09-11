# Native evolutionary covariance implementation validation (2026-09-11)

This report concerns the actual NWKIT implementation described in
[SHIFT_COVARIANCE.md](../../SHIFT_COVARIANCE.md). No commit, push or release was
performed. GeneGalleon forwards opt-in `full` covariance and `shared` alpha; its
defaults remain diagonal covariance and trait-specific alpha.

## Runtime and scope

Apple M2 Max host; Linux aarch64 GeneGalleon Docker runtime, Python 3.12.14,
NumPy 1.26.4 and SciPy 1.17.1. BLAS/OpenMP threads were held at one.
Final local image: `local/genegalleon:nwkit-shift-covariance-dev`, image ID
`sha256:c8c69075b5d591499bf3195ec912ef5b5e1aa1cd163f03864639ed6a669e9429`.
It installs the local NWKIT source on the existing GeneGalleon development
runtime. This is Docker validation, not Apptainer/SIF validation.

The implementation fingerprint includes all package source files, numerical
library versions, platform and thread environment. Read the exact fingerprints
from each JSON result; calibration preceded the final lazy CLI registration and
compatibility-alias corrections. Those corrections do not change the fitting or
simulation algorithms. Timing uses the installed image without a source mount,
so concurrent repository work cannot alter that run.

## Correctness and integration

- Independent dense Gaussian oracle versus vector pruning; diagonal nesting;
  separable profile versus general likelihood and numerical ML; transformed units,
  trait permutations, covariance parameter counts and unbounded-layout exclusion.
- Correlated simulation: 30,000 independent draws agree with independently
  reconstructed means/covariance, for fixed and stationary roots; known and
  correlated extra observation error, missing masks and reproducible seeds.
- Joint candidate screening and profiled scores, actual CLI model/simulation
  round trips, complete-search bootstrap and output-transaction checks.
- Final installed NWKIT: **82 passed** across CLI contracts and the three new test
  modules (`final-installed-tests.log`). Two marker warnings arose because the
  isolated test mount omitted project pytest configuration.
- GeneGalleon installed-package adapter/config tests: **46 passed**
  (`gg-installed-tests.log`); after final CLI registration changes the affected
  native adapter suite again passed **12/12** (`gg-final-tests.log`). No source
  `PYTHONPATH` override was used for these installed-package checks.
- 100-tip, two-trait general-path smoke with missing coordinates, known sampling
  errors, estimated trait-specific alpha, full process covariance and estimated
  additional measurement variance converged with finite likelihood
  (`general-smoke.json`, runnable `general_smoke.py`). This is a convergence smoke,
  not parameter-recovery or selection-calibration evidence.
- Broad non-slow suite: **4086 passed, 57 skipped, 13 deselected**, initially six
  failures. Four CLI lazy-import failures were fixed and covered by the final
  installed tests. The other two failures reproduce with the pre-change package
  snapshot (`baseline-failures.log`): root `candidates_out` compatibility alias,
  and archived null-contract seeded replay. They remain unresolved here.
- Ruff lint and mypy (226 source files) passed; maintainability hard limits passed.
  Global formatting still reports only the unrelated
  `tests/test_reconciliation_exports.py`. Changed files are formatted. GeneGalleon
  shell syntax and its eight-entrypoint config schema validation passed.
- Wheel/sdist reproducibility and contents checks passed in an isolated copy
  (`dist-check.log`), preserving the working repository's build artifacts.
  The new covariance guide is included in the source distribution.

This was not a full slow suite, coverage/security audit, SIF run or broad
statistical adoption study. No such certification is implied.

## Fitted-null bootstrap pilot

`calibration.json` records 90 paired jobs (180 method results), with no failed
jobs. Each method fits its own null means/covariance, generates from that fitted
null and repeats all 199 zero/one-shift candidate fits for each of 19 bootstrap
draws. Alpha is fixed at 1 throughout. This is **not** the earlier oracle-null
prototype and does **not** validate estimated-alpha selection.

The exactly ultrametric, balanced tree has 100 tips. Generating evolutionary
correlation is 0.8, with 2, 5 or 10 traits. Data are complete with no observation
errors. Non-null clades contain 10–30 tips and have tip displacement norm 2,
aligned across traits or opposed in the first two. There are only 10 evaluation
replicates per cell. The minimum bootstrap p-value is 1/20 = 0.05. Counts below
are descriptive pilot observations, not precise rates or error-control guarantees.

| Traits | Null: any shift, diagonal / full | Aligned: exact branch, diagonal / full | Opposed: exact branch, diagonal / full |
|---:|---:|---:|---:|
| 2 | 2/10 / 0/10 | 3/10 / 0/10 | 4/10 / 10/10 |
| 5 | 9/10 / 0/10 | 2/10 / 0/10 | 2/10 / 10/10 |
| 10 | 10/10 / 1/10 | 2/10 / 0/10 | 0/10 / 10/10 |

The diagonal fitted null omits the generating cross-trait covariance, so its
bootstrap is misspecified in this experiment. Its frequent null detections must
not be interpreted as useful power. The full model's opposed-direction detections
are promising, while its aligned-direction detections here are weak. More
replicates, larger bootstrap budgets, estimated alpha, other trees/correlations,
multiple shifts and observation-error scenarios are needed for adoption claims.

Reproduction from the NWKIT root in the stated runtime:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
python tools/benchmark_shift_covariance.py --part calibration \
  --output reviews/trait-covariance-2026-09-11/calibration.json \
  --tips 100 --traits 2,5,10 --replicates 10 \
  --calibration-replicates 19 --workers 4
```

## Timing and numerical equivalence

Warmup plus three sequential measurements per method; medians in seconds. No calibration or test jobs ran concurrently. Both search methods use shared alpha fixed at 1, exact zero/one-shift enumeration (199 candidates), and refit covariance and means. These are complete, error-free 100-tip data; timings exclude simulation and imports.

| Traits | Diagonal search | Full search | Full / diagonal | Full fixed-layout ML, exact profile |
|---:|---:|---:|---:|---:|
| 2 | 0.1750 | 0.1782 | 1.018 | 0.001126 |
| 5 | 0.2955 | 0.3049 | 1.032 | 0.001586 |
| 10 | 0.5422 | 0.5330 | 0.983 | 0.002673 |

Small differences around a ratio of one are timing noise, not evidence that full covariance is intrinsically faster.

At two traits, general numerical ML took 2.219 s versus 0.001126 s for the exact profile, with absolute likelihood difference 4.55e-13 and maximum covariance-entry difference 1.03e-07. This comparison includes removing numerical covariance optimization; it is not a pure pruning-kernel speed comparison.

A first repeated numerical-ML run was manually stopped while processing five traits; its unfinished log is `timing-general-unfinished.log`. No completed five/ten-trait numerical-optimization time is claimed. The final tool defaults to `--pruning-max-traits 2`. At five/ten traits it evaluates the exact profile solution using general vector pruning and checks likelihood agreement, without reoptimizing covariance. These likelihood errors are recorded in `timing.json`.

The JSON also separately records fixed-layout fits estimating shared alpha. It does not measure exhaustive estimated-alpha search, observation-error fits or bootstrap wall time. Peak RSS is whole-process high-water usage, not a paired per-method memory benchmark.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
python tools/benchmark_shift_covariance.py --part timing \
  --output reviews/trait-covariance-2026-09-11/timing.json \
  --tips 100 --traits 2,5,10 --repeats 3 --pruning-max-traits 2
```

The matching package Python sources are preserved in `calibration-source.tar.gz` and `timing-source.tar.gz`. The calibration image `sha256:18ad8304f1341c5ff67393d54c1452d6a4914388761a423161391ab8ee560c14`, with all three thread variables set to one, was verified to reproduce the exact saved implementation fingerprint.

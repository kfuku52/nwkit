# Multivariate ASR pruning measurements

The vector-pruning backend removes the dense fitter's 1,000-observed-coordinate
limit for missing/noisy data. It trades Python computation time for bounded
per-node matrix storage; it is not a general speed optimization. Small supported
fits retain the dense backend. Complete exact MV-BM retains its specialized path.

## Reproducible workload

Run from the repository root, using the same environment for both backends:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONPATH=. \
  python tools/benchmark_vector_asr.py --backend dense --model BM --tips 256
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONPATH=. \
  python tools/benchmark_vector_asr.py --backend pruning --model BM --tips 256
```

Repeat with `--model OU`. The fixed seed is 311; the balanced tree has two traits,
every fifth second-trait value missing and known SEs 0.1/0.2. Every run includes
all-node reconstruction. OU fixes alpha at 0.4. The harness performs one warmup
and three measured fits, reports the median and process peak RSS (native units;
bytes on macOS), and prints fitted covariance and likelihood for equivalence
checks. Peak RSS includes imports and warmup; it is not isolated allocation size.

## Exploratory results

Recorded during development on macOS 26.6.2 x86_64, Python 3.10.14, NumPy 1.26.4,
SciPy 1.15.2, with one BLAS/OpenMP thread. Other development work was active, so
timings are exploratory and should not be treated as stable performance ratios.

| Workload | Dense median seconds | Pruning median seconds | Dense peak MiB | Pruning peak MiB |
| --- | ---: | ---: | ---: | ---: |
| BM, 256 tips / 460 observed coordinates | 4.303 | 93.262 | 202.7 | 145.7 |
| OU, 256 tips / 460 observed coordinates | 4.829 | 144.934 | 210.1 | 146.3 |

The pruning likelihood differences were below `3e-11` (BM) and `1e-10` (OU).
Covariance entries agreed within `3e-6` absolute error. These measurements show
the substantial time cost of the generic Python pruning path, despite lower
observed RSS. Re-measure in an otherwise idle environment before making speed
or memory guarantees for another workload or platform.

A separate real MV-BM fit at 560 tips / 1,008 observed coordinates (same fixture)
completed with `fit_status=ok`, restricted log likelihood `-1795.0328955294872`
and finite mean/covariance summaries at all 1,119 nodes. This checks behavior
beyond the old dense limit; it is not a timed speed comparison. Missing-data,
measurement-error, root-prior and dense numerical-oracle checks live in
`tests/test_vector_gaussian.py`, `tests/test_vector_fit.py`, and
`tests/test_vector_simulation.py`.

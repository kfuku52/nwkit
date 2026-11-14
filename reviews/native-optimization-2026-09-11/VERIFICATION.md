# Verification at measurement time

See [the subsequent audit](AUDIT.md) for the later full checks and fixes.

- `OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 python tools/check.py quick -- tests/test_shift_screen_blocks.py tests/test_shift_native_*.py tests/test_gaussian_whitening.py`: **218 passed**, 34.32 seconds. The quick runner excludes `slow`; none of these files declares a slow marker.
- The same command passed repository-wide Ruff lint, format (481 files) and incremental mypy (220 source files).
- `python tools/check_maintainability.py`: all hard limits respected; 3,290 functions, mean 6.78, maximum 50 (existing exception). Other pre-existing baseline warnings were reported.
- Added 11 cases: eight blocked/full-matrix comparisons across alpha 0/0.7/1000/infinity and standardized/unstandardized designs, one selected-column order check, and two dense-covariance tests that change numerical parameters and observation order while reusing topology.
- Existing independent dense GLS, OU boundary, missingness/measurement-error, shared-regime, calibration, support, CLI and replay tests are included.
- The Python environment emits a Requests dependency warning about urllib3/chardet/charset_normalizer. Tests passed with that warning; dependencies were not changed.
- The entire repository test suite, distribution builds, other Python/OS/BLAS combinations, and statistical production-adoption studies were not run. This is a native SHIFT optimization, not a push/release or a production-adoption claim.

Performance measurements run after these tests, sequentially in fresh subprocesses, using frozen source copies. Numerical comparisons and per-run records are retained in `final-measurements/`.

# Integrated local validation — 2026-09-10

The reviewed NWKIT changes were combined on `master` for one local commit. They include experimental native SHIFT inference and its fixed-covariance scaling improvements, ASR output/input safeguards, iterative PCA tree copying, the RADTE boundary-initialization repair, and the associated scientific validation tools, data, plots, and documentation. Native SHIFT remains experimental; this validation does not establish its production statistical gates.

## Verified source

The integration copy matched all 1,053 indexed files byte-for-byte before distribution checks (index tree `fc350895934fe3cbb958a96747754c81f2465951`, based on `e2465b0`). [Python source hashes](source-hashes.json) record the tested files. This report and its logs were added afterward. The RADTE initializer helper refactor was incorporated before pytest imported the source; its Ruff, format, mypy and Bandit checks also passed separately. The documentation-only RADTE validation update and `.gitattributes` were synchronized before the complete file comparison.

## Checks

- `OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 python tools/check.py full`: **passed**. Ruff lint and format (473 files), mypy (217 source files), dependency consistency, Bandit and dependency audit passed. The full test suite had **3,902 passed, 79 skipped, 6 warnings in 910.99 seconds**. Coverage was **85%**. Maintainability hard limits passed (3,276 functions; mean complexity 6.77, maximum 50). See [full log](full-check.log) and [supplemental RADTE type check](radte-type.log).
- `SOURCE_DATE_EPOCH=1789010333 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 python tools/check.py dist`: **passed**. Direct and sdist-built wheels, package contents and metadata, and distribution reproducibility were checked. See [distribution log](dist-check.log).
- `git diff --cached --check`: passed. Whitespace exceptions preserve raw logs, benchmark patch context, and the source FASTA fixture.

Checks ran locally with Python 3.10.14 on the same host used for the recorded performance work. The dependency audit emitted requests/cache warnings but completed successfully with no known vulnerabilities. Test skips and existing complexity warnings are retained in the logs. This is not the hosted CI matrix or a new performance measurement. Builds ran in the integration copy, preserving the working repository's build directories.

The earlier [native scaling report](../native-scale-2026-09-10/REPORT.md) describes checks at that stage. Its then-outstanding formatting issues are resolved by this integrated validation. Its performance measurements and stated statistical limitations remain unchanged.

## Next work

Measure 1,000-tip searches with covariance estimation enabled, increasing shifts from 10 to 30 to 100 and recording elapsed/CPU time, peak memory, convergence and nonfinite failures. Optimize verified refitting bottlenecks while preserving numerical results. Then assess shared-regime searches and the complete calibration/support-bootstrap pipeline. The existing 3.44× full-search improvement applies to the documented fixed-covariance benchmark only.

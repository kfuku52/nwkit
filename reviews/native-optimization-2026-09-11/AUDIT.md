# Native SHIFT optimization audit — 2026-09-11

## Fixes

- Include the optimization report, measurements and harness in the sdist, and require their presence in distribution validation.
- Move the global Python cache/bytecode exclusion after every graft. A setuptools FileList reproduction confirmed that the previous order admitted report `__pycache__` files. Reject cache/bytecode members during archive validation.
- Reject reused timing records when Python, NumPy, SciPy, platform, thread settings, harness hashes, or baseline sources differ. Require real source packages and compare complete search configuration except the output path.
- Use explicit comparison failures instead of assertions that disappear under `python -O`. Reject matching NaN/infinite screening arrays and different array keys.
- Reject colliding JSON/NPZ output names and invalid fixture sizes in the screen harness; record SciPy and thread settings.

## Numerical and evidence review

The three numerical implementation files are unchanged from the final measured snapshot. Additional independent dense-covariance tests cover 20 randomized topologies, unary nodes, zero/negative transition slopes, root variances, arbitrary observation subsets/order, and more structures than the 16-entry cache can retain. No numerical mismatch was observed.

The original measured harnesses are retained byte-for-byte in `measured-source/`, with hashes matching the historical `source-integrity.json`. They are historical evidence: to replay the original layout, place them at their original paths in a separate source copy. Use the current harness for new measurements. Historical manifests without the newly required provenance fields are deliberately rejected by the strengthened reuse check; their recorded measurements have not been rewritten.

The previously reported timings remain measurements of the same numerical implementation. Only the validation harness and packaging were changed in this audit.

## Completed validation

- Frozen audit source: `master` at `080c069` plus this optimization/audit change set. Concurrent CLI/root/reconcile work was excluded from the audit copy and the commit scope.
- `OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 MPLBACKEND=Agg python tools/check.py full`: **4037 passed, 79 skipped, 6 warnings**, 963.58 seconds. Combined line/branch coverage **86%** exceeds the repository gate. Repository-wide lint/format, uncached mypy (220 files), pip consistency, Bandit, dependency vulnerability audit and maintainability passed. [Full log](audit-full.log)
- Distribution-order regression was reproduced directly with setuptools FileList. The final archive check rejects bytecode/cache entries; this avoids adding setuptools as a new dependency for the plain pytest suite.
- `SOURCE_DATE_EPOCH=1789090000 python tools/check.py dist` in a separate source copy: independent wheel/sdist builds, byte reproducibility, metadata and archive-content checks passed. Existing workspace build artifacts were preserved. [Build log](audit-dist.log)
- All local links from the archived optimization report and verification document resolve inside the sdist; no Python bytecode/cache members remain. [Artifact hashes and checks](distribution-audit.json)
- Final edited Python files passed scoped Ruff lint/format again after packaging/test portability adjustments. The numerical files remain identical to the measured snapshot.
- The resource-based screening benchmark is POSIX-only; its three argument-rejection integration cases are explicitly skipped on Windows. Numerical tests are not subject to that skip. Other Python/OS/BLAS combinations and production-adoption simulation studies were not run here.

The new full run covers the prior 218 native checks and the additional numerical/harness audit cases. Skipped optional tests remain skips, not evidence that their external backends passed. The logs include a Requests dependency warning in the local environment; pip consistency and the fresh dependency vulnerability audit nevertheless passed.

## Next priority

Profile 1,000-tip workloads with 10, then 100 shifts, estimated covariance, known/estimated measurement error, and distinct/shared regimes. Measure calibration and support separately from single-search time. Use the resulting cost breakdown to choose between QR updates, covariance optimizer improvements, and active-set/matrix-free screening; do not assume the small-K speedup extrapolates to these workloads.

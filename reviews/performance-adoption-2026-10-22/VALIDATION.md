Delivery validation for NWKIT 0.43.55

The complete release entrypoint exited **0** on implementation revision
`f3d1bd3f`:

```sh
DYLD_FALLBACK_LIBRARY_PATH=/Users/kf/miniforge3-arm64/lib \
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 MKL_NUM_THREADS=1 \
NWKIT_TEST_WORKERS=3 NWKIT_CI_DIAGNOSTICS=0 \
/tmp/nwkit-perf-20261022-9fYwTD/venv/bin/python tools/check.py release
```

- Ruff lint and formatting: passed (589 files).
- Uncached mypy: passed (260 source files).
- `pip check`: passed; no broken requirements.
- Bandit: passed at the repository's required severity/confidence settings.
- `pip_audit .`: passed; no known vulnerabilities found.
- Complete regression suite, including slow regressions: **5,220 passed,
  127 skipped**, 366.36 seconds with three workers.
- Combined line/branch coverage: **85.93%**, above the unchanged 80% gate.
- Maintainability: all applicable complexity ceilings passed. Existing
  comparison-baseline warnings and the existing draw exception were retained;
  no ceiling or test tolerance was raised.
- Independent direct/sdist wheel builds, independent sdists, metadata/archive
  checks, Twine and byte-for-byte reproducibility: all passed.

The release command does not print individual skip reasons. Source and runtime
inspection identifies unavailable optional R/kfl1ou, IQ-TREE/library worker,
PAML/MCMCTree, GeneRax and independent research-simulator integrations, explicit
opt-in scientific studies, and platform-specific cases among the skipped
coverage. These are not claimed as passes or live independent-backend
validation. No optional scientific executable was installed. CairoSVG's real
native PNG rendering passed using the existing ARM64 Cairo library.

Replicated scientific studies were not enabled: the model, uncertainty method,
seeds, tolerances and search/work bounds were retained. Existing small real
bootstrap, independent-reference, numerical-invariance, serialization, CLI,
concurrency and rollback regressions were part of the release run. Local
validation uses one macOS ARM64 Python/NumPy/SciPy environment; GitHub's normal
push workflow provides its selected additional CI environments.

The full local log is preserved at
`/tmp/nwkit-performance-adoption-20261022/release.log`, SHA-256
`3331408b9fc6b5c44bc8ca6b44cf2d06c9fa0afe1b31bd77e9e7c24846817b36`.
The checked implementation has no subsequent runtime/test changes; the final
record/control files are documentation and experimental evidence only. The
documentation-only delivery revision is followed by `tools/check.py dist`,
as required by DEVELOPMENT.md, before pushing. Original build artifacts remain
in the backup directory recorded in REPORT.md.

The fixed refreshed remote base is
`3da1df5724b432e59b3b1f04d60c649852cb55b3`. The intended publication is a normal
push of the existing `master` branch to `origin/master`, with version 0.43.55.
No manual tag, release or branch-protection change is part of this task.

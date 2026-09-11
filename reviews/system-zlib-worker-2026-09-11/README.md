# External IQ-TREE worker: system zlib linking

GeneGalleon [run 34595419636](https://github.com/kfuku52/genegalleon/actions/runs/34595419636)
failed while linking its AMD64 worker, before SIF conversion or runtime tests.
IQ-TREE built successfully, but `libiqtree.a` retained references to `gzread`,
`gzopen`, and the other system zlib entry points. NWKIT's adapter setup copied
compiler and global linker flags but omitted the external zlib selected by CMake.
The original errors and command are retained in `amd64-link-failure.log`.

NWKIT 0.43.17 adds the matching configuration's `ZLIB_LIBRARY_RELEASE` or
`ZLIB_LIBRARY_DEBUG` after the static archive in the link command. Builds without a CMake
system-zlib record do not add an external zlib link argument. The regression test
constructs a real static archive that calls zlib: linking without the selected
library fails, while the adapter command links and runs successfully. Paths with
spaces, Debug selection, absent libraries, and bundled builds are also checked.

The focused suite passed 20 tests in a GeneGalleon ARM64 container with a C++
toolchain and zlib development files. Ruff, formatting, mypy, pip check,
Bandit, dependency auditing, and the complexity limit passed.

Before this change, NWKIT [run 34574344999](https://github.com/kfuku52/nwkit/actions/runs/34574344999)
already failed its archived shift-engine audit on Linux, an errors-in-variables
regression case on Python 3.13, and Windows collection because `resource` is
unavailable. Those independent failures are outside this build-link correction.


An actual IQ-TREE 3.1.4 library and the corrected worker also built and passed
capability checks on ARM64 (`arm64-library-build.log`). This architecture's
IQ-TREE build additionally merges `liblinux_arm/libz.a` into its archive even
when CMake finds system zlib, so the old adapter can link there too.

For a direct replay of the failing dependency boundary, a validation-only copy
of that IQ-TREE archive had the members from `liblinux_arm/libz.a` removed with
`ar d`. The library code and adapter were otherwise unchanged. The old adapter
then failed on the same gzip symbols, while the corrected adapter linked and
passed capability checks (`system-zlib-replay.log`). This is an ARM64 reproduction
of the missing external-library condition, not an AMD64 or SIF execution claim.


The full local suite completed with **4,172 passed, 54 skipped, one failed**.
The sole failure was `test_archived_engine_requires_explicit_scope_and_intact_snapshot`
with `Seeded complete-search replay disagrees`, the same unrelated test/error
already present in the preceding hosted run (`full-checks.log`). It is retained
and is not skipped or weakened by this change. The full/release gate is therefore
not reported as green.

The independent wheel/sdist builds passed metadata, archive-content, and
byte-reproducibility checks (`distributions.log`).

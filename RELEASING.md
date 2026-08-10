# Releasing NWKIT

1. Update `nwkit.__version__` for every change merged to `master` and move the
   relevant entries in `CHANGELOG.md` from `Unreleased` into a dated version
   section.
2. Run the local release checks:

   ```sh
   python -m venv .venv
   . .venv/bin/activate
   python -m pip install -U pip
   python -m pip install -c constraints-dev.txt -e ".[dev,image]"
   python tools/check.py release
   ```

   The release runner derives `SOURCE_DATE_EPOCH` from the current commit when
   it is not already set, runs all quality and security gates, builds both
   wheel paths, and compares their bytes and contents.

3. Commit the version and changelog changes and merge them into `master`. The
   `Tests` workflow validates the push, after which the release-tag workflow
   checks the version from the exact tested commit.
4. Patch-only versions whose patch component is nonzero (for example,
   `0.34.2`) do not receive a tag or GitHub Release. Major and minor versions
   whose patch component is zero (for example, `0.35.0` or `1.0.0`) receive an
   annotated `v<version>` tag automatically.
5. The tag starts the `Release` workflow, which verifies the source and tests,
   builds byte-for-byte reproducible distributions, and creates the GitHub
   Release.
6. Bioconda discovers tagged upstream releases, so its `nwkit` recipe is
   updated for major and minor releases only. Patch-only versions are
   intentionally not autobumped.

Do not create release tags manually unless recovering the automated workflow.
If recovery is necessary, point the annotated tag at the commit that passed
`Tests` and preserve the existing `v<version>` tag format.

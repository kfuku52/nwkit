# Releasing NWKIT

1. Update `nwkit.__version__` and move the relevant entries in `CHANGELOG.md` from `Unreleased` into a dated version section.
2. Run the local release checks:

   ```sh
   ruff check nwkit tests
   pytest tests/ -q
   python -m build
   python -m twine check dist/*
   check-wheel-contents dist/*.whl
   ```

3. Commit the version and changelog changes, merge them into the default branch, and confirm the Tests workflow is green.
4. Create and push a matching annotated tag, for example:

   ```sh
   git tag -a v0.30.0 -m "NWKIT 0.30.0"
   git push origin v0.30.0
   ```

5. The Release workflow verifies the tag, rebuilds and checks the distributions, and creates the GitHub Release.
6. After the GitHub Release succeeds, update the Bioconda recipe version, source hash, Python requirement, dependencies, and license metadata in a separate change.

# Documentation and implementation audit — 2026-09-22

## Scope and source state

Started on `master` at `a5ffa952ad3b1775cfbe6c8f1e7ed7f60d0a3688`.
The initial worktree contained changes to workflows, agent/development/release
instructions, README, changelog, version metadata and a validation skill.
They were preserved and became commit `46a614d` during this audit.
This audit adds a separate patch version, 0.43.28. Historical release entries
and research results were not rewritten. Runtime behavior was not changed.

Prioritized installation, first execution, shared CLI/TSV contracts, tree/table
conversion, scalar and branch-specific ASR, signal, PCA and DTT outputs.
Read parser definitions, handlers and relevant tests, not just help text.
Checked literal `nwkit` shell commands in `docs/guides/` and example READMEs
against registered command options; checked local Markdown link destinations
across README, development/release guides, docs and example READMEs. The final scan covered 64 documents, 278 local links and 58 literal commands.
All checked paths exist and all checked literal long option names are registered. This
check does not validate every option value, shell expansion or Markdown anchor.

## Findings and disposition

| Class | Location and problem | Evidence and action |
| --- | --- | --- |
| A | `docs/guides/CONVERT.md`: opening conversion chain used only `--age-ci drop` for plain Newick. It fails when `age` or other normalized properties remain. | Reproduced exit 2 with a two-tip NHX tree containing `age` and an interval. `tree_formats._render_attributes` rejects remaining properties; `tests/test_convert.py` covers deliberate property removal and refusal of implicit loss. Changed the example to `--properties drop`, explained the data loss, retained the annotated output, and identified the initial MCMCtree file as user-supplied. |
| A | `docs/README.md`: the combined “DTT, PCA, signal, and stochastic maps” link led only to DTT. | All four separate guides exist and describe separate registered commands/features. Replaced the combined entry with direct links and added conversion and quick-start routes. |
| A | `README.md`: Git installation lacked environment isolation and the Git prerequisite; there was no local first analysis. | `pyproject.toml` requires Python >=3.10 and defines the console entry point; runtime imports and local source installation were tested. Added fresh-environment commands, checkout-install alternative, Windows activation/workaround link and a quick-start link. Git/Bioconda installation itself was not executed. |
| A | Shared CLI guide omitted environment precedence/cache details and described every boolean as taking a value. | `cli.main` parses explicit arguments/defaults without a global config loader; `--debug` is `store_true`. `iqtree_library.find_worker`, `image.resolve_*worker_count`, `image.resolve_*cache_dir` and `util.resolve_download_dir` define the documented exceptions. Added concise descriptions and linked the existing worker precedence guide. Cache/worker rules were inspected in code; no image service was contacted. |
| A | Shared output guide did not explain the default `nwk2table` columns, the meaning of `sister`, or whether `table2nwk` consumes `age`. | `nwk2table._node_row`, `_sister_branch_id`, `table2nwk_main` and their tests establish the schema and reconstruction fields. Added column/age-unit descriptions and warned that one sister is not a complete polytomy description. |
| A | No self-contained installed-user example connecting validation, missing traits and output interpretation. | Added `docs/guides/QUICK_START.md`; executed its shell blocks, checked seven-node TSVs and a table round-trip, and reran ASR to verify deterministic replacement. Interpretation follows the existing ASR guide; no new scientific assumptions were introduced. |
| B | `info -i -` reports an invented input-file path despite reading stdin. | The shared CLI guide and help define `-` as stdin. `info_main` unconditionally displays `os.path.realpath(args.infile)`. Reproduction below exits 0 but prints `Tree file PATH: <cwd>/-`; that file does not exist. Tree statistics are correct. Left the implementation and stdin contract unchanged; a separate fix should distinguish stream/inline input from real paths and test its reporting. `tests/test_info.py` covers statistics but not this stdin path label. |

The audit host was macOS 27.0 (26A428), ARM64.

A second **B** installation/runtime issue was reproduced in a fresh Python
3.10.21 environment: `python -m pip install .` installed 0.43.28 and SciPy
1.15.3, and `pip check` plus top-level imports passed. The quick start stopped
with exit 2 on `nwk2table --age yes`: SciPy's `_spropack.cpython-310-darwin.so`
could not load (`section '__DATA/__thread_bss' has a zero-fill section type,
but offset field is not zero`). `python -c 'import scipy.sparse.linalg'`
reproduces the failure without NWKIT. This isolates the observed failure to
loading the dependency binary in this host/environment; it does not establish
that all Python 3.10 or macOS installations fail. No dependency pins, supported
Python declaration or runtime code were changed. README/quick start now
explain why help/metadata checks alone are insufficient. `DEVELOPMENT.md` now
imports the compiled SciPy submodules in preflight (category A: incomplete
verification instructions); that strengthened command passes on Python 3.14.7
and exposes the Python 3.10.21 failure. Python 3.10 sample
execution is **not verified**; diagnosing/rebuilding the dependency is deferred.

No category C scientific/specification conflict was established within this
scope. That is not a certification of the unreviewed models or backends.

### B reproduction (small substitute input)

```sh
printf '%s\n' '(A:1,B:1);' | nwkit info -i -
```

Actual first line: `Tree file PATH: <working-directory>/-`. The input came from
the pipe, and no file named `-` was created. Related implementation:
[`nwkit/info.py`](../nwkit/info.py), `info_main`, construction of `lines`.

### Conversion reproduction (small substitute input)

```sh
printf '%s\n' '(A:1,B:1)[&&NHX:age=1:age_ci_low=0.5:age_ci_high=1.5:age_ci_kind=HPD:age_ci_level=0.95];' > dated.nhx
nwkit convert -i dated.nhx --to newick --age-ci drop -o dated.nwk
nwkit convert -i dated.nhx --to newick --properties drop -o dated.nwk
```

The first conversion exits 2 (`Plain Newick cannot retain these properties:
age`); the corrected conversion exits 0 and writes `(A:1,B:1);`.
This is a reduced synthetic input, not a run of the full MCMCtree workflow.

## Execution evidence

All trial outputs used fresh temporary directories. The existing `.venv` was preserved. Release validation ran in a separate
source copy so its build cleanup did not clear the checkout's build, dist
or direct-dist directories. The existing `.venv/bin/python`
failed with `bad CPU type in executable`; it was not repaired or replaced.
A separate macOS ARM64 Python 3.14.7 environment was created using
`python -m venv`, then installed with:

```sh
python -m pip install -c constraints-dev.txt '.[dev,image]'
python -m pip check
python -c 'from ete4 import Tree; import numpy, scipy, pandas, matplotlib, PIL; assert len(list(Tree("(A:1,B:1);", parser=1).leaves())) == 2'
```

Installation, metadata checks and runtime imports passed (ETE4 4.4.0).
Fetched dependencies were confined to isolated environments; this does not test
Bioconda resolution or the Git URL installation on other platforms.

The following published shell blocks were executed verbatim. A temporary
working directory contained an `examples` symlink to the checkout, so input
paths and numerical settings were unchanged while outputs were isolated.

| Published command location | Result |
| --- | --- |
| `examples/signal/README.md`, first block | Exit 0; four rows. K = 2.164516 and 0.476939 within 1e-6; lambda = 1 and 0 with boundary status. All 26 documented columns present. |
| `examples/branch_gaussian/README.md`, all three blocks | Exit 0; ASR 5 rows, normalized model 4 rows, likelihood 1 row, prior simulation 15 rows. Log likelihood -5.487990617573634 matches the published value. Both process JSON files exist. |
| `docs/guides/PCA.md`, first block | Exit 0; scores 8 rows, loadings 9, eigenvalues 3, ancestors 45; JSON and nonempty PNG created. Output headers match the guide. |
| `docs/guides/DTT.md`, first block | Exit 0 with the published 999 simulations and two workers; DTT 8 rows and summary 1 row, JSON and nonempty PNG. All 11 primary columns match. |
| `docs/guides/QUICK_START.md`, both blocks | Exit 0 for each command; validation status `ok`, 7 node/ASR rows, root age 2, tip D marked imputed; restored tree validates. Repeating ASR replaces its TSV with identical bytes. |

Also ran `nwkit COMMAND --help` for validate, nwk2table, table2nwk, convert,
asr, signal, pca, dtt, image and radte. Output schemas were cross-checked with
handlers and tests. Figure files were checked for existence/size; this audit
does not claim a new visual/layout review.

Focused verification:

```sh
python tools/check.py quick -- tests/test_cli.py tests/test_cli_contracts.py tests/test_interface_conventions.py tests/test_convert.py tests/test_nwk2table.py tests/test_table2nwk.py tests/test_validate.py tests/test_branch_gaussian_cli.py tests/test_signal.py tests/test_pca.py tests/test_dtt.py tests/test_image_helpers.py -rs
```

Ruff lint/format and mypy passed; **422 tests passed**, with no skips.
No new test duplicates those existing contracts; the new quick-start commands
were checked directly from their Markdown source. No dedicated documentation
renderer or link checker is configured; the local link/option scans above
supplement the distribution/document inclusion check.

Delivery verification ran `python tools/check.py release` in a separate source
copy with the final package documentation and version. It exited **0**:
**3976 passed, 83 skipped** in 392.28 seconds, **85%** branch-aware coverage;
Ruff, mypy, dependency metadata, Bandit, pip-audit and maintainability gates
passed. The existing complexity-baseline warnings concern unchanged runtime
functions; the baseline was not relaxed. Wheel/sdist metadata, contents and
byte reproducibility passed. `git diff --check` passed as well.

The host has no `iqtree3`, `nwkit-iqtree-worker`, `mcmctree` or `Rscript`, and
`NWKIT_TEST_RSCRIPT` is unset. These optional integration paths remain
unverified; skipped tests are not passes. The release runner did not request
pytest's per-case skip-reason summary, so the aggregate skip count does not
provide an individual reason for every skipped case.

## Limits

No live image/taxonomy service, R/kfl1ou, IQ-TREE library build, PAML execution,
large calibration study, paid service or external write was used for sample
validation. Optional-backend unit tests with fixtures are not live integrations.
Git/Bioconda installs, Windows/Linux installations, external web/wiki links and
remote Markdown anchors were not verified. The Windows ETE workaround is
linked to the existing development guide rather than claimed as tested here.

Regression families, SHIFT/RADTE statistical interpretation, all ASR model
variants, every docstring and historical validation-study commands were not
exhaustively audited. Local links in historical documents were inspected but
historical results and changelog descriptions were not modernized. This audit
establishes the listed cases, not repository-wide absence of discrepancies.

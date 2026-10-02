# Native MUL-Reconciliation Validation

Executed 2026-10-02. This validates the exact fixed-tree D+L objective and
GeneGalleon replacement, not biological error rates for inferring polyploidy.
See the [method and output contracts](../guides/MUL_RECONCILE.md).

## Independent and Original-Program Evidence

| Check | Executed result |
| --- | --- |
| Independent enumeration | Five fixed gene topologies over all 20 three-species hypotheses: minimum D+L, every optimal duplication/loss pair, and tie counts matched. |
| Seeded topology/copy-count enumeration | Seed 20261002, 24 rooted gene trees with 2-7 tips over all 40 four-species hypotheses: all 960 gene/candidate combinations matched complete leaf-assignment enumeration. |
| Positive controls | Known auto- and allopolyploid placements scored zero and improved upon the singly labeled baseline. |
| Negative/missing-species controls | Matching single-copy topology preferred the baseline; two X copies on `((A,X),B)` scored one duplication plus two ancestral losses. Root losses reconstruct from JSON node rows. |
| Original author's manual data | All candidate scores matched for 25 trees and 10 hypotheses with H1=`x,y,z`. Only absent terminal semicolons were added to the author's tree lines. No genes filtered. |
| Small original-program comparison | Five null/auto/allo/missing/duplicate gene cases: 19/20 candidate scores matched; the whole-root auto candidate was 15 in GRAMPA versus exact 13. |
| Difference adjudication | Every leaf assignment was scored with the original `reconLCA` function, independently of native DP. Per-gene minima `[2,1,4,2,4]` sum to 13, matching the separate ancestor-path enumeration and NWKIT. The grouping stage excludes a better assignment; no legacy heuristic was copied into NWKIT. |
| Runtime/CLI contracts | Sequential and two-process scoring agreed; report/score transactional preservation, output caps, stdout ownership, empty/unrooted/polytomous/unknown-species inputs and input-overwrite rejection passed. |

Independent enumeration uses explicit ancestor paths and all leaf-occurrence
products, not NWKIT's LCA index or DP recurrence. DP ties are compared as a
multiset of duplication/loss totals, including all optimal leaf assignments.
The whole-root difference prevents claiming universal numerical equivalence
to the original grouped algorithm; it does not validate a polyploidy P-value.

Reference source: [author's GRAMPA repository](https://github.com/gwct/grampa),
commit `4c01a982b11b9adc7738b81cbff53ad6a361b7ea`. This is test-run provenance,
not a default upstream pin. The objective is described by
[Thomas et al. (2017)](https://pubmed.ncbi.nlm.nih.gov/28419377/) and the
[official documentation](https://gwct.bio/grampa/readme.html).

Reference data SHA-256:

```text
manual_species_tree.tre  9d2d6343b1d59bbef8750fbf05d35c653824bcc30460ec5bc73456ba92dacb27
manual_gene_trees.txt    deb4553706828128fe8f488a4a7d6b1e9874e0823200f3e43ca4e73676a4cc94
```

Original-program tests are opt-in via `NWKIT_GRAMPA_REFERENCE`; they were
explicitly enabled in both host and Docker runs here, not skipped.

## Execution and Workflow Evidence

NWKIT's existing checker was invoked from its repository root:

```bash
NWKIT_GRAMPA_REFERENCE=/tmp/nwkit-grampa-reference-20261002/grampa.py \
  /tmp/nwkit-wgd-dev-20261001/bin/python tools/check.py quick -- \
  tests/test_mul_reconcile.py tests/test_mul_reconcile_reference.py \
  tests/test_cli_contracts.py tests/test_cli.py -x
```

Host CPython 3.12.14: 125 tests passed; whole-repository Ruff, formatting and
mypy passed. The same four test files passed all 125 tests under the fresh
ARM64 GeneGalleon Docker runtime with owned NWKIT source mounted and
`PYTHONPATH=/Users/kf/repos/nwkit`. This is a source-overlaid development
runtime, not a published or newly built deployment image.
An attempted cross-repository co-collection failed before executing tests
because the repositories' pytest marker configurations conflict. The suites
were rerun separately with their own configurations.

GeneGalleon checks used its existing `workflow/tests/run_in_runtime.sh` and
`dev` entrypoints, from the GeneGalleon repository root. Environment:

```bash
export GG_TEST_RUNTIME=docker
export GG_CONTAINER_DOCKER_IMAGE=local/genegalleon:subgenome-20261001
export GENEGALLEON_DOCKER_EXTRA_BINDS=/Users/kf/repos/nwkit
```

Executed commands and results:

- `bash workflow/tests/run_in_runtime.sh env PYTHONPATH=/Users/kf/repos/nwkit python -m pytest -q workflow/tests/test_native_mul_reconcile_runtime.py workflow/tests/test_parse_grampa.py -x`: 11 passed on final source. The real shell wrapper/native CLI/summary parser ran with `grampa.py` replaced by an executable that fails if invoked. Legacy output paths, species/gene summaries, nonempty best tree, and metadata passed. Legacy/taxonomic/TSV mapping, branch-length preservation and invalid-preparation output preservation passed.
- The same wrapper with `test_native_mul_reconcile_runtime.py`, `test_parse_grampa.py` and `test_genome_evolution_protein_mode.py`: 63 passed, including broader protein/config/cache workflows. This broader run preceded the final serialization refinements; the affected focused 11 checks were rerun afterwards. Counts overlap and are not additive.
- `bash ./dev check static`: 293 passed in Docker. These are static assertions, not scientific runtime evidence.
- `bash ./dev config-check`: passed on host (8 entrypoints, 23 common parameters).
- Whole GeneGalleon `ruff check workflow container/scripts`: passed on host. `git diff --check` passed in both repositories.
- `bash ./dev lint`: host aggregate blocked by Bash 3.2; Docker aggregate completed all tracked shell syntax checks under Bash 5 but then stopped because the runtime lacks Ruff. Docker Git ownership was scoped to this repository using per-command `GIT_CONFIG_*`, not a global configuration change. Separate host Ruff/config checks passed. The aggregate command itself is not reported as passed.
- NWKIT maintainability check: all hard limits passed; existing unrelated baseline complexity warnings remain.

Docker image ID:
`sha256:8c03f1aec42d761e44ad2077212fb97e5e2a495d5b7729fe324fe56665e01568`.
Runtime freshness checks passed without bypasses. No SIF/Apptainer, native
amd64, full production analysis, or fresh multi-platform build was executed.
No commit or push was made. Legacy container GRAMPA installation checks were
not removed; the three replaced stages no longer call the executable.

Final native source SHA-256:

```text
mul_reconcile.py        68bdf282967f636b77fa73042b4ec0f40ea7c4c209cf06883894978b7aeba54d
mul_reconcile_model.py  5722204e5d63cc707538d731a207c99779d55569a2898217098d8a1a066014ac
mul_reconcile_cli.py    3339d8dff8daae73dd4a6c5fed6fd9003a1e5b5298ad29c1e5526b383a66254f
```

## Follow-Up Audit on 2026-10-02

The initial execution record above is retained. A second audit reproduced and
fixed these issues without changing the objective or scientific thresholds:

- Recursive traceback failed on a valid 1,201-tip comb gene tree. Iterative
  backpointer unranking now exports all 2,401 node rows: 1,200 duplications,
  one ancestral loss, score 1,201 and exactly one optimal assignment.
- Detailed mapping rows were accumulated before writing. They now stream into
  the existing multi-output transaction; mapping-cap failures still preserve
  every prior output. Unrequested check rows are no longer collected.
- MUL CLI options lacked underscore aliases required by NWKIT's existing
  interface conventions. Canonical kebab-case and compatibility aliases pass
  the repository-wide interface tests.
- GeneGalleon's whitespace normalization merged multiple H1 selectors.
  Whitespace delimiters are now preserved while species underscores normalize.
- A shared `grampa_out` scratch directory and sequential output publication
  could collide with other work or partially replace results. Unique scratch
  directories and the existing recoverable nine-file bundle publisher preserve
  results on native, invalid-input and publication failures. Stage callers no
  longer delete the previous summary before running.
- Three stage contracts tracked only the summary. They now track all related
  results, including native metadata and the best tree. Empty input with any
  prior result fails without recording a stale analysis as current.
- Preparation accepted duplicate output targets before collapsing them into a
  dictionary. Distinct-path validation now precedes that dictionary. Directory
  input also avoids one command-line argument per gene tree, retaining the
  existing exclusion of hidden `.nwk` files.
- Summary reading treated literal `NA` and embedded `#` as missing data or
  comments, and filename inventory as CSV. Literal labels and quoted filenames
  are preserved; gene collections use the shared Newick stream parser instead
  of splitting by physical line. Multiline trees remain one gene tree.
- Legacy GRAMPA container requirements remained. Both architecture manifests
  now remove them, and the native MUL score/map/model contract is mandatory.

Additional scientific checks executed:

- Seed 20261003: 30 supplied-MUL tests with two A occurrences and three X
  occurrences matched independent enumeration of **every node assignment**,
  not just score pairs; no duplicate traceback assignments were emitted.
- Seed 20261004: 12 random 2-6-tip genes over 20 hypotheses (240 pairs)
  matched exhaustive assignments scored by the original GRAMPA `reconLCA`,
  including minimum scores and exact optimal-assignment counts.
- Child order and branch-length perturbations left candidate scores unchanged.
  The earlier 960-pair independent check and all four original-source tests
  also passed, with original tests explicitly enabled and no skips.

Final NWKIT host command used `tools/check.py quick --` with
`test_mul_reconcile.py`, `test_mul_reconcile_reference.py`,
`test_cli_contracts.py`, `test_interface_conventions.py`, `test_cli.py`,
`test_output_transaction.py`, and `test_rooting_state.py`: **247 passed**;
whole-repository Ruff, formatting and mypy passed. The same 247 tests passed
in the initial fresh ARM64 container with mounted source, before packaging
changes. Native installed file hashes below match that tested source.

An incremental development image was built from
`local/genegalleon:subgenome-20261001` using the temporary recipe
`/tmp/gg-mul-runtime-20261002.qcopE8/Dockerfile`. It consumes the updated
Conda/tool declarations and validation scripts, removes only GRAMPA through a
normal Conda transaction, and installs the current local NWKIT wheel. Other
source layers and their recorded revisions are inherited from the base.
Build-input hashes were computed with the existing repository helper and the
base's source-revision manifest; no defaults or upstream pins were changed.
The image records `nwkit_local_source.sha256`. `pip check` and all **62 required
runtime checks passed**, including real native MUL execution. `grampa.py` is
absent from PATH. This is not a fresh full production or multi-platform build.

Image: `local/genegalleon:mul-reconcile-20261002`,
ID `sha256:258fa6200c4049c6b689db12ab27d6500175a9d05aacbe7c6c5c223f18b3e5c0`.
Its standard daily freshness check matches current container inputs without
bypass. GeneGalleon validation used this installed NWKIT with no `PYTHONPATH`
overlay, through its existing `run_in_runtime.sh` and `dev` commands:

- `python -m pytest -q -rs -p no:cacheprovider` with native wrapper, legacy
  parser and protein-mode files: **70 passed** before the final summary-reader
  refinements. The affected wrapper/parser checks were rerun: **19 passed**.
- Earlier expanded wrapper/parser/protein-mode/workspace-helper coverage:
  **194 passed** in the base/source-overlay runtime, including real bundle
  rollback and helper lock tests. Counts overlap and are not additive.
- `bash ./dev check static`: **294 passed** on the final container recipe.
- `bash ./dev config-check`: **8 entrypoints and 23 common parameters** passed
  in Docker and on host. The native container reconciliation probe passed again.
- Parser and optional-output provenance checks: **6 passed, 52 deselected**;
  literal labels, multiline trees and optional-output present/absent states passed.
- All tracked shell syntax passed under Docker Bash 5. The same aggregate
  syntax attempt under host Bash 3.2 failed in the unrelated transcriptome core;
  it is not runtime evidence. The earlier aggregate lint limitation remains:
  the image lacks Ruff, so separate whole-repository host Ruff checks were used.
- `git diff --check` passed in both repositories. SIF, native amd64, full
  production workflows and remote-source deployment builds remain unverified.

Audited native source and installed-wheel SHA-256:

```text
mul_reconcile.py        f6501e2c2e1bb28104bfbc599a5fcfdb1e839bdf82e6b8f27181981f9d7f4a1e
mul_reconcile_model.py  81de463c92fbe9381b97e404712fce192dfd92dde4e85af0086cebcc2a5cb85d
mul_reconcile_cli.py    834daf08df1becef5d12b8545a92db5da5b009953036b9cf4573b1c0f7622665
```

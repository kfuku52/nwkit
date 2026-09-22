# First tree and trait analysis

Install NWKIT using the [installation instructions](../../README.md#installation),
then check `nwkit --version` and `nwkit asr --help` in that environment.
`python -m nwkit` is equivalent to the `nwkit` entry point. This example needs
only the Python runtime dependencies; it does not call R, IQ-TREE, PAML or image
services. Native Cairo plus the image extra is needed for SVG rasterization,
not for these TSV/Newick outputs.

If installation succeeds but a command reports a native-library import error,
check `python -m pip check` and
`python -c 'import scipy.linalg, scipy.sparse.linalg'` in the same environment.
Metadata checks and `--help` alone do not establish that compiled dependencies
can load. Keep the traceback and use the
[isolated environment troubleshooting guidance](../../DEVELOPMENT.md) before
changing an existing environment. A failed example is not a validated install.
For the specific SciPy 1.15.3 `_spropack` / `__thread_bss` error on macOS 27
ARM64 with Python 3.10, use the
[verified same-version wheel recovery](../../DEVELOPMENT.md#scipy-1153-wheel-on-macos-27-arm64--python-310).

## Create small inputs

Run these POSIX shell commands in a new temporary directory. `printf` writes
actual tabs in the trait table, including an empty value for tip D.

```sh
NWKIT_EXAMPLE_DIR=$(mktemp -d)
cd "$NWKIT_EXAMPLE_DIR"
printf '%s\n' '((A:1,B:1):1,(C:1,D:1):1);' > tree.nwk
printf 'leaf_name\tx\nA\t1.0\nB\t1.4\nC\t2.0\nD\t\n' > traits.tsv
nwkit validate -i tree.nwk --require-rooted yes --require-all-lengths yes \
  --fail-on-issue yes -o validation.tsv
nwkit nwk2table -i tree.nwk --age yes -o nodes.tsv
nwkit asr -i tree.nwk --trait traits.tsv --state-column x \
  --trait-type continuous --model BM --output summary -o asr.tsv
```

The binary root is treated as rooted in the absence of a declaration. This
tree has seven nodes and root-to-tip length two. The trait key `leaf_name`
must match the tip labels exactly. Blank trait cells are missing observations;
they are not zero. Identifier strings such as `NA` remain literal names.
See [the shared TSV contract](CLI_TSV_CONVENTIONS.md#tip-keyed-tsv-files) for
duplicate keys, other missing tokens and unmatched-row policies.

## Read the results

- `validation.tsv` has one row for the tree, including `status`,
  `is_rooted`, `rooting_state` and `issues`. Without `--fail-on-issue yes`,
  reported invalid trees do not by themselves produce a failing exit status.
- `nodes.tsv` has seven rows with `branch_id`, `parent`, `name`, `dist`,
  `support`, `sister` and `age`. IDs are level-order, with root 0 and root
  parent -1. Ages range from zero at tips to two at the root, in the same
  units as the branch lengths. Absent support is empty. See
  [the output vocabulary](CLI_TSV_CONVENTIONS.md#output-tsv-vocabulary).
- `asr.tsv` has seven rows for trait `x`, including `mean`, `variance`, `sd`,
  `ci_lower`, `ci_upper`, `ci_level`, `observed_value`, `observed_se` and
  `is_imputed`. Tip D has an inferred value rather than an observed zero.
  Means, SDs and interval endpoints are in trait units; variance is in squared
  trait units. These are marginal latent-trait intervals conditional on the
  model parameters and tree, not intervals for new noisy measurements. See
  [ASR](ASR.md) for the default parameter estimation and root treatment.

The commands explicitly name their output files in the current directory;
they do not create a separate result directory or choose a prefix automatically.
Repeating them replaces those files, without an overwrite prompt or automatic
resume. Use another temporary directory to retain a previous run. The
[output transaction rules](CLI_TSV_CONVENTIONS.md#related-outputs-and-node-editing)
are command-specific: do not assume every command or every ASR sidecar is one
atomic transaction. `--audit run.jsonl`, when supplied, appends provenance
records rather than replacing the log.

For a table round-trip:

```sh
nwkit table2nwk -i nodes.tsv -o restored.nwk
nwkit validate -i restored.nwk --require-rooted yes --require-all-lengths yes \
  --fail-on-issue yes -o restored-validation.tsv
```

This restores the table's topology, names, lengths and support. It is not a
general NHX annotation round-trip; use [convert](CONVERT.md) to retain those
properties. For examples with committed input files, continue with
[signal](../../examples/signal/README.md), [PCA](../../examples/pca/README.md)
or [branch-specific Gaussian ASR](../../examples/branch_gaussian/README.md),
running their commands from the repository root.

# CLI and TSV conventions

This document defines the shared interface conventions for NWKIT commands.
Command-specific help remains authoritative for options and columns that are
unique to one command.

## CLI option names

- Canonical long options use kebab-case: `--species-map-tsv`,
  `--quoted-node-names`, and `--retain-per-clade`.
- Existing snake_case spellings remain accepted as unadvertised compatibility
  aliases. Using one writes a deprecation warning to standard error naming its
  canonical replacement. New documentation and scripts should use the
  canonical kebab-case names.
- Boolean options use explicit `yes|no` values. `mcmctree --add-header`
  continues to accept the historical value-less form as well as `yes|no`.
- `-i/--infile`, `-o/--outfile`, and `-of/--outformat` are reserved for the
  primary input, primary output, and Newick output format. Commands producing
  tables, text, images, or MCMCtree syntax do not expose `--outformat`.
- Secondary tabular audits use `--report` where one file is produced. Other
  metadata outputs use a noun followed by `-out`, such as `--tree-out`,
  `--model-out`, `--manifest-out`, and `--attribution-out`.
- `-` means standard input for a supported input or standard output for the
  primary output. Only one input may own standard input in a command. Auxiliary
  metadata and report outputs require file paths so that result streams cannot
  be mixed. `intersection --seqout` is a second primary result and may use
  standard output only when `--outfile` is a file path.

## Tree input containers

Shared tree readers accept standard Newick/NHX plus the PAML/MCMCtree
containers commonly encountered in dating workflows:

- a PAML treefile beginning with `nTips nTrees`;
- a direct-label NEXUS `TREE` or `UTREE` statement, including MCMCtree
  `FigTree.tre` (NEXUS `TRANSLATE` tables are not applied); and
- the annotated species-tree block in MCMCtree's main text output.

MCMCtree/FigTree `95%HPD` and legacy `95%` node-age comments are normalized to
`age_ci_low`, `age_ci_high`, `age_ci_kind`, and `age_ci_level` NHX properties.
This normalization happens before command-specific processing, including for
STDIN, so a header-bearing `nwkit mcmctree --add-header` result remains valid
input to other NWKIT commands.

### Rootedness and root polytomies

Tree inputs accept `--input-rooted auto|yes|no` (default `auto`). This controls
how the current top-level node is interpreted; it does not move the root or
resolve any polytomy.

- `auto` respects an explicit leading `[&R]` (rooted) or `[&U]` (unrooted),
  including tokens in supported NEXUS tree statements. It also accepts the
  root-only NHX property `nwkit_rooted=yes|no|unknown`.
- Without a declaration, roots with at most two children retain the legacy
  rooted interpretation. Roots with three or more children are `unknown`, not
  necessarily unrooted. A root name, numeric support (including `0`, `1`, or
  `100`), or stem length is not evidence of rootedness.
- `yes` treats the current top-level node as the root, including a multifurcating
  root. `no` treats the tree as unrooted, even when that node has two children.
  Overriding an explicit opposite declaration emits a warning. Malformed or
  mutually conflicting declarations still fail, as do other input-validation
  errors.

Auxiliary trees have independent options: `--infile2-rooted`,
`--reference-rooted`, `--species-tree-rooted`, `--reconciliation-tree-rooted`,
`--root-source-rooted`, `--name-source-rooted`, `--support-source-rooted`,
`--length-source-rooted`, `--property-source-rooted`, and
`--densitree-trees-rooted`, wherever the corresponding input exists. In
`regress`, `--input-rooted` applies to the primary `--tree`, `--gene-tree`, or
every member of `--gene-tree-ensemble`; it does not override the species tree.
Collection inputs apply their selected mode to every tree, while `auto` reads
each tree's declaration independently.

For example, either input below enables ASR at a multifurcating root:

```sh
nwkit asr -i '[&R](A:1,B:1,C:1);' --trait traits.tsv --state-column state
nwkit asr -i '(A:1,B:1,C:1);' --input-rooted yes --trait traits.tsv --state-column state
```

ASR supports both internal and root polytomies. Rooted RF distance compares the
displayed nontrivial descendant clades without inventing binary resolutions;
normalized RF divides by the sum of the two displayed-clade counts (zero when
both counts are zero). Contrasts and reconciliation retain their strictly
bifurcating-tree requirement. Conventional phylogenetic regression, root transfer,
and MCMCtree output still require exactly two root children. A rooting override
does not waive any of these structural requirements.
`consensus --comparison rooted|unrooted` continues to select the comparison for
unmarked inputs; an explicit unrooted declaration is incompatible with rooted
comparison unless the corresponding input override is used.

Shared Newick writers preserve explicit or forced states, and any state that
would otherwise be lost through an operation such as collapsing or pruning,
using `[&&NHX:nwkit_rooted=yes]` (or `no`/`unknown`) on the root. This is readable
by ETE's Newick parser and independent of the selected annotation properties.
Legacy unmarked binary output is unchanged. Rooting operations replace the
input interpretation with rooted status; label and support changes do not.
The `nwkit_rooted` property and private `_nwkit_rooting_*` bookkeeping properties
are reserved and cannot be reassigned through annotation/property transfer.
The Python `read_tree` and `read_trees` helpers expose the same policy through
the keyword-only `rooted="auto"|"yes"|"no"` argument.

`validate` retains the boolean `is_rooted` column and adds `rooting_state`
(`rooted`, `unrooted`, or `unknown`) and `rooting_source` (for input trees:
`marker`, `nhx`, `override`, `topology`, or `unknown`).
`--require-same-rooting yes` compares all immediate root-child clades, independent
of their order, rather than assuming two root children, and distinguishes unknown
from explicitly unrooted trees. `diff --comparison rooted` likewise compares the
full root partition when both trees are rooted; its existing `root_split` summary
row uses `reason=projected_root_partition` for multifurcating roots and lists all
sides separated by `|`. Binary roots retain `reason=projected_root_bipartition`.
An explicitly unrooted input instead produces an unresolved root summary with
`reason=input_explicitly_unrooted`. Audit JSON's
`primary_input` additionally records `first_tree_rooting_state`,
`first_tree_rooting_source`, and `first_tree_rooting_declared`; the latter retains
the original `yes`/`no`/`unknown` declaration, or an empty string if absent, even
when a CLI override changes its interpretation.

## Tip-keyed TSV files

`annotate --table` and the `--trait` files used by `asr`, `asrcompare`,
`contrast`, `draw`, `monophyly`, `sample`, and `skim` follow one shared
contract:

- The first-class key is `leaf_name`.
- `leaf_name` values are exact tree tip labels and must be non-empty. They are
  normally unique. `contrast --biological-id` instead accepts multiple rows
  per tip, and `draw --tip-image-manifest` accepts multiple ranked image rows
  and uses the first. Values such as `001`, `NA`, `NaN`, and `null` are
  preserved as literal labels rather than inferred as numbers or missing data.
- Required command-specific data columns must be present. Extra columns are
  retained where the command can use or report them.
- `--missing-values CSV` controls missing markers in non-key columns. The
  default is the empty string, `NA`, `NaN`, `nan`, `?`, `missing`, and
  `unknown`.
- `--unmatched warn|error|ignore` controls table-only rows and tree-only tips.
  The default is `error` for `contrast` and `warn` for the other tip-metadata
  commands listed above. Table-only rows are not applied to the tree; commands
  that can operate on absent metadata treat tree-only tips as missing. Choosing
  `warn` or `ignore` does not waive a command's requirement for estimable traits.

## Specialized input TSV schemas

| Option or command | Required columns | Additional rules |
|---|---|---|
| `--species-map-tsv` | `leaf_name`, and at least one of `species_label` or `taxonomy_query` | Every row must define at least one mapping value. |
| `image --species-name-tsv` | `leaf_name`, `species_name` | Legacy image-only mapping; prefer `--species-map-tsv` for new workflows. |
| `draw --tip-image-manifest` | `leaf_name`, `local_path` | Multiple rows per tip are allowed and the first is used. Relative paths use the manifest directory or `--tip-image-root`; broken paths are rejected. |
| `rename --name-tsv` | `old_name`, `new_name` | Both values are non-empty and `old_name` is unique. |
| `--taxid-tsv` | `leaf_name`, `taxid` | `taxid` must be a non-missing integer. |
| `--weight-tsv` | `weight` | Positive finite weights. Optional `tree_id` values are unique, 1-based input-tree numbers; without it, row order is used and row count must equal tree count. |
| `contrast --reconciliation` | Output columns from `nwkit reconcile` | Stable `gene_clade_id` values must uniquely and exactly cover the contrast tree. Event/status enums, orientation clades, coverage states, and `tree_id` consistency are validated. |
| `contrast --trait` with `--biological-id` | `leaf_name`, biological-ID column, and selected numeric traits | Multiple rows per tip are biological observations. Optional technical IDs are nested within leaf and biological ID; optional batch values are categorical and complete. Missing trait cells are allowed, but every tip needs an estimable mean. Under pooled variance, a trait with exactly one observation at every fitted tip remains exact tip data even when another selected trait is replicated. |
| `contrast --trait` with `--within-variance known-se` | `leaf_name`, selected means, and one SE column per trait | Exactly one row per tip. SEs are finite and non-negative. Optional sample-size columns contain positive integers. |
| `regress --data` | `leaf_name`, selected `--responses` and `--predictors`, plus optional response- and predictor-replicate columns | Conventional regression requires every tree tip. Either role can use biological/technical replicates, batch adjustment, or role-specific known-SE columns. A non-replicated role must have one value per tip, or the same value on all rows created by the other role's replicate design. |
| `regress --evolution-covariance` | First column `leaf_name`; remaining columns named for tree tips | Wide square covariance matrix. Row and column names must each exactly match all tree tips; values are reordered by name and must be numeric, finite, symmetric, and positive definite. |
| `regress --expression` | `leaf_name`, selected `--responses`, and optional replicate columns | End-to-end mode applies the same replicate, technical-replicate, batch, missing-value, and known-SE rules as `contrast --trait`. Gene-tip labels must match `--gene-tree`. |
| `regress --species-traits` | `leaf_name`, selected `--predictors`, and optional predictor-replicate columns | End-to-end mode accepts continuous and categorical predictors. Every `--species-tree` tip needs an estimable value or category distribution. Continuous predictor replicates use the `--predictor-*` biological/technical ID, batch, variance, known-SE, and sample-size options; categorical replicates follow the policies in the regression guide. |
| `regress` input | Reconciled output columns from `nwkit contrast` | Selected responses must be speciation contrasts. Eligible rows require stable gene, lineage, species-event, and orientation IDs plus positive finite contrast variance. |
| `regress --predictor-contrasts` | Species-tree output columns from `nwkit contrast` | Selected predictor traits must uniquely identify every joined `branch_clade_id`; event taxa and signed orientation must exactly match the response table. |
| `regress --response-sampling-covariance` | `tree_id`, `trait`, `contrast_id_1`, `contrast_id_2`, `sampling_covariance` | Complete explicit covariance or sparse factor loadings, as described below. Required when response contrasts contain replicate metadata. |
| `regress --predictor-sampling-covariance` | `tree_id`, `trait`, `contrast_id_1`, `contrast_id_2`, `sampling_covariance` | Complete explicit covariance or sparse factor loadings for every selected predictor. The predictor contrast table must also contain positive `contrast_variance`; together they define the phylogenetic prior and observation-error covariance of the latent predictor. |
| `table2nwk` input | `branch_id`, `parent` | Optional `name`, `dist`, `support`, and `rooted`; exactly one root has `parent = -1`. Only its row may set `rooted` to `yes`, `no`, or `unknown`; all other rows must leave it empty. |

All supported input TSV options may use `-` for standard input, subject to the
single-standard-input rule.

### Sampling-covariance representations

The two sampling-covariance inputs above describe scalar contrasts, grouped by
`tree_id` and `trait`. The optional `covariance_representation` column selects
the meaning of each row. Omitting the column means `covariance`; when supplied,
each row must name one of the representations below.

- `covariance`: both IDs identify contrasts and `sampling_covariance` is their
  covariance. Every selected unordered pair, including each diagonal, must
  occur exactly once. NWKIT writes the upper triangle; the reader reconstructs
  its symmetric counterpart and requires finite, positive-semidefinite matrices.
- `factor-loading`: `contrast_id_1` identifies a contrast, `contrast_id_2`
  identifies a latent error column, and `sampling_covariance` is a loading in
  a factor matrix `U`. The covariance is `M = U U'`, not the table of loadings
  itself. Omitted loadings are zero. Every selected contrast must occur at
  least once, including an explicit zero row for a contrast with no sampling
  error. Contrast/latent-column pairs must be unique and loadings finite.

One `tree_id`/`trait` group cannot mix the two representations. Reconciled
scalar analyses can emit factor loadings above 500 contrasts, avoiding a dense
matrix while retaining off-diagonal covariance. See the
[regression guide](PHYLOGENETIC_REGRESSION.md#conventional-tip-level-regression)
for the propagation and numerical precision rules.

These scalar re-input schemas are distinct from tip-level regression audit
tables and categorical-predictor audits containing `trait_2` cross-column
covariances. Do not interpret those cross-column blocks as separate scalar
sidecars; use the raw-input regression mode to propagate categorical predictor
uncertainty.

## Output TSV vocabulary

- `branch_id` identifies nodes with 0-based level-order numbering. The root is
  `0`; `parent = -1` denotes no parent. Branch IDs are local to one exact tree
  serialization and are not cross-file join keys.
- `branch_clade_id`, `gene_clade_id`, `species_event_id`, and
  `lineage_clade_id` use a `clade-sha256:` digest of sorted descendant tip
  labels. They are invariant to Newick child order and are the supported keys
  for matching corresponding clades across files.
- Qualified IDs retain the same suffix, for example `target_branch_id` and
  `source_branch_id`.
- `node_class` values are `root`, `intnode`, and `leaf`.
- Columns containing comma-separated taxon labels end in `_taxa`, such as
  `descendant_taxa`, `target_taxa`, and `intruder_taxa`. Their values are sorted
  for deterministic output. Counts use `num_*_taxa`.
- Empty result sets still contain their documented header row.
- Missing output cells are empty. Literal names that resemble missing markers
  remain intact in identifier columns and in `nwk2table`/`table2nwk`
  round-trips.
- `nwk2table` adds the optional root-only `rooted` column when needed to preserve
  a rooting interpretation. Without this column, `table2nwk` uses the legacy
  root-degree inference described above.
- `reconcile` qualifies gene- and species-tree IDs as `gene_branch_id` and
  `species_branch_id`. It retains repeated species-tree IDs when distinct
  paralog lineages map to the same speciation event. `tree_id` plus
  `lineage_clade_id` identifies lineage segments separated by observed or
  pruning-collapsed non-speciation events. Stable cross-file species-event joins
  use `species_event_id`, not `species_branch_id`.
- `contrast` uses `numerator_*` and `denominator_*` columns to make signed
  contrasts auditable. With a reconciliation table, the order follows the
  canonical sorted order of the mapped species-tree daughter clades.
- Replicate-aware `contrast` output adds `sampling_variance`,
  `evolutionary_variance`, `replicate_model`, and `min_n_biological`. Its
  sampling-covariance sidecar stores the full propagated `L D L'` upper
  triangle rather than only marginal variances. The optional tip summary
  records mean estimation and replicate counts. Every contrast row also records
  `evolution_model`, `evolution_parameter_name`, `evolution_parameter`, and
  `branch_length_mode` for the transformed-tree calculation.
- Conventional `regress --tree/--data` reports one row per response coefficient,
  including `(intercept)` by default. It records the number of species,
  evolutionary model and generic shape parameter, whether the parameter was
  estimated or fixed, branch-length mode, evolutionary rate, sampling-variance
  contribution, likelihood, and optimization diagnostics. Its optional
  sampling-covariance sidecar uses `tree_id`, `trait`, `leaf_name_1`,
  `leaf_name_2`, and `sampling_covariance`.
  Predictor replicates have separate `--predictor-sampling-covariance-out` and
  `--predictor-tip-summary-out` files. Result rows distinguish the input
  predictor sampling variance from the posterior latent-predictor variance and
  report the fitted predictor evolutionary rate and optimization diagnostics.
- Conventional `regress --model-comparison-out` reports one row per response and
  candidate evolutionary model. It contains the ML log likelihood, AIC, AICc,
  BIC, likelihood parameter count, estimated shape parameter, convergence and
  boundary diagnostics, delta AIC and AIC weight, and delta AICc and AICc
  weight. An individual AICc is empty when its sample size is insufficient;
  delta AICc and AICc weights are empty for the response if any candidate lacks
  AICc. AIC-derived cells remain available.
- Reconciled `regress` reports one coefficient row per `tree_id`, response, and predictor.
  `n_gene_contrasts` is the physical input-row count, whereas
  `n_species_events` and `degrees_of_freedom` expose the effective independent
  species-event count. Variance-component diagnostics state whether event and
  lineage effects were fitted. The optional random-effects table contains
  conditional modes, their grouping IDs, fitted variance component, and group
  observation count. Precomputed response and predictor inputs must retain the
  `tree_id` and evolutionary-transform metadata emitted by `contrast`; mixed
  transforms or multiple predictor-tree IDs are rejected. Each coefficient row
  records `predictor_tree_id` and both transforms. Raw mode additionally records
  whether each shape parameter was fixed, estimated, or not applicable;
  optimizer convergence and boundary diagnostics; predictor marginal
  log-likelihood; and whether bootstrap inference re-estimated the gene
  parameter.
- End-to-end `regress --out-prefix PREFIX` additionally writes deterministic
  `.reconciliation.tsv`, `.gene-contrasts.tsv`, `.species-contrasts.tsv`,
  `.random-effects.tsv`, and `.regression.tsv` files. Replicate-aware input also
  writes `.response-sampling-covariance.tsv` and `.response-tip-summary.tsv`
  for responses and `.predictor-sampling-covariance.tsv` and
  `.predictor-tip-summary.tsv` for predictors, whenever the corresponding role
  has replicates or known SEs.
  Applicable files are staged and committed as one recoverable transaction,
  with concurrent writers serialized by prefix/output locks; an ordinary write
  or commit failure does not leave a partial mixed bundle.
  Without `--out-prefix`, intermediate tables remain in memory and only the
  requested primary and random-effect outputs are written.
  Gene and species contrast rows, and every final reconciled result, retain the
  independently selected response and predictor evolutionary models, parameter
  names and fixed or automatically estimated values, and branch-length modes.

## ASR defaults and Newick annotations

`asr --trait-type auto` and `asrcompare --trait-type auto` are the defaults.
Only non-missing values in the selected column for tips present in the tree
determine the type: numeric values select continuous BM; otherwise the trait is
discrete Mk. Non-finite numeric inputs are rejected on the continuous path, not
reclassified as categories. An all-missing column cannot be auto-classified.
Use `--trait-type discrete` for numeric category codes, which retain their
original spelling, or explicitly select `continuous` to reject nonnumeric data.
Detection is reported on STDERR and in optional model metadata. Incompatible
mode-specific options are errors; they do not override automatic type detection.

`asrcompare` routes model-specific options only to applicable candidates that
consume them. Its TSV retains one row per requested candidate, including
not-applicable, not-fitted, failed, nonregular, and statistically equivalent
rows. Compatibility groups include trait dimensionality, likelihood convention,
and root prior. Count and rank fields are integral when present; unavailable
values are empty. Shared input-preparation and per-candidate fit durations are
reported separately; shared preparation includes tree/trait parsing and cached
auxiliary inputs. Structurally equivalent rows identify their fitted
representative in `equivalent_to` and do not duplicate criterion weights. The
equivalence contract distinguishes binary GTR from directly bounded ARD/F81,
and collapses exact neutral/fixed or one-regime reductions only when every
relevant option and root contract matches. MV-BM and MV-OU sample sizes count
distinct observed phylogenetic positions rather than scalar coordinates.

Discrete ASR defaults to ER, equal root-state probabilities, and probability
output. Likelihoods, marginals, and conditional node-state sampling use log-space
messages to avoid child-product underflow at large polytomies; zero-length
branches retain their exact no-transition constraint.
Stochastic maps distinguish exactly zero rates from small positive rates: scaling branch lengths
and inversely scaling a fixed rate matrix preserves the process and seeded map
counts. Bounds for estimated rates remain expressed in inverse branch-length
units; adjust `--rate-bounds` when changing those units.

Continuous ASR defaults to BM, a flat prior on the root value, REML rate
estimation, and mean/variance/interval summaries. `--sigma2` fixes the rate;
`--standard-error-column` optionally provides finite, non-negative known SEs for
observed tips. Absent SE input means exact observations. Intervals include root
uncertainty but condition on the fixed/fitted rate and input tree. Zero-length
edges enforce identical latent values; identical exact observations at one such
position count once. A singular zero-rate boundary has an empty residual
log-likelihood field and an explicit fit status, never a substituted finite
density. See [ASR.md](ASR.md) for schemas, boundary rules, and examples.

Both modes use the existing rootedness contract, including declared/forced root
polytomies. Rootedness is separate from the prior on the root's trait value.
Annotated Newick uses the shared writer, preserving quoted numeric internal
names, root declarations, and missing-support handling. Discrete NHX fields use
`asr_state`/`asr_probability`; continuous fields use `asr_mean`, `asr_variance`,
`asr_sd`, and interval properties.

## Protecting related outputs

`draw`, `asrcompare`, `rootcompare`, `regress`, and `intersection` share
recoverable output installation. All file destinations are validated before
writing, staged beside their targets, and installed after generation succeeds.
Handled write/rename failures restore existing outputs and remove newly
installed files. Recovery backups are retained if restoration itself fails.
Existing file modes and the process creation mask are respected.

`draw` rejects outputs that alias each other or any input, including species maps
and referenced tip images. Checks cover symbolic links, hard links, and equivalent
names on case-insensitive filesystems. Image and layout-report failures preserve
the previous complete pair. As with any separate filesystem renames, this is not
a crash-atomic multi-file transaction. Output locks coordinate NWKIT writers, not
arbitrary external programs. A stdout write can trigger file rollback on failure,
but already emitted stdout bytes cannot be retracted.

## Compatibility names

The following older names remain accepted, but canonical output and new
examples use the replacements below. Option aliases emit a deprecation
warning; legacy ASR target values are normalized without an option-name warning.

| Older interface | Canonical interface |
|---|---|
| snake_case long options such as `--species_map_tsv` | kebab-case, for example `--species-map-tsv` |
| `dist -d RF` / `dist --dist RF` | `dist --metric rf` |
| `image --name-tsv` | `image --species-name-tsv` |
| `image --manifest` / `--attribution` | `--manifest-out` / `--attribution-out` |
| `sample --output-table` | `sample --report` |
| `skim --output-groupfile yes` | `skim --group-table-prefix PATH` |
| ASR `node_type` with value `tip` | `node_class` with value `leaf` |
| ASR `--target tip` / `missing_tip` values | `leaf` / `missing-leaf` |
| clade frequency `leaf_set`, `num_leaves` | `descendant_taxa`, `num_taxa` |
| monophyly `*_leaves`, `num_*_leaves` | `*_taxa`, `num_*_taxa` |
| transfer/compose `*_taxon_count` | `num_*_taxa` |
| report `node_id`, `target_node_id`, `source_node_id` | `branch_id`, `target_branch_id`, `source_branch_id` |

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

## Tip-keyed TSV files

`annotate --table` and the `--trait` files used by `asr`, `contrast`, `draw`,
`monophyly`, `sample`, and `skim` follow one shared contract:

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
  The default is `warn`. Table-only rows are not applied to the tree; commands
  that can operate on absent metadata treat tree-only tips as missing.

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
| `pgls --data` | `leaf_name`, selected `--responses` and `--predictors`, plus optional response- and predictor-replicate columns | Conventional PGLS requires every tree tip. Either role can use biological/technical replicates, batch adjustment, or role-specific known-SE columns. A non-replicated role must have one value per tip, or the same value on all rows created by the other role's replicate design. |
| `pgls --evolution-covariance` | First column `leaf_name`; remaining columns named for tree tips | Wide square covariance matrix. Row and column names must each exactly match all tree tips; values are reordered by name and must be numeric, finite, symmetric, and positive definite. |
| `pgls --expression` | `leaf_name`, selected `--responses`, and optional replicate columns | End-to-end mode applies the same replicate, technical-replicate, batch, missing-value, and known-SE rules as `contrast --trait`. Gene-tip labels must match `--gene-tree`. |
| `pgls --species-traits` | `leaf_name`, selected numeric `--predictors`, and optional predictor-replicate columns | End-to-end mode requires an estimable mean for every `--species-tree` tip and predictor. Predictor replicates use the `--predictor-*` biological/technical ID, batch, variance, known-SE, and sample-size options. |
| `pgls` input | Reconciled output columns from `nwkit contrast` | Selected responses must be speciation contrasts. Eligible rows require stable gene, lineage, species-event, and orientation IDs plus positive finite contrast variance. |
| `pgls --predictor-contrasts` | Species-tree output columns from `nwkit contrast` | Selected predictor traits must uniquely identify every joined `branch_clade_id`; event taxa and signed orientation must exactly match the response table. |
| `pgls --response-sampling-covariance` | `tree_id`, `trait`, `contrast_id_1`, `contrast_id_2`, `sampling_covariance` | Contains one upper-triangle row for every selected response-contrast pair. Each model matrix must be complete, finite, symmetric by construction, and positive semidefinite. It is required when response contrasts contain replicate metadata. |
| `pgls --predictor-sampling-covariance` | `tree_id`, `trait`, `contrast_id_1`, `contrast_id_2`, `sampling_covariance` | Contains one complete upper triangle per predictor from replicate-aware species-tree contrasts. The predictor contrast table must also contain positive `contrast_variance`; together they define the phylogenetic prior and observation-error covariance of the latent predictor. |
| `table2nwk` input | `branch_id`, `parent` | Optional `name`, `dist`, and `support`; exactly one root has `parent = -1`. |

All supported input TSV options may use `-` for standard input, subject to the
single-standard-input rule.

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
- Conventional `pgls --tree/--data` reports one row per response coefficient,
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
- Conventional `pgls --model-comparison-out` reports one row per response and
  candidate evolutionary model. It contains the ML log likelihood, AIC, AICc,
  BIC, likelihood parameter count, estimated shape parameter, convergence and
  boundary diagnostics, delta AIC and AIC weight, and delta AICc and AICc
  weight. An individual AICc is empty when its sample size is insufficient;
  delta AICc and AICc weights are empty for the response if any candidate lacks
  AICc. AIC-derived cells remain available.
- Reconciled `pgls` reports one coefficient row per `tree_id`, response, and predictor.
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
- End-to-end `pgls --out-prefix PREFIX` additionally writes deterministic
  `.reconciliation.tsv`, `.gene-contrasts.tsv`, `.species-contrasts.tsv`,
  `.random-effects.tsv`, and `.pgls.tsv` files. Replicate-aware input also
  writes `.response-sampling-covariance.tsv` and `.response-tip-summary.tsv`
  for responses and `.predictor-sampling-covariance.tsv` and
  `.predictor-tip-summary.tsv` for predictors, whenever the corresponding role
  has replicates or known SEs.
  Applicable files are staged and committed as one recoverable transaction,
  with concurrent writers serialized by a per-prefix lock; an ordinary write
  or commit failure does not leave a partial mixed bundle.
  Without `--out-prefix`, intermediate tables remain in memory and only the
  requested primary and random-effect outputs are written.
  Gene and species contrast rows, and every final reconciled result, retain the
  independently selected response and predictor evolutionary models, parameter
  names and fixed or automatically estimated values, and branch-length modes.

## Compatibility names

The following older names remain accepted with a deprecation warning, but
canonical output and new examples use the replacements below.

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

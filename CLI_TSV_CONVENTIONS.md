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

## Tip-keyed TSV files

`annotate --table` and the `--trait` files used by `asr`, `draw`, `monophyly`,
`sample`, and `skim` follow one shared contract:

- The first-class key is `leaf_name`.
- `leaf_name` values are exact tree tip labels, must be non-empty, and must be
  unique. Values such as `001`, `NA`, `NaN`, and `null` are preserved as literal
  labels rather than inferred as numbers or missing data.
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
| `rename --name-tsv` | `old_name`, `new_name` | Both values are non-empty and `old_name` is unique. |
| `--taxid-tsv` | `leaf_name`, `taxid` | `taxid` must be a non-missing integer. |
| `--weight-tsv` | `weight` | Positive finite weights. Optional `tree_id` values are unique, 1-based input-tree numbers; without it, row order is used and row count must equal tree count. |
| `table2nwk` input | `branch_id`, `parent` | Optional `name`, `dist`, and `support`; exactly one root has `parent = -1`. |

All supported input TSV options may use `-` for standard input, subject to the
single-standard-input rule.

## Output TSV vocabulary

- `branch_id` identifies nodes with 0-based level-order numbering. The root is
  `0`; `parent = -1` denotes no parent.
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

# Reusable workflow APIs

`SpeciesParser("extended-taxonomic")` and `--species-parser extended-taxonomic`
retain cultivar groups, subgroups and hybrid binomials. `legacy` and `taxonomic`
keep their existing defaults and behavior. Explicit species maps override all
presets; species identity and the taxonomy lookup query remain separate.

`ParentIndex(parents, root=None, children=None)` accepts caller-owned node IDs.
It validates parent/child connectivity iteratively, offers lazy ancestor chains
and `common_ancestor(left, right)` through the shared binary-lifting LCA index.
It never replaces IDs with traversal indices. `numbered_tree` adapts ape edge
tables and builds all internal clade keys in one operation.

`branch_features` provides clade-rank IDs, branch tables, species coverage and
species-overlap statistics. `map_internal_support_by_split` preserves strict
topology/range/partial-label validation and delegates transfer to the native
kernel. Equal support on repeated unary splits requires explicit opt-in in
the general transfer API; conflicting values remain ambiguous.

`tool_trees` exports topology-only and validated plain GeneRax trees without
rounding branch lengths. `convert --node-label q1` understands ASTRAL support
annotations as well as NHX, and `--comments drop` removes ordinary comments.
Support fields such as quartet frequency and local posterior are not rescaled
or reinterpreted. `calibration_scale` changes calibration ages separately from
branch lengths and leaves tail probabilities unchanged.
GeneRax input and output paths must be distinct; target serialization is checked
before atomically replacing an existing output.

`recphyloxml_records.read_event_records` preserves raw event order, names and
attributes, including incomplete transfer evidence. `expected_clades` checks
the reference topology before event mapping. Use `reconciliation_import` and
`reconciliation_document` for resolved histories with paired transfers instead.

`tree_features` exposes maximal foreground clades, projected clade candidates,
constraint split inclusion and patristic distances. `root_candidates` composes
native MAD/midpoint rooting with an explicit eligible-root collection and
supports a pre-publication input-proof callback.
Constraint verification also checks the complete result tip set and rejects
duplicate labels. Foreground extraction uses linear auxiliary memory, including
for deeply unbalanced trees.

`dated_presentation.render_dated_tree` provides dated-tree layout, geological
intervals, annotation lanes and count panels. Side-panel interpretation is
provided by a loader callback; category labels and colours are configurable.
`r_render` preserves a traditional basic `ggtree` style as an optional R backend.
`posterior_diagnostics.diagnose_chains` uses R `posterior` for multi-chain
rank-normalized split R-hat, bulk/tail ESS and MCSE. R packages are optional
backend dependencies; no diagnostic threshold or model is selected by this API.
Chains must identify distinct files, including through links. The `chain`
profile requires unique `Gen` and parameter columns and finite numeric samples.

Public trait-preparation functions are `contrast.read_mixed_replicate_traits`,
`ordinary_regression.ordinary_predictor_values` and
`ordinary_regression.sampling_covariance_table`. Private spellings remain for
compatibility. Scientific aggregation and multiplicity policies belong to callers.

R clients can submit one JSON array to `python -m nwkit.species_parser --infile labels.json --outfile species.json --species-parser extended-taxonomic`; the batch uses the same preset and preserves input order. Dated presentation ages and geological boundaries are expressed in Ma.

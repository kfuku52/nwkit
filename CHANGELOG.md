# Changelog

All notable changes made after the `v0.21.1` tagged release are tracked here.

## [Unreleased]

## [0.40.8] - 2026-08-31

### Fixed

- Profile feasible zero-variance boundaries explicitly in conditional EIV fits,
  then calculate coefficient information over the remaining free parameters.
  This avoids vanishing log-variance gradients and replaces the globally tighter
  stopping rules from 0.40.7, which caused line-search failures in sparse EIV
  models on Linux. Boundary likelihood and coefficient covariance are checked
  against independent calculations.

## [0.40.7] - 2026-08-31

### Fixed

- Improved conditional EIV gradient accuracy and convergence near a zero
  variance component. Equivalent dense and structured covariance fits now
  agree with an independently profiled likelihood reference, including the
  Windows case exposed by the 0.40.6 validation run.

## [0.40.6] - 2026-08-31

### Fixed

- Reject oversized dense Gaussian searches before allocating trial covariances
  instead of silently reporting lambda=0 when nonzero candidates exceed the
  resource limit. Explicit dense opt-in still explores the complete model.
- Stabilized Gaussian, latent-predictor, and conditional EIV calculations under
  response/predictor offsets, using residual and sampling scales for Gaussian
  variance bounds and reporting exact fits as variance-boundary cases.
- Restored ASR's default equal root prior, preserved numeric internal-node names
  and missing support through the common Newick writer, and made stochastic maps
  invariant to inverse rate/time scaling. Impossible state bridges now fail
  explicitly instead of returning zero transitions.
- Protected drawing inputs and paired image/report outputs against hard-link,
  symbolic-link, and case aliases. Drawing, root comparison, regression, and
  intersection now share recoverable file installation with failure rollback,
  permission preservation, and retained backups when recovery itself fails.
- Removed hurdle-likelihood overflow warnings with stable zero-truncation
  derivatives, and made categorical sampling factors explicitly use sparse
  matrices across SciPy's sparse-array migration.

### Changed

- Moved shared model metadata out of numerical imports, keeping CLI help/version
  startup lightweight. Scalar searches retain at most two full fits, and ordinary
  automatic lambda fitting reuses one validated Brownian covariance.
- Split the drawing renderer into typed setup, rendering, quality, and output
  stages while preserving deterministic SVG/report output.
- Added real parser-to-handler smoke coverage for every subcommand and numerical,
  filesystem-failure, and cache-lifetime regression tests.
- Made `quick` use incremental mypy and skip explicitly marked slow tests, with
  pytest target/option forwarding. Full checks and source CI retain all tests.
- Replaced the average-complexity gate with per-function ceilings and a checked
  baseline updater. Added a reproducible checkout-comparison benchmark runner and
  a dedicated development guide.
- Scoped CI to changed paths while retaining full platform/version coverage for
  dependency changes, weekly/manual runs and release versions; cached the verified
  ETE Windows wheel and validated documentation-only distributions without the
  numerical runtime dependencies.

## [0.40.5] - 2026-08-29

### Changed

- Updated the pinned development tools to diff-cover 10.5.1, Hypothesis
  6.165.10, mypy 2.3.1, and Ruff 0.16.4.

## [0.40.4] - 2026-08-28

### Fixed

- Made MV rooting treat the two displayed root branches of a two-tip tree as
  one physical edge. The root is now placed at that edge's midpoint without
  shortening the tip-to-tip distance, including in `rootcompare` evaluations.
- Counted a direct one-tip tree's branch length in de novo consensus length
  aggregation, making one-tip `mean` and `median` consensus output idempotent.
- Clamped the RF maximum reported by `shuffle` to zero for a one-tip tree, and
  skipped RF calculation for whitespace-only tip labels as well as missing
  labels.
- Treated whitespace-only tip labels as empty wherever unique named tips are
  required. `intersection` and `subtree` now use the same validation as other
  tip-identity-dependent commands.
- Made singleton-node removal also promote one-tip roots, preserving the tip's
  metadata and the complete finite root-to-tip branch length. Consequently,
  `intersection`, `prune`, `sample`, `skim`, de novo `consensus`, `collapse`,
  and `sanitize` now write a direct leaf such as `A:4;` instead of `(A:4);`.
- Collapsed singleton nodes after `constrain --collapse` removes duplicate gene
  tips, so that gene collapsing cannot leave unary internal nodes.
- Normalized singleton root wrappers before midpoint, outgroup, MAD, MV, and
  taxonomy rooting, and before validating a `root --method transfer` source.
  This prevents rerooting from turning a wrapped root into an empty artificial
  tip while preserving pairwise distances and the root stem length.

## [0.40.3] - 2026-08-28

### Fixed

- Kept the source-distribution mode assertion POSIX-only so the reproducibility
  regression test also passes on Windows, where `chmod(0644)` is not represented
  as POSIX permission bits.

## [0.40.2] - 2026-08-28

### Changed

- Updated Ruff to 0.16.2, GitHub Actions checkout to 7.0.1, and setup-python
  to 7.0.0 after validating the selected Dependabot updates against the
  current codebase.

### Fixed

- Normalized source-distribution archive metadata and added a byte-for-byte
  comparison of two independent sdist builds to the release gate.

## [0.40.1] - 2026-08-26

### Changed

- Reduced the README to installation and subcommand navigation, linked
  `rootcompare` to its Wiki guide, and moved non-duplicated workflow details to
  the corresponding Wiki pages.

## [0.40.0] - 2026-08-25

### Added

- Added `rootcompare`, which runs applicable rooting methods by default, writes
  every best position (including ties and automatic-method failures) to TSV,
  and reuses the unrooted `draw` renderer to mark the selected branches in PDF.
  Taxonomy sources are evaluated independently and are enabled automatically
  when species names can be parsed or NCBI taxids are supplied. Physically
  identical endpoint ties are merged while their equivalent split IDs remain
  available in TSV, reference root-edge ratios are retained when usable, large
  trees automatically use equal-angle layout with density-aware labels and
  canvas sizing, and the output files are staged before replacement.

### Fixed

- Midpoint rooting now places the root at the exact tree-diameter midpoint,
  including midpoints that coincide with a branch endpoint.
- Avoided retaining per-edge taxon sets during ordinary MAD and MV rooting,
  restoring linear-memory behavior when tie details are not requested.

## [0.39.0] - 2026-08-22

### Added

- Added `root --method reconciliation` to select the minimum-cost gene-tree
  root over all physical edges using weighted LCA-reconciliation duplication
  and implied-loss counts, without requiring Notung.

### Fixed

- Prevented reconciliation rooting from creating an unnamed tip when the gene
  tree has a singleton root, and made large finite total scores printable
  without floating-point overflow.
- Made LCA duplication/loss scoring safe for one-pass child iterators and
  replaced quadratic reroot-annotation split construction with compact
  bitmask indexing. Reconciliation rooting now avoids a redundant analysis-tree
  copy and per-edge rational-number allocation, and skips split indexing when
  branch annotations or missing lengths do not require it.

## [0.38.0] - 2026-08-18

### Changed

- Reorganized `nwkit regress` help around its three input workflows and the
  response, predictor, reconciled-model, inference, and diagnostic roles. The
  command now selects a workflow only from primary input paths and reports
  mixed or incomplete inputs in that workflow's terms.
- Replaced ambiguous regression options with role-specific names, including
  `--response-contrasts`, `--reconciled-model`, `--predictor-reference`,
  response-prefixed replicate and audit options, and the `event|contrast`
  event-weighting choices. Removed the replaced spellings and hidden
  underscore aliases from `regress`.
- Renamed the earlier reconciled cluster-robust estimator from `legacy` to
  `cluster-hc1` in the CLI, Python API, and output metadata.

### Fixed

- Rejected every explicitly supplied precomputed-mode option that requires
  original tip observations or trees instead of silently ignoring factor
  coding, categorical-replicate, or coefficient-regularization settings.

## [0.37.0] - 2026-08-18

### Changed

- Renamed the phylogenetic regression command from `pgls` to `regress`
  without a compatibility alias, reflecting that it fits Gaussian PGLS as
  well as non-Gaussian phylogenetic GLMMs.
- Renamed the command-facing implementation modules and end-to-end output
  bundle terminology. The primary bundle output now ends in
  `.regression.tsv`, and its transaction lock ends in
  `.regression-bundle.lock`.

## [0.36.5] - 2026-08-15

### Changed

- Rejected SVG document type declarations before XML parsing instead of
  stripping them with an unsafe backtracking expression.
- Corrected the cyclomatic-complexity ratchet to record non-increasable legacy
  renderer ceilings while retaining the strict budget for all other blocks.

### Fixed

- Removed an SVG validation denial-of-service path caused by exponential
  regular-expression backtracking on malformed document type declarations.
- Restored clean formatting and the complete quality and distribution gates.

## [0.36.4] - 2026-08-15

### Changed

- Made sparse reconciled errors-in-variables PGLS use exact, bounded-memory
  precision diagonals at every dataset size, removing a hidden stochastic
  objective change above 512 observations while retaining practical scaling.
- Cached repeated grouped marginal-variance calculations outside optimization
  and used deterministic small posterior factorizations for reproducible
  cross-platform inference.

### Fixed

- Fixed exact posterior-variance diagnostics for sparse factor-loading models
  and strengthened threshold, dense-equivalence, and cache-reuse regression
  coverage.
- Fixed drawing type annotations and made image tests robust to legitimate
  font-rendering differences and unavailable optional Cairo support.

## [0.36.3] - 2026-08-15

### Fixed

- Allowed one shared biological-replicate table to contain both replicated
  traits and traits with exactly one observation per tip. The latter now keep
  ordinary exact-tip semantics instead of failing pooled residual-variance
  estimation.

## [0.36.2] - 2026-08-14

### Changed

- Added versioned variants of the NWKIT logo and updated the README to use the
  logo describing phylogenetic tree processing, visualization, and analysis.

## [0.36.1] - 2026-08-14

### Added

- Added rectangular, slanted, cladogram, circular/fan, radial, unrooted,
  spiral, and fractal tree layouts, with optional tidy subtree packing,
  label-aware spacing, automatic display-only clade collapsing, annotation
  tracks, branch-property styling, scale/depth guides, collision auditing, and
  reproducible layout reports.
- Added common readers for PAML treefile headers, MCMCtree `FigTree.tre` and
  annotated main-output trees, posterior node-age conversion to dated NHX, and
  fixed- or varying-topology DensiTree overlays and empirical path envelopes.

### Changed

- Updated the README, package metadata, and CLI description to present NWKIT as
  a toolkit for phylogenetic tree processing, visualization, and comparative
  analysis rather than only Newick tree processing.

### Fixed

- Fixed batch-adjusted response replicates in conventional PGLS, rejected
  duplicated raw TSV headers, serialized overlapping explicit output sets, and
  retained sparse factor-loading and posterior-precision representations
  throughout reconciled errors-in-variables PGLS.
- Made equal-event Gaussian pseudo-likelihood invariant to uneven identical
  paralog copies, including shared predictor uncertainty; separated and
  whitened lineage-slope variance components by predictor source; propagated
  predictor-dependent covariance into lineage intervals; tightened sampling
  covariance and grouped-uncertainty validation; and bounded categorical
  covariance audits to marginal sparse products for large datasets.

## [0.36.0] - 2026-08-12

### Added

- Added `reconcile` and `contrast` commands for auditable reconciled
  speciation contrasts, stable clade identifiers, explicit sampling coverage,
  GeneRax transfer metadata, and multi-family lineage clustering.
- Added `pgls` for regression on reconciled phylogenetic contrasts, with equal
  species-event weighting and species-event-clustered uncertainty so additional
  paralogs cannot inflate the effective sample size.
- Added opt-in preservation of custom NHX properties during pruning, dropping,
  and sanitization, including lineage-boundary markers for collapsed
  duplication and transfer events.
- Added biological/technical replicate handling, batch adjustment, known
  standard errors, and exact propagation of tip-mean uncertainty through PIC
  transforms and conventional PGLS.
- Added REML/ML reconciled-contrast models with identifiable species-event
  random effects, partially pooled lineage-specific slopes, conditional modes,
  boundary diagnostics, and variance-component-refitted bootstrap inference.
- Added an end-to-end raw-input PGLS mode with optional separate
  annotation-bearing and dated gene trees, transactional inspectable output
  bundles, and in-memory stage reuse.
- Added conventional tip-level PGLS with an intercept, ML/REML fitting, Wald or
  parameter-refitted bootstrap inference, sampling covariance, and per-tip
  estimation audits.
- Added one shared evolutionary-model layer for Brownian, Pagel lambda,
  fixed-root OU, Pagel kappa and delta, early-burst/ACDC, independent, and
  named custom covariance models.
- Added fixed or estimated generic shape parameters to conventional PGLS,
  including parameter-refitted bootstrap inference and generic model metadata.
- Added conventional ML evolutionary-model comparison with log likelihood,
  AIC, AICc, BIC, both AIC and AICc weights, parameter counts, convergence,
  and boundary diagnostics.
- Added model-aware transformed-tree contrasts and independently selectable
  gene-expression and species-trait models in end-to-end reconciled PGLS.
- Added automatic reconciled-PGLS shape-parameter estimation: response-specific
  gene parameters maximize the reconciled likelihood, predictor-specific
  parameters use species-trait marginal ML, and parametric bootstrap refits
  automatically estimated gene parameters from simulated tip values.
- Added a mathematical guide to reconciled speciation contrasts, including the
  reconciliation criterion, PIC recursion, covariance propagation, repeated-
  paralog weighting, and a hand-calculated minimal example.
- Added biological/technical replicates, batch adjustment, and known standard
  errors for PGLS predictors as well as responses. Conventional and reconciled
  PGLS now condition on a phylogenetic latent predictor and propagate its
  posterior covariance through the regression instead of treating noisy
  species means as exact.
- Added typed categorical and ordered predictors and responses, including
  replicate-aware factor uncertainty and shared species effects that prevent
  paralog pseudoreplication in reconciled tip-level PGLMMs.
- Added Poisson, negative-binomial, zero-inflated, hurdle, Gamma, lognormal,
  beta, beta-binomial, and censored-Gaussian phylogenetic likelihoods, with
  offsets/trials/censor bounds, likelihood-level biological replicates,
  coefficient regularization, likelihood/profile inference, and
  family-specific parametric bootstrap.
- Added equally weighted gene-tree/reconciliation ensembles with Rubin
  within/between-tree coefficient uncertainty and explicit support metadata.
- Added ML/REML multivariate Gaussian PGLS for conventional and reconciled
  tip-level models, full response covariance components, and observed-data
  likelihoods for partially missing response vectors.
- Added lineage-slope intervals and reliability, likelihood-ratio/bootstrap
  tests for heterogeneous or jointly nonzero lineage effects, lineage
  leave-one-out refits, and Mk stochastic-map trait-origin sensitivity
  diagnostics for categorical predictors.

### Changed

- Made complete daughter-clade sampling the safe default for reconciled
  speciation contrasts while retaining an explicit partial-coverage sensitivity
  mode, and use child-order-independent clade IDs for cross-file joins.
- Made the hierarchical Gaussian REML model the PGLS default, with the earlier
  cluster-HC1 estimator retained as an explicit sensitivity model.
- Require non-empty gene-family IDs and complete replicate covariance inputs,
  and record replicate, covariance, random-effect, planned-path, and raw-tree
  provenance in the audit trail.
- Made PGLS a strict three-mode command: conventional tip-level data,
  end-to-end reconciled raw input, or precomputed reconciled contrasts.
- Retained `reconcile` and `contrast` as inspectable low-level commands while
  allowing the raw-input PGLS command to run the complete workflow.
- Replaced model-specific conventional PGLS correlation options with
  `--evolution-model` and `--evolution-parameter` without compatibility aliases.
- Record the exact evolutionary model, parameter, and branch-length mode in
  contrast tables and both model choices in all reconciled PGLS results;
  precomputed mode now requires and validates this metadata and one predictor
  tree ID instead of silently accepting mixed transforms.

### Fixed

- Parse and validate GeneRax `S`, `D`, and `H=Y@source@destination`
  annotations, preserving transfer placement instead of silently replacing it
  with descendant-species LCA.
- Reject non-finite contrasts, invalid reconciliation joins, incomplete or
  non-positive-semidefinite sampling covariance, confounded replicate designs,
  and non-identifiable requested random effects.
- Keep technical replicates out of biological sample-size counts and reject
  technical observations spanning batches.
- Stage and commit complete output bundles under a per-prefix lock, restoring
  pre-existing outputs after a failed commit instead of leaving mixed results.
- Retain unexpected numerical warnings while suppressing only the known SciPy
  finite-difference warning from non-finite variance-component trial points.
- Validate custom covariance names, dimensions, numeric values, symmetry, and
  positive definiteness before fitting, and prevent any PGLS output from
  replacing that input.
- Reject non-ultrametric trees for delta and contrast-based OU transforms,
  invalid fixed parameter domains, and
  information-criterion calculations with insufficient sample size.
- Search multiple deterministic parameter basins before local refinement when
  estimating conventional PGLS evolutionary parameters, and reject raw-only
  options in conventional mode even when their explicitly supplied value equals
  the reconciled default.
- Reject ambiguous replicate-column roles, colliding or incomplete gene-tip
  species mappings, mixed precomputed transform metadata, and output/input path
  aliasing before analysis begins.

## [0.35.2] - 2026-08-10

### Added

- Added one local/CI validation entry point with quick, full, distribution,
  and release modes, plus an exact development-tool constraints file.
- Added Hypothesis checks for FASTA preservation, chunked Newick parsing, and
  portable filenames, and an 80% changed-line coverage gate for pull requests.
- Added grouped monthly Dependabot updates for Python and GitHub Actions.

### Changed

- Enabled Ruff import and bugbear rules and type-checked function bodies across
  the package, with full annotation enforcement starting in the FASTA and CLI
  convention modules.
- Removed every grade-F cyclomatic-complexity block by splitting MAD rooting,
  rendering, annotation, validation, and property transfer into focused helpers;
  the maintainability gate now prevents grade-F regressions.
- Split pure image metadata, license, search, and media-type logic from the image
  orchestration module while preserving its public imports.
- Deduplicated CI work, added cancellation and job timeouts, and made local,
  test, build, and release validation use the same constrained commands.

### Fixed

- Explicitly dispatch the release workflow after an automated tag is created,
  because tags pushed with the GitHub Actions token do not emit another
  workflow-triggering push event.
- Removed process-global input-format state from tree serialization so one read
  cannot change the inferred output format of another tree.
- Expanded ignored generated artifacts and made distribution checks reject
  stale modules while verifying reproducible wheel contents.

## [0.35.0] - 2026-08-09

### Added

- Added branch-coverage, whole-package type, dependency-vulnerability,
  static-security, and maintainability-baseline checks to CI.
- Added current-Python test coverage on macOS and Windows alongside the full
  Python 3.10-3.14 Linux matrix.

### Changed

- Hardened SVG parsing with `defusedxml` and documented the native Cairo setup
  required by optional CairoSVG processing on common platforms.
- Made release wheels reproducible by deriving archive timestamps from the
  release commit and comparing direct and source-distribution-built wheels.

### Fixed

- Replaced recursive ETE deep copies with a property-preserving iterative tree
  copy, allowing rooting, root transfer, and shuffling to handle deeply
  unbalanced valid trees beyond Python's recursion limit.
- Converted missing native Cairo-library failures into actionable NWKIT errors
  instead of exposing an import-time `OSError`.

## [0.34.4] - 2026-08-09

### Changed

- Reorganized the largest root, image, and utility test modules by concern and
  centralized shared test factories and fakes.
- Reduced the full test-suite runtime by replacing unnecessary taxonomy-cache
  validation, subprocess startup, file-lock polling, and process-pool startup
  with focused deterministic test doubles.
- Added strict pytest marker configuration and slow-test duration reporting,
  while limiting clean-source-distribution reruns to artifact-focused CLI and
  documented-example smoke tests.

## [0.34.3] - 2026-08-09

### Changed

- Reduced large-tree runtime and memory use in annotation aggregation, ASR,
  property transfer and composition, tree-collection processing, path
  distances, taxonomy constraint construction, automatic parsing, and MAD
  rooting.
- Deferred optional table and HTTP imports, streamed consensus and clade
  frequency inputs in one pass, and avoided process-pool startup for small
  collections.

### Fixed

- Preserved stable path distances and exact MAD root positions when float
  acceleration would lose significance, while retaining fast paths for
  well-conditioned inputs.
- Preserved large-integer annotation sums, skipped unused numeric validation
  on single-tip trees, and retained the former nonredundant taxonomy-helper
  topology.

## [0.34.2] - 2026-08-06

### Changed

- Restricted automated release tags, GitHub Releases, and downstream Bioconda
  updates to major and minor versions whose patch component is zero.
- Made release tagging wait for the full test workflow and use its exact tested
  commit.

## [0.34.1] - 2026-07-31

### Added

- `consensus --comparison rooted|unrooted` makes clade-versus-split semantics
  explicit and supports root-position-invariant unrooted consensus trees.
- Release validation now tests the source distribution and installed wheel,
  including CLI round trips and embedded tip-image rendering.

### Changed

- MAD rooting now uses an independently implemented tree-DP algorithm that
  preserves the established result while reducing large-tree runtime from
  cubic behavior to approximately quadratic scaling.
- Ranked and maximum-PD sampling use exact internal arithmetic so very small
  gains remain distinguishable beneath extremely long shared branches.
- Pillow is a core dependency for safe raster validation; CairoSVG remains an
  optional dependency for SVG rasterization.
- Source distributions include the test suite, while generated caches and the
  removed legacy MAD implementation are excluded from release artifacts.

### Fixed

- Consensus and clade-frequency calculations now handle extreme finite weights
  and branch lengths without overflowing, preserve arbitrary root positions,
  and reject unrooted input where rooted clade frequencies are required.
- Midpoint, minimum-variance, MAD, outgroup, and transferred rooting preserve
  root stems, annotations, missing lengths, and finite scale across extreme
  values; duplicate tip labels are rejected before split-keyed restoration.
- Tree comparison treats the two halves of a bifurcating root as one physical
  edge and reports length, annotation, rooting, and exit-status differences
  consistently; branch, path, and delta metrics remain stable at the finite
  floating-point limits.
- ASR preserves custom ambiguous-state separators, emits collision-free state
  identifiers, and rejects failed or non-finite optimizer results.
- Unrooted monophyly reports use ETE's actual foreign leaves; validation,
  rescaling, table conversion, sanitization, and composition now preserve
  finite, parseable, and format-compatible output at boundary values.
- Newick collection parsing recognizes single- and double-quoted labels,
  escaped quotes, nested comments, and semicolons inside labels or comments;
  failed serialization restores the caller's tree unchanged.
- Root-edge property transfer preserves tiny target ratios, distributes huge
  representable half-edges without overflowing their physical total, and
  rejects only derived values that cannot be represented finitely.
- Phylogenetic-diversity sampling rejects negative lengths and stops before
  writing either output when a gain, total, or pruned edge would overflow.
- Intersection writes are transactional and descriptor-backed, reject special
  output files and aliases, retain output modes where supported, and roll back
  partial multi-file commits.
- Image candidate normalization, license classification, attribution, filename
  generation, cache refreshes, media-type detection, SVG/CSS sanitization, and
  raster processing now handle malformed, adversarial, and concurrent inputs
  deterministically.
- TimeTree and taxonomy downloads stream bounded responses, validate encodings
  and content lengths, and close resources on every success and failure path.
- Draw rejects impractical tip-image sizes and always releases Matplotlib
  figures, including when rendering or output fails.

### Security

- Provider media requests enforce HTTPS host allowlists, DNS/IP checks,
  redirect limits, cross-origin header stripping, response-size limits, and
  media-appropriate `Accept` headers to reduce SSRF and content-confusion risk.
- Provenance and cache locking now resist symlink, inode-replacement, FIFO,
  path-collision, oversized-input, and stale-owner races; Windows process
  liveness checks no longer use the terminating `os.kill(pid, 0)` behavior.
- ETE taxonomy redirects remain on the expected NCBI HTTPS host, traversal
  caches use a restricted unpickler, and taxonomy archives are bounded and
  checksum-verified before atomic installation.
- Direct wheel builds remove stale deleted modules from prior `build/` caches,
  and CI inspects wheel and source archives to prevent legacy code or notices
  from re-entering release artifacts.
- Requests and urllib3 minimum versions exclude known-vulnerable legacy
  releases.

## [0.34.0] - 2026-07-23

### Added

- `draw --tip-image-manifest` renders images and silhouettes from an `image`
  manifest in an aligned tip column, with configurable size and padding,
  offline path resolution, unmatched-tip policy, and image-asset provenance.
- `draw` supports configurable figure geometry, typography, branch colors,
  aligned or branch-end tip labels, root markers, categorical tip badges,
  filtered node-probability pies, property labels, legends, and transparent
  output.

### Changed

- PhyloPic selection prefers the exact taxon's curated `primaryImage`, then
  ranks same-rank fallbacks by license openness, vector availability, drawable
  aspect ratio, and resolution before considering genus or family fallbacks.
- Image manifests record `is_primary`, `is_vector`, and `selection_reason` so
  representative-image choices are auditable.

## [0.33.0] - 2026-07-22

### Added

- Taxonomy caches now support configurable freshness checks through
  `--taxonomy-cache-max-age-days` and explicit refreshes through
  `--refresh-taxonomy-cache`.

### Changed

- `image` now validates downloaded media, avoids filename collisions, bounds
  and refreshes query caches, indexes NCBI image metadata in SQLite, and
  records complete attribution for reused local files.
- Consensus and clade-frequency collection stream tree files, ASR stochastic
  maps aggregate simulations incrementally, and descendant-taxon caches use
  compact representations for better behavior on large or unbalanced trees.
- Provenance records include composition-manifest dependencies and every CLI
  role associated with a physical input or output file.

### Fixed

- ETE4 taxonomy downloads and database builds are checksum-verified, bounded,
  atomic, integrity-checked, and isolated from the caller's working directory.
- Multi-output commands reject paths that resolve to the same file instead of
  silently overwriting an earlier output.
- MCMCtree calibration values and TimeTree responses reject non-numeric,
  non-finite, inverted, or out-of-range values.
- OpenTree and TimeTree requests validate response shapes and sizes, retry
  transient failures, and preserve configured taxonomy-source fallbacks.
- Cache and audit locks now tolerate incomplete lock metadata, verify ownership
  before release, and serialize concurrent JSONL audit appends.
- Audit stderr capture is bounded while retaining warning records.

## [0.32.0] - 2026-07-21

### Added

- `dist` now calculates RF, normalized RF, weighted RF, branch-score,
  topological path, and branch-length path distances. `--metric all` is the
  default, and `--comparison rooted|unrooted` controls clade- or split-based
  metrics.

### Changed

- `dist` uses repeatable or comma-separated `--metric` values and writes one
  metric per row in a stable long-form TSV. The former `-d/--dist` spelling
  remains accepted with a deprecation warning; `-d/--dist RF` retains its
  historical two-column output for pipeline compatibility.

## [0.31.0] - 2026-07-21

### Added

- A documented shared contract for CLI spelling, standard input, tip-keyed TSV
  validation, missing values, unmatched rows, branch IDs, node classes, taxon
  sets, and empty-result headers.

### Changed

- Canonical multiword long options now use kebab-case throughout. Historical
  snake_case and other replaced spellings remain compatibility aliases and
  emit a warning naming the canonical replacement when used.
- Tip-keyed tables used by annotation, ASR, drawing, monophyly checks, and
  sampling now share `leaf_name`, duplicate, missing-value, and unmatched-row
  handling through `--missing-values` and `--unmatched`.
- Node reports now use 0-based level-order `branch_id` values, `node_class`, and
  the `*_taxa` / `num_*_taxa` vocabulary. ASR uses `leaf` rather than `tip`,
  including for canonical `--target` values.
- Image-specific mappings are named `--species-name-tsv`; image output paths
  are `--manifest-out` and `--attribution-out`; sample tables use `--report`;
  and skim group tables use `--group-table-prefix`.
- Auxiliary outputs require file paths, and at most one input may read from
  standard input in a command.

### Fixed

- `image` now honors `--species-parser` and `--species-map-tsv`.
- `info` and `printlabel` now honor `--outfile`; non-Newick commands no longer
  advertise an unused `--outformat`.
- `table2nwk` preserves node names such as `NA`, `NaN`, and `null`.
- Empty TSV results retain stable headers.

## [0.30.0] - 2026-07-17

### Added

- `transfer` and `compose` now accept repeatable
  `--root-edge-policy TARGET_PROPERTY=POLICY` options for deterministic
  bifurcating-root ambiguity handling. Available policies are `auto`, `skip`,
  `equal-only`, `matching-side`, `mean`, `min`, `max`, and `edge-total`.
  Composition manifests support both top-level policy mappings and per-property
  overrides.
- Transfer reports now include the selected root-edge policy, resolution,
  target/source candidate counts, and a deterministic JSON record of every
  source candidate value.

### Changed

- Conflicting support, names, and NHX values on the two halves of a source root
  edge now follow the half with matching projected descendant taxa by default.
  Unique and equal values retain their edge-wide behavior. Numeric reducers and
  conservative rejection remain explicitly selectable.
- Ambiguous root-edge branch lengths now use the source edge total while
  retaining the target root-position ratio. If the target edge total is zero,
  the length is divided equally. An explicit `matching-side` policy can instead
  retain source-side lengths.

### Fixed

- Root transfer and the root alignment used by `compose` and `transfer` now
  preserve internal-node support, names, and NHX properties on their canonical
  unrooted branch splits instead of leaving them attached to nodes whose
  descendant taxa change during rerooting.
- `compose` annotation and value sources no longer change the target rooting
  when `--root-source` is absent.
- Failed root transfers no longer leave the target tree partially rerooted.
- All `root` methods now preserve internal support, names, and NHX properties
  on canonical branch splits, preserve root and tip metadata, and leave their
  input tree unchanged. `root` output also retains custom NHX properties.
- Collapsing a singleton root no longer loses names or custom properties.
- `validate --require-same-rooting yes` now compares the actual root
  bipartition for trees with identical leaf sets, rather than only comparing
  rooted/unrooted status.
- When branch lengths are transferred from a differently rooted tree, the
  source now supplies the total length of the aligned root edge while the
  target/root source retains its root-position ratio. This avoids replacing an
  intentional root position with the rerooting library's temporary 1:1 split.
- Split-based length transfer now recognizes the leaf and internal children of
  a root placed on a pendant edge as two halves of one physical edge. This
  prevents half-edge lengths from being lost when target and source rootings
  differ.
- `diff --comparison unrooted` now counts each canonical root edge once and
  selects terminal/internal splits independently of the displayed root, so a
  pure reroot no longer appears as an unrooted topology difference.
- Split-based transfer now treats equal support, names, and NHX values on the
  two children of a bifurcating root as one root-edge value.

## [0.29.0] - 2026-07-17

### Added

- Rooted-clade or root-independent canonical-split matching for `transfer` and
  `compose`, selected with `--match-basis` and recorded in TSV reports.

### Changed

- Partial-taxon correspondences are now reported as projected rather than
  exact matches. Projected support and branch-length transfer is blocked by
  default and requires `--allow-projected-values yes`; strict policy rejects
  all projected matches.
- `intersection` now uses a built-in FASTA reader and writer that preserves
  retained records verbatim; non-FASTA sequence formats are no longer accepted,
  and Biopython is no longer a runtime dependency.

## [0.28.0] - 2026-07-17

### Added

- `compose` for assembling a topology, root, node names, support values, branch
  lengths, and arbitrary NHX properties from compatible source trees, with an
  optional per-clade provenance report and JSON manifest.
- `diff` for rooted-clade or unrooted-split comparisons that report leaf-set,
  root, topology, value, and arbitrary-property differences as TSV.
- `annotate` for joining tip metadata to trees and aggregating categorical or
  numeric values onto internal nodes.
- A shared `--audit` option that appends command arguments, input and output
  hashes, detected input semantics, random seeds, warnings, and runtime status
  as JSON Lines.
- Auditable clade mapping across partially overlapping tip sets through
  `--taxon-mode intersection` for `transfer`, `compose`, `diff`, and transferred
  roots.
- Arbitrary NHX-property transfer and renaming with `--property` and
  `--property-map`, plus strict or compatible-only transfer policies and TSV
  reports.

### Changed

- Root alignment used internally for annotation and branch-value mapping now
  preserves source branch lengths instead of redistributing an already matching
  root edge.
- Named-tree parsing preserves support values explicitly stored as NHX
  properties, allowing names and support to coexist in composed output.

## [0.27.0] - 2026-07-15

### Added

- A reproducible `--seed` option for the `shuffle` and `skim` commands.
- Wheel-content, installed-wheel smoke, lint, and current-Python checks in CI.
- A tag-driven GitHub Release workflow that validates the tag against the package version.
- Standard `nwkit` and `python -m nwkit` entry points with `--version` support.

### Changed

- The minimum supported Python version is now 3.10; CI covers Python 3.10 through 3.14.
- Packaging metadata and command entry points are declared in `pyproject.toml`.
- Runtime modules use explicit imports instead of wildcard imports from `nwkit.util`.
- TimeTree requests and ETE documentation links use HTTPS.

### Fixed

- Test and legacy script packages are no longer included in wheels.
- Tests no longer silently skip when untracked legacy fixture data is absent.
- `mark` output format selection no longer depends on process-wide `sys.argv` state.
- The image-provider User-Agent now follows the package version.

[Unreleased]: https://github.com/kfuku52/nwkit/compare/v0.35.0...HEAD
[0.35.0]: https://github.com/kfuku52/nwkit/compare/8d1250764cc5dfc6620655bf2fa3082d6acb59d9...v0.35.0
[0.34.4]: https://github.com/kfuku52/nwkit/compare/557c461e0062631b13b97f9cc4ad66396e6363c1...8d1250764cc5dfc6620655bf2fa3082d6acb59d9
[0.34.3]: https://github.com/kfuku52/nwkit/compare/d44ba92b697039ad83759c62070e45bfd5de1367...557c461e0062631b13b97f9cc4ad66396e6363c1
[0.34.2]: https://github.com/kfuku52/nwkit/compare/c737e58c10cf67933b1ad0d8c0e231bbdc4b5ed9...d44ba92b697039ad83759c62070e45bfd5de1367
[0.34.1]: https://github.com/kfuku52/nwkit/compare/5f7e41f1342216726cb3595f80861cec9020ac18...c737e58c10cf67933b1ad0d8c0e231bbdc4b5ed9
[0.34.0]: https://github.com/kfuku52/nwkit/compare/846a53825e6536c61aee80547ad0d4befa0501a5...5f7e41f1342216726cb3595f80861cec9020ac18
[0.33.0]: https://github.com/kfuku52/nwkit/compare/9b43aef933bda115fbb844e8f824d819f0408c2c...846a53825e6536c61aee80547ad0d4befa0501a5
[0.32.0]: https://github.com/kfuku52/nwkit/compare/bb0af9f322948f608dd566a12cf97be5dfc150dc...9b43aef933bda115fbb844e8f824d819f0408c2c
[0.31.0]: https://github.com/kfuku52/nwkit/compare/ffa752a904f1f2ca2c0310fff45aa8943bea6fa2...bb0af9f322948f608dd566a12cf97be5dfc150dc
[0.30.0]: https://github.com/kfuku52/nwkit/compare/1bcd1ca7be22f8f9a54bd60031bc8f47b0f35743...ffa752a904f1f2ca2c0310fff45aa8943bea6fa2
[0.29.0]: https://github.com/kfuku52/nwkit/compare/9d7353df9c4328c68c51a896c1ffa60b41339338...1bcd1ca7be22f8f9a54bd60031bc8f47b0f35743
[0.28.0]: https://github.com/kfuku52/nwkit/compare/v0.27.0...9d7353df9c4328c68c51a896c1ffa60b41339338
[0.27.0]: https://github.com/kfuku52/nwkit/compare/v0.21.1...v0.27.0

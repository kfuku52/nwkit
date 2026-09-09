# Reconciled gene-tree dating

`nwkit radte` dates a rooted gene tree while **sharing one age parameter for
every species speciation event**. All corresponding paralog nodes receive
exactly that age, including in uncertainty samples. The default native backend
needs neither R, Notung, nor MCMCTree.

This is a new, experimental estimator, not a reproduction of `ape::chronos` or
MCMCTree posteriors. Numerical tests check likelihoods, derivatives, integration,
constraints, and serialization. Broad simulation-based accuracy and interval
coverage validation is still needed before treating it as an established
replacement for a Bayesian dating analysis.

The objective functions and rate integration are derived in [RADTE_MATH.md](RADTE_MATH.md).
Measured workloads and their limitations are recorded in [RADTE_VALIDATION.md](RADTE_VALIDATION.md).

External species-age intervals, fixed/bounded/ensemble comparisons, and separate
input versus conditional uncertainty summaries are described in
[Species-age uncertainty](RADTE_SPECIES_UNCERTAINTY.md).

## Inputs and quick start

```sh
nwkit radte --gene-tree examples/radte/gene.nwk \
  --species-tree examples/radte/species.nwk \
  --species-map-tsv examples/radte/species-map.tsv \
  --reconcile lca --max-age 30 --out-prefix results/family
```

The species tree must be rooted, binary, ultrametric, and measured in time units.
The rooted gene tree must have positive substitution branch lengths. Tip ages
are zero. Species node ages are fixed to the species tree unless supplied in
`--species-node-bounds-tsv` with columns `node`, `age_min`, `age_max`.
Internal species names used in that file must be unambiguous. Time units are
preserved; estimated rates are substitutions per site per time unit.

Choose one reconciliation source:

* `--generax-nhx family.nhx`: retain GeneRax `S` and `D` annotations; transfers
  are rejected because a bifurcating dating model cannot represent them. If its
  species-node names differ from the dated species tree, supply
  `--reconciliation-species-tree generax-species.nwk`. This reads names (including
  numeric labels) in Newick format 1 and maps identical rooted descendant clades;
  only the dated `--species-tree` supplies ages. Incompatible topologies fail.
* `--gene-tree family.nwk --notung-parsable family.parsable.txt`: read existing
  Notung output; no Notung executable is needed.
* `--gene-tree family.nwk --reconciliation events.tsv`: reuse a `nwkit
  reconcile` table or a RADTE `.events.tsv` table.
* `--gene-tree family.nwk --reconcile lca`: run the existing native LCA
  reconciliation directly. Use the standard species-map options. LCA assumes
  duplication/loss and a correct rooted topology; it is not a GeneRax likelihood
  search or a model of incomplete lineage sorting.

Duplications are constrained to their reconciled species branch; duplications
above the species root require an explicit `--max-age`. All time-order and
calibration constraints are retained. Inconsistent or infeasible constraints
fail with an error instead of being silently removed.

## Native statistical models

Without an alignment, branch lengths are treated as exact observations:
`log(b_e / duration_e) = mu + u_e`. Gaussian log rates are marginalized, their
mean and variance are profiled, and a constrained optimizer estimates ages.
The default has independent log rates. `--rate-correlation rho` uses a stationary
Gaussian AR(1) process per gene-tree edge, integrating the root rate. It is not
a continuous-time Brownian clock. This mode conditions on the supplied root
split, branch lengths, topology, and calibration domain; it does not measure
sequence sampling uncertainty.

Add `--alignment aligned.fasta` to use sequence information. The alignment
must contain exactly the gene tips. DNA defaults to GTR and proteins to LG,
with four gamma categories. Explicit choices include JC69, HKY, F81, Poisson,
and LG-F. HKY kappa, GTR relative exchangeabilities, and gamma shape are fitted
on the unclocked tree unless fixed with their corresponding options. Frequencies
and other fitted substitution parameters are conditional in the dating step.
Use an explicit model when a protein alphabet can also be interpreted as DNA.

```sh
nwkit radte --gene-tree family.nwk --species-tree species.nwk \
  --species-map-tsv species-map.tsv --reconcile lca --max-age 100 \
  --alignment family.fasta --substitution-model hky \
  --uncertainty bootstrap --bootstrap-replicates 100 --out-prefix results/family
```

### Codon alignments

Select `--substitution-model gy94`, `ecmk07`, or `ecmrest` explicitly for CDS.
DNA/codon content cannot be distinguished reliably by alphabet or sequence length,
so the standalone command does not infer a coding frame from DNA input.

| Model | Exchangeabilities | Default frequencies |
|---|---|---|
| `gy94` | One-base changes, transition multiplier κ and nonsynonymous multiplier ω | `f3x4` |
| `ecmk07` | Published unrestricted Kosiol et al. (2007) matrix, including multi-base changes | Published matrix frequencies |
| `ecmrest` | Published restricted matrix; direct changes differ at one base | Published matrix frequencies |

GY94 estimates one shared κ and ω in the unclocked prefit; `--kappa` and `--omega`
fix them. Dating and profile intervals condition on those fitted values. Plain ECM
has no fitted κ or ω multiplier; specifying these controls for ECM is an error.

`--codon-frequencies f` replaces frequencies with observed unambiguous codon counts
plus a 0.5 pseudocount per sense codon (the analogue of `+F`, with explicit smoothing).
`f3x4` uses position-specific base counts, `f1x4` pools positions, and `fq` uses equal
sense-codon frequencies. Base-count modes add 0.5 per base per position, then
normalize products over sense codons. Ambiguous codons contribute a base only when
that position is resolved by all compatible sense codons. `model` selects published
frequencies for ECM and is invalid for GY94.

The first implementation supports **standard genetic code 1 only**, for all three
models. `--genetic-code` rejects other codes. Sequence lengths must be multiples of
three. In-frame `---` and ambiguous codons are integrated over compatible sense
states; partial gaps and codons compatible only with stops are errors, including
terminal stops. Stops must be dealt with explicitly during upstream alignment
preparation. Bootstrap resamples complete codon columns.

All codon branch lengths and rates use **expected nucleotide changes per codon
site**, or that quantity per time unit. Multi-base ECM transitions contribute their
Hamming distance to the normalization. This differs from unit codon-event rate:
for a fixed matrix, convert to event units by multiplying by `-sum(pi * diag(Q))`.
Do not reuse DNA-per-nucleotide branch lengths or unconverted external codon
likelihood summaries as if their units were identical. The unclocked sequence
prefit estimates branch lengths from the alignment in the selected model's units.
The manifest records the unit, state order, genetic code, frequencies, and fitted
parameters. Independent IQ-TREE tests compare fixed-tree ECM likelihoods after
converting its codon-event branch units.

```sh
nwkit radte --gene-tree family.nwk --species-tree species.nwk \
  --species-map-tsv species-map.tsv --reconcile lca --max-age 100 \
  --alignment family.cds.fasta --substitution-model gy94 \
  --uncertainty profile --interval-level 0.95 --out-prefix results/family \
  --figure-out results/family.pdf
# Use --substitution-model ecmk07 --codon-frequencies f for empirical exchangeabilities
# with frequencies estimated from this alignment.
```

The original ECM numerical data are packaged with their source recorded in each
file: [Goldman laboratory supplementary material](https://www.ebi.ac.uk/research/goldman/empirical-codon-models/),
Kosiol, Holmes & Goldman (2007), *Molecular Biology and Evolution* 24:1464–1479.

The sequence likelihood uses scaled pruning and analytic branch derivatives.
The two edges adjacent to the root are combined into their unrooted length
before constructing a full-covariance quadratic approximation in log lengths.
Native marginal inference integrates Gaussian log rates analytically conditional
on one root-rate contrast, then integrates that contrast by adaptive-order
Gauss-Hermite quadrature. It estimates the log-rate mean and SD together with
ages. This construction avoids treating a reversible likelihood's root split
as an independently observed branch length.

The default `--inference auto --likelihood auto` uses this marginal method when
the quadratic approximation is usable and passes exact likelihood and
standardized-score checks at fitted rate states. If the approximation is
unavailable or fails those checks (including checks during profile-interval
exploration), it refits the point estimate and requested profile intervals using exact sequence likelihood
and conditional joint MAP. That fallback estimates rate SD from the non-root
branches and holds it fixed; it is a different estimator, explicitly recorded
in diagnostics and the manifest. Numerical quadrature failures are reported.
`--inference marginal` forbids fallback. `--likelihood exact` selects exact
conditional joint MAP. No native output is labeled an MCMC posterior.

`--rate-sd` fixes log-rate SD. Zero specifies the strict-clock limit for
sequence inference. Two-tip sequence dating requires this option because rate
variance cannot be estimated from one unrooted length.
In tree-only mode, zero SD requires branch lengths consistent with a strict
clock under the chronology constraints; inconsistent inputs are rejected.
Rate-bootstrap intervals are unavailable at the strict-clock limit.

## Uncertainty and diagnostics

`--uncertainty none` is the default. Optional methods are:

* `laplace`: conditional curvature intervals, refused at active bounds,
  singular information, or the strict-clock limit. This is an unadjusted normal
  approximation; estimating rate variance from a few branches can give severe
  undercoverage even with a long alignment.
* `studentized`: small-sample adjusted curvature intervals. When rate variance
  is estimated, this uses residual branch degrees of freedom, an `n/df` variance
  correction, and a Student t critical value. Intervals use a logit transformation
  of the feasible age domain, preserving all hard bounds without clipping.
  Supplied `--rate-sd` uses a normal critical value without variance inflation.
  This is an approximate alternative, not a posterior or a guarantee of 95%
  coverage. Point ages and rates are unchanged. Active bounds, singular
  information, the strict-clock limit, and absent residual degrees of freedom
  still make intervals unavailable.
* `profile`: conditional likelihood-ratio intervals; endpoints limited by
  calibrations are explicitly marked. These use asymptotic reference thresholds.
* `bootstrap`: site resampling with an alignment, or parametric Gaussian-rate
  simulation for branch-only input. It refits each replicate. Marginal sequence
  replicates must pass exact approximation checks; fewer than 90% successful
  fits (and fewer than 20) makes intervals unavailable.
* `input-ensemble`: refit supplied gene-tree and/or species-chronogram samples
  as described below. Percentiles describe variation between conditional fits.

Laplace, studentized, profile, and site/rate bootstrap condition on the supplied topology,
reconciliation, and species calibration domain. Shared ages reduce parameter
count but do not ensure identifiability. Inspect nonunique/local-optimum,
boundary, approximation, and interval diagnostics. Fixed species calibrations
are assumptions, not independent evidence of estimation accuracy. `--starts`,
`--maxiter`, and `--seed` control reproducible optimization, including each
profile fit. Profile searches move outward from the point estimate and reuse
nearby fitted ages and nuisance parameters. Initial ages are interpolated
within the feasible chronology domain before optimization; constraints and
convergence requirements remain unchanged. Numerical solver failures trigger
intermediate half-steps, with at most 32 failed attempts per requested age.
Approximation-validation failures are propagated immediately. A constrained
profile objective improving on the reference by more than
a scale-aware tolerance aborts interval publication and asks for a point-estimate
refit with more starts; only smaller numerical differences are clamped to zero.
For branch-only residual objectives, this tolerance is `1e-7` times the larger
of the absolute reference objective and the specified rate variance (zero when
unspecified). Sequence objectives use `1e-7 * max(1, abs(reference objective))`.
This prevents a unit-sized tolerance from hiding a materially better fit when
the log-rate variance is very small.

For example, add `--uncertainty studentized --interval-level 0.95` to a native
dating command to request the adjusted intervals. The manifest records
`conditional-studentized-curvature` (or `conditional-bounded-normal-curvature`
with a supplied SD), observation count, residual degrees of freedom, variance
factor, and age transformation in its uncertainty status and diagnostics.
See [the derivation](RADTE_MATH.md#small-sample-curvature-adjustment) and
[the independent-family coverage checks](RADTE_VALIDATION.md#small-sample-interval-validation).

The manifest's `optimizer_attempts` records profile fits with `phase` (`profile`
or `profile-quadrature`), `profile_group`, `profile_age` in input time units,
`successful_starts`, and `objective_spread`, alongside each attempt's success,
objective, iteration count, and solver message. `diagnostics` also flags profile
local optima and nonunique age solutions. Failed runs retain diagnostics on the
in-memory fit and report an error; the existing output bundle is preserved.

When all positive ages can be multiplied by a common factor inside their hard
bounds, dividing the mean rate by that factor gives exactly the same likelihood.
Such runs explicitly report `absolute-scale-unidentified`; their point ages are
representatives of a likelihood ridge, not uniquely estimated dates. Anchor at
least one positive age, or use an appropriate probabilistic calibration analysis
such as the reference backend. Extra optimizer starts cannot resolve this lack
of absolute-time information.

### Input-tree samples

Use `--species-tree-ensemble species-samples.nwk --uncertainty input-ensemble`
to propagate complete species chronograms, preserving correlations among all
node ages. Each Newick must retain the same species topology and tips. The full
chronogram is held fixed during each conditional refit. If a bounds TSV is also
supplied, samples must lie inside its resulting calibration domain.

`--gene-tree-ensemble gene-samples.nwk` supports sampled branch lengths and
topologies with the same gene tips. Use internal LCA for changing topologies,
or fully annotated GeneRax samples when the primary source is GeneRax.
Reconciliation files remain specific to their original clades. When both
ensembles are provided, they must have equal counts and are paired **by Newick
order**; nodes are never sampled independently. No tree files are modified.

The main dated tree keeps the reference-input point estimate. `.age-samples.tsv`
retains input sample indices, and the manifest records failed samples and
inference methods. Intervals require at least 20 successful fits, at least 90%
overall success, and at least 90% representation of a reference age parameter.
Nodes report `sample_event_presence` and `sample_clade_presence`: a shared
species age can remain defined even when a particular gene clade is absent.
Duplication percentiles are conditional on the same reference duplication
clade/event being represented. Missing intervals are not replaced by zero width.

These percentiles propagate the uncertainty supplied by the input ensemble;
they do not additionally integrate within-fit uncertainty or automatically
combine different posterior distributions. If the tree samples and alignment
derive from the same data, interpreting their combination requires care about
double use of information. A single-tree likelihood summary cannot be reused
for a gene-tree ensemble; supply the original alignment or use branch-only mode.

## Files and reuse

Use `--figure-out report.pdf` (or `.svg` / `.png`) to publish a dated-tree and
interval report with the results. Saved bundles can be redrawn using
`nwkit draw --radte-prefix family --species-tree species.nwk --outfile report.svg`
without running inference. See [reconciliation and dating figures](RECONCILIATION_PLOTS.md).

`--out-prefix family` writes `.dated.nwk`, `.nodes.tsv`, `.species.tsv`,
`.events.tsv`, `.shared-ages.tsv`, `.age-samples.tsv`,
`.conditional-intervals.tsv`, `.uncertainty-components.tsv`, `.likelihood.json`,
`.mcmctree-trace.tsv`, and `.manifest.json`. Empty sample/trace tables are normal
when their methods were not requested. The manifest records inputs and hashes,
options, estimator, diagnostics, and output hashes. Outputs are staged together
and rolled back on write failure. Gene names and NHX annotations survive dating;
time branch lengths retain double precision.

The likelihood JSON contains a reusable quadratic summary when available and
an explicit unavailable marker otherwise. Use `--likelihood-summary` with the
same gene topology and identifiers to reuse it. Without the original alignment,
only a local log-length trust-region check is possible; exact validation and
site bootstrap are unavailable. A saved summary is conditional on the fitted
substitution model and alignment used to create it.
Joint-MAP rate-variance estimation uses the summary's unclocked non-root branch
lengths, so reusing a summary does not substitute the input Newick lengths.

## Optional MCMCTree reference

`--backend mcmctree --alignment family.fasta --substitution-model hky` invokes
an installed executable directly, without R. `--mcmctree-bin` selects its path.
Mirrored node labels enforce the same speciation-age groups in every sample.
The backend supports JC69/HKY, separate chains, and explicit burn-in, thinning,
sample-count, rate-prior, and variance-prior options; see `nwkit radte --help`.
`--mcmctree-likelihood approximate` runs PAML's own `usedata=3`/BASEML
precomputation once and then `usedata=2` for each chain. BASEML must be beside
the selected MCMCTree executable or on PATH. The default uses direct sequence
likelihood. Both routes retain their control files and approximation artifacts.
PAML uses soft calibration priors, unlike the native hard intervals, and not
every duplication bound can be selected. Tables and the manifest distinguish
these policies. Unrepresented species events are not posterior estimates.

Scratch directories retain control files and logs. Basic split R-hat and
autocorrelation ESS are reported; convergence requires both their thresholds,
and is not established by a short smoke run. These diagnostics are not
rank-normalized R-hat. Native MAP/marginal estimates and PAML posterior means
need not agree, particularly under weak data or different prior support.

## Reproducible validation workloads

```sh
python tools/benchmark_radte.py --outdir /tmp/radte-benchmark \
  --species 2 8 32 --families 20 --sites 2000 --warmups 1 --repeats 3
# Optional direct and approximate PAML references on the same simulations:
python tools/benchmark_radte.py --outdir /tmp/radte-paml-benchmark \
  --species 2 --families 20 --methods branch sequence paml paml-approximate
```

Use a new output directory, and avoid concurrent CPU-intensive work while
timing. Each family has two copies of a balanced species tree, a true root
duplication age of 20, species-root age 10, and independently simulated lognormal
rates (default SD 0.3). JC69 alignment simulation is independent of likelihood
pruning. Vary `--sites` and `--rate-sd` to assess sensitivity. Every repeat
includes process startup, input parsing, model prefit, inference, and output.
The JSON records commands, environment, failures, root-age errors, interval
availability/coverage, equality checks, diagnostics, wall time, and peak RSS
of the largest process (not total simultaneous process-tree memory).

Use one result per simulated family for bias/RMSE/coverage; timing repeats are
not independent replicates. Report interval availability alongside coverage,
and separate failed PAML convergence diagnostics from converged comparisons.
These balanced root-duplication cases do not establish robustness to losses,
topology errors, internal duplications, or correlated calibration uncertainty.
Additional workloads use `--scenario nested` (at least four species),
`--scenario loss`, `--calibration-width 0.2` (non-root species ages ±20%, root
fixed), `--copy-rate-ratio 3` (a clade-wide rate shift), or `--mapping-errors 1`
(a deliberately wrong tip-to-species assignment). Mapping-error results measure
sensitivity to violated assumptions, not recovery under the stated model.

## References

* Original [RADTE](https://github.com/kfuku52/RADTE).
* dos Reis and Yang (2011), [approximate likelihood dating](https://academic.oup.com/mbe/article/28/7/2161/1051613).
* Le and Gascuel (2008), [LG model](https://doi.org/10.1093/molbev/msn067).
  The bundled numerical LG parameters follow the published model as tabulated
  in IQ-TREE; no IQ-TREE implementation code is included.

The native marginal root-contrast integration above is this implementation's
construction; the references do not establish its accuracy or coverage.

## Optional IQ-TREE sequence engine

The native sequence implementation remains available and is the default.
`--backend native --sequence-engine iqtree` delegates substitution-model fitting,
sequence log likelihoods and branch scores to the `iqtree` executable. NWKIT still
handles reconciliation, shared speciation ages, the clock model, conditional
intervals and report generation. MCMCTree is never executed by this adapter.

```sh
nwkit radte --gene-tree gene.nwk --species-tree dated_species.nwk \
  --alignment cds.fasta --substitution-model gy94 \
  --sequence-engine iqtree --uncertainty profile --out-prefix iqtree_dates

# Complete IQ-TREE model syntax replaces separate frequency/gamma controls.
nwkit radte --gene-tree gene.nwk --species-tree dated_species.nwk \
  --alignment cds.fasta --sequence-engine iqtree \
  --iqtree-model 'GY+F3X4+R4' --uncertainty profile --out-prefix freerate_dates
```

Use the reconciliation/species-mapping arguments appropriate to the input tree,
as in the examples above. `--iqtree-executable` selects an executable path and
`--iqtree-threads` sets worker threads (default 1). `--iqtree-mode persistent`
is the default and requires an IQ-TREE build with `--likelihood-session` support,
plus IQ2MC export support for the initial unclocked prefit. This session mode is
a local IQ-TREE source extension; ordinary released binaries are not assumed
to contain it. The adapter validates the protocol at runtime.
`--iqtree-mode subprocess` explicitly selects the original per-evaluation CLI
route, which requires IQ2MC export support but no session extension.
No fixed upstream version is embedded in NWKIT.

Supported base models are JC, HKY, GTR, F81, Poisson, LG, WAG, JTT, GY, MG,
ECMK07 and ECMrest. Supported modifiers are IQ-TREE's compatible combinations of
`+F`, `+FQ`, `+F1X4`, `+F3X4`, `+I`, `+Gk` and `+Rk` (explicit category counts,
minimum 2). Substitution and rate parameters can be supplied in braces using IQ-TREE syntax;
explicit frequency vectors are not exposed. GY braces
contain omega followed by kappa. A model without gamma/FreeRate has homogeneous
sites. ModelFinder, partitions, ascertainment corrections, nonreversible models
and mixture models are not exposed by this adapter. Codon inputs currently require
standard genetic code 1, and IQ-TREE requires at least three sequences. Unsupported
models or failed evaluations raise errors; they do not select the native engine
implicitly.

Without `--iqtree-model`, the existing substitution/frequency/gamma options build
the IQ-TREE model; GY94 still defaults to F3x4 and four gamma categories. For that
interface, fix both GY94 kappa and omega or estimate both. A complete
`--iqtree-model` cannot be combined with separate frequency/gamma/parameter
controls. It takes precedence over `--substitution-model`.

IQ-TREE frequency estimation follows IQ-TREE's definitions, including its treatment
of absent states; native frequency pseudocounts are not imposed. IQ-TREE branch
rates use its substitution-per-site normalization. In particular ECMK07 counts
codon replacement events, whereas native ECMK07 counts nucleotide changes per
codon; absolute rates from those two engines therefore need a unit conversion.
The model string, frozen parameters, executable version/hash, normalization label
and alignment hash are recorded in the manifest and likelihood metadata.

The adapter first fits the supplied reconciled topology without a clock, freezes
its model parameters and checks that the fixed-model likelihood reproduces the
fit. Tip aliases preserve identifiers and identical sequences. Subsequent calls
fix topology and branch lengths; bipartitions map IQ-TREE's unrooted edges back
to NWKIT, summing the two root branches. Explicit tree-output precision preserves
very short positive branches. A persistent worker returns double-precision branch
scores directly through a first-derivative-only `SCORE` request, without writing
a full Hessian. Explicit diagnostic evaluations can request diagonal curvature
through `EVAL`. The
subprocess route uses IQ2MC scores; its approximate cross derivatives are
**not** treated as the observed Hessian. NWKIT forms the
full log-length curvature by differencing IQ-TREE scores, retaining the existing
exact likelihood/gradient checks and exact profile refit when auto's quadratic
approximation fails. Substitution parameters remain conditional throughout dating.
Bootstrap replicates resample complete codons (or NT/AA columns) and refit the
requested IQ-TREE model.

Marginal profile fits share immutable covariance and likelihood matrices across
age constraints; a changed model, Hessian, mapping or correlation invalidates
that shared structure. This does not alter the profile search or its validation.

After the separate unclocked prefit, the persistent worker loads the frozen model
and alignment once. Repeated evaluations send only unrooted branch lengths over
pipes; alignment patterns, model eigensystems and allocated buffers remain in
memory. Partial likelihoods are invalidated when branch lengths change. Each
bootstrap replicate owns a separately fitted session. A bounded evaluation cache
avoids duplicate requests. Workers are reaped when their likelihood object is
released or explicitly closed, including initialization and protocol failures.
A timeout, malformed response or numerical failure stops evaluation; it never
silently restarts the worker or switches engines. The connection mode and
protocol version are included in the provenance metadata.

To compare the two IQ-TREE connection modes with the same binary in the target
runtime, run `python tools/benchmark_radte_iqtree.py --output benchmark.json`.
The benchmark uses simulated codon alignments with 4/16/64 tips, three repetitions,
alternating mode order and a warmup. It checks likelihood/score equivalence and
records setup time, 20 uncached evaluations and largest-child peak RSS. These
kernel timings do not establish complete dating/profile throughput; measure
that separately on representative inputs.

The subprocess exporter in the tested IQ-TREE 3.1.3 runtime can return an
infinite branch score for extremely short sibling codon edges (observed around
`8e-11` and `6e-11`). The local persistent extension detects spectral
cancellation and recomputes the same fitted model within IQ-TREE using stable
matrix exponentiation and original-state pruning. Regression tests compare
these boundary values and scores against an independent matrix exponential,
including a subsequent return to ordinary lengths. NWKIT still rejects
nonfinite results; it does not floor branches or switch engines. These checks
do not establish accuracy for every possible numerical boundary.

# Preliminary RADTE validation, 2026-09-09

The implementation passes its numerical and input/output contract tests.
Small simulation comparisons support useful speed and point-estimation
behavior, but **do not establish general high accuracy or nominal interval
coverage**. The native estimator remains experimental. In particular,
conditional MAP fallback was common and curvature intervals were often
unavailable in the wider-calibration cases.

Raw measurements and settings are in
[`examples/radte/validation-results.json`](examples/radte/validation-results.json).
The scripts retain individual command lines, inputs, process logs, manifests,
and complete results in their output directories. See [RADTE.md](RADTE.md) for
reproduction commands and [RADTE_MATH.md](RADTE_MATH.md) for the estimator.

## Numerical and behavioral checks

* Closed-form two-tip JC69 likelihood versus scaled pruning; branch derivatives
  versus finite differences for DNA and protein models; transition probabilities
  versus an independent matrix exponential.
* Full root-pair reduction, quadratic gradient/Hessian checks, and local
  approximation validation against exact likelihood.
* Marginal root-rate integration versus independent two-dimensional Gaussian
  quadrature, including posterior rate means; finite differences of age, mean,
  and variance gradients with independent and correlated rates.
* Invariance to the arbitrary split of reversible root-edge lengths; fixed-model
  sequence refits also remove dependence of estimated rate variance on supplied
  starting branch lengths.
* Shared species ages in point estimates, bootstrap refits, paired input
  chronograms, and real MCMCTree mirror samples. Input-ensemble tests check
  correlated species dates and changing gene-clade/event presence.
* Infeasible chronology, transfers, malformed annotations, missing leaves,
  nonpositive lengths, input/output collisions, and rollback of failed writes.
* Actual MCMCTree direct and approximate (`usedata=3` then `2`) runs; posterior
  sample counts, mirror equality, soft-bound labeling, and convergence diagnostics.

## Same-alignment reference comparison

Ten independent simulated four-tip gene families each contain two copies of
a two-species tree. The true duplication age is 20, species age is fixed at 10,
independent log-rate SD is 0.3, and each alignment has 2,000 JC69 sites. Seeds are
1701–1710. This is a deliberately small assessment, not a representative power
study. Branch-only inference receives the simulated true substitution lengths,
so it has different information from the sequence methods.

| Method | Median full-process time | Max process RSS | Duplication RMSE | Interval availability | Coverage among available |
|---|---:|---:|---:|---:|---:|
| Native, exact input branches | 1.25 s | 134 MiB | 3.54 | 10/10 | 9/10 |
| Native, sequence auto | 3.12 s | 152 MiB | 3.08 | 9/10 | 9/9 |
| MCMCTree, direct likelihood | 19.98 s | 252 MiB | 3.98 | 10/10 | 9/10 |
| MCMCTree, approximate likelihood | 4.08 s | 240 MiB | 3.98 | 10/10 | 9/10 |

Times include interpreter startup, reconciliation, all sequence precomputation,
inference, diagnostics, and output. Native measurements have one excluded warmup
and three measured runs per family. PAML reference accuracy uses one run per
family, two chains per run, 2,000 burn-in iterations, 20,000 saved samples per
chain, and thinning by 10. A separate initial comparison repeated the first
families three times after warmup and produced identical seeded posterior means;
the definitive table uses the complete ten-family batch. Peak RSS is the largest
process, not the sum of concurrently resident processes.

All ten direct and approximate PAML runs passed the implemented basic split
R-hat (<1.05) and ESS (≥200) thresholds. All native runs satisfied their hard
constraints. All methods had zero discrepancy between shared ages. Nine native
sequence fits rejected the quadratic approximation and used exact conditional
MAP. The remaining marginal fit reached the numerical variance boundary and
correctly declined a curvature interval.

On this workload native sequence inference was approximately 6.4 times faster
than direct PAML and 1.3 times faster than approximate PAML. These are observed
workflow timing ratios for different estimators, **not equivalent-posterior
speedups**. PAML uses soft priors and posterior means; native inference uses
hard calibrations and conditional MAP or marginal likelihood. The apparent RMSE
difference has substantial Monte Carlo uncertainty with only ten families.
Availability must be reported alongside coverage; nine available intervals do
not establish that a nominal 95% method is calibrated.

Environment: Python 3.10.14 (x86_64), NumPy 1.26.4, SciPy 1.15.2, macOS 26.6.2,
PAML 4.10.10. The JSON records the processor and executable hash. These numbers
should not be transferred directly to another architecture or Python/BLAS build.

## Scaling and sensitivity checks

With 1,000 sites, one family per size and three measured repetitions after
warmup:

| Gene tips | Branch-only time | Sequence time | Sequence max RSS |
|---|---:|---:|---:|
| 16 | 1.28 s | 4.26 s | 163 MiB |
| 64 | 1.27 s | 14.04 s | 246 MiB |

Both sequence cases used exact fallback. The 64-tip likelihood prefit contained
a boundary branch, so a regular quadratic approximation was unavailable. Its
curvature age interval was also unavailable. These two runs are scaling probes,
not an accuracy assessment at those sizes.

Three-family probes additionally covered internal duplications with non-root
species calibrations widened by 20%; losses with a threefold copy-wide rate
shift; and deliberately incorrect tip-to-species mappings. All runs completed
with valid shared ages and retained native chronology constraints. Sequence
duplication RMSEs were 0.38 (true age 7.5, internal duplication), 2.21 (true age
20, losses/rate shift), and 0.79 (true age 20, one wrong mapping).

Curvature intervals were unavailable for all internal-duplication and wrong-map
cases. A separate branch-only profile run on the internal-duplication case
returned calibration-limited intervals. The wrong-map probe evaluates only the
original root duplication; it does not validate newly inferred events or show
that erroneous reconciliation is harmless or automatically detectable.

## Original R RADTE workflow

`tools/benchmark_radte_reference.py` uses the original repository's 13-tip
GeneRax example and species bounds, with `chronos_model=discrete`, lambda 1,
maximum age 1000, and seed 1. Three measured runs after warmup gave median
6.08 s / 148 MiB for original R RADTE and 1.41 s / 133 MiB for native RADTE.
The rooted topology was identical, and both outputs were ultrametric to their
serialization precision. Root ages differed (452.44 versus 485.82).

The objective functions, selected constraints, and output bundles differ;
the R workflow also produces figures. The interval-only calibrations in this
example admit a scale ridge in the native model. Its reported age is therefore
not a uniquely identified estimate. This comparison checks input compatibility
and provides a workflow baseline; it is not an accuracy or numerical-equivalence
claim. Native inference now explicitly diagnoses the feasible common-scale
ridge even with a single optimizer start.

## Remaining scientific limits

Broader, prospectively specified simulation studies are needed across family
sizes, sequence lengths, calibration dependence, rates, and topology errors.
The present evidence is insufficient to promote curvature intervals as general
95% uncertainty estimates or claim that the marginal approximation always
replaces MCMCTree. Existing fast dating programs were not benchmarked as if they
implemented the same duplication and shared-event constraints.

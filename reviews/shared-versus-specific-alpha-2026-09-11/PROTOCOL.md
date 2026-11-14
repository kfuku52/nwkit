# Shared versus trait-specific alpha: prespecified comparison

Both fitted models estimate alpha from the observations and estimate full
cross-trait evolutionary covariance. No fixed generating alpha is passed to
fitting. Use the installed NWKIT in GeneGalleon image
`sha256:c8c69075b5d591499bf3195ec912ef5b5e1aa1cd163f03864639ed6a669e9429`.
The two models are different statistical models; output equality is not expected.
At the same layout, the trait-specific model contains the shared model, so a
materially smaller maximized likelihood is a numerical warning to report.

## Generating models

- 100-tip balanced, exactly ultrametric binary tree, height 1, fixed OU root.
- Alpha shared: all rates 1. Alpha different: geometric spacing from 0.25 to 4.
- Diffusion correlation 0.8. Diagonal diffusion entries are set to give marginal
  process tip variance 1. Off-diagonal **tip** correlations consequently differ
  between the equal/unequal-alpha generating models. This distinguishes a valid
  diffusion covariance from an arbitrary proposed tip covariance.
- Complete observations without sampling or extra measurement error. This exposes
  the shared-alpha separable fast path and is not a benchmark of noisy/missing data.
- Null, aligned shift, opposed shift. Shift clade has 10–30 tips; displacement at
  its tips has norm 2. Aligned shifts spread this over all traits; opposed shifts
  use opposite signs in the first two. Regime optima compensate for attenuation
  separately for each generating alpha, so displacement is held constant.
- Exact tree, parameters, trait/tip order, seed and observed-data digest are
  recoverable from the script and each saved outcome. Both fitted models receive
  the same generated observations for every paired comparison.

## Search pilot

Two traits; ten independent datasets per generating-alpha/scenario cell, hence
60 datasets and 120 model runs. Native group-lasso/beam search: at most one shift,
four candidate branches, five layout refits, beam width four, 200 lasso iterations
and 10,000 screening evaluations. Both models use the same budgets; each produces
its own covariance-aware candidate set. The true branch is never injected into
the search. These finite budgets make discovery approximate.

The same retained fits are scored using AIC and BIC. With at most one shift,
all one-shift layouts have the same parameter count, so retaining the
maximum-likelihood fit per complexity suffices to rerank these criteria. Report
any-shift, exact-branch and false-branch counts, failed runs, timeouts and actual
candidate counts. These criteria are research scores, not calibrated significance
tests. This experiment does not establish type-I error control or bootstrap power.

Default numerical settings are retained: three optimizer starts and 300
iterations. Each search has a 180-second wall-time limit; timed-out/failed runs
are retained, not redrawn or counted as non-detections. Search jobs may run in
four processes. Their wall times are operational diagnostics, not speed estimates.

## Timing

Fixed true layout, opposed shift, 2/5/10 traits, both generating-alpha settings.
For each fitted model: one warmup followed by three measurements on the identical
dataset. Run serially, with one BLAS/OpenMP thread and no concurrent test or
simulation jobs. Include alpha/covariance/mean optimization in the timed region;
exclude imports and data generation. Each individual fit has a 120-second limit.
If warmup or a measured fit times out, report that stage and do not invent a
median. Report likelihood/alpha diagnostics, within-model repeat consistency,
shared-submodel likelihood nesting and whole-worker peak RSS. RSS includes imports
and warmup; it is not allocation volume attributable solely to optimization.

The initial `preflight/` was an execution check, not the experiment: its unequal
alpha fixture incorrectly supplied a tip covariance that implied an invalid
diffusion matrix. Those six fixture failures are retained as development evidence
and excluded from formal comparisons. The corrected fixture supplies valid
**diffusion** covariance directly. No NWKIT inference code is changed for this
comparison. Scripts and outcomes are kept separate from the earlier fixed-alpha
covariance pilot.

## Supplement: observation errors (declared before measurement)

The separable fast path is unavailable when observation errors are present.
After the main search pilot finishes, run a separate serial fixed-layout timing
comparison with two traits, both generating-alpha settings and both fitted-alpha
models. Use the same opposed-shift fixture and seed rule, known sampling SE 0.1
for each tip/trait and generating additional measurement variance 0.04 per trait.
Both fitted models estimate the additional diagonal measurement variance along
with alpha and full evolutionary covariance. No missing coordinates are added.
Use warmup plus three repeats and the same 120-second per-fit timeout. This
supplement measures computation, not noisy-data detection power. It is stored in
its own directory/manifest and does not replace the predeclared error-free pilot.

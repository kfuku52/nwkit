# THRESHOLD posterior diagnostics

The `rank_split_v1` diagnostics assess computation of the identified threshold
posterior. Passing them is not evidence that Brownian evolution, the tree, or
its branch lengths are scientifically correct. An artificial drifting trace
can expose a diagnostic defect without proving a real fitted chain failed.

## Posterior and initialization

The root is N(0,1); each child liability conditional on its parent is normal
with mean equal to the parent liability and variance equal to branch length.
Binary thresholds are zero unless fixed explicitly. Estimated ordinal thresholds
have first threshold zero and **flat joint prior density over the positive
ordered remaining thresholds**. The uniform full-conditional threshold update
uses this prior; it is not a normal or independent-gap prior.

Every state must have an unambiguous observation when thresholds are estimated.
In particular, a highest-category observation bounds the largest threshold above
by that tip's liability. Integrating over K-2 free thresholds is bounded by the
positive part of that Gaussian liability to power K-2, divided by (K-2)!.
All such Gaussian moments are finite. The likelihood is positive on a region
of nonzero volume, so the posterior normalizer is finite and positive for
finite positive branch lengths. The flat threshold prior itself is improper;
it must not be used as a prior-predictive generator for simulation-based
calibration. Fixed-threshold models have a proper Gaussian prior and support
prior-predictive calibration directly.

Chains use independent SeedSequence children and dispersed feasible initial
values. Unconstrained initial liabilities use N(0,4); constrained values are
uniform inside finite intervals or exponential distances into half-lines.
Ambiguous tips select among allowed initial categories. Estimated threshold gaps
are jittered by log-normal factors. These distributions initialize the chain;
they do not replace its target prior. A seed is reproducible within this method
version; results need not match the previous initialization algorithm.

Disjoint observation intervals are selected with log-domain Gaussian interval
masses. Subtracting saturated ordinary CDFs and then choosing a nearest interval
can bias even symmetric tails. Unresolvable masses raise an error rather than
silently choosing a different conditional law.

## Definitions and checks

For a continuous quantity, R-hat is the maximum of rank-normalized split R-hat
and rank-normalized folded-split R-hat. Folding is absolute deviation from the
pooled median. Odd-length chains omit the central draw for splitting, and ranks
use average ranks for ties and the Blom normal-score transform.

Bulk ESS uses rank-normalized split draws. Tail ESS is the minimum ESS of the
indicators below the pooled 5% and 95% quantiles. Autocovariances use FFTs and
multiple lags, pooled within/between-chain variance, and Geyer's initial positive
and monotone sequences. The antithetic correction allows ESS above draw count,
with a stability cap of total split draws times log10(total split draws).

Mean ESS uses untransformed split draws. Mean MCSE is pooled sample standard
deviation divided by sqrt(mean ESS). The second-moment rows diagnose the other
quantity needed for the reported liability variance; their MCSE is for E[X²],
not directly for the variance or SD.

For category k, each draw uses its own thresholds to form I(category=k). Its
mean is the reported probability. Rank split R-hat, bulk ESS, mean ESS and mean
MCSE apply to that indicator. Folded R-hat and tail ESS are inapplicable here:
folding a balanced binary series can produce a constant, and category-probability
precision is a mean-estimation question.

A continuous row passes when R-hat <= 1.01 and bulk/tail ESS >= 400 and all
required diagnostics are finite. A category row also needs mean ESS >= 400
and absolute MCSE <= 0.01. The 0.01 threshold does not ensure small relative
error for rare probabilities. Four chains are recommended. Two are supported;
one chain or fewer than eight retained draws per chain is insufficient for this
contract. More draws do not increase the fixed ESS cutoff. Thinning does not
solve autocorrelation and is one by default.

## Status and TSV schema

`--liability-diagnostics-out` is a file-only auxiliary output with the same
input-overwrite and output-collision protections as other ASR auxiliary files.
The primary category and liability-moment table columns are unchanged.

| Column | Meaning |
|---|---|
| `branch_id`, `name` | Existing branch ID/name; threshold rows use -1/empty |
| `variable` | `liability`, `liability_second_moment`, `category`, or `threshold` |
| `state_or_threshold` | State label or one-based threshold index; otherwise empty |
| `rhat` | Continuous rank/split/fold maximum, or category rank/split R-hat |
| `ess_bulk`, `ess_tail` | Multi-lag bulk/tail ESS; tail is missing for categories |
| `ess_mean`, `mcse_mean` | Mean-estimation ESS and Monte Carlo standard error |
| `status` | `ok` or explicit reason(s), joined with `+` |

`structural_constant` rows include fixed thresholds, exact observed tip
categories and prohibited categories at ambiguous tips. They are excluded from
aggregate checks. `constant_trace` identifies a constant unknown continuous
quantity in any chain or half-chain. `unresolved_rare_category` identifies an
unknown category indicator that is constant in any chain or half-chain; zero
visits cannot establish probability zero. Tail indicators for continuous
quantities may be constant in individual chains and still have pooled ESS;
this does not bypass the checks on the original continuous trace.

Other reasons include `nonfinite_trace`, `insufficient_draws`,
`mcmc_rhat_unavailable`, `diagnostic_unavailable`, `mcmc_rhat`, `mcmc_low_ess`,
and `mcmc_probability_precision`. Missing numbers are not replaced with a
successful R-hat or total-draw ESS. Known structural constants cannot make the
whole fit fail. Unknown unavailable variables cannot make the fit pass.

`--model-out` includes `mcmc_diagnostic_version`, `mcmc_ess_bulk_min`,
`mcmc_ess_tail_min`, `mcmc_probability_mcse_max`, `mcmc_monitored_variables`,
`mcmc_unavailable_variables`, and `mcmc_problem_variables` (branch:variable:label).
`mcmc_rhat_max` and `mcmc_ess_min` retain their column names, but the latter now
means the minimum of available bulk and tail ESS. These extrema must be read
alongside `fit_status` and unavailable counts. They are not comparable without
qualification to old lag-one ESS. No old statistic is used for a pass decision.

## Storage and computation

Liability/threshold traces require about 8 × chains × retained draws ×
(nodes + thresholds) bytes. Up to 64 MiB is stored in memory; larger arrays use
a temporary disk-backed mapping, removed on normal completion or exceptions.
This is an allocation threshold, not an RSS ceiling: the operating system may
keep mapped pages resident. Diagnostics additionally need
workspace proportional to the draws for one variable, and the output table is
proportional to nodes × categories. Category indicators are reconstructed one
node/category at a time instead of storing a full one-hot trace array.

Rank sorting and FFT autocovariances cost approximately O(MS log(MS)) per
variable; every node and nonstructural category is assessed. A large tree may
therefore require substantial scratch disk space and diagnostic time. Sampling
and the diagnostic target set do not silently change to reduce resource use.
Full traces are not persistently exported by this option.

## Independent checks

- `tests/data/threshold_diagnostics.json` contains reference results from R
  `posterior` 1.7.0 for independently generated iid, drift, location/scale
  mismatch, heavy-tail, lag-two, antithetic, tied and odd-length cases. Exact
  input arrays are stored in `threshold_diagnostic_draws.npz` with a recorded
  SHA-256, so future RNG changes do not change a reference test's input.
  Regenerate with `PYTHONPATH=. python tools/validate_threshold_diagnostics.py
  --output tests/data/threshold_diagnostics.json` in an R environment with
  `posterior`. NWKIT diagnostics are not called by the exporter. Tests require
  agreement to relative 1e-6/absolute 1e-8.
- `tests/threshold_posterior_support.py` integrates independently derived
  binary root/internal-node posteriors and a free ordinal threshold posterior,
  including a disjoint ambiguous observation. Integration tolerances are
  tightened independently before comparison with MCMC.
- A [small prior-predictive pilot](examples/threshold/README.md) compares short
  and longer computations on balanced, pectinate, short-branch, polytomy, rare
  ordinal and missing-tip cases. It retains failed diagnostics and root interval
  coverage denominators; it is not a precise calibration experiment.
- Tests also cover stuck nonroot nodes, fixed quantities, missing categories,
  MCSE, seed initialization, scratch-file cleanup and CLI output safety.

The main algorithm reference is [Vehtari et al. (2021)](https://arxiv.org/abs/1903.08008).
Reference diagnostics are provided by [R posterior](https://mc-stan.org/posterior/).

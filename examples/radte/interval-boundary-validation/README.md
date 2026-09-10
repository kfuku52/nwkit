# RADTE variance boundaries and exact contrasts, 2026-09-10

This work adds a true zero-variance boundary to marginal fitting, an exact
branch-only contrast method, and an expanded validation runner. It does **not**
establish a general sequence 95% interval. The default uncertainty remains
`none`; the studentized formula and its historical evidence are retained.

The working checkout started at `8d3a4a487bc57b1363c957e1e2bc1e2bb7fab5d6`.
Each evidence directory records the source hashes actually used. Later
diagnostic classification, CLI naming and reporting edits are not retroactively
assigned those hashes. In particular, the first low-SD pilot calls the zero
boundary an active-bound failure; the final code labels it
`unavailable-estimated-zero-rate-variance`. The numerical comparison predates
the explicit strict-clock guard added to the exact-contrast method; its
positive-SD cells do not exercise that guard, which has a separate unit test.

## What changed

- The marginal model estimates `tau=sigma² >= 0`. Its analytic right derivative
  at zero includes both the Gaussian covariance contribution and curvature of
  the combined-root log length. This removes the old positive `exp(-9)` SD floor.
  Estimated zero variance still makes regular curvature intervals unavailable.
- `nwkit radte --uncertainty exact-log-duration` checks whether a branch-only
  chronology reduces to a linear Gaussian contrast in `log(age-offset)`. In
  that restricted model it uses an exact t pivot with estimated variance, or
  a normal pivot with supplied positive SD and fixed known rate correlation.
  It preserves point estimates and intersects the confidence set with the
  hard domain. Other chronologies and sequence input are rejected explicitly.
- The general calibrated-profile prototype regenerates rates and sites and
  fits marginal likelihoods under each candidate age. Failed refits are not
  discarded. Its finite-grid hull and nuisance plug-in remain approximations.
  It is available to the research runner, **not** the public uncertainty CLI.

The local derivation and restrictions are in [RADTE_MATH.md](../../../RADTE_MATH.md).
MCMCTree was not used as a truth standard; no MCMCTree accuracy/coverage study
was run. The available PAML command-contract tests did execute their small
external inference cases during regression testing.

## Independent branch-only reference experiment

The [protocol](protocol.json) was written before running these new seeds.
It fixed three cells, exactly 1,000 families per cell, a six-edge chronology,
true duplication age 20, species age 10 and maximum age 100. An explicit dense
genealogical covariance and whitened GLS contrast form the independent reference.
There was no coefficient tuning or optional extension after observing outcomes.

| Generating SD / known rho | Exact returned / all | Exact covered / returned | Laplace covered / returned | Studentized covered / returned |
| --- | ---: | ---: | ---: | ---: |
| .1 / 0 | 1000/1000 | 960/1000 | 809/1000 | 961/1000 |
| .3 / .5 | 1000/1000 | 946/1000 | 803/1000 | 945/1000 |
| .6 / .9 | 1000/1000 | 942/1000 | 809/1000 | 943/1000 |

Exact endpoints agreed with the independent reference to at most `3.11e-15`
in normalized time. All methods returned intervals for every family in these
cells, so conditional coverage equals the correct-interval return fraction.
This is a branch-only experiment and does not supply sequence-rate variance
or sequence calibration evidence.

**The prespecified statistical gate did not pass in all cells.** The three
simultaneous one-sided lower confidence bounds for the exact method's correct
return fraction were .9447, .9287 and .9242, compared with the required .93.
The last two fail. The numerical reference and algebraic Gaussian identity are
distinct evidence; neither permits relabeling that empirical gate as passed.
No extra families were generated to make it pass. The 95% procedure is exact
only under its stated mathematical special-case assumptions, not a claim that
all measured cells or broader applications met the adoption criteria.

[Family-level paired CSV](exact-reference/cases.csv),
[summary](exact-reference/summary.json), [tested source hashes](exact-reference/metadata.json),
and [environment](environment.json) are retained. Reproduce into a new directory:

```sh
PYTHONPATH=. OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
python tools/validate_radte_exact_intervals.py \
  --protocol examples/radte/interval-boundary-validation/protocol.json \
  --output /tmp/radte-exact-reference
```

## Sequence and topology development probes

The independent generator uses dense Gaussian draws for branch rates and
`scipy.linalg.expm` for sequence transitions, rather than native rate precision
or pruning simulation. A 20-family low-SD marginal-only pilot used new seeds
92430000–92430019, four gene tips, 2,000 JC69 sites, SD .1 and rho 0.
Only **3/20** point fits passed the exact quadratic-approximation check; all
three estimated variance at zero and had no curvature interval. The remaining
**17/20** are retained as failures. Thus R/N=0 and C/N=0; C/R is undefined.
This is deliberately marginal-only, not the auto fallback workflow of the
historical coverage study.

On family 9, a single-age null probe at the true age 20 generated 19 independent
rate-plus-sequence replicates. **12/19 failed** approximation checks. They counted
as exceedances in the conservative diagnostic Monte Carlo tail, and the rate
of failure exceeds the 10% interval-availability threshold. This probe is not
a calibrated P-value claim, an interval or a coverage experiment. Its role is
to identify why the general sequence prototype is not ready for adoption.

The two-family pectinate smoke probe additionally exercised rho=.5, four
species/eight gene tips, 250 sites and 20% internal species calibration width.
It completed both point fits; studentized returned one interval and Laplace
none. Two cases establish execution, not statistical calibration.

The [low-SD cases](low-sd-marginal-pilot/cases.jsonl),
[null probe](low-sd-marginal-pilot/sequence-null-probe.json), and
[pectinate cases](pectinate-smoke/cases.jsonl) preserve all failures and the
source/argument metadata alongside them. The probe input hashes can be checked
after regenerating the specified seeds; sequences were not copied into this
evidence directory.

```sh
PYTHONPATH=. OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
python tools/validate_radte_intervals.py --generator independent \
  --species 2 --sites 2000 --families 20 --rate-sd .1 \
  --seed 92430000 --inference marginal --output /tmp/radte-low-sd-pilot

PYTHONPATH=. OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
python tools/validate_radte_intervals.py --generator independent \
  --shape pectinate --species 4 --sites 250 --families 2 \
  --generate-rho .5 --fit-rho .5 --calibration-width .2 \
  --seed 92420000 --output /tmp/radte-pectinate-smoke
```

## Expanded runner contract and next acceptance gate

The interval runner now supports independent/legacy generators, balanced or
pectinate trees, root/nested/loss scenarios, separate generating and fitted
correlation/model, copy-wide rate shifts, calibration widths, branch-only input,
explicit inference/likelihood selection, multiple named targets, and separate
optimizer seeds. `truth.json` may contain `target_ages` for externally supplied
multi-event cases. Missing or ambiguous target matches remain in the denominator.
Fixed targets are excluded from inferred-age summaries. A failure in one
interval calculation does not erase the other methods' results.

JSONL records include numerical boundary indices and gradients, the right
variance score at zero with other parameters held fixed, input hashes and
calibration-profile details. A positive fixed-parameter score alone is not a
profiled KKT certificate. Summaries separate R/N, C/R and C/N, include binomial
intervals, unavailable reasons and full-domain return frequency. Validation
runs require `--study-role validation --protocol PATH`; development is the
default. A supplied protocol is hashed, not automatically judged scientifically
sufficient. Seed separation and frozen nuisance/selection rules remain the
study author's responsibility.

Before general sequence adoption, the approximation failures must be addressed
with a validated marginal-likelihood calculation, not by removing checks or
silently switching to a conditional MAP penalty. Subsequent work must retain
the original plan's independent rate-correlation, internal/multiple duplication,
calibration-boundary, alignment-length and model-misspecification cells,
nuisance sensitivity, grid-resolution checks and availability criteria.
HKY generating/JC69 fitting is now runnable; gamma/codon misspecification,
multiple-duplication scenario generation and topology-error calibration remain
unvalidated. This work does not claim to have completed those studies.

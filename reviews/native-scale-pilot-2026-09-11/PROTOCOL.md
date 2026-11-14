# Scale pilot protocol — frozen before results

Baseline: bb70467. Single-thread Python/BLAS, sequential fresh processes.
Seed 2026091107. Independent OU trait covariance, fixed root, alpha-height 1+j,
process tip variance 1, strong terminal shifts with alternating signed effects.
Shared truth uses up to five nonbackground regimes. Known errors vary 0.01–0.09;
estimated-error strata also simulate additional variance 0.04. Complete observations.
These fixtures stress computation and do not establish general recovery or false-positive rates.

First measure known-layout covariance refits for K=10/100, then searches at K=10/100.
Use balanced 1-trait no-error first; include known/estimated errors, a pectinate tree,
and four traits as diagnostics. Preserve all failed fits. Timeouts are censored results.
Full-search budgets: pool 128, refits 256, screening 100000, beam 2.
A K=100 cap does not prove K=100 was reached; record largest fitted K and exhausted budgets.
Profile separately from uninstrumented timings. Repeat comparable before/after runs
three times only if a bottleneck-driven implementation change is selected.

Workflow pilot: 32 tips, one true shift, pool 8, refits 8, screening 200,
shared regimes enabled, 19 calibration draws at level 0.05, two support draws.
Each support draw repeats calibrated selection. This tests execution and cost;
it cannot satisfy the statistical adoption gates in NATIVE_SHIFT_VALIDATION.md.

## Profile-driven follow-up, declared before variant timings

The baseline K=10 search profile spent 35.9/58.9 seconds in refits and
10.7 seconds in repeated node-group propagation (overlapping caller times).
The known-layout K=100 profile spent 1.36/2.51 fit seconds in QR.
Compare an exact integer-label interval implementation against parent propagation.
Separately evaluate SciPy economic QR against NumPy reduced QR; retain it only
if repeated same-input fits improve and independent numerical tests agree.
No optimizer bounds, starts, failure checks or search budgets may be relaxed.

The known-error fixed-layout run required 15,455 likelihood evaluations versus
146 without errors. Before timing variants, extend v2 to batch independent
numeric whitening rotations by level when average level width is at least eight;
retain scalar traversal for narrow/deep trees and the original determinant
summation order. This targets the measured whitening-construction cost without
changing variance optimization. v1 changes only integer regime propagation;
v2 combines that with economic QR and level-wise numeric construction.

## Final candidate decision

The SciPy economic QR prototype regressed in all three K=100 fixed-layout trials
and all three pectinate trials; stop that comparison rather than adopting it.
Retain NumPy QR. The final comparison repeats the baseline and the interval +
wide-level whitening candidate three times each, alternating order, for K=100
fixed-layout, pectinate K=10 fixed-layout, K=10 search, and known-error K=10
fixed-layout workloads. The 48 initial dense/interval tests passed on this candidate.
A subsequent 128-tip nested-workflow check will exercise the wide-level branch
inside calibration and support, in addition to the baseline 32-tip pilot.

## Rounding compatibility correction

The NumPy-hypot prototype failed the predeclared 1e-6 parameter comparison in the
estimated-error K=100 fit (alpha absolute difference 1.59e-6; alpha-height about
1.59e-5). Do not widen the tolerance. Retain Python math.hypot for both norms and
divisors inside vectorized level processing. A 32-case comparison found 4,892
rotation/scale entries differing in the NumPy-hypot prototype, but bitwise equality
for every checked rotation, scale, determinant and whitened value with the corrected
implementation. The estimated-error fit then matched the complete baseline result
exactly, including optimizer diagnostics.

With scalar-compatible hypot, choose an average level width of at least 32 before
batching; keep the original scalar path for smaller/narrower structures. Repeat the
final after trials on this exact source and reuse only the verified twelve baseline
trials. Thus the final source comparison is not fully interleaved. Baseline source,
harness, environment and thread checks are enforced by run_final_after.py. Preserve
all rejected prototype records separately. The independent wide-tree dense tests
now use 512 tips, including partial observation masks, to cover the final cutoff.

# Known-error calibration diagnostic

## Results

All 100 prespecified pairs (200 datasets) completed without fitting failures.

| True alpha-height | No-shift: plug-in / oracle | One-shift: plug-in / oracle |
| --- | ---: | ---: |
| 0 | 0/20 / 0/20 | 0/20 / 0/20 |
| 0.01 | 3/20 / 3/20 | 2/20 / 2/20 |
| 2 | 0/20 / 0/20 | 0/20 / 0/20 |
| 100 | 0/20 / 0/20 | 0/20 / 0/20 |
| infinity | 1/20 / 1/20 | 1/20 / 2/20 |

At the prespecified 5% level, no-shift rejection was **4/100 for both methods**
(Wilson 95% interval 1.6–9.8%). The unconditional one-shift-family test rejected
**3/100 with plug-in calibration** (1.0–8.5%) and **4/100 with the oracle**
(1.6–9.8%). All no-shift decisions agreed. In the one-shift test, one case was
oracle-only and none was plug-in-only.

Thus the previously observed 4/20 rate was **not reproduced** in this new study.
These results do not justify replacing the production calibration on the basis
of that small subgroup. They also do not prove uniform 5% control: there are only
20 trees per alpha cell, process variance is fixed at 1, errors follow one scale
pattern, and the later-stage test was unconditional with one specified effect size.

Matching decisions do not mean matching p-values. Probabilities differed in
90/100 pairs for each stage; median absolute differences were 0.015 and 0.020,
and the largest difference was 0.47 in each stage. For example, no-shift case 6
gave 0.525 versus 0.995. Both probabilities use finite Monte Carlo samples.
The oracle test has a valid rank-test construction under its known generator;
its probability is not an exact numerical tail area, and the true generator is
unavailable for real data.

The next priority is a prespecified grid-sensitivity and signal-to-noise study,
including zero/small process variance and a broader range of error scales.
The research envelope and tail-certified backend below remain experimental;
their conservative construction and synthetic tests alone are insufficient
grounds for changing production selection.

This is an independently seeded diagnostic of the known-error plug-in approximation,
not a replacement selection method. The frozen protocol specifies 100 random
ultrametric 16-tip trees, 20 for each true alpha-height value 0, 0.01, 2, 100 and
infinity. Process tip variance is 1 and tip observation variances range from
0.01 to 0.25. Each tree yields a no-shift dataset and a paired dataset with a
single clade's effective tip mean increased by 4.

For each dataset, the complete-search likelihood improvement is calibrated using
199 simulations from (a) fitted nuisance parameters and (b) the true generating
mean and covariance. Both use the same standard-normal bootstrap draws. The
second calibration is an **oracle diagnostic** unavailable for real data.

The no-shift test uses the no-shift family; the later-stage test uses the entire
family with at most one shift. Both are executed unconditionally. In particular,
these are not error rates conditional on passing an earlier selection test,
and they do not characterize convergence selection. Tree topology and shift
location are generated independently of the trait noise. Null and one-shift
results within a tree are paired, so they must not be pooled as independent rows.
Wilson intervals describe each cell's marginal rejection rate. With 20 trees per
cell, these are exploratory estimates, not evidence of exact 5% control.

All source snapshots, inputs, fitted null/alternative models and probabilities
are saved. The audit regenerates all inputs and refits all observed candidate
models; it replays both bootstrap generators for the first prespecified replicate
of each alpha cell. It does not independently replay the remaining probabilities.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
python tools/diagnose_shift_known_error.py \
  --output /tmp/shift-known-error --per-alpha 20 --replicates 199 --workers 4
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
python tools/verify_shift_known_error.py /tmp/shift-known-error --workers 4
```

See `protocol.json`, `records.jsonl.gz`, `summary.json` and `audit.json` for the
frozen design and completed results. Generation must finish before auditing.

## Research envelope prototype

`tools/shift_known_error_envelope.py` provides a separate no-shift prototype.
For every covariance on the search's alpha/process-variance grid it can repeat
the full model search on null draws and take the largest Monte Carlo p-value.
Duplicate zero-process-variance covariances are evaluated once. It checks the
fitted covariance first, then the remaining alpha values at that variance, then
the rest of the fixed grid. This ordering cannot justify early rejection.

An evaluated probability above the level establishes non-rejection for the full
grid; the unresolved maximum is bounded above by 1. A computation budget reached
while all evaluated probabilities are small yields `reject: null` (unresolved),
not a significance claim. Rejection requires evaluating every grid point.

For a true no-shift covariance belonging to this finite grid, the complete
maximum is at least the Monte Carlo probability under that true covariance.
This argument does not cover continuous nuisance values outside the grid, nor
the unknown location and effect of an existing shift in a later-stage test.
The prototype is not integrated into the CLI and does not change production
selection. Its output bounds are computational bounds, not confidence intervals.

`envelope-probe.json` explores two previously rejected datasets, deliberately
selected from the earlier stress study. These probes are not an independent
false-positive-rate evaluation. Their eight-evaluation budget may leave the
result unresolved.

The two exploratory probes (historical case IDs 93 and 105) both remained
unresolved after eight evaluations. Their lower bounds were 0.05 and 0.025,
respectively, with upper bound 1; each complete grid contains 1,324 distinct
covariances. This does not establish either significance or non-significance.
To replay the probes with the current prototype, run
`python examples/shift/known-error-diagnostic/envelope_probe.py --output /tmp/envelope-probe.json`
with a new output path. The exact prototype source used for the saved probes is preserved in
`research-source/shift_known_error_envelope.py` and matches the hash in their JSON.

The current prototype additionally records a failed bootstrap evaluation and
returns unresolved if the production profiler rejects a numerical boundary.
In particular, draws generated at the largest variance grid point can cause
fitted variance to reach the upper grid boundary. The production safeguard is
retained; a failed calculation is not converted to a p-value or a non-rejection.
A four-tip test reproduces this boundary case. Extending or certifying the
variance search consistently for observation and bootstrap fits is needed
before this becomes a practical full-grid method.

## Tail-certified research search

`tools/shift_variance_tail.py` now supplies `TailCertifiedSearch`, a research-only
backend for **strictly positive** known errors. It extends the same geometric
variance grid instead of accepting a maximum at a fixed artificial upper limit.
For contrast covariance `C(v) = v A + D`, with positive-semidefinite `A` and
positive-definite `D`, log-determinant is nondecreasing in `v`. At any candidate
mean, the log likelihood is at most

`-0.5 * (d * log(2*pi) + logdet(C(v)))`.

Once this bound at the grid endpoint is below every candidate's best likelihood
for every data column and every alpha, all larger variances are excluded. The
backend extends by the existing geometric ratio until this condition holds;
non-representable or uncertified tails still fail explicitly. The result is a
maximum on the geometric variance grid with its upper tail ruled out. This
does not certify interpolation between grid points or continuous alpha.

Tests compare the selected likelihoods with independent dense GLS, probe the
excluded tail, verify that expanding the cache for another batch preserves old
fits, and complete a full nuisance-grid calculation on a four-tip synthetic
example. This backend can be supplied to `known_error_envelope`; its class name
is recorded in the result. It is not used by the paired study above and is not
installed as the production CLI backend. Larger-sample calibration and grid
sensitivity must be checked before adoption.

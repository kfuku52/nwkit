# Final bounded calibration experiment (2026-09-10)

The user authorized an improvement attempt and a main-checkout commit, followed
by stopping this investigation if no acceptable improvement is demonstrated.
This protocol is written before the new development simulations.

The estimand remains the existing no-shift contrast test over the declared
finite mean/alpha search, including explicitly labelled scaled-effect drift and
independent-tip limits. It does not turn a shared mean into a test proving
convergence. This response will not silently remove limits or weaken error
thresholds to recover power.

## Candidate and mathematical scope

Compare the current full-alpha envelope with confidence-set restricted envelopes
of Berger and Boos (1994), with fixed beta=0.005 and 0.01. On the no-error null,
mean and process scale drop out. For each alpha in the existing 27-point grid,
construct a (1-beta) confidence set by inverting a Monte Carlo test of the
null-only profile likelihood ratio. Use `(1+exceedances)/(B+1)` with B=999.
For each retained alpha compute the original full-search shift-test probability;
return `min(1, beta + max(retained probabilities))`. An empty set has supremum
zero. Sharing simulation noise between these two tests does not invalidate the
union-bound argument; each pointwise probability must itself be valid.

For a true alpha belonging to the finite grid, rejection implies either exclusion
of the true alpha (probability at most beta) or a true-alpha shift probability
at most level-beta (probability at most level-beta). Hence level control follows
for that finite-grid model. This is NOT a continuous-alpha theorem, and does not
cover known-error variance or later-stage unknown mean parameters.

References: [Berger and Boos, JASA 1994](https://doi.org/10.1080/01621459.1994.10476836)
and [Dufour, Journal of Econometrics 2006, author manuscript](https://jeanmariedufour.github.io/Dufour_1995_MCT_W.pdf).
The OU-specific confidence statistic and implementation are derived here.

## Frozen development design and stopping rule

Use 100 independent replicate blocks on the balanced eight-tip tree from the
existing independent branch-innovation generator, true alpha*H=2.1,
sigma2*H=0.25, finite optimum effect=2 and no observation error. Each block
contains null, single, distinct and shared alternatives with common generating
innovations; a separate Monte Carlo seed is independent across blocks. Conditions
inside a block are paired and never counted as independent replicates.
Use master seed 20260930. The tested full search has max_shifts=2,
convergence=False, the unchanged 27-point alpha grid and level=0.05.
The shared family is contained in the unrestricted alternative, so first-stage
convergence-on/off equivalence is separately checked, not treated as more data.

Both beta values are fixed in advance. Continue to a fresh validation design
ONLY if a candidate improves detection by at least five percentage points in
the critical distinct-shift cell, loses at most five points in each other
alternative, has no calculation failures and has null rejection at most 7.5%
in development. These are screening criteria, not a claim of calibrated power
or Type-I error. If neither candidate passes, reject production adoption,
record uncertainty (including exact paired-gain bounds), commit the guards,
audits and negative result to the main checkout, and stop this investigation.
No post-result beta adjustment, effect-size changes or sample extension.

A timing-only pilot of two blocks at seed 20260928 may precede the development
run. It is excluded from all scientific acceptance counts. Its purpose is only
to verify execution and cost. Sources and design are frozen in each output
bundle before simulations. Keep every attempted block and failure.

If screening passes, select the candidate using distinct-shift net gain, with
smaller beta breaking ties, then freeze a separate new-seed validation protocol
before looking at its data. It must include cellwise null control (simultaneous
95% upper <=7.5%), power/noninferiority, tree/alpha/known-error sensitivity and
matched observation/bootstrap optimization. Positive development results alone
never replace the production method. Continuous-parameter or later-stage claims
require an additional valid construction and verification; a denser grid alone
is insufficient. A failure of these acceptance gates also ends this bounded
attempt without production adoption.

The current main checkout has advanced since the previous response. Its extra
mixed-zero-error boundary guards affect neither the no-error development data
nor the stored prior no-error evidence. Preserve these guards and all unrelated
changes during integration; integrate only the actual response delta.

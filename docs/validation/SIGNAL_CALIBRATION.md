# Pagel lambda calibration study

Use `tools/validate_signal_calibration.py` to measure null rejection,
alternative detection, fit availability, and optional profile-interval coverage
on independently generated Gaussian traits. This is a research study, not a
change to the default `nwkit signal` test.
The [200-dataset pilot](../../examples/signal/calibration-pilot-20260928/README.md)
is retained with individual records; it is not the confirmation run.
The known-SE pilot used only 19 inner draws because its nested rate and lambda
optimizations remain costly; its 0.05-step P-values are exploratory.
The [independent 20,000-dataset no-SE null confirmation](../../examples/signal/calibration-confirm-20260928/README.md)
found bootstrap rejection of 4.42–5.48% across the four conditions, with all
P-values available. Two of four prespecified Wilson-containment gates passed;
the other two overlapped a boundary. The conventional chi-square tail was
strongly conservative. This does not validate known-SE plug-in calibration
or the existing profile interval.
The [separate 5,000-dataset interior profile study](../../examples/signal/profile-ci-confirm-20260928/README.md)
returned 5,000/5,000 intervals covering true lambda 0.6, well outside the
prespecified 93–97% criterion. The profile interval is conservative under
that weak-signal eight-tip condition. Every one of those intervals was the
entire feasible range `[0, 1]`; the coverage result represents a lack of
resolution, not evidence of a precise interval procedure.

The generator constructs balanced and pectinate positive-definite tree
covariances from clade contributions, without calling NWKIT's covariance
builder. It draws from a direct Cholesky factor of
`0.7 * (diag(C) + lambda * (C-diag(C))) + diag(SE^2)` with root mean 1.2.
The cases are balanced and pectinate eight-tip nulls, an eight-tip null with a
fixed six-tip observed pattern, a balanced eight-tip null with heterogeneous
known SEs, a balanced 32-tip null, and an eight-tip alternative with lambda 0.6.
The missing pattern and SEs are held fixed across independent outer datasets.
The same simulated input is tested by the conventional chi-square tail and the
new parametric bootstrap. Each replicate gets its own seed derived from the
case name, master seed, replicate index, and stream, so case order and worker
count cannot change its input or bootstrap draws.

For each method, report rejections out of **all generated** datasets and out
of those with an available P value. A failed or unidentifiable fit is never
counted as a non-rejection. The threshold comparison is `p <= level`; with
`B=199` the minimum Monte Carlo P value is 0.005. Report Wilson Monte Carlo
intervals and fit status counts. Optional `--profile-ci` measures coverage of
the existing chi-square profile interval; it does not change that interval.
Use `--inner 0 --profile-ci` to measure profile coverage without bootstrap
draws; the summary then reports zero bootstrap P-value availability.
No conclusion about calibration should be drawn from a small smoke run or
from conditional rejection alone.

Before seeing confirmation outcomes, use distinct master seeds for exploratory
and confirmation runs. For a 5% null test, require 5,000 independent datasets
per case and a Wilson 95% Monte Carlo interval contained in 0.04–0.06 before
calling a supported case calibrated. Check alternative detection separately.
This engineering band is a study criterion, not a universal statistical
theorem. Known SEs make the bootstrap a nuisance plug-in procedure, so a good
no-SE result cannot validate that case. Trees, topology errors, estimated SEs,
model misspecification, and correlated traits outside the listed cases remain
untested.

For the existing profile interval, separately assess at least 5,000 new
datasets with true lambda 0.6. Call that condition calibrated only if all
intervals are available and the Wilson 95% Monte Carlo interval for coverage
lies within 0.93–0.97. Report boundary coverage at lambda 0 separately; it
cannot validate interior values. This criterion is set before the interior
confirmation run and does not change the production cutoff.

Write trial output to a fresh temporary directory:

```sh
PYTHONPATH=. OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  python tools/validate_signal_calibration.py \
  --output /tmp/nwkit-signal-study --outer 200 --inner 199 \
  --scenarios balanced-8-null,pectinate-8-null --workers 4 --seed 20260928
```

The output contains a protocol with exact case settings and source hashes,
one JSONL record per generated dataset, and a summary. `completed=true`
is written only after source hashes are rechecked and all records are saved.
Use a new directory for every run; do not overwrite previous evidence.

With **all SEs zero**, the likelihood ratio is invariant to
`x -> a + b*x` for `b != 0`: profiling the mean and diffusion rate changes
both the null and free-lambda log likelihood by the same `-n*log(abs(b))`.
Consequently the simulated null likelihood-ratio distribution does not depend
on the fitted mean or rate. Under the stated Gaussian model, the plus-one
Monte Carlo test has finite-simulation level control. With nonzero known SEs,
this invariance no longer holds and the fitted-null bootstrap needs empirical
calibration. These statements condition on the supplied tree, observed-tip
pattern, and correctness of the Gaussian model.

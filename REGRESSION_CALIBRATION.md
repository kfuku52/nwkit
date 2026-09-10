# Fixed-model regression calibration

`tools/validate_regression_calibration.py` evaluates RSC Gaussian contrast
models and unpenalized phylogenetic binomial, Poisson and NB2 GLMMs. It also
has raw-tip RSC cases that run reconciliation, contrast construction and
biological/technical replicate handling through `nwkit regress`.

This is a validation harness, not a new production inference mode. The
scientific null is a prespecified coefficient equal to zero. Neither testing
selected variables nor calibrating nested-CV prediction scores is in scope.
Selection remains a separate exploratory procedure without post-selection
p-values. Multi-family error control must be evaluated separately from the
calibration of individual tests.

## Running

Use a GeneGalleon container, with this NWKIT source mounted read-only and on
`PYTHONPATH`. A frozen source copy is recommended during concurrent work.
The tool records Python/package versions, thread settings, every Python source
hash, and a source archive, and refuses to finish if source files change.
Use the container image identity in `--runtime-label`; this is provenance, not
a pinned upstream source default.

```bash
python tools/validate_regression_calibration.py --list-cases
python tools/validate_regression_calibration.py \
  --output /results/new-wald-pilot --replicates 200 --workers 2 \
  --methods wald,oracle --runtime-label 'GeneGalleon image identity and source overlay'
python tools/validate_regression_calibration.py \
  --output /results/new-bootstrap-pilot --cases rsc-e2,rsc-e5,rsc-e20 \
  --replicates 200 --bootstrap-replicates 999 --workers 2 \
  --methods wald,parametric-bootstrap,oracle,null-bootstrap \
  --runtime-label 'GeneGalleon image identity and source overlay'
```

Output directories must be new. Use `--replicate-start` for independent shards;
seeds depend on the case name, master seed, replicate and stream, not worker
count, case ordering or method ordering. The same case/master seed/replicate
produces the same input for all methods and for before/after comparisons.
Use an independent master seed for a final confirmation study.

`--case-file` accepts a JSON array of `Case` fields; omitted fields use the
defaults in `tools/regression_calibration_design.py`. This supports additional
sample sizes, effect sizes, copy imbalance, predictor correlations, sampling
variances and missingness without changing analysis code. For example:

```json
[
  {"name": "rare-binary-n120", "engine": "glmm", "size": 120,
   "family": "binomial", "baseline": 0.01, "phylogenetic_variance": 0.5},
  {"name": "rsc-e10-correlated-error", "engine": "rsc", "size": 10,
   "copies": 5, "predictors": 2, "correlation": 0.9,
   "predictor_variance": 0.5, "sampling_variance": 1,
   "biological_replicates": 5}
]
```

The binomial `baseline` is the probability at zero linear predictor before
random effects, not the realized or marginal prevalence. Count predictors
are generated as integer copy numbers and enter the fit as `log1p(count)`.
The generator's `dispersion` field is NB **size** `r` (variance `mu+mu^2/r`);
NWKIT's fitted `response_dispersion` is **alpha** `1/r`.
MAR/MNAR missingness uses a logistic intercept setting, so the actual missing
fraction is recorded rather than forced. Raw-tip missingness currently applies
to biological expression observations; complete loss at a leaf is ineligible.

## Methods and interpretation

* `wald` and `parametric-bootstrap` call the existing production estimators.
  RSC Wald uses the existing event-based t degrees of freedom. The coefficient
  bootstrap uses centered coefficient samples and percentile intervals.
* `oracle` is an independently calculated Gaussian GLS result with the true
  covariance known. It is a generator/linear-algebra control, not an available
  practical estimator. Under outcome-dependent missingness its unconditional
  Gaussian model is also misspecified.
* `null-bootstrap` is a validation-only constrained-null comparison. It
  refits null and alternative variance parameters with ML. For RSC, covariance
  component definitions remain fixed between the two fits, and the statistic
  is an event-balanced **composite objective** difference, not a chi-squared
  likelihood ratio. The reference is currently limited to exact predictors.
  It produces a p-value, not a confidence interval.
* GLMM `profile-likelihood` and `likelihood-ratio` use NWKIT's existing paths.
  A finite likelihood or successful optimizer alone is not evidence of
  calibrated inference.

RSC physical-generation cases preserve shared species-event and paralog-lineage
effects. Cases ending in `-working` instead generate from the event-inflated
working covariance. The latter does not turn the composite objective into an
ordinary normalized Gaussian likelihood. Both are sensitivity experiments;
neither should be mislabeled as a generic phylogenetic validation.
Contrast-scale predictor errors are generated from a conditional Gaussian
distribution, shared across paralogs in the same event. Raw biological
predictor replicates are a separate case. Estimated expression SE cases
regenerate sample-variance uncertainty in each outer dataset; production
bootstrap refits condition on the supplied sampling covariance.

The raw-tip technical-replicate case deliberately duplicates each biological
observation exactly. It tests invariance to technical duplication, not the
adequacy of a technical-error model. Raw cases use a known, correctly reconciled
tree and Brownian generation; they do not validate uncertain dating,
reconciliation, tissue/batch comparability or annotation errors.

## Evidence and acceptance

Each record contains generated data, truth, input hash, seed, fit diagnostics
and elapsed time. The summary reports both all-generated and available-only
rejection rates; interval availability, conditional coverage and joint
delivery/coverage; error reasons; and bootstrap attempted/successful fits.
An invariant binary sample is never regenerated to conceal unavailability.

Production bootstrap's discarded refits are observed without changing its
calls. The reference null bootstrap uses a fixed number of attempts: any
failures make the point p-value unavailable and yield lower/upper p-value
bounds. Raw-pipeline bootstrap currently explicitly reports that its internal
attempt accounting is not instrumented; use the direct-engine experiment for
refit-failure comparisons. Neither missing tests nor missing intervals mean
non-significance or successful coverage.

The proposed confirmation criteria for correctly specified, supported cases
are Wilson 95% Monte Carlo intervals within 0.04–0.06 for a nominal 0.05 test
and 0.93–0.97 for 95% coverage, plus a numerical-failure upper bound at 1%.
These are prespecified engineering bands, not universal statistical constants.
Pilot studies with 200 datasets generally cannot establish these bands; 5,000
datasets give a Monte Carlo SE about 0.0031 for a rate of 0.05. Report all
diagnostics and unsupported conditions instead of selecting favorable cells.

The runtime estimate for a cell is outer replicates times the measured base
fit plus bootstrap replicates times measured refit cost, including failed
attempts. Reference null comparisons usually require two refits per draw.
The reported worker peak RSS is cumulative, not a per-fit memory measurement.
Set BLAS threads to one when using multiple workers. Check `protocol.json`
for completion; partial records are not a completed experiment.

`tools/regression_calibration_laplace.py` independently integrates selected
fitted GLMM likelihoods by importance sampling with randomized Sobol points.
It reconstructs the Laplace expression from the prior, conditional likelihood
and posterior Hessian, then compares increasing sample sizes and independent
scrambles. This measures integration error **at fitted parameters**; it does
not establish an accurate reference optimum, SE or frequentist coverage.

See the [recorded pilot evidence](examples/regression_calibration/README.md).
The validation approach follows [Morris et al.](https://doi.org/10.1002/sim.8086);
the need to distinguish Laplace approximation error from statistical calibration
is consistent with [Ogden](https://arxiv.org/abs/1808.06341).

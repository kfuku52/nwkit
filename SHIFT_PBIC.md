# pBIC coordinate correction

The mismatch identified by the joint-search reference is in **kfl1ou**, not
NWKIT's adapter. A local, unreleased correction has been implemented in the
kfl1ou checkout based on version 3.0.9. The installed release of 3.0.9 does not
contain it. NWKIT's production interface and default criterion are unchanged.

The correction restores score consistency; it does not establish statistical
calibration or make convergence-model weights validated posterior probabilities.
Previous pBIC scores, selected models and weights must be recomputed with the
corrected backend. Historical evidence bundles are retained unchanged.

## Cause and correction

For an unconstrained OU model, let `m_b` be a fitted mean displacement and
`delta_b` its corresponding optimum change. On an ultrametric tree,

```text
m_b = w_b * delta_b
w_b = 1 - exp(-alpha * parent_age_b)
```

The free-fit solver uses `m_b`, while the pBIC coefficient determinant is defined
using optimum changes and their OU design matrix. This is the convention in
[equation 5 of Khabbazian et al. (2016)](https://doi.org/10.1111/2041-210X.12534).
Transforming a covariance determinant therefore requires

```text
log det Cov(intercept, delta) = log det Cov(intercept, m) - 2 sum_b log(w_b)
```

The old free-fit score omitted this term. All-singleton convergence groups use
absolute optima; their conversion to an intercept and parent-to-child optimum
changes has determinant magnitude one. They consequently describe the same
coefficient-volume penalty as the corrected free fit.

Two additional inconsistencies were corrected:

- The convergent score used coefficient covariance evaluated at the initial,
  unconstrained alpha, even after refitting the constrained model at another
  alpha. It now uses covariance from the actual constrained fit.
- The free score always counted alpha as estimated. It now counts the covariance
  parameters actually estimated, excluding fixed alpha and including estimated
  observation error where applicable.

At exactly zero alpha, a shifted mean displacement does not identify a finite
OU optimum change. Such a shifted model receives infinite pBIC rather than
negative-infinite evidence. The no-shift Brownian case remains eligible. Small
positive alpha remains sensitive to bounds; the correction is not a remedy for
weakly identified optima or the limitations of the convergence approximation.

Changes apply to diagonal-trait pBIC. BIC, AICc, pBICess, and kfl1ou's separate
full trait-covariance pBIC extension are unchanged.

## Minimal reproduction

The eight-tip fixture uses two terminal shifts and fixed observations from the
held-out pilot. Its public-API script is
`tools/reproduce-pbic-coordinates.R` in the kfl1ou checkout. No backend score
formula is copied into NWKIT.

| Alpha | Representation | Released 3.0.9 | Corrected checkout |
|---|---|---:|---:|
| Estimated | Free | 32.77624199 | -36.09058668 |
| Estimated | Singleton groups | -36.09058668 | -36.09058682 |
| Fixed at 0.4 | Free | 33.17362872 | 26.65565545 |
| Fixed at 0.4 | Singleton groups | 26.65565545 | 26.65565545 |

The estimated-alpha difference is explained by `2 sum log(w_b) = -68.86682867`.
For fixed alpha, the coordinate correction is `-4.43853173` and removing the
extra alpha penalty contributes `-log(8) = -2.07944154`. The likelihoods are
unchanged. The corrected estimated-alpha scores differ by `1.42e-7` because
separate optimizations return slightly different alpha estimates; fixed-alpha
scores agree exactly in this fixture.

[Before](examples/shift/pbic-correction/before.csv) and
[after](examples/shift/pbic-correction/after.csv) CSVs retain full precision.
[Backend provenance](examples/shift/pbic-correction/backend-provenance.json)
records the base commit and modified-source hashes. Both installations report
3.0.9 because this work has not been committed or released; use the source
hashes and isolated library identity to distinguish them.

## Regression and replay evidence

The R regressions independently reconstruct the OU mean design by propagating
regime weights along branches. They compute dense phylogenetic covariance and
its information determinant, rather than calling the backend score helper.
Coverage includes fixed/random roots, estimated/fixed alpha, disconnected and
nested shifts, known observation variances, missing-trait mappings, shared
optima, group ordering and the exact Brownian boundary. The original backend
failed 18 assertions in the initial reproduction suite; the corrected expanded
suite passes all 39.

The package check passed all 1,097 assertions without skipped tests, errors,
warnings or notes. This check excluded vignette building and the PDF manual. A subsequent full
`--as-cran` check also passed the tests and vignette rebuilding, but failed the
PDF-manual step because `pdflatex` is unavailable and warned about missing
`qpdf`. Its environment notes also include compiler flags and unavailable HTML
validation/math-rendering tools. The full check is therefore not clean.
NWKIT's 179 shift-related tests also pass against the corrected library.

The [replay bundle](examples/shift/joint-validation-pbic-fixed/) uses the exact
same 12 datasets as the earlier joint pilot. It is a software regression replay,
not independent statistical replication. All 4,138 candidate attempts completed.
All 12 BIC summary records are unchanged. Both additional pBIC cases pass the
representation audit, with maximum information-penalty gaps below `2e-6` across
210 checks; the four previous candidate failures did not recur.

For the two pBIC cases, the fixed-root comparison is tied within `2e-6`; the
random-root joint fit improves the score by 0.350323. Neither method recovers
the true shared partition in either case. The subsequent
[independent simulation and alpha-bound sensitivity study](SHIFT_ALPHA.md)
compares accuracy on new data before considering model averaging.

## Reproduce with an isolated backend

From the corrected kfl1ou checkout, choose a new library directory:

```bash
mkdir -p /tmp/kfl1ou-pbic-library
R CMD INSTALL --library=/tmp/kfl1ou-pbic-library .
R_LIBS=/tmp/kfl1ou-pbic-library Rscript tools/reproduce-pbic-coordinates.R /tmp/pbic-after.csv
```

Then from the NWKIT checkout:

```bash
R_LIBS=/tmp/kfl1ou-pbic-library PYTHONPATH=. python tools/validate_shift_joint_pilot.py \
  --output /tmp/shift-joint-pbic-fixed --rscript Rscript --seed 20260910
PYTHONPATH=. python tools/summarize_shift_joint.py \
  --run /tmp/shift-joint-pbic-fixed --output /tmp/shift-joint-pbic-fixed-evidence
```

Compare backend source hashes with the recorded provenance. Installing an
unchanged released 3.0.9 will reproduce the original mismatch, not the correction.

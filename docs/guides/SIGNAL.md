# Phylogenetic signal

`nwkit signal` estimates Blomberg's K and Pagel's lambda for one or more
continuous traits. It reports signal in the supplied traits, not residual
signal from a regression. Categorical labels must not be encoded as numbers
and interpreted as continuous observations.

```sh
nwkit signal -i examples/signal/tree.nwk \
  --trait examples/signal/traits.tsv --columns clustered,mixed \
  --n-sim 999 --seed 1 -o signal.tsv

# Known sampling SEs, paired with traits in the same order:
nwkit signal -i examples/signal/tree.nwk \
  --trait examples/signal/traits.tsv --columns clustered,mixed \
  --standard-error-column clustered_se,mixed_se -o signal_se.tsv
```

Use `--method K|lambda|both` (default `both`) and `--test yes|no` (default
`yes`). Lambda profile intervals remain available with `--test no`.
`--ci-level` defaults to 0.95. Tree or trait input may use stdin, but not both;
`-o -` writes only the result TSV to stdout. File output is staged before
replacement and cannot overwrite either input (including file aliases).

## Inputs and missing observations

The trait TSV uses `leaf_name` and the shared missing-value/unmatched policies.
`--columns` selects distinct numeric trait columns. `--standard-error-column`
selects one distinct SE column per trait, in matching order. SEs are known,
independent sampling standard errors in trait units, not standard deviations
of individual observations. Observed values require finite nonnegative SEs.
Missing SEs are allowed only where the corresponding trait is missing.
No SE columns means exact observations; explicit zero SEs give the same fit.
There is no estimation of unknown species-level sampling variances here.

Each trait uses its own observed tips; missing values are omitted, never
imputed for the signal test. At least three tips and a nonconstant trait are
required. The analysis covariance is rooted at the retained tips' MRCA:
the shared stem above that MRCA is excluded, as when pruning the tree.
Root stem lengths are ignored. Non-ultrametric trees and rooted polytomies
are supported under `--input-rooted auto|yes|no`. Non-root branch lengths
must be explicit, finite and nonnegative. A singular retained BM covariance
is reported without adding jitter or resolving zero-length branches.

## Statistics and inference

Let `C` be the BM covariance (shared root-to-MRCA distances), `x` the observed
trait vector, `n` its length, and `a = (1' C^-1 x)/(1' C^-1 1)`.

```
K = [(x-a)'(x-a) / ((x-a)' C^-1 (x-a))]
    / [(trace(C) - n/(1' C^-1 1))/(n-1)]
```

K compares the observed pattern with the BM expectation; it is not a
correlation coefficient or a fraction of variance explained. For error-aware
K, first fit the diffusion rate by ML under `V = sigma2*C + diag(SE^2)`,
profiling the root mean. Apply the finite-sample correction `n/(n-1)` to
that rate and substitute the resulting `V` for `C` in the K formula. This is
the Ives/phytools-style error-aware statistic; it is labelled by
`measurement_error=yes`, and its `sigma2` column is the corrected rate.
An all-zero SE vector uses the exact-observation calculation.

The K test is upper-tailed tip-label randomization: shuffle value/SE pairs
together and refit the rate for each error-aware permutation. With `B`
permutations, `p=(1 + number(K_perm >= K_observed))/(B+1)`, including numerical
ties. Default `B=999`; the minimum attainable p-value is `1/(B+1)`. This tests
exchangeability across tip labels, not `K=1`. Seeded results are stable under
trait selection/order and input row/tree child order because tips are sorted
and each trait has a separate seed derived from its name.

For lambda, the observation covariance is

```
V(lambda) = sigma2 * [diag(C) + lambda*(C-diag(C))] + diag(SE^2)
```

The root mean and diffusion rate are profiled by **ML**, including at lambda=0;
this is not the default REML fit of `asr`. Lambda is constrained to **[0,1]**,
matching NWKIT's existing lambda model. Both endpoints and multiple bounded
search intervals are evaluated. The reported `sigma2` is the uncorrected ML
rate for lambda. Fits and likelihoods are returned in original trait and
branch-length units.

The lambda test compares the ML fit to lambda=0 and reports the conventional
`chi-square(1)` likelihood-ratio tail probability. This is the convention
used by phytools, not a finite-sample exact calibration: lambda=0 is a boundary,
and small samples or weakly identifiable variance components can invalidate
ordinary asymptotic approximations. No parametric-bootstrap calibration is
implemented by this command. A zero LR gives p=1.

The confidence limits are the connected profile-likelihood interval containing
the selected maximum, with cutoff `logL_max - chi2_quantile(ci_level,1)/2`.
The mean and rate are reoptimized at each lambda. Limits may terminate at 0 or
1; they are approximate likelihood intervals and do not include tree or
measurement-error-estimation uncertainty. They do not claim to enumerate
separate confidence regions in a multimodal likelihood.

`--p-adjust bh` (default) applies Benjamini–Hochberg across tested traits,
**separately for K and lambda**. `num_tests` records each family's size;
unestimable traits are excluded. `--p-adjust none` retains raw p-values.
The usual BH assumptions apply; arbitrary dependence among traits is not
covered by a universal FDR guarantee.

## Output

There is one row per trait/method, in requested trait order and K/lambda order.
Unavailable fields are empty, never substituted zeros.

| Columns | Meaning |
|---|---|
| `trait`, `method` | Input column and K/lambda |
| `num_taxa`, `num_missing_taxa` | Observed and omitted tree tips |
| `estimate`, `status`, `message` | Statistic, estimability/boundary status, numerical failure explanation |
| `measurement_error` | Whether SE columns were supplied |
| `sigma2`, `root_mean` | Error-aware K corrected rate, or lambda ML rate and mean |
| `log_likelihood`, `null_log_likelihood`, `likelihood_ratio` | Lambda ML fit and null comparison |
| `test_method`, `p_value`, `p_adjusted`, `p_adjust`, `num_tests` | Inference and multiple-testing family |
| `num_simulations`, `seed` | K randomization controls |
| `ci_level`, `ci_lower`, `ci_upper` | Lambda connected profile interval |
| `lambda_lower`, `lambda_upper` | Search limits, currently 0 and 1 |

Statuses are `ok`, `boundary`, `constant_trait`, `insufficient_taxa`,
`singular_covariance`, `unidentifiable_lambda`, and `fit_failed`.
Star trees have no shared off-diagonal covariance and cannot identify lambda.
An identifiable lambda is also unavailable when the fitted diffusion rate is
zero or the likelihood profile is numerically flat. Such rows have no
lambda estimate, interval, or p-value. Constant and insufficient traits retain
both requested rows. Invalid tables, lengths and options fail the command;
numerical fit failures retain a row with their reason. An unbounded singular
zero-rate likelihood is not reported as a successful tiny positive rate.

The diffusion-rate search uses log(rate), with an explicit zero-rate candidate.
Exact observations are separated before whitening the noisy observations; this
avoids numerical rank loss when SEs are heterogeneous. Unit conversions use
power-of-two scaling, avoiding overflow from squaring a large trait scale before
dividing by the time scale. Non-finite or unrepresentable results are reported
as failed fits instead of receiving p-values.

The implementation uses dense tip covariances: storage is quadratic in tip
count, and matrix factorizations are cubic. Error-aware permutations refit
sampling-error-adjusted rates and can be costly. There is no performance
claim of equivalence to linear-time signal algorithms such as `fast.SSC`.

## References and validation

- Blomberg, Garland & Ives (2003), *Testing for phylogenetic signal in comparative
  data: behavioral traits are more labile*, Evolution 57:717–745.
- Pagel (1999), *Inferring the historical patterns of biological evolution*,
  Nature 401:877–884.
- Ives, Midford & Garland (2007), *Within-species variation and measurement error
  in phylogenetic comparative biology*, Systematic Biology 56:252–270.
- [phytools phylosig documentation](https://search.r-project.org/CRAN/refmans/phytools/html/phylosig.html).

Tests compare K with phytools 2.3.0 and Gaussian likelihoods against an
independent direct multivariate-normal optimization. Numerical SE fits use
tighter optimization tolerances than the default R calculation. Phytools can
search lambda above 1 up to a tree-specific maximum, so unconstrained phytools
lambda estimates need not equal NWKIT's. Permutation p-values also need not
match R because the RNG and finite-simulation correction differ.

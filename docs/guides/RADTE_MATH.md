# Native RADTE likelihood and shared chronology

This note specifies the implemented estimator. It is not a claim of universal
accuracy, consistent model selection, or calibrated finite-sample intervals.
See [RADTE.md](RADTE.md) for inputs, diagnostics, and executable examples.

## Chronology

Let `a` contain one age for each species-tree node and one for each gene
duplication. A gene speciation node indexes its corresponding species age;
it does not introduce another parameter. The incidence matrix `D` therefore
gives all gene-branch durations as `d = D a`. Another incidence matrix `C`
contains gene ancestry, species ancestry, and the lower/upper species events
containing each duplication. The feasible domain is

\[
  l\leq a\leq h,\qquad Ca\geq\epsilon,
\]

where fixed calibrations have equal bounds. The implementation normalizes
time by the species-tree height and uses a positive numerical duration floor
of `1e-10` in normalized units. It detects cycles before solving and propagates
age bounds through the directed constraint graph. All species events remain
in this graph even if no gene speciation observes them. Ages supported only
by this graph may be nonidentifiable within their intervals.

## Gaussian log rates

For each gene edge, write its instantaneous constant branch rate as

\[
 r_e=\exp(\mu+u_e),\qquad u\sim N(0,\sigma^2 K_\rho).
\]

`rho=0` gives independent branch log rates. For `0<rho<1`, introduce a stationary
Gaussian root state and specify each child-edge state as `rho` times its parent
state plus independent noise of variance `sigma²(1-rho²)`. Marginalize the
unobserved root state. Thus correlations decay per edge along the genealogy;
this is not a time-scaled diffusion. Tests compare the sparse precision with
the inverse of the independently constructed dense covariance.

## Branch-only observations

Given exact positive substitution lengths `b`, define

\[
 y_e(a)=\log b_e-\log d_e(a),\qquad Q=K_\rho^{-1}.
\]

The likelihood of `b` after integrating the rates is lognormal. Its Jacobian
is `product(1/b_e)`, constant with respect to ages. With fixed `rho`, profiling
the mean gives

\[
 \widehat\mu(a)=\frac{\mathbf1^TQy}{\mathbf1^TQ\mathbf1},\qquad
 S(a)=(y-\widehat\mu\mathbf1)^TQ(y-\widehat\mu\mathbf1).
\]

If variance is estimated, `sigma² = S/n`. Minimizing `S/2` therefore finds the
same age optimum as the profiled likelihood when `S>0`. The gradient is

\[
 \nabla_a(S/2)=-D^T\operatorname{diag}(1/d)Q(y-\widehat\mu\mathbf1).
\]

At zero residual variance the profiled likelihood reaches its strict-clock
limit; regular curvature and likelihood-ratio arguments no longer apply.
With a user-fixed nonzero SD, the age optimum is unchanged, but the
likelihood-ratio statistic is `(S(a)-S(best))/sigma²`, rather than
`n log(S(a)/S(best))`. The code distinguishes these for profile intervals.

This input model treats both sides of the root as observed. A reversible
sequence model generally identifies only their sum, so branch-only inference
is explicitly conditional on the supplied root split.

## Sequence likelihood

Scaled Felsenstein pruning computes a reversible sequence likelihood using
compressed site patterns. Discrete gamma uses conditional-mean category rates,
normalized to mean one. The sequence model and fitted nuisance parameters are
held fixed during dating. The ML prefit optimizes **unrooted** branch lengths:
for the two root edges `p,q`, their single identifiable length is `b_p+b_q`.

Let `v` be the vector of log unrooted lengths. At the unclocked fit `v0`, the
quadratic approximation to negative log likelihood is

\[
 \widetilde f(v)=f_0+g^T(v-v_0)+\tfrac12(v-v_0)^TH(v-v_0).
\]

The full Hessian includes branch covariances. It must be positive definite and
numerically identifiable; the prefit must be interior. Define

\[
 V=H^{-1},\quad y=v_0-Vg,\quad f_*=f_0-\tfrac12g^TVg.
\]

Then the approximated likelihood is `exp(-f*)` times an unnormalized Gaussian
kernel with center `y` and covariance `V`.

## Integrating rates without an identifiable root split

Put `c=e_p-e_q`, `z=c^T u`, and `k=c^T K c`. Gaussian conditioning gives

\[
 z\sim N(0,\sigma^2k),\qquad
 u\mid z\sim N(qz,\sigma^2K_c),
\]

\[
 q=Kc/k,\qquad K_c=K-Kcc^TK/k.
\]

For every non-root edge `e`, its log length conditional on `z` is
`mu + log(d_e) + u_e`. For the combined root edge it is exactly

\[
 v_{root}=\mu+u_q+\log(d_p\exp z+d_q).
\]

This identity handles the sum of two lognormal branch lengths without
approximating that sum as another lognormal random variable. Let `M` select
each non-root edge and edge `q` for the combined root row. Conditional on `z`,
`v` is Gaussian with mean `m(a,mu,z)` defined above and covariance
`sigma² M K_c M^T`. Set

\[
 W=V+\sigma^2 M K_c M^T.
\]

Integrating the conditional Gaussian rates gives the likelihood

\[
 \widetilde L(a,\mu,\sigma)=
 e^{-f_*}\frac{|V|^{1/2}}{|W|^{1/2}}
 E_{z\sim N(0,\sigma^2k)}
 \left[\exp\{-\tfrac12(m-y)^TW^{-1}(m-y)\}\right].
\]

The remaining expectation is one-dimensional. The implementation uses
Gauss-Hermite quadrature with log-sum-exp arithmetic and compares successive
orders for both values and gradients. It includes derivatives of `W`, its
log determinant, and the scaled quadrature nodes when estimating `sigma`.
Constrained optimization estimates the ages, mean log rate, and variance by
marginal likelihood, not by maximizing over every rate realization.

The optimization coordinate for an estimated rate variance is now
`tau = sigma² >= 0`; the upper numerical protection remains `tau <= exp(4)`.
This includes the exact strict-clock boundary instead of the former
`log(sigma) >= -9` floor. For positive tau the variance score is obtained by
the chain rule. At zero, write the conditional mean derivatives with respect
to root contrast z as `b=m'(0)` and `h=m''(0)`, with `h_root=w(1-w)` and other
entries zero. With `v=V^-1(m(0)-y)`, `P=M K_c M^T`, and `k=Var(z)/tau`, the
right derivative of the negative log likelihood is

```text
d f / d tau | 0 = 1/2 tr(V^-1 P) - 1/2 v' P v
                 + k/2 [b' V^-1 b + h' v - (b' v)^2].
```

It follows by differentiating the Gaussian expectation at zero variance.
The root log-sum second derivative is retained. A nonnegative right score at
a zero-variance optimum satisfies the one-sided variance KKT condition; it
does not justify an inverse-Hessian interval for an estimated boundary variance.

Numerical checks include independent two-dimensional Gaussian integration on
a two-tip example, central finite differences for age/mean/SD gradients with
independent and correlated clocks, and root-split invariance. These establish
the implemented identities locally; they do not establish statistical coverage
or global optimization for arbitrary gene families.

## Exact conditional joint MAP and approximation checks

Exact sequence inference minimizes

\[
 -\log L_{seq}\{d(a)\exp(\mu+u)\}+\frac{u^TQu}{2\sigma^2},
\]

holding SD fixed. Without a supplied SD, the branch-only fit estimates it from
non-root branch observations, so arbitrary allocation of the root length does
not determine the variance estimate. This is an empirical-Bayes conditional
joint MAP in log-rate coordinates. It is not the marginal estimator above.

The default sequence route checks the quadratic kernel against exact sequence
likelihood at high-weight conditional rate means and displacements along its
largest conditional covariance direction. The score discrepancy is standardized
by `H^-1`; raw derivatives in short-branch coordinates are not comparable.
Failed checks trigger exact conditional inference only when fallback is allowed.
This finite set of checks cannot bound error everywhere in the integral.

Conditional Laplace intervals require interior ages and nuisance parameters,
positive information, and intervals within the hard calibration domain.
Profiles recheck approximation validity along the constrained solutions.
Bootstrap replicates preserve shared events and refit the specified estimator;
they report failures instead of silently changing marginal replicates into
conditional MAP estimates.

## Small-sample curvature adjustment

`--uncertainty laplace` retains its unadjusted Gaussian definition. It can
undercover when the number of independent branch-rate observations is small.
For exact conditional sequence MAP, the fitted rate SD is held fixed inside
the information calculation; its estimation uncertainty is absent. Marginal
inference includes an SD coordinate, but a local Gaussian approximation to
the joint ML fit does not remove small-sample variance bias or provide an
accurate reference distribution near a zero-variance boundary. Increasing
alignment length does not increase the number of independently realized rates.

`--uncertainty studentized` is a separate, approximate finite-sample correction.
Let `C` be the free-age block of the inverse **full** observed-information matrix;
nuisance rate/mean/SD directions are retained before inversion. Write

\[
  X=[\mathbf 1,\ J],\quad \nu=n-\operatorname{rank}(X),\quad
  \widetilde C=(n/\nu)C.
\]

Here `J` is the local Jacobian of the observed log durations with respect to
free ages. Tree-only inference uses its observed branches. Exact conditional
sequence inference uses the non-root branches from which its SD was estimated.
Marginal sequence inference uses identifiable unrooted lengths, combining the
root pair before taking the log. Columns with no influence do not consume a
degree of freedom. Fitted latent rates are random effects and are not subtracted
as fixed-effect parameters. A four-tip root-duplication example has `(n, nu)`
equal to `(4, 3)` for conditional sequence MAP and `(5, 3)` for marginal inference.
It does not have 2,000 rate observations merely because it has 2,000 sites.

The `n/nu` factor converts an ML residual-variance estimate to its residual-df
version in the local homoscedastic Gaussian regression. The critical value is
`t_(1-alpha/2, nu)`, not `z_(1-alpha/2)`. If the user supplies `--rate-sd`, there
is no variance correction and the normal critical value is used. The correction
does not alter the SD or point estimates saved by the dating fit.

Propagate **all** hard age and positive-duration constraints through the shared
event graph to obtain each age's feasible range `(L_i, U_i)`. Use

\[
 g_i(a)=\log\frac{a-L_i}{U_i-a},\qquad
 s_i=g_i'(\hat a_i)\sqrt{\widetilde C_{ii}},
\]

and transform the two limits `g_i(a_hat) +/- critical * s_i` back with the
inverse logit. These are marginal intervals: their Cartesian product is not
a simultaneous feasible chronology or a joint confidence region. The method
does not clip a Gaussian interval at a calibration, remove a constraint, or
replace a point estimate. It refuses active-bound, singular, strict-clock, and
zero-residual-df cases instead of turning these into narrow intervals.

The t/variance adjustment follows the familiar unknown-variance Gaussian
regression calculation ([NIST confidence intervals](https://www.itl.nist.gov/div898/handbook/prc/section1/prc14.htm)).
Its extension to this nonlinear model, an unobserved root-rate split, and
sequence measurement error is a **local approximation**, not an exact pivot,
REML implementation, or universally calibrated confidence procedure. In
particular, it neither estimates missing topology/reconciliation uncertainty
nor integrates substitution-model parameter uncertainty. Low rate variation
can still put the marginal variance estimate at zero and make
intervals unavailable. The validation report records interval availability,
coverage conditional on availability, and correct intervals returned across
all families separately.

## Exact log-duration contrasts

For branch-only input, suppose there is exactly one free age a, every affected
duration is `d_e=s_e(a-a0)` with positive s_e and common a0, and unaffected
durations are fixed. Divide observed branch lengths by s_e for affected edges
and by their fixed durations otherwise, then take logs. The resulting model is

```text
y = 1 mu + I_affected theta + error,   theta=log(a-a0),
error ~ N(0, sigma² K_rho).
```

For fixed positive-definite K, GLS gives beta_hat and residual sum of squares S.
The contrast for theta has a t distribution with `n-2` residual degrees of
freedom when using `S/(n-2)`, or a normal distribution with supplied positive
sigma. Exponentiating its endpoints and adding a0 is an exact monotone
transformation. Intersecting this confidence set with the known feasible age
domain does not change coverage of a true age in that domain. The unrestricted
GLS contrast is used even when the constrained point estimate is at a boundary.
An empty intersection is reported explicitly; zero residual variance is not
used to manufacture a narrow interval. Exactness is under the stated Gaussian
branch-length model and known K, not robustness to estimated lengths or rho.

`exact-log-duration` checks eligibility algebraically. Internal free ages that
both lengthen child branches and shorten parent branches do not generally meet
the condition. The independent dense-whitening reference lives in
`tools/radte_interval_reference.py`; the production path uses the normal
equations and does not call the reference.

## Research-only calibrated profile

`tools/validate_radte_intervals.py --methods calibrated-profile` develops a
constrained nuisance plug-in parametric-bootstrap LR test. It regenerates
Gaussian rates on the rooted genealogy and then native CTMC alignment sites;
site resampling is a different observation experiment. Sequence fits must use
marginal inference with the native alignment and retain exact approximation
checks. Informative IUPAC ambiguities require an explicit coarsening model and
are rejected by this generator; fully missing entries retain their mask.

Failed refits count as exceedances for a conservative Monte Carlo tail bound,
and more than 10% failures makes the candidate interval unavailable. A one-sided
99% binomial upper bound accounts for simulation uncertainty at a fixed
candidate age, but not nuisance plug-in error or simultaneous grid error.
All grid points are inspected, and the hull of accepted points is padded by
adjacent cells. A finite grid can miss components and is not exact inversion;
the output is labeled `experimental-calibrated-profile-grid-hull`. The exact
branch-only special case bypasses this approximation. This prototype has not
passed the criteria for a general-purpose 95% confidence method and is not a
public `--uncertainty` choice.

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

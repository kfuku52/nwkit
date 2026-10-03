# Bounded Conditional MUL/MSC Estimation

`mul-reconcile --score-model msc --msc-fit age|ne|joint` fits population
parameters separately for each second-parent candidate. It retains the
[conditional direct-parent model](MUL_MSC.md): one disomic allotetraploid event
is assumed, species topology/ages and H1 are supplied, and observed copies are
conditioned on. It is not WGD detection, full AlloppNET or DL+ILS.
Default `--msc-fit fixed` (also when omitted) preserves fixed-MSC outputs;
default `--score-model dl` preserves D+L. Fit-only options in either other
mode fail rather than being ignored.

## Inputs and Example

```bash
nwkit mul-reconcile -i genes.nwk --species-tree species.nwk \
  --species-regex '.*_([^_]+)$' --score-model msc --msc-fit joint \
  --h1 X --species-time-unit generations \
  --hybridization-age-bounds 0.1 1.9 --population-size-bounds 0.1 10 \
  -o fits.tsv --report genes.tsv --msc-profile-out profiles.tsv \
  --model-out fits.json
```

These bounds are illustrative for the tiny synthetic chronogram in the pilot,
not biological defaults. Supply justified bounds and actual generations for
real analyses. Years, substitutions/site, or RADTE ages are not automatically
converted to generations. Ne is diploid/subgenome effective population size;
branch duration is divided by `2*Ne`. All rooted/binary, ultrametric, sampling,
root-stem and species-mapping requirements of fixed MSC still apply.

| Fit | Required | Omit |
| --- | --- | --- |
| `age` | `--hybridization-age-bounds`; coalescent units, or generations with fixed `--effective-population-size` | fixed `--hybridization-age`, Ne bounds |
| `ne` | generations, `--population-size-bounds`, fixed `--hybridization-age` | fixed Ne, age bounds |
| `joint` | generations, both bounds | fixed Ne and fixed age |

Bounds must be finite, positive and strictly ordered. Ne is not fitted again
to lengths already expressed in coalescent units. No species node age, branch
Ne, rate heterogeneity, founder jump or ghost-parent age is estimated.
For age fitting each candidate receives the intersection of supplied bounds
with the open H1 and H2 stems. A candidate is not excluded just because it
cannot exist at another candidate's initial age. Empty intersections, baseline
and autopolyploid candidates remain explicit exclusions. Strict stem endpoints
are moved one floating-point step inward; source ages stay fixed.

## Computation and Diagnostics

The fit uses a full coarse grid in normalized age and log-Ne coordinates,
then bounded L-BFGS-B from the best grid starts. `--msc-grid-points` defaults
to 5 per coordinate (minimum 3), `--msc-fit-starts` to 3, and `--msc-maxiter`
to 200. All attempted starts, their convergence, and the best coarse-grid score
are recorded. At least one start must converge; a converged result worse than
the coarse grid by more than `1e-7` log-likelihood units fails. This does not
certify a global optimum. Increase grid/starts for sensitivity checks rather
than treating one optimizer result as exact.

`--msc-max-evaluations` defaults to 5000 unique parameter vectors per candidate,
including diagnostics and profiles. Exceeding it aborts; there is no grid-only,
fixed-parameter or family-filter fallback. Coalescent work/assignment limits
still apply to each distinct family calculation at each parameter vector.
Species-colored topology duplicates are grouped with their original family
multiplicity; every input family remains represented in outputs.

Statuses distinguish numerical estimates from evidence:

- `flat`: log likelihood varies by at most `1e-8` on the three-point-per-axis
  diagnostic grid; fitted values are withheld.
- `boundary`: an optimum is within `1e-5` of a normalized search boundary;
  fitted values are withheld, with the numerical solution retained separately.
- `locally-unidentified`: observed-pattern log-probability sensitivity has
  deficient rank; fitted values are withheld.
- `locally-distinguishable`: the above local sensitivity has full rank at
  both finite-difference steps. This is **not global identifiability**, precise
  parameter estimation, adequacy of the biological model or calibrated support.
- `excluded`: the candidate has no supported event/time interpretation.

The sensitivity matrix uses centered steps `1e-4` and `2e-4` in normalized
coordinates, weighted by square-root family multiplicity. Singular values
above `max(1e-9, largest*1e-7)` count toward rank; both step ranks are reported.
Near-boundary step lengths shrink to remain interior. These are numerical
diagnostics, not scientific confidence thresholds.

For example, with only species A and two X homoeologs in the synthetic sister
candidate, probabilities depend on `(2-age)/(2*Ne)`. Different age/Ne pairs
with that ratio cannot be distinguished. Likewise, scaling every species age,
attachment age and Ne by the same factor leaves probabilities unchanged;
absolute time information must come from outside the gene topology model.
If no polyploid copies are observed, attachment age can be completely
uninformative.

The direct-parent model equates donor attachment and hybridization age. A
ghost-parent divergence plus a distinct hybridization time is not identified
merely by optimizing this field: the inferred attachment must not be relabeled
as a general biological WGD age.

## Profiles and Output Contract

Every fit records a finite-grid nuisance-refitted profile for each fitted
parameter, including the fitted coordinate and search endpoints. At each knot
the other coordinate is refitted with the same grid/start checks. A profile
above the fitted maximum by more than `1e-6` fails. Profiles are likelihood
diagnostics, **not confidence intervals**, P-values, or Bayesian posteriors.
No likelihood-drop cutoff is converted into a claimed coverage interval.

The new primary schema is `nwkit-mul-msc-fit-v1`:
`mul.tree, h1.node, h2.node, log_likelihood, delta_log_likelihood,
hybridization_age, effective_population_size, dated.tree, status, reason`.
Unknown fitted points and excluded scores are blank, never fabricated zeros.
Fixed supplied parameters are distinguished from fitted parameters in JSON.
Likelihoods for bounded/unidentified fits remain available as conditional
candidate scores, but dated point trees are withheld for those rows.

JSON method is `conditional-direct-parent-MUL-MSC-bounded-fit-v1`. It records
settings, attempted versus reported parameter names, each fit's numerical
solution, status/rank, all optimizer attempts, exclusions and full profiles.
`--report`/`--check-out` add `family_weight` (1 for CLI families) to fixed MSC's
gene columns; weighted contributions reconstruct the total. `--msc-profile-out`
exports parameter, value, profile score/drop, nuisance solution and attempts.
Only the main table may use stdout. All six possible files form one staged
bundle, protected against each other and every declared input.

Candidates within `1e-7` of the maximum log likelihood are recorded as
numerical ties, not support sets. Requesting `--tree-out` fails without
replacing any output if the best parent is tied or its parameters are withheld.
A unique maximum is not proof that its parent is correct: the true donor may
be absent from the candidate list or the model may be misspecified.

## Evidence

See [fitting validation](../validation/MUL_MSC_FIT_VALIDATION.md) and the
[independent sequence/tree pilot](../../examples/mul-msc/README.md). Analytic
three-tip estimators and confounding, a full five-tip forest oracle, CLI
serialization/rollback and unchanged fixed/D+L cases provide distinct evidence.
The pilot separately records true and estimated gene trees, bounded versus
reported points, and deliberate model violations. Its tiny replicate count
cannot establish general biological accuracy, event-test error rates or
uncertainty coverage. Estimated-tree inputs are still treated as fixed trees
by the likelihood; their uncertainty is not integrated by this command.

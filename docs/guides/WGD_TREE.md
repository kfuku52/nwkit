# Fixed Gene-Topology DL/WGD Likelihood

`nwkit wgd-tree` evaluates a rooted binary species-colored gene topology under
the background and single-event fits produced by `wgd-count`. It does not fit
new parameters to that gene tree or use its branch lengths/sequence likelihood.

```sh
nwkit wgd-tree -i gene.nwk --input-rooted yes --species-tree species.nwk \
  --count-model count-model.json --tree-id family01 \
  -o origins.tsv --model-out likelihood.json
```

The count JSON must contain full species-tree identity. Clades, parentage,
branch lengths and tips must match; harmless child reordering is accepted.
Gene tips must be uniquely named and map to species-tree tips through the
shared species parser/regex/map options. Transfers, unresolved polytomies,
unmatched species, and non-doubling multiplication candidates are unsupported.
The gene and species rooting declarations are independent: `--input-rooted`
controls the gene tree, and `--species-tree-rooted` controls the species tree.
Both default to `auto`; explicitly unrooted trees fail unless the user explicitly
overrides that declaration. An override does not reroot a tree.
The default size limit is 128 observed gene tips. Unsupported inputs or
numerical failures stop the command, with no weaker fallback.

## Model

The dynamic program integrates linear duplication/loss and species splitting,
tip detection, gamma family-rate categories and a retained-copy WGD pulse.
Extinction/survival probabilities and the positive-geometric ancestral-count
prior match the native count model. A unit-length Yule root stem with rate
`log(root_mean)` defines ancestral topology probabilities; family-rate
categories do not rescale this root prior. Colored unordered child topologies
use their appropriate orientation multiplicity.

Node-specific markers propagate the WGD split contribution through this
likelihood. Mixture assignments are ratios of weighted likelihood sums, not
unweighted averages of per-category probabilities. Conditioning is applied
after averaging categories. `--ascertainment observed` requires at least one
observed gene; `root-clades` additionally requires genes in each immediate
species-root child clade. This is gene-tree selection, not arbitrary count
filter correction.

Branch flow is analytic. With extinction parameter `a`, positive-count
geometric parameter `b`, survival `s`, geometric success `d`, and downstream
survival `S`, the shared integrating factor is `H=s*d/(d+b*S)^2` and the
exposure is `z=b/(d+b*S)`. Writing a topology component as `L_g=H*F_g(z)` gives
`F'_g=c_g*F_left*F_right`, so each component is a finite positive polynomial.
Its coefficients and node markers are evaluated with log-space convolutions.
Hidden extinct lineages remain integrated by the branching-process PGF; this
is not a latent-copy truncation or an approximation to the time integral.

Double and platform-dependent extended arithmetic must agree within
`--tolerance` for log likelihood and `--origin-tolerance` for assignments.
These differences diagnose floating-point sensitivity, not a certified
roundoff bound; on platforms without extra precision this comparison is less
informative. The report records the arithmetic precision and differences.
Analytic and independent count-likelihood tests provide separate validation.

## Interpretation

`conditional_wgd_probability` is the probability of a latent node origin
**given the observed topology, supplied parameters, root prior, sampling and
one specified candidate doubling**. It is not a probability that WGD occurred,
a parameter/tree posterior, a Bayes factor, or a new calibrated event test.
The likelihood difference also uses separately supplied count-fit parameters;
it is not an independent significance test. All internal gene nodes are
reported, including nodes not classified as duplications by LCA reconciliation.

The default evaluates all supplied count candidates; `--event-ids` restricts
stable species-clade identifiers. Tables include gene-clade and branch IDs,
per-candidate likelihoods and assignment meaning. JSON records the background,
candidate results, numerical errors and limitations. Related files are staged
together; declared inputs cannot be overwritten. Only the main table supports
standard output. The gene tree, species tree, or count-model JSON may be the
single `-` standard-input source; two simultaneous standard-input owners fail.

Later family origins, ILS, HGT, allopolyploid donor trees, triplication
polytomy resolution, tree ensembles and joint multiple-event histories are not
modeled. This is an experimental independent implementation, not a claim of
Whale equivalence.

## Verification

`tests/test_wgd_tree_model.py` checks analytic Yule/pure-loss limits, distinct
colored topology sums against independent count likelihoods, mixture weights,
ascertainment, extreme loss and WGD-after-SSD node attribution.
`tests/test_wgd_tree_numerical_regressions.py` checks tiny branch identity,
rare extinction, underflow-scale ascertainment and exact 32/64/128-tip combs.
`tests/test_wgd_tree.py` checks CLI output, full-tree identity and rejected
inputs. `tests/test_wgd_tree_flow_reference.py` checks marked branch transport
against independent ODEs and pure-loss/geometric-root marker formulas.
These bounded implementation checks do not establish real-genome
false-positive rates or robustness to topology/parameter uncertainty.

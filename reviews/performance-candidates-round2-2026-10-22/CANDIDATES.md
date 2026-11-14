# Ten new hypotheses after the 0.43.55 integrations

These are hypotheses, not adoption decisions. They differ from the first
round's K-permutation cache, GLMM factor solves, cluster grouping, grouped
projection, diagonal test, RHS batching, ASR log cache, gamma reuse, DTT batching
and observed-only multivariate covariance proposals.

| ID | Target | Proposed improvement | Essential boundary |
| --- | --- | --- | --- |
| 01 | `signal_stats.lambda_fit` | One diagonally scaled eigendecomposition for the zero-error lambda profile | Preserve the optimizer, intervals, null fit and error-aware refits |
| 02 | `GaussianTreeProcess.covariance` | Fill covariance by descendant/LCA blocks instead of querying LCA for every pair | Arbitrary requested nodes, duplicates, signed/zero/extreme slopes and order |
| 03 | `CountModel._prune_log` | Tile logsumexp across rare families rather than call once per family | Preserve unbounded process, all-missing subtrees, tails and underflow recovery |
| 04 | `_log_matrix_product` in WGD | Memory-bounded batches for log-semiring pulse compositions | Impossible states, extreme negative logs, same truncation and event semantics |
| 05 | Threshold Gibbs sampling | Prepare invariant conditional precision/transition data once | Preserve node permutation, draw order, thresholds and posterior diagnostics |
| 06 | Vector Gaussian `_propagate` | Validate branch SPD covariance without computing an unused inverse and logdet | Preserve validation/failure behavior and all partial/exact observations |
| 07 | `phylogenetic_pca.fit_pca` | Reuse final factor/whitening for reported profile likelihood | Preserve GLS centering, scale/rank/SVD checks, signs and full output |
| 08 | Stochastic map preparation | Reuse uniformization parameter work between guards and context construction | Run work/memory guards before allocation; preserve rare endpoint bounds and RNG |
| 09 | ER discrete pruning | Exploit diagonal-plus-constant transitions for O(k) log messages | Preserve zeros, ambiguous/missing tips, branch regimes and posterior fields |
| 10 | Full-OU moment identifiability | Reuse selected rows/products instead of full covariance products per design entry | Preserve moment design, Jacobian/rank tolerance, parameterization and diagnostics |

Each numbered directory will contain the independently owned prototype,
primary references, comparable timings and numerical/test evidence.

## Primary-source starting points

These sources motivate algebra/model boundaries. The proposed NWKIT-specific
speedups are inferences from source inspection and still require measurement.

- [Pagel (1999)](https://www.nature.com/articles/44766) and
  [SciPy `eigh`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.eigh.html):
  lambda-model context and the symmetric spectral solver, candidate 01.
- [Ho & Ané (2014)](https://academic.oup.com/sysbio/article/63/3/397/1649891):
  tree-structured Gaussian inference background, candidate 02. This proposal
  targets the dense covariance consumer rather than replacing existing pruning.
- [Csűrös (2010), Count](https://academic.oup.com/bioinformatics/article/26/15/1910/189891),
  [SciPy `logsumexp`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.logsumexp.html),
  and [Blanchard, Higham & Higham](https://arxiv.org/abs/1909.03469):
  gene-family likelihood context and floating-point log-domain stability,
  candidates 03–04. Algebraic equivalence alone does not prove numerical equality.
- [Felsenstein (2012)](https://pubmed.ncbi.nlm.nih.gov/22218305/):
  liability/threshold-model context, candidate 05. Retain the current Gibbs
  transition kernel and random schedule rather than substituting a sampler.
- [SciPy `cho_factor`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.cho_factor.html):
  SPD factorization and validation, candidate 06. The unused inverse is a
  separate operation from checking that factorization succeeds.
- [Revell (2009)](https://onlinelibrary.wiley.com/doi/abs/10.1111/j.1558-5646.2009.00804.x):
  phylogenetic centering and PCA context, candidate 07.
- [Hobolth & Stone (2009)](https://arxiv.org/abs/0910.1683):
  endpoint-conditioned CTMC simulation and uniformization, candidate 08.
  Cache only repeated preparation; retain the conditional tail bound.
- [Felsenstein (1981)](https://pubmed.ncbi.nlm.nih.gov/7288891/):
  likelihood pruning, candidate 09. The ER diagonal-plus-constant simplification
  is an algebraic inference for the already supported model.
- [Bartoszek et al. (2012)](https://pubmed.ncbi.nlm.nih.gov/22940235/)
  and [NumPy `einsum`](https://numpy.org/doc/stable/reference/generated/numpy.einsum.html):
  multivariate OU context and selected contractions, candidate 10. The current
  stable full-OU parameterization and rank diagnostic remain authoritative.

Current web manuals can describe newer versions than the installed numerical
libraries. Actual experiments use the environment recorded in PROTOCOL.md.

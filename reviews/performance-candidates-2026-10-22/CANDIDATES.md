# Candidate selection and primary-source rationale

The source inspection found existing bitmask clades, binary-lifting LCA,
low-rank/Woodbury covariance solvers, sparse tree precision, bounded fit caches,
and symmetric-generator spectral transition reuse. Those broad suggestions
would duplicate existing implementation. The ten candidates below target
specific remaining work. Each candidate is assigned to a distinct subagent;
at most three agents are active concurrently and numerical runs use one lock.

| ID | Target | Hypothesis | Correctness boundary |
| --- | --- | --- | --- |
| 01 | `signal_stats.permutation_test` | Reuse fixed Cholesky and GLS constants across K permutations | Same permutations and tail counts; positive sampling errors need their existing refits |
| 02 | Dense GLMM mode and coefficient covariance | Use `cho_solve` for already computed Cholesky factors | Retain Newton/backtracking/optimizer and positive-definiteness checks |
| 03 | `regress._fit_model` sandwich meat | Aggregate event scores once and multiply the score matrix | Same event weights, CR1 correction and inferential outputs |
| 04 | GLMM grouped predictor uncertainty | Contract covariance and slopes in bulk | Same observation/trait index order and cross-observation blocks |
| 05 | `gaussian.is_diagonal` | Replace row-by-row reductions with bulk exact-zero counting | Preserve exact zeros, NaN/Inf behavior and early-exit tradeoffs |
| 06 | `regress._profile_covariance_fit` | Solve design and response together for each factor | Preserve residual quadratic stability and ML/REML optimization |
| 07 | Discrete ASR inside likelihood | Cache log probabilities of shared transition arrays within one call | Preserve zero/impossible states, regime-specific rates and mutation lifetime |
| 08 | WGD log transition kernel | Hoist rate-independent log-factorial/combinatorial terms | Same unbounded branching process and unnormalized omitted tail |
| 09 | DTT Brownian simulations | Batch tree propagation across nodes/replicates | Same seed streams, stable subtree moments and bounded memory |
| 10 | Partial-response multivariate PGLS | Construct observed covariance directly instead of full Kronecker then slice | Preserve Fortran response ordering, all cross-trait terms and returned fit schema |

The hypotheses are inferences from source inspection and the following
primary literature/official APIs. These sources do not establish speedups in
NWKIT: only each candidate's comparable local measurements do that.

- [Ho & Ané (2014), linear-time Gaussian and non-Gaussian trait models](https://academic.oup.com/sysbio/article/63/3/397/1649891): structured likelihood and repeated quadratic forms; broad tree-structured optimizations already exist here.
- [Bastide et al. (2020), efficient general Gaussian phylogenetic inference](https://arxiv.org/abs/2003.10336): tree structure for likelihood/gradient computation; background on preserving model structure.
- [Blomberg, Garland & Ives (2003), phylogenetic signal](https://onlinelibrary.wiley.com/doi/abs/10.1111/j.0014-3820.2003.tb00285.x): K and tip-label randomization, candidate 01.
- [Cameron & Miller (2015), cluster-robust inference](https://cameron.econ.ucdavis.edu/research/Cameron_Miller_JHR_2015.pdf): cluster scores and sandwich covariance, candidate 03.
- [Felsenstein (1981), likelihood pruning](https://pubmed.ncbi.nlm.nih.gov/7288891/): discrete likelihood recursion, candidate 07.
- [Harmon et al. (2003), disparity through time](https://pubmed.ncbi.nlm.nih.gov/12920297/): DTT and Brownian simulation background, candidate 09. NWKIT's documented squared-distance statistic remains authoritative.
- [Csűrös (2010), Count](https://academic.oup.com/bioinformatics/article/26/15/1910/189891): gene-family birth/death likelihood context, candidate 08. No replacement with a killed finite CTMC is permitted.
- [SciPy `cho_solve`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.cho_solve.html) and [`solve_triangular`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.solve_triangular.html): factor reuse and matrix right-hand sides, candidates 01, 02 and 06. Finite checks must not be silently weakened.
- [NumPy `bincount`](https://numpy.org/doc/stable/reference/generated/numpy.bincount.html): weighted grouped sums, candidate 03.
- [NumPy `einsum`](https://numpy.org/doc/stable/reference/generated/numpy.einsum.html) and [`einsum_path`](https://numpy.org/doc/stable/reference/generated/numpy.einsum_path.html): contractions and memory/performance tradeoffs, candidate 04.
- [NumPy `count_nonzero`](https://numpy.org/doc/2.0/reference/generated/numpy.count_nonzero.html): exact nonzero counts, candidate 05.
- [SciPy `logsumexp`](https://cache-docs.scipy.org/doc/scipy/reference/generated/scipy.special.logsumexp.html): stable log-domain summation, candidates 07 and 08.
- [SciPy `splu`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.linalg.splu.html): sparse ordering options. NWKIT already uses symmetric minimum-degree ordering with explicit pivot validation, so a generic ordering suggestion was excluded.
- [Clavel, Escarguel & Merceron (2015), mvMORPH](https://besjournals.onlinelibrary.wiley.com/doi/10.1111/2041-210X.12420) and [NumPy `kron`](https://numpy.org/doc/stable/reference/generated/numpy.kron.html): multivariate covariance structure, measurement error and missing data, candidate 10.
- [NumPy `SeedSequence`](https://numpy.org/doc/stable/reference/random/bit_generators/generated/numpy.random.SeedSequence.html): reproducible spawned simulation streams, candidate 09.
- [pyperf measurement protocol](https://pyperf.readthedocs.io/en/latest/run_benchmark.html) and [result analysis](https://pyperf.readthedocs.io/en/latest/analyze.html): warmups, sequential runs, calibrated loops, raw distributions and host-noise limitations. Our scripts use paired timing directly; they do not claim to have run pyperf's full multi-process protocol.

Current web manuals can describe newer releases than this experiment's
SciPy 1.17.1 and NumPy 2.5.4. The scripts and numerical checks use the installed
versions; references explain the API/algorithm rather than certify binaries.

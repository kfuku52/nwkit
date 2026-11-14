# Candidate 05: prepare invariant threshold-Gibbs conditionals

**Verdict: conditional; defer adoption of this prototype.** Its sampled paths,
all reported fields and final chain RNG states were exactly equal, and the
existing independent-reference tests passed. However, complete fits with all
tips observed showed no useful improvement. A 128-tip/two-observation control
gave only **1.05x**, and the other apparent gains were noisy. The ordinary
workload is dominated by truncated draws and posterior diagnostics rather than
the conditional constants targeted here. Production code was not changed.

## Hypothesis and preserved model

`_conditional_parameters` repeatedly looks up fixed branch transitions,
accumulates the same conditional precision and computes the same conditional
variance/SD at every node update. The prototype prepares these invariants and
ordered parent/child tuples once per complete fit. It retains the original
numerator's multiplication, subtraction, addition and division order. It uses
the same node permutations, normal/truncated-normal calls, ordinal threshold
updates, initial values and diagnostics. Prepared data has per-fit lifetime;
there is no global mutable tree cache or simultaneous Gibbs update.

The background is [Felsenstein's threshold model (2012)](https://www.journals.uchicago.edu/doi/abs/10.1086/663681):
discrete observations constrain latent continuous liabilities, which are
sampled conditional on the phylogeny. The repository's
[ASR guide](../../../docs/guides/ASR.md) and
[diagnostic contract](../../../docs/guides/THRESHOLD_DIAGNOSTICS.md) remain
authoritative for its identified N(0,1) root/unit diffusion and ordinal prior.
[SciPy `truncnorm`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.truncnorm.html)
documents the standardized interval parameters retained by the sampler;
[NumPy SeedSequence](https://numpy.org/doc/stable/reference/random/bit_generators/generated/numpy.random.SeedSequence.html)
documents the independently spawned streams. These sources motivate correctness
boundaries, not a speedup claim. The precomputation is an inference from the
current implementation's fixed process, not a different statistical method.

## Complete-operation results

Each call includes per-fit preparation, initialization, all chains and sweeps,
retained trace storage, liability/category/threshold summaries, rank/ESS/MCSE
diagnostics and the returned fit. Fixture generation, imports and the
comparison instrumentation are excluded from timings. Both implementations
have the same settings and seeds. These short runs test computation, not
posterior convergence; their diagnostic warnings are retained and compared.

| Workload | Median baseline (s) | Median prototype (s) | Ratio | Evidence |
| --- | ---: | ---: | ---: | --- |
| Binary 4 tips, all observed | 0.06685 | 0.06813 | 0.98x | 3/7 pairs faster; ranges overlap |
| Binary 32 tips, all observed | 0.55103 | 0.54399 | 1.01x | 4/7 pairs faster; ranges overlap |
| Binary 128 tips, only 2 observed | 0.63816 | 0.60722 | 1.05x | 7/7 faster; ranges narrowly separated |
| Ordinal 24 tips, estimated thresholds, ambiguous/missing tips | 0.62345 | 0.54290 | 1.15x | 4/7 faster; broad overlapping ranges |
| Same ordinal fixture, fixed thresholds | 0.47656 | 0.45228 | 1.05x | 5/7 faster; ranges overlap |

The sparse case used 32 retained draws, 64 burn-in sweeps, thinning 8 and two
chains. Its raw ranges were 0.62274–0.65054 s and 0.59039–0.62114 s. The
separation is only about 1.6 ms, and other user CPU jobs were active. This is
weak evidence for a small benefit in that restricted setting, not a reliable
general THRESHOLD speedup. Raw seven alternating pairs and exact settings are
in [results.json](results.json).

An amortized sweep of just the conditional-parameter helper was 1.30x/1.23x/
0.99x at 4/32/128 tips. This excludes RNG, draws, initialization and diagnostics,
and does not establish a consistent scaling benefit. Preparation itself was
about 1.8/16.5/148 microseconds at those sizes; the JSON retains an identical-
function timing control for that preparation. It is included in every complete
fit comparison above.

## Memory and remaining bottleneck

Whole-fit peak traced allocations were about 0.080→0.076 MB at 4 observed tips,
0.247→0.261 MB at 32 observed tips, 0.879→0.954 MB in the sparse 128-tip case,
and 0.218–0.223→0.232–0.236 MB for ordinal cases. The small apparent decrease
is not claimed as a memory improvement. Constructing the prepared data alone
peaked at 2.4/20.1/80.8 kB for 4/32/128 tips. Storage grows with the tree's
nodes/edges and is released after the fit. Tracemalloc is not total RSS and
does not observe every native allocation.

The explanatory [single-call profile](profile.txt), outside timed comparisons,
used the 32-tip all-observed workload. In 0.993 profiled seconds,
`_conditional_parameters` accounted for 0.0073 cumulative seconds over 6,048
updates, while `_truncated_draw` accounted for 0.884 s and diagnostics for
0.088 s. Profiling changes absolute costs, but this supports the observed
small complete-fit effect. Optimizing the constants cannot remove the larger
truncated-normal/log-domain work. No sampler, tail-stability check or diagnostic
coverage was replaced to obtain a larger result.

## Correctness and tests

- All five complete-operation comparisons had exactly equal retained traces,
  posterior arrays, every `ThresholdFit` field (including the complete
  diagnostics DataFrame), and final RNG state of every chain.
- **952 direct helper comparisons** matched exactly on balanced, pectinate and
  star trees, with both the supported Brownian process and general signed/
  zero-slope affine transitions used as helper stress controls.
- Eight failure controls matched exception type and message: zero/negative/
  non-finite branches, threshold count/order, too few observations, missing
  ordinal categories and an empty allowed-state set.
- Changing a branch between fits recomputed the metadata correctly; the next
  fit matched the original, and repeated calls also matched exactly.
- Prototype-injected existing tests: **21 passed in 40.55 s**, no skips.
  They cover binary/ordinal/ambiguous sampling, CLI exports, disk/memory trace
  equality, independently integrated binary and ordinal posteriors and fixed
  ordinal references. See [tests.log](tests.log). Existing scientific
  tolerances were unchanged.
- Ruff lint and formatting checks passed for all three experimental Python
  files. The full release suite, other platforms, quiet-host replication,
  default-size long MCMC runs and replicated scientific calibration were not
  run. This proposal changes deterministic bookkeeping only, not the target
  posterior or uncertainty method.

Integration is deferred pending a convincing complete-workload gain on larger
representative datasets under quieter conditions, with enough benefit to
justify the extra per-fit metadata and code. This is not a claim that the
proposal is numerically invalid; the present performance evidence is weak.

## Reproduction

Experiments used the immutable `86c8f2c537769d852c17e2c789340f6c9de4ef54`
archive, NWKIT 0.43.55, isolated from unrelated live-worktree edits. Results
record `nwkit.__file__` and the baseline target-module SHA-256. The environment
was Python 3.12.14, NumPy 2.5.4, SciPy 1.17.1, pandas 3.0.6, macOS ARM64.
Current documentation can describe newer versions than those measured.
All substantial runs held the common exclusive lock and used one numerical
thread. The source scripts set these variables before numerical imports.

```sh
cd /tmp/nwkit-perf-round2-20261022-rkzoGe/baseline
PY=/tmp/nwkit-perf-20261022-9fYwTD/venv/bin/python
DIR=/Users/kf/repos/nwkit/reviews/performance-candidates-round2-2026-10-22/05-threshold-conditionals
$PY "$DIR/benchmark.py"
$PY "$DIR/benchmark.py" --tests > "$DIR/tests.log" 2>&1
$PY "$DIR/additional_checks.py"
```

The test subprocess invokes the repository entrypoint exactly as
`python tools/check.py test -- -p candidate05_plugin tests/test_threshold_asr.py tests/test_threshold_reference.py -rs`.
The prototype wrapper refreshes runtime module globals so existing monkeypatches
of trace storage remain effective; production integration would not require
that experimental wrapper.

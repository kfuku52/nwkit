# Candidate 03: bounded batching of rare-family log pruning

**Verdict: adoptable within the tested many-rare-family scope.** The final
1 MiB workspace-estimate prototype accelerates the complete fixed-parameter
family-likelihood operation by **1.24x at 128 rare families** and **1.31x at
512**, both at maximum count 86. Every one of seven paired samples improved,
the raw ranges were separated, and all compared outputs were exactly equal.
This does not establish a speedup for routine datasets, fitted searches or
bootstraps. Production code was not changed.

## What is new

`CountLikelihood._prune_log` currently calls `scipy.special.logsumexp` once per
family for every child edge. The prototype batches only the independent family
axis. The latent-state sum remains the contiguous last-axis reduction, and
the branch and root accumulation order is retained. It does not modify
transitions, truncation bounds, omitted tails, detection, root priors, family
mixtures or ascertainment. This is separate from the previous round's
log-gamma reuse and this round's pulse-matrix-composition candidate.

The inference is supported by SciPy's documented array-axis reduction and
NumPy broadcasting, with explicit tiling because broadcasting can increase
working storage. See [SciPy logsumexp](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.logsumexp.html)
and [NumPy broadcasting](https://numpy.org/doc/stable/user/basics.broadcasting.html).
The scientific context is independent-family birth/death pruning and retained
genome-multiplication likelihoods: [Csűrös (2010)](https://academic.oup.com/bioinformatics/article/26/15/1910/189891)
and [Rabier, Ta & Ané (2014)](https://pubmed.ncbi.nlm.nih.gov/24361993/).
These references motivate the implementation hypothesis; they do not establish
the NWKIT timings. The [repository model guide](../../../docs/guides/WGD_COUNT.md)
remains authoritative for NWKIT's narrower assumptions.

## Whole-operation measurements

Final results are in [results.json](results.json). Times below are median
seconds per `family_log_likelihoods` call on a three-tip tree. Inputs are
constructed once, but each timed call constructs emissions and probability
transitions, runs scaled pruning, identifies rare rows, constructs log
transitions when necessary, recovers the rare likelihoods and applies
ascertainment. Model construction, imports and CLI parsing are excluded.

| Workload | Baseline | Prototype | Ratio | Interpretation |
| --- | ---: | ---: | ---: | --- |
| Normal 129-family control, bound 20, no log fallback | 0.000681 | 0.000673 | 1.01x | Same executed path; no improvement claim |
| 1 rare family, bound 86 | 0.03396 | 0.03422 | 0.99x | Retains original one-row reduction; no measurable change |
| 16 rare families, bound 86 | 0.03951 | 0.03762 | 1.05x | Large overlapping ranges; inconclusive |
| 128 rare families, bound 86 | 0.08857 | 0.07170 | **1.24x** | Seven improvements; ranges 0.08185–0.10333 vs 0.06650–0.07544 |
| 512 rare families, bound 86 | 0.23548 | 0.18007 | **1.31x** | Seven improvements; ranges 0.22736–0.29453 vs 0.17012–0.19355 |
| 128 rare families, bound 160 | 0.25705 | 0.23955 | 1.07x | Batch=1 retains original reduction; overlapping ranges, no gain claim |
| 64 rare families with pulse and detection, bound 86 | 0.09791 | 0.09198 | 1.06x | Transition composition dominates; modest overlapping ranges |

Rare rows in the timed rare cases have 69–84 observed copies (or 135–150 in
the larger-state case) with tiny positive duplication rates. These intentionally
exercise the existing underflow-recovery path and resemble the independent
high-copy regression fixtures. No claim is made that most empirical families
occupy this regime. Every rare case recorded the number of actual fallback
rows; the normal control recorded zero.

The transition-prepared complete `_prune_log` operation, including its
binomial emissions and all tree messages, was 3.84x/6.43x faster at 16/128
families with bound 20, and 1.53x/1.46x at 128/512 with bound 86. Those gains
are kernel-level evidence, not complete likelihood or CLI gains.

## Memory and tile selection

The final prototype estimates four float-sized temporaries for each
family-by-parent-state-by-child-state tile and uses a **1 MiB estimate**, with
at most 64 families per tile. A one-family call or a state dimension too large
for a multi-family tile uses the original list reduction. This estimate is
not a hard cap on total allocations or RSS: logsumexp has its own temporaries,
and existing input, tree-message and output storage also contributes.

Final whole-operation traced peaks were 1.69→2.72 MB at 128 families and
5.226→5.225 MB at 512 families/bound 86. The latter operation's unchanged
probability-path allocations dominate its overall peak. Pulse/detection was
1.33→2.38 MB. These are tracemalloc peaks and exclude some native allocations.

The initial 4 MiB estimate gave a larger 1.53x whole-operation gain at 512
families but raised traced peak from 5.23 to 9.72 MB. The deliberate
[workspace refinement](workspace-results.json) compared 256 KiB, 1 MiB and
4 MiB with seven alternating pairs each. At 512/bound 86, 1 MiB still gave
1.33x with separated ranges and no overall peak increase. At bound 160,
1 MiB selects the original loop; its ratios around 1.00 are correctly treated
as noise. The [initial 4 MiB results](results-4mib.json) are retained rather
than overwritten. The smaller final tile trades some maximum speed for much
less additional storage.

## Correctness and checks

- All six final pruning benchmarks, seven complete-operation benchmarks,
  **212 numerical edge/helper comparisons**, and three state-doubling
  comparisons were exactly equal, including infinity/NaN masks.
- Checks include pure birth/death/zero-rate models, root mean 1 and 1.7,
  retentions 0/0.3/1, detection probabilities, family-specific missingness,
  entirely missing descendant subtrees, shuffled/repeated requested rows,
  impossible states, negative logs down to order 1e6, last short tiles,
  tiny workspace settings and a 513-state batch=1 fallback.
- Prototype-injected existing tests passed through the repository entrypoint:
  **102 passed in 8.46 s**, no skips. These include independent Yule and pulse
  sums, brute-force and CTMC references, missing-subtree integration,
  ascertainment, adaptive state doubling, unrenormalized tails, CLI outputs
  and the existing seeded small bootstrap. See [tests.log](tests.log).
- The same selection also passed with the initial 4 MiB prototype
  ([tests-4mib.log](tests-4mib.log)). Repeating it was justified by the changed
  tile dispatch. Ruff lint and formatting checks passed for these three
  experimental Python files.

The complete release suite, other platforms, fitted-search speed, opt-in
replicated scientific studies and additional empirical datasets were not run.
No scientific model, tolerance, seed, optimizer bound or inference procedure
was changed. Integration would still need the normal affected numerical/CLI
checks and delivery gates. Do not claim quiet-host or universal performance.

## Reproduction and source identity

All numerical/test runs used the common exclusive measurement lock, one
BLAS/OpenMP thread, warmups and seven alternating measurements. Other user
CPU jobs remained active. The initial unchanged-path control falsely appeared
1.31x faster with overlapping ranges, illustrating why small differences were
not used for adoption; its final control was 1.01x. Raw samples, loops and
spreads are retained. Python 3.12.14, NumPy 2.5.4, SciPy 1.17.1, pandas 3.0.6,
ETE 4.4.0, macOS ARM64; full metadata is in JSON. Current web manuals describe
SciPy 1.18, but the measured installed implementation was 1.17.1.

The immutable archive was `86c8f2c537769d852c17e2c789340f6c9de4ef54`,
NWKIT 0.43.55. `nwkit.__file__` and the target source SHA-256 are recorded
in results.json; imports resolve to
`/tmp/nwkit-perf-round2-20261022-rkzoGe/baseline/nwkit`, isolating this work
from another task's unrelated live-worktree edits.

```sh
cd /tmp/nwkit-perf-round2-20261022-rkzoGe/baseline
PY=/tmp/nwkit-perf-20261022-9fYwTD/venv/bin/python
DIR=/Users/kf/repos/nwkit/reviews/performance-candidates-round2-2026-10-22/03-wgd-log-pruning
$PY "$DIR/benchmark.py"
$PY "$DIR/benchmark.py" --tests
$PY "$DIR/refine_workspace.py"
```

The script sets thread variables before numerical imports, obtains the lock,
and injects the prototype only in its own process. The test command is
`python tools/check.py test -- -p candidate03_plugin tests/test_wgd_count_model.py tests/test_wgd_count.py -rs`
from the immutable source with the prototype directory on `PYTHONPATH`.

# Candidate 08: reuse stochastic-map uniformization preparation

**Verdict: adoptable within the tested low-draw, preparation-bound scope.**
For 32 tips with distinct branch lengths, a complete one-draw count operation
was **1.30x** faster and a complete one-draw history sample **1.58x** faster.
Their seven raw timing ranges were separated, and all compared output fields
were exactly equal. A regime/hidden-projection history case gave 1.25x.
At the documented default of 100 draws, both operations were about 1.01x with
overlapping ranges: **no default-workload speedup is established**.
Production sources, existing tests and Git history were not changed.

## New work and preserved semantics

The existing implementation already caches uniformization contexts by
generator identity and exact branch length, and caches backward endpoint
messages inside those contexts. This proposal does not repeat those caches.
It removes duplicated *parameter preparation*: `_validate_stochastic_map_work`
calculates `(omega, lambda, max_n)` for each distinct generator/length pair,
then context construction recalculates them. Histories calculate them again
for the segment guard. Each calculation includes a Poisson quantile and, where
applicable, a branch transition and the rare-endpoint-relative tail bound.

The prototype makes the work guard optionally return those already validated
tuples, passes them explicitly to context construction, and reuses them for
the history segment guard. All data is local to one operation, keys retain
the existing `(id(matrix), float(length))` convention, and the fit retains the
matrix objects while preparation is in progress. The tuple dictionary is
released immediately after contexts are constructed. No cross-call or global
data cache is introduced. Default guard calls still return `None`.

Work, dense-array, record and segment limits still run before contexts and
their large arrays are built. `_build_uniformization_context` retains its
own backward-message memory check before allocating arrays. The original
Poisson quantile, conditional tail target, cutoff extension, context formulas,
rate/length validation, root/node sampling and bridge algorithm are unchanged.
The optional tuple handoff is internal; it must originate from the same
operation's validated guard, rather than a saved cache or user-supplied cutoff.

## Complete-operation measurements

The [final results](results.json) use fixed fitted inputs, fresh preparation
and sampling on every call, warmups, seven alternating pairs and one numerical
thread under the study lock. Imports, input construction and rate fitting are
excluded; the count operation includes final table construction, while the
history operation returns all conditional draws and branch segments. These
are API-operation measurements, not whole-CLI or fitting speedups.

| Complete operation | Baseline median ms | Candidate median ms | Ratio |
| --- | ---: | ---: | ---: |
| 4 tips, repeated length, 1 count draw | 0.959 | 0.881 | 1.09x |
| 32 tips, unique lengths, 1 count draw | 15.574 | 11.983 | **1.30x** |
| 32 tips, unique lengths, 20 count draws | 84.765 | 81.309 | 1.04x |
| 32 tips, repeated length, 20 count draws | 63.914 | 63.974 | 1.00x |
| 32 tips, unique lengths, 1 history draw | 19.292 | 12.226 | **1.58x** |
| 32 tips, unique lengths, 20 history draws | 92.604 | 86.434 | 1.07x |
| 32 tips, repeated length, 20 history draws | 68.705 | 68.585 | 1.00x |
| 16 tips, unique lengths, 2 regimes/4 hidden states, 4 history draws | 18.208 | 14.551 | **1.25x** |
| 32 tips, unique lengths, default 100 count draws | 335.814 | 332.333 | 1.01x |
| 32 tips, unique lengths, default 100 history draws | 361.250 | 357.025 | 1.01x |

The one-draw 32-tip history ranges were 19.124–19.649 ms versus
12.154–12.579 ms; count ranges were 15.122–18.757 versus 11.548–12.344 ms.
Every one of their seven pairs improved. The smaller 20-draw unique-length
effects also had separated ranges in this run, but preparation is an ever
smaller fraction as draws increase. Repeated-length controls provide no
meaningful gain. The 100-draw ranges overlap and only four/five pairs improve.

The non-timed call audit confirmed that, for 62 distinct branches,
`_uniformization_parameters` calls drop **124→62** for counts and **186→62**
for histories. Separate setup-only measurements, including the same guards,
contexts and map specification, gave 1.74x/2.40x. These explain the mechanism
but are not substituted for the complete-operation results above.

Peak traced allocations were 424,743→415,376 bytes for the one-draw count
case and 219,798→209,879 bytes for one-draw histories. Repeated-length
20-draw controls added about 1.1 KB/0.1 KB. Default 100-draw histories were
2,003,766→2,003,511 bytes. These are `tracemalloc` peaks, not total RSS or
all native allocations; no general memory improvement is claimed. The added
metadata is linear in the number of distinct generator/length pairs.

The [initial run](initial-results.json) is retained. It preceded a private
guard-return compatibility refinement and formatting; the public formulas
were identical. All principal results were rerun with the final prototype.

## Equivalence and affected checks

- All ten measured complete operations compared every count-table field or
  every history field: states, identifiers, depths, lengths, node states,
  simulation indices, branch segments, counts and the exported history table.
  Every comparison was exactly equal.
- Four additional count/history audits covered homogeneous and regime/hidden
  cases. All final states of every created RNG were exactly equal (four
  generators for counts; eight, including event-time streams, for histories).
- Zero rates with zero-length branches, and a 65-state rare-endpoint chain,
  matched complete outputs and RNG states. The latter uses a fixed root state
  and a rare conditioned endpoint, retaining the same conditional cutoff.
- Five fault-injection checks preserved exception types/messages and observed
  zero context constructions before rejection: count work, history work,
  segment bound, record bound and dense-array bytes. The small injected limit
  does not alter the production limit; it proves check ordering in the prototype.
- Six preparation-only controls repeated calls after in-place generator and
  branch-length changes, comparing every map-specification field. These test
  cache lifetime, not a claim to refit those mutated scientific inputs.
- The final prototype was injected before collection of existing tests:

```text
python tools/check.py test -- -p candidate08_plugin tests/test_asr.py \
  tests/test_stochastic_history.py tests/test_stochastic_history_audit.py \
  -k 'stochastic or uniformization or bridge or history or segment' -rs
57 passed, 24 deselected in 7.24s
```

Those checks include independent two-state and asymmetric bridge integrals,
the rare-endpoint cutoff, real serial/two-worker reproducibility, legacy-count
equality for extended models, CLI tables/figures, rollback and memory/work
guards. The direct cutoff helper is intentionally unchanged; the additional
rare-chain whole-operation comparison also exercises the new handoff.
See [tests.log](tests.log) and [additional-results.json](additional-results.json).
Ruff lint/format and strict JSON checks pass for the experimental artifacts.

No full release gate, cross-platform run, new external backend, or replicated
population-calibration study was run. No seed, inference method, scientific
tolerance, model or search/work limit was relaxed. Integration still needs
affected entrypoint/CLI coverage and normal delivery gates. The additional
private parameter plumbing should be weighed against the narrow benefit:
this is most relevant to few requested histories and many unique branches.

## Reproduction and primary sources

Every substantial run used immutable snapshot
`86c8f2c537769d852c17e2c789340f6c9de4ef54` (0.43.55), isolating unrelated live
checkout changes. JSON records imported source paths and target SHA-256s.
The actual environment was CPython 3.12.14, NumPy 2.5.4, SciPy 1.17.1,
pandas 3.0.6, macOS ARM64. Unrelated host CPU jobs remained active, so these
are exploratory local results rather than quiet-host guarantees. Online SciPy
documentation may describe a later version than the one measured.

```sh
cd /tmp/nwkit-perf-round2-20261022-rkzoGe/baseline
PY=/tmp/nwkit-perf-20261022-9fYwTD/venv/bin/python
DIR=/Users/kf/repos/nwkit/reviews/performance-candidates-round2-2026-10-22/08-stochastic-preparation
$PY "$DIR/benchmark.py"
$PY "$DIR/additional_checks.py"
$PY "$DIR/benchmark.py" --tests
```

The scripts set numerical thread limits before importing NumPy, take the
common lock, and patch only their own processes. See [the study protocol](../PROTOCOL.md)
for recreating the exact source archive or environment.

- [Hobolth & Stone (2009)](https://arxiv.org/abs/0910.1683),
  [published DOI](https://doi.org/10.1214/09-AOAS247), provides primary context
  for endpoint-conditioned CTMC sampling and uniformization. This candidate
  retains that sampler and reuses only repeated deterministic preparation.
- [SciPy Poisson](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.poisson.html)
  documents the quantile/log-PMF operations whose repeated preparation is
  removed. The specific local speedup is established by measurements here.
- [NWKIT stochastic-map guide](../../../docs/guides/STOCHASTIC_MAPS.md)
  defines the authoritative RNG, conditional-tail and preflight contracts.

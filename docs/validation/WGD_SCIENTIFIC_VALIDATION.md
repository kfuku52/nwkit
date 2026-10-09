# Native WGD and Ks Scientific Validation

Audit date: 2026-10-01 UTC / 2026-10-02 JST. Status: experimental, conditional inference. This report
does not certify genome-wide WGD/SSD error rates, equivalence to Whale, or SIF
compatibility. No scientific classification threshold was relaxed to pass a
test. Input validation, numerical optimization, uncertainty calculations, and
evidence interpretation are separate audit targets.

## Reproducibility and scope

The executable specifications are:

- `tests/test_wgd_tree_history_validation.py`: independent forward gene
  histories, colored topology probabilities, and WGD-origin attribution.
- `tests/test_ksrate_scientific_validation.py`: independent additive-distance
  references, exact median ranks, and repeated-data confidence-interval coverage.
- `tests/test_wgd_count_scientific_validation.py`: independent forward counts
  and refitted, search-wide count-inference experiments.
- `tools/prepare_wgd_empirical_counts.py` and
  `tests/test_wgd_empirical_data.py`: audited public-data preparation, not an
  alternative test runner or a fixture-regeneration operation.

Use the repository checker for the full scientific studies, including slow
tests:

```bash
python tools/check.py test -- \
  --run-studies \
  tests/test_wgd_tree_history_validation.py \
  tests/test_ksrate_scientific_validation.py \
  tests/test_wgd_count_scientific_validation.py -s
```

Host studies used Python 3.12.14, NumPy 2.5.3, SciPy 1.18.1, NWKIT 0.43.36,
and a separate development
environment because this checkout's existing `.venv` executable was not usable
on this host. The original environment was not modified. The new history/Ks/data
studies also ran in GeneGalleon's source-bound Docker runtime
`local/genegalleon:subgenome-20261001`; 56 tests passed in 85.32 seconds before
the additional int64-overflow adapter regression. A read-only pytest-cache
warning was harmless. Docker evidence is not SIF evidence.

Docker used Python 3.12.14, NumPy 1.26.4, SciPy 1.17.1, and the source-bound
NWKIT 0.43.36. Its arm64/Linux image ID was
`sha256:8c03f1aec42d761e44ad2077212fb97e5e2a495d5b7729fe324fe56665e01568`.
The runtime wrapper confirmed the daily owned-upstream snapshot and current
container-input freshness. After the Ks arithmetic correction, 69 Docker
Ks/science/data-adapter tests passed in 56.28 seconds. Both numerical environments
reproduced 4,716/5,000 percentile and 5,000/5,000 conservative coverage.

## Independent gene-history oracle

The oracle simulates individual exponential duplication/loss waiting times,
species splits, per-copy pulse retention, and terminal detection. It prunes
unobserved descendants and unary nodes, and records the actual pulse origins.
It does not call the production transition kernel, likelihood, or simulator.
The Yule stem has duration proportional to `log(root_mean)`, without multiplying
that duration by a family rate category. Observational selection is applied
after category mixing, exactly as in the stated model.

Each of seven predeclared regimes used 60,000 independent histories. Comparison
included common colored topologies with at most six tips and at least 120
observations. Conditional origin probabilities were compared within each
topology, summing equivalent origin markers only where identical subtrees make
orientation unidentifiable. This is stronger than comparing total copy counts
or total origin mass alone.

| Regime | Seed | Selected histories | Topologies checked | Largest topology Z | Largest bounded origin Z |
| --- | ---: | ---: | ---: | ---: | ---: |
| Background null | 611 | 58,183 | 56 | 2.288 | 0 |
| Internal WGD, rate mixture, partial detection | 612 | 50,013 | 29 | 1.552 | 2.029 |
| Terminal WGD, strong loss | 613 | 42,421 | 37 | 2.884 | 2.475 |
| SSD after WGD | 614 | 51,202 | 24 | 2.112 | 2.303 |
| Heterogeneous branch/family rates | 617 | 44,416 | 34 | 2.153 | 0 |
| Critical duplication = loss | 618 | 57,857 | 50 | 2.671 | 2.057 |
| Retained pulse without background SSD/loss | 619 | 58,514 | 15 | 1.799 | 1.245 |

All fixed six-standard-error checks passed. The largest absolute topology
probability difference was 0.003969. For a bounded origin count `X` in `[0,K]`,
the standard-error upper bound uses the model mean `mu` and
`Var(X) <= mu * (K - mu)`. This avoids an invalid zero empirical variance when
rare genuine origins happen not to occur. With `K=1` this is the binomial
variance. Predicted impossible origins must never occur. Neither this numerical
check nor the conditional origin probabilities are calibrated WGD-occurrence
tests or posterior probabilities.

The host run had eight passing tests in 19.28 seconds, including a perfect-WGD
trace sanity check. Docker reproduced the study outcomes.

## Ks interval coverage and correction

Twenty independent additive-tree references cover asymmetric 20-species trees.
The interval study uses three species, additive true distances, 80 independent
families, and a shared lognormal family rate multiplier whose population median
is one. Pair dependence within a family is retained. All pairs are observed.
The percentile bootstrap is checked against direct family-row resampling, so
the coverage finding is not an implementation mismatch.

The fixed repeated-data study used seed 616, 5,000 datasets, and 199 percentile
bootstrap draws per dataset. Exact binomial Monte Carlo confidence intervals
describe study coverage, not a confidence interval for any one phylogeny.

| Primary interval | Covered | Coverage | Exact 95% Monte Carlo interval | Median width |
| --- | ---: | ---: | --- | ---: |
| Family bootstrap percentile | 4,716 / 5,000 | 0.9432 | [0.936421, 0.949455] | 0.0626324 |
| Pair-median Bonferroni | 5,000 / 5,000 | 1.0000 | [0.999262, 1.000000] | 1.311467 |

At nominal 95%, the percentile interval undercovers in this particular regime.
The new `--ci-method pair-median-bonferroni` uses noninterpolated binomial
order-statistic intervals for pair population medians, divides the error budget
over every used pair, and propagates simultaneous endpoints through
`FS + FO - SO` and the median across complete trios. Pairwise independence is
not required for the Bonferroni guarantee; IID families within each pair are
required. Too few pair families yield an unavailable primary interval, not a
discarded trio or a spuriously narrow interval. Twenty-one rank/reference cases
check the finite-sample binomial construction.

The standalone CLI keeps `family-bootstrap-percentile` as its compatible
default, explicitly described as approximate. GeneGalleon's WGD workflow now
requests `pair-median-bonferroni` and reads the primary endpoints. Bootstrap
diagnostics remain in separate columns. The conservative interval was about
21 times wider in this experiment: improved coverage has a substantial power
cost. This study does not validate hidden paralogy, count-dependent ortholog
selection, non-IID gene families, gene conversion, or saturated dS estimates.
See the [Ks guide](../guides/KSRATE.md) and
[NIST median-interval reference](https://itl.nist.gov/div898/software/dataplot/refman1/auxillar/mediancl.htm).

Three regression cases additionally reproduced float overflow in a finite pair
median, a cancellable focal correction, and simultaneous endpoint propagation.
Stable midpoints and cancellation-first arithmetic now preserve finite results;
a genuinely unrepresentable correction is rejected. Subnormal values and
negative corrected distances remain supported. Nonfinite family weights and
integer accumulation overflow are rejected before median calculation.

## Public yeast count-data protocol

The public SPIMAP resource distributes count families derived from the fungi
dataset of Butler et al. (2009). Downloads were not substituted with a synthetic
fixture or edited to resemble a known event:

- [Raw 16-species count matrix](https://compbio.mit.edu/spimap/pub/data/real-fungi-butler2009.counts.txt)
  SHA256 `e771f3f67cb0c156a505529c91d01727bc60ba3af95d5328f237a7daff084ac2`.
- [Dated species tree](https://compbio.mit.edu/spimap/pub/config/fungi.stree)
  SHA256 `f40484af23f9eda39fc7907bb6da07d37df2aae74838514c1c0a180e519b9560`.
- [Resource description](https://compbio.mit.edu/spimap/index.shtml).

The actual raw download has 9,209 rows. Presence in every immediate root-child
clade retains 3,923 rows. These are the observed data counts, not a forced
published sample-size match. No family is excluded for large copy counts; the
resource's separate filtered file is not used because its copy-count cap is not
represented in this likelihood's ascertainment model. Uniform sampling without
replacement, seed 620, selects 500 families and retains a family with 84 copies.
The stable family labels are one-based source-row labels because the raw file
does not contain original family identifiers. The selection audit records all
9,209 rows.

Prepare a fresh study directory with:

```bash
python tools/prepare_wgd_empirical_counts.py \
  --counts RAW_COUNTS --species-tree SPECIES_TREE \
  --output NEW_STUDY_DIRECTORY --max-families 500 --seed 620
```

The exploratory protocol scans every positive non-root branch, using
root-clade ascertainment, terminal/internal background regimes, fixed gamma
shape 1 with four rate categories, multiplicity 2, and fractions
0.25/0.5/0.75. Numerical state convergence is required, not waived for the
large-copy family. Known WGD labels are used only for independent post-fit
comparison, not family selection or candidate restrictions. A scan
without bootstrap has no significance claim. A single positive-control dataset
cannot establish false-positive rates, sensitivity across genomes, or robustness
to arbitrary lineage-specific SSD and annotation errors.

The independent positive-control lookup is the incoming branch of the MRCA of
`scer`, `scas`, and `cgla` (six descendant tips in this supplied
tree). The shared event is established by
[Scannell et al. (2006)](https://www.nature.com/articles/nature04562); the count
method's benchmark use is described by
[Rabier et al. (2014)](https://pubmed.ncbi.nlm.nih.gov/24361993/).
This lookup is not supplied as a candidate restriction or fit parameter.

An additional budget/sampling check uses 100 uniformly selected families with
the same selection seed 620, without a copy-number filter. Its largest observed
copy count is five. This run was specified and started before observing either
scan's fitted scores; it is not a subset selected for agreement with a known
event. Both selection audits are preserved. The exact-integer data adapter also
preserves values beyond float's exact-integer range; its corrected preparation
produced byte-identical 500-family counts and selection audits for this dataset.

### Executed public-data outcomes

The 100-family scan completed all 30 non-root branches and all three event
fractions. The background log likelihood was `-602.1347451516194`. Doubling
the state limit from 8 to 16 changed any background family log likelihood by
at most `3.56124e-11`; the largest event and burst errors across candidates
were `3.76502e-8` and `2.19874e-8`, below the unchanged `1e-7` tolerance.
No bootstrap was requested (`--bootstrap 0`), so every candidate is explicitly
`not_calibrated`, with unavailable p-values.

| Rank | Branch/clade | LR versus background | Burst AIC minus WGD AIC | Burst nuisance bound |
| --- | --- | ---: | ---: | --- |
| 1 | `spar` terminal branch 30 | 12.047515 | 0.437395 | Reached |
| 2 | `sbay` terminal branch 26 | 10.684742 | -0.384924 | Not reached |
| 3 | Five-species branch 7, excluding `scas` | 7.529554 | 2.840751 | Reached |
| 4 | Known six-species WGD stem, branch 3 | 5.040923 | -9.950497 | Not reached |

A positive AIC difference favors WGD; a negative difference favors the burst
alternative. The known WGD stem ranked fourth and its burst alternative was
preferred by about 9.95 AIC units. This is **not a successful positive-control
detection**. Thirteen of the 30 burst fits reached a nuisance bound, including
the highest-scoring branch. Event-only convergence flags would conceal this
qualification even though every event fit was away from nuisance bounds.

The 500-family scan initially failed because a very small, mathematically
nonzero family probability underflowed. The repaired log-domain calculation
gave finite log likelihood `-24932.124092089856` at the same failing parameters,
including `-2817.6238453491324` for the affected family. A separate 500-family
background fit completed with log likelihood `-3747.7995181105543`; doubling
its state limit from 86 to 172 gave zero represented family-log-likelihood
change. No large-copy family was discarded or censored.

The ensuing 500-family full scan stopped again, this time at strict optimizer
stationarity certification: the three starts reported function-value stagnation
but independently evaluated projected gradients were `4.78744e-5`,
`5.47346e-4`, and `9.45098e-5`, exceeding `1e-5`. No candidate table or model
was published from this failed scan. A validated background fit is not a
completed all-branch scan. The follow-up numerical diagnosis and any completed
rerun are recorded separately below rather than replacing this failed attempt.

The completed 100-family snapshot used count model SHA256
`f29271372a2e38f4916d75581b42af3b989d6fe6e3d00f402adc42956fd38480`,
fit source `3e536676769d93e64cff905c9877c5cb412d2f0013a726a2142414b42d0f6911`,
and producer `9db354d83eea4ffd8eed5cb952aa45646856dd165bef47f711b45118b2772e44`.
These are run provenance, not pinned dependency defaults. Subsequent additive
producer diagnostics expose both background and branch-burst nuisance bounds.

The executed exploratory CLI configuration was:

```bash
python -m nwkit wgd-count --infile fungi.nwk --counts prepared100/counts.tsv \
  --ascertainment root-clades --rate-model terminal-internal \
  --family-gamma-shape 1 --family-rate-categories 4 \
  --event-fractions 0.25,0.5,0.75 --multiplicity 2 --bootstrap 0 --seed 621 \
  --max-states 512 --max-iterations 400 \
  --outfile candidates100.tsv --model-out model100.json
```

The failed 500-family scan used the same configuration, with `prepared/counts.tsv`
and separate output paths. Existing archived outputs must not be overwritten
when repeating either analysis.

### Numerical follow-up

The stationary-fit failure was isolated to the first branch-burst comparator
(node 1), after that branch's three WGD fractions had passed. Three-point
finite differences had an `h^2` truncation bias in the stiff log-root-mean
coordinate near `0.00146`. Paired per-family differencing alone agreed with
scalar differencing and did not remove this bias. At the last failed point,
the three-point root gradient was `-9.45098e-5`, whereas a separate complex-step
root derivative was `2.44121e-4`: the point really was nonstationary.

Richardson derivatives cancel the leading stencil error. A higher-order
polishing phase is used only if no initial start passes stationarity, within
each start's remaining iteration budget. Verification uses twice the polishing
step. The projected-gradient threshold `1e-5`, state tolerance `1e-7`,
parameter bounds, statistical models, and all 500 families are unchanged.

The isolated comparator then completed in 316.56 seconds with log likelihood
`-3632.4422771034297`, state-doubling error `2.84217e-14`, and no nuisance
boundary. Its verified projected gradient was `4.77820e-6`; the complex-step
root derivative was `4.74496e-6`, with consistent Richardson checks at step
factors 0.5, 1, 2, and 4. Complex-step verification checks the root derivative,
not an independent oracle for the entire pruning model.

The **full 500-family grid remains uncompleted**. The repaired isolated
comparator is not a WGD-detection result. Detailed first-failure parameters,
derivative checks, commands, final source hashes, and limitations are preserved
in [the follow-up JSON](wgd_count_followup_20261002.json) and the local archive.

## Replicated count inference

Independent Gillespie histories produced 64 families per dataset on three
species/four non-root branches. Every branch was searched at its midpoint.
Each of 36 dataset runs refitted 19 null bootstrap draws (684 total), with
unchanged alpha 0.05, state tolerance `1e-7`, state cap 128, and 200 optimizer
iterations. No final dataset or bootstrap draw failed or was discarded. These
are bounded pilot experiments: the smallest nonzero p-value is 0.05, and
replicate counts are inadequate for reliable population error-rate certification.

| Regime | Conditional support | Exact 95% Monte Carlo interval | Elapsed seconds |
| --- | ---: | --- | ---: |
| Homogeneous null | 0 / 12 | [0, 0.264648] | 413.32 |
| Local SSD, misspecified homogeneous fitted null | 1 / 12 | [0.002108, 0.384796] | 439.97 |
| Pure WGD | 4 / 4 | [0.397635, 1] | 172.56 |
| Strong WGD | 4 / 4 | [0.397635, 1] | 162.24 |
| Known branch/family rate heterogeneity | 0 / 2 | [0, 0.841886] | 521.13 |
| Known detection and fixed missing annotation masks | 0 / 2 | [0, 0.841886] | 135.08 |

Data seeds are `170000..170011`, `171000..171003`, `172000..172003`,
`173000..173011`, `174000..174001`, and `175000..175001` respectively.
Each bootstrap seed is its data seed plus 740000. Pure-WGD copy data are
deterministic under this truth; only bootstrap seeds vary. The heterogeneous
experiment fits the supplied regimes and discrete family mixture, and detection
probabilities are known, not estimated. It does not certify inference under
unknown annotation/rate heterogeneity or continuous gamma rates.

In the deliberately misspecified local-SSD regime, the raw WGD-versus-background
test rejected in **12/12** datasets. The branch-burst AIC gate prevented support
in 11 cases but still admitted false WGD support for data seed `173002`.
Consequently the AIC gate is not calibrated protection against an SSD null.
Both WGD regimes recovered and supported the true branch in 4/4 cases, with the
wide intervals above. These results do not measure the full GeneGalleon workflow's
error rates or prove reliable sensitivity for ancient low-retention events.

After the stationarity repair, the existing 59-draw null calibration passed
again. Three existing local-SSD datasets (`173000..173002`) were replayed as
regressions: counts, support decisions, p-values, and every bootstrap statistic
were identical. These are reused datasets, not three new independent replicates.
The false-support dataset `173002` had no nuisance bound in any of the three
fits, so the new boundary diagnostics do not remove this model-misspecification
false positive.

The recorded count-study snapshot's input counts, candidate outputs, and all bootstrap
statistics were identical to its preceding pre-underflow-fix study snapshot.
The arithmetic repair therefore did not change those finite-probability study
outcomes. Docker count tests for that snapshot passed 100 cases in 91.06 seconds;
the separate existing 59-draw slow calibration passed in 100.66 seconds on the
host. Timings include concurrent work and are not performance comparisons.

The exact invocations are the existing checker with
`NWKIT_WGD_COUNT_STUDY_REPLICATES=12`, `4`, or `2` and `-m slow -k`
selecting respectively `study and (null or local-ssd)`,
`study and (pure-wgd or strong-wgd)`, or
`study and (heterogeneous or missing-annotation)`. All use
`tests/test_wgd_count_scientific_validation.py -s`, with
`OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1`.
Study JSON retains the full protocol, all candidates, failures, and bootstrap
statistics. The archived run manifest is
`/tmp/nwkit-count-validation-manifest-20261001.json`; its six result files and
full hashes identify the executed artifacts, rather than dependency defaults.
To replay these protocols on the current checkout, also pass `--run-studies`;
the original invocations above describe the archived runs.

## Evidence interpretation and remaining limits

Branch-rate heterogeneity can mimic WGD signals. The explicit branch-burst
comparison is useful but does not test against arbitrary SSD histories. The
importance of those assumptions is discussed in the primary
[Whale paper](https://academic.oup.com/mbe/article/36/7/1384/5475503).
Fixed-tree results condition on the supplied topology and rate/observation model;
they do not average over gene-tree, species-tree, rooting, or allopolyploid donor
uncertainty. Families originate at the species root in the current count model.
There is no model for de novo later family origins or horizontal transfer.

The GeneGalleon audit checks malformed/mismatched counts and gene identifiers,
isoform inflation, full-locus coordinates, overlapping loci, dependent block
replication, event/species matching, cross-child duplicate pairs, dS validity and
coverage, unrooted trees, inherited NHX evidence, and evidence/output mutation.
Unsupported or inconsistent evidence remains unresolved; it is not converted
to terminal SSD. WGD-supported and SSD-supported classifications remain
evidence-supported hypotheses rather than probabilities of historical truth.

Count evidence now requires explicit false values for event, background, and
branch-burst nuisance-bound diagnostics before it can support classification.
A missing/unknown diagnostic or a bounded comparator remains unresolved;
the event-only flag is insufficient. Eight new consumer regressions reproduced
the previous acceptance of bounded or incompletely described comparator fits.
This stricter evidence gate is not the pilot study's conditional-support metric,
which includes the bootstrap and AIC gates but not the entire downstream
genomic evidence workflow.

Explicit `[&U]` and NHX no/unknown declarations are now honored. Counts accept
explicitly rooted root polytomies; reconciliation retains its bifurcation
requirement. A user-provided rooted declaration is not independently established
as the biological root. Block replication now means connected components of
blocks that do not share qualifying anchor loci, not proven statistical or
spatial independence. Validating more complete annotation bounds and conservative
components can reduce apparent evidence; this is not a changed scientific
threshold.

Source provenance now fingerprints all NWKIT/KF package Python sources. Existing
genome evidence and subsequent classifications must be regenerated through the
normal reviewed rebuild procedure. Do not manually rehash an old result to make
it appear current. Original OrthoFinder inputs are not converted or silently
collapsed; known same-locus isoforms are rejected at preflight.

## Final affected checks

Final count source hashes match the follow-up JSON. The production hashes were
unchanged across the final parent checks. Slow tests were deliberately selected
separately from the quick lane; deselection is not reported as a pass.

| Check | Runtime | Result |
| --- | --- | --- |
| `tools/check.py quick -- tests/test_wgd_*.py tests/test_ksrate*.py tests/test_cli_contracts.py` | Host development Python | 381 passed, 15 slow deselected, 85.62 s; Ruff lint/format and mypy passed |
| Four `workflow/tests/test_wgd_*.py` files through `run_in_runtime.sh`, with source-bound NWKIT/KF | GeneGalleon Docker | 150 passed, 81.87 s, including four real runtime integrations |
| Five count/scenario/scientific files, `tools/check.py test`, `-m 'not slow'` | Docker | 116 passed, 7 slow deselected, 103.19 s |
| Existing 59-draw null calibration | Host | 1 passed, 6 deselected, 100.85 s |
| Three-dataset SSD regression replay | Host | 1 study test passed, 13 deselected, 122.50 s |
| `tools/check_maintainability.py` | Host | All hard limits passed; unrelated baseline warnings retained |
| Ruff on `workflow` and `container`; `dev config-check` | Host static / Docker config | Passed; eight entrypoints and 23 shared parameters |
| `git diff --check` | All three checkouts | Passed |

The exact final GeneGalleon command, run from its repository root, was:

```bash
GG_TEST_RUNTIME=docker \
GG_CONTAINER_DOCKER_IMAGE=local/genegalleon:subgenome-20261001 \
GENEGALLEON_DOCKER_EXTRA_BINDS=$'/Users/kf/repos/nwkit\n/Users/kf/repos/kfFractBias' \
bash workflow/tests/run_in_runtime.sh \
  env PYTHONPATH=/Users/kf/repos/nwkit:/Users/kf/repos/kfFractBias/src \
  python -m pytest -q -rs -p no:cacheprovider \
  workflow/tests/test_wgd_evidence.py workflow/tests/test_wgd_ssd_contracts.py \
  workflow/tests/test_wgd_ssd_gene_stage.py workflow/tests/test_wgd_ssd_runtime.py -x
```

`dev lint` did not complete as an aggregate command: host Bash 3.2 lacks the
required `mapfile`, and the Docker runtime lacks Ruff. Separate shell syntax,
Ruff, and configuration checks passed, but they are not labeled a successful
`dev lint` invocation. SIF compatibility, a full release lane, realistic
genome-wide classification error rates, and the complete 500-family yeast grid
remain unverified. No commit, push, or release was performed.

## Preserved run artifacts

The local research archive is
`/Users/kf/repos/genegalleon/workspace/output/validation/wgd-science-20261001.HC6sei`.
It contains the six complete replicated-count JSON results, their run manifest,
the first underflow reproduction and repaired 500-family background fit, and
the full `public-yeast/` study directory, including unmodified downloads,
both sampling audits, the 100-family candidate table, and model diagnostics.
Original temporary artifacts were preserved as well. This directory contains
research evidence, not disposable workflow output or a new memory store.
The `count-followup/` subdirectory additionally preserves the first stationary-fit
failure, prefix log, numerical proof scripts/results, frozen follow-up JSON,
and the post-fix SSD replay. The scripts are diagnostic research artifacts,
not a second repository test runner.

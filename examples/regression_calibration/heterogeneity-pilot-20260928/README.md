# RSC heterogeneity and error pilot

This is a separate **200-dataset-per-case pilot** with master seed 20261020,
199 bootstrap refits per dataset, and the frozen source at commit
`5cf3df31f1a9547ca6e5dfa5fec7611ff9e492e2`. It uses four prespecified
20-event RSC null cases: independent missing observations, known response
sampling variance from biological replicates, estimated response SE from
replicates, and predictor measurement error. Wald, original coefficient
bootstrap, and opt-in studentized bootstrap use the same generated data.
The [audit](audit.json) verified all 800 datasets and retained source/input
records; no failed fit or interval is removed from the denominator.

| Condition | Studentized null rejections / 200 | Studentized 95% coverage / returned | Intervals returned / 200 |
| --- | ---: | ---: | ---: |
| MCAR missingness | 9 (4.5%) | 188/200 (94.0%) | 200 |
| Known response sampling variance | 8 (4.0%) | 190/200 (95.0%) | 200 |
| Estimated response SE | 8 (4.0%) | 190/200 (95.0%) | 200 |
| Predictor measurement error | 13 (6.5%) | 174/196 (88.8%) | 196 |

The predictor-error case has four fit failures with non-positive-definite
errors-in-variables coefficient information. Its conditional interval coverage
is below 95%, and correct-interval delivery is 174/200 (87.0%). For comparison
in that same case, Wald covered 182/196 and the original coefficient bootstrap
178/196; this small pilot does not establish a winner. It does show that the
studentized improvement seen in the two-event and five-event confirmation
cannot be generalized to predictor measurement error. The 200-dataset Wilson
intervals are broad; inspect [summary.json](summary.json) for all methods,
availability, failures, widths, and Monte Carlo intervals. No default or
acceptance criterion was changed from these results.

Regenerate into a fresh directory with the archived source and the cases,
methods, seed, replicate counts, and environment in [protocol.json](protocol.json).
The retained `records.jsonl.gz` includes the generated data and per-dataset
diagnostics; `source.tar.gz` preserves exact implementation inputs.

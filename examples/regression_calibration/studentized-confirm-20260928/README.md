# Independent RSC studentized-bootstrap confirmation

This run used the frozen source at commit `5cf3df31f1a9547ca6e5dfa5fec7611ff9e492e2`,
master seed `20260929`, 5,000 new independent datasets per case, and 199
bootstrap draws per dataset. The two prespecified raw-tip RSC null cases differ
in the number of species events (`rsc-e2` and `rsc-e5`); both generate from the
physical species-event and paralog-lineage model. All 10,000 point fits and
intervals were available. All 1,990,000 studentized bootstrap refits succeeded.
The source, exact protocol, per-dataset records, summary, and independent
integrity audit are retained here.

| Case and method | Null rejections / 5,000 | 95% intervals covering truth / 5,000 |
| --- | ---: | ---: |
| rsc-e2 oracle | 244 (4.88%) | 4,756 (95.12%) |
| rsc-e2 Wald | 1,489 (29.78%) | 3,511 (70.22%) |
| rsc-e2 coefficient bootstrap | 1,449 (28.98%) | 3,464 (69.28%) |
| rsc-e2 studentized bootstrap | 221 (4.42%) | 4,692 (93.84%) |
| rsc-e5 oracle | 271 (5.42%) | 4,729 (94.58%) |
| rsc-e5 Wald | 626 (12.52%) | 4,374 (87.48%) |
| rsc-e5 coefficient bootstrap | 608 (12.16%) | 4,309 (86.18%) |
| rsc-e5 studentized bootstrap | 239 (4.78%) | 4,678 (93.56%) |

The prespecified Wilson-interval criterion for null rejection was met only for
`rsc-e5` studentized bootstrap. Its coverage Wilson interval overlaps the
lower 93% acceptance boundary; for `rsc-e2`, the rejection Wilson interval
overlaps the lower 4% boundary. Thus this confirmation supports a substantial
improvement over Wald and the original coefficient bootstrap in these cases,
but **does not meet every promotion criterion**. The opt-in method is not made
the default. These cases do not cover heterogeneous known or estimated
sampling errors, missingness, other tree shapes, or variable selection.

To verify the retained evidence from the repository root:

```sh
PYTHONPATH=. python tools/verify_regression_calibration.py \
  examples/regression_calibration/studentized-confirm-20260928
```

To regenerate it into a fresh directory using the archived source rather than
the current checkout, extract `source.tar.gz` and run its
`tools/validate_regression_calibration.py` with `--cases rsc-e2,rsc-e5`,
`--replicates 5000`, `--bootstrap-replicates 199`, `--seed 20260929`, and
`--methods wald,parametric-bootstrap,studentized-bootstrap,oracle`. The
remaining environment and worker settings are in `protocol.json`.

The archived source is the execution source, not the current checkout. The
protocol and audit carry SHA-256 hashes of the source and evidence files.

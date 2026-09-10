# Independent calibration validation

Seed 20260912; B=199; stagewise nominal level 0.05. Plug-in tests and finite nuisance grids do not guarantee uniform 5% error.

| Family | Truth | Root | Completed | Any shift | True partition |
|---|---|---|---:|---:|---:|
| primary | null | OUfixedRoot | 50/50 | 5/50 | 45/50 |
| primary | null | OUrandomRoot | 50/50 | 1/50 | 49/50 |
| primary | single | OUfixedRoot | 50/50 | 27/50 | 21/50 |
| primary | single | OUrandomRoot | 50/50 | 25/50 | 22/50 |
| primary | distinct | OUfixedRoot | 50/50 | 31/50 | 30/50 |
| primary | distinct | OUrandomRoot | 50/50 | 31/50 | 30/50 |
| primary | convergent | OUfixedRoot | 50/50 | 18/50 | 17/50 |
| primary | convergent | OUrandomRoot | 50/50 | 24/50 | 24/50 |
| weak_pull | null | OUfixedRoot | 25/25 | 0/25 | 25/25 |
| weak_pull | null | OUrandomRoot | 25/25 | 0/25 | 25/25 |
| weak_pull | convergent | OUfixedRoot | 25/25 | 0/25 | 0/25 |
| weak_pull | convergent | OUrandomRoot | 25/25 | 2/25 | 0/25 |
| known_error | null | OUfixedRoot | 25/25 | 0/25 | 25/25 |
| known_error | null | OUrandomRoot | 25/25 | 1/25 | 24/25 |
| known_error | convergent | OUfixedRoot | 25/25 | 10/25 | 9/25 |
| known_error | convergent | OUrandomRoot | 25/25 | 8/25 | 8/25 |
| sixteen_tips | null | OUfixedRoot | 25/25 | 3/25 | 22/25 |
| sixteen_tips | null | OUrandomRoot | 25/25 | 0/25 | 25/25 |
| sixteen_tips | convergent | OUfixedRoot | 25/25 | 25/25 | 22/25 |
| sixteen_tips | convergent | OUrandomRoot | 25/25 | 24/25 | 23/25 |

The source snapshot and `protocol.json` record the exact implementation used.
`audit.json` verifies regenerated inputs and source hashes and adds a separate
250-dataset null check of the default mode without convergence candidates
(10/250 false selections). `default-null-records.jsonl.gz` retains those fits.
See [method and interpretation](../../../SHIFT_CALIBRATION.md).

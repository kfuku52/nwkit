# Independent calibration validation

Seed 20260915; B=199; stagewise nominal level 0.05. Plug-in tests and finite nuisance grids do not guarantee uniform 5% error.

| Family | Truth | Root | Completed | Any shift | True partition |
|---|---|---|---:|---:|---:|
| primary | null | OUfixedRoot | 50/50 | 1/50 | 49/50 |
| primary | null | OUrandomRoot | 50/50 | 1/50 | 49/50 |
| primary | single | OUfixedRoot | 50/50 | 23/50 | 22/50 |
| primary | single | OUrandomRoot | 50/50 | 29/50 | 25/50 |
| primary | distinct | OUfixedRoot | 50/50 | 29/50 | 27/50 |
| primary | distinct | OUrandomRoot | 50/50 | 26/50 | 24/50 |
| primary | convergent | OUfixedRoot | 50/50 | 24/50 | 24/50 |
| primary | convergent | OUrandomRoot | 50/50 | 28/50 | 27/50 |
| weak_pull | null | OUfixedRoot | 25/25 | 3/25 | 22/25 |
| weak_pull | null | OUrandomRoot | 25/25 | 4/25 | 21/25 |
| weak_pull | convergent | OUfixedRoot | 25/25 | 1/25 | 0/25 |
| weak_pull | convergent | OUrandomRoot | 25/25 | 1/25 | 0/25 |
| known_error | null | OUfixedRoot | 25/25 | 0/25 | 25/25 |
| known_error | null | OUrandomRoot | 25/25 | 1/25 | 24/25 |
| known_error | convergent | OUfixedRoot | 25/25 | 6/25 | 6/25 |
| known_error | convergent | OUrandomRoot | 25/25 | 9/25 | 8/25 |
| sixteen_tips | null | OUfixedRoot | 25/25 | 1/25 | 24/25 |
| sixteen_tips | null | OUrandomRoot | 25/25 | 0/25 | 25/25 |
| sixteen_tips | convergent | OUfixedRoot | 25/25 | 23/25 | 21/25 |
| sixteen_tips | convergent | OUrandomRoot | 25/25 | 22/25 | 20/25 |

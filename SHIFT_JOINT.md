# Small-tree joint shift and convergence reference

The research tools enumerate shift locations and shared-optimum groups together,
then fit each declared model through the public kfl1ou `fit_OU` API. They compare
this reference with NWKIT's exhaustive unconstrained search followed by backward
convergence. The production `shift` command and its defaults are unchanged.

The held-out pilot found lower BIC scores in 2 of 12 datasets, but both methods
recovered the true shared-optimum tip partition in 4 of 12. This demonstrates
search omissions in this sample, while also showing that wider discrete search
alone does not resolve its recovery failures. It is not an error-rate estimate.

## Scope and comparability

The reference accepts positive, fully bifurcating ultrametric trees with at most
16 tips and at most two shifts. For each branch configuration it enumerates all
set partitions of the background and shift identities. It retains disconnected
convergence and a nested return to the background optimum. Equalities making a
shift identical to its parent are removed. Configurations with an empty ancestral
regime are excluded as rank deficient; at this size limit this is the pair of
root-child edges. Thus coverage refers to this identifiable candidate space.
Alternative branch configurations with the same tip partition are retained:
their OU histories and resulting means need not be equivalent.

Eight-tip trees have 195 candidates; sixteen-tip trees have 899. Each candidate
refits continuous parameters with the backend, retaining the baseline's root
model, known observation variances, alpha bounds and criterion. There is no
independent certification of the continuous global optimum. The result is the
best returned fit in the enumerated space, not a proof of globally optimal OU
inference. Failed candidates remain in the ledger and invalidate complete
coverage. The selected model's mean and Gaussian likelihood are independently
recomputed in Python, and shared-optimum equalities are checked.

A baseline refit checks that the comparison reproduces the current two-stage
score. Each all-singleton grouping is additionally fitted without convergence
constraints. These represent the same model, so their information penalties
must agree even if numerical optimization returns slightly different likelihoods.
The audit compares `score + 2 * log_likelihood` between representations. An audit
failure sets `score_improvement` to JSON `null`; raw scores remain diagnostic.

## Important pBIC finding

The findings below describe the original released-3.0.9 run, preserved unchanged.
The subsequent [pBIC correction and replay](SHIFT_PBIC.md) identifies and fixes
the causes in kfl1ou; both replayed pBIC cases pass the representation audit.


With kfl1ou 3.0.9, the free and all-singleton constrained representations can have
the same likelihood but different pBIC scores. The two held-out pBIC fits had
maximum penalty discrepancies of 68.8668 and 7.4692. Both fail the equivalence
audit; neither is evidence that joint search improves pBIC. One also had four
singular candidate fits out of 195. Every failure is retained.

Inspection of the installed backend indicates that the two pBIC paths use
different coefficient representations in their determinant penalties. This is
a suspected explanation, not a validated correction. The reference does not
replace the backend's formula. Comparisons across these representations, and
model weights derived from them, require resolving this discrepancy first.
All 1,980 BIC equivalence checks across the 12 cases passed (maximum discrepancy
below 2.1e-13). This checks representation consistency; it does not establish
statistical calibration of BIC or pBIC.

## Prespecified held-out pilot

The master seed is `20260910`, distinct from the preceding simulation pilot.
Each dataset uses `SeedSequence([master_seed, case_index])`. The grid contains
one dataset per cell: all four scenarios on eight-tip trees under both root
models, plus convergence with SE 0.2 on eight tips and convergence on sixteen
tips without observation error. Nonzero generating shifts have effect 2.
BIC is evaluated in every cell; pBIC is an additional diagnostic on the two
eight-tip convergent cases without observation error. There are no bootstrap
replicates or tuned criterion constants.

All 3,748 BIC candidate fits completed. The true shared partition was represented
among the fitted candidates in every case. In all eight cases where the selected
partition was wrong, its BIC was lower than the best candidate with the true
partition. Continuous fitting, the criterion and finite data therefore remain
possible explanations; this experiment does not separate those causes.

| Case | Tips / truth | Root | SE | BIC decrease | Shared truth recovered, two-stage → joint | Best true-partition BIC minus joint BIC |
|---|---|---|---:|---:|---|---:|
| 00 | 8 / null | fixed | 0 | 0 | no → no | 2.6424 |
| 01 | 8 / null | random | 0 | 0 | no → no | 3.8387 |
| 02 | 8 / single | fixed | 0 | 0.1701 | yes → yes | 0 |
| 03 | 8 / single | random | 0 | 0 | no → no | 2.8110 |
| 04 | 8 / distinct | fixed | 0 | 1.8244 | no → no | 3.1356 |
| 05 | 8 / distinct | random | 0 | 0 | yes → yes | 0 |
| 06 | 8 / convergent | fixed | 0 | 0 | no → no | 3.0884 |
| 07 | 8 / convergent | random | 0 | 0 | no → no | 8.8415 |
| 08 | 8 / convergent | fixed | 0.2 | 0 | no → no | 1.3384 |
| 09 | 8 / convergent | random | 0.2 | 0 | no → no | 3.1269 |
| 10 | 16 / convergent | fixed | 0 | 0 | yes → yes | 0 |
| 11 | 16 / convergent | random | 0 | 0 | yes → yes | 0 |

Partition metrics ignore label permutations. A selected point estimate is
reported; the table does not claim uniqueness among tied models. These twelve
datasets are diagnostic examples, not enough replication to compare accuracy
or generalize the observed proportions.

## Reproduce and inspect

From a checkout with NWKIT's Python dependencies and kfl1ou >= 3.0.9 installed:

```bash
PYTHONPATH=. python tools/validate_shift_joint_pilot.py \
  --output /tmp/shift-joint --rscript Rscript --seed 20260910
PYTHONPATH=. python tools/summarize_shift_joint.py \
  --run /tmp/shift-joint --output /tmp/shift-joint-evidence
```

A single input directory needs `tree.nwk` and `traits.tsv` with columns
`leaf_name`, `value`, and `se`. Optional `truth.json` uses the simulation
harness schema. Existing output directories are refused.

```bash
PYTHONPATH=. python tools/validate_shift_joint.py \
  --case /tmp/shift-joint/case-06 --output /tmp/shift-joint-single \
  --criterion BIC --root-model OUfixedRoot --max-shifts 2 --rscript Rscript
```

The raw run retains every candidate score, error and warning, the generated R
script, baseline and winning RDS objects, tip predictions, branch effects,
inputs and execution source snapshots. The portable
[evidence bundle](examples/shift/joint-validation/) includes all candidate
records, model tables, inputs, seeds and execution source text without binary
R objects. Its exporter checks input hashes, winner likelihood, score audit,
truth metrics and source snapshots before publishing a complete bundle.
These are consistency checks, not authentication or independent optimization.

Candidate enumeration tests compare against an independent label-assignment
construction. Real R integration covers both root models with and without
known SEs. Regression tests check score comparability and preserve quoted,
multiline warnings and candidate failures.

Before model averaging, resolve the pBIC representation discrepancy and expand
joint-search validation with more independent draws, tree shapes and effect
sizes. Any later candidate weights must preserve shared-optimum constraints.

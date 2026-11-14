# RADTE exact-sequence interval comparison

The [protocol](protocol.json) was fixed before these 600 new independent JC69
families were generated. Each cell uses four gene tips, 2,000 sites, a true
root-duplication age of 20, exact sequence likelihood, and conditional joint
MAP with estimated rate SD. The same point fit was used for Laplace,
studentized-curvature, and ordinary profile intervals. Each method targets
nominal 95% coverage. These results describe this specific estimator and
generator; they do not validate marginal RADTE or external GY94 data.

| Generating log-rate SD | Method | Returned / 200 | Covered / returned | Correct return / 200 | Median width |
| --- | --- | ---: | ---: | ---: | ---: |
| 0.3 | Laplace | 198 | 169 (85.4%) | 169 | 9.51 |
| 0.3 | Studentized | 200 | 196 (98.0%) | 196 | 18.94 |
| 0.3 | Profile | 200 | 167 (83.5%) | 167 | 9.92 |
| 0.6 | Laplace | 125 | 97 (77.6%) | 97 | 14.47 |
| 0.6 | Studentized | 200 | 193 (96.5%) | 193 | 40.09 |
| 0.6 | Profile | 200 | 162 (81.0%) | 162 | 21.25 |
| 0.1 | Laplace | 200 | 185 (92.5%) | 185 | 4.63 |
| 0.1 | Studentized | 200 | 200 (100.0%) | 200 | 8.79 |
| 0.1 | Profile | 200 | 187 (93.5%) | 187 | 4.66 |

The ordinary profile's 95% label is not calibrated in the base and high-SD
cells. Studentized curvature improves correct-return coverage in these cells,
but its intervals are substantially wider and are conservative in the low-SD
cell. The exact joint-MAP fit did not reproduce the marginal model's known
estimated-zero variance failures; those remain unresolved. No production
method or default was changed based on these comparisons.

A follow-up stratification of these same families points to rate-SD estimation
as a contributor. In the base cell, profile covered 97/129 families whose
fitted SD was below the generating 0.3, versus 70/71 at or above it. In the
high-SD cell the corresponding counts were 96/133 and 66/67; in the low-SD
cell they were 90/103 and 97/97. Profile misses occurred on both sides of the
true age (base: 15 below and 18 above; high-SD: 17 below and 21 above).
This stratification uses the generating SD, which is unavailable in real data,
and was chosen after viewing results. It diagnoses these cases; it does not
calibrate a new interval or justify an SD-based selection rule.

Each case directory retains per-family rows, source and input hashes, summary,
and compressed original inputs. An independent check of all 1,800 rows verified
the 200 complete families per cell, identical point fits across methods, all
input hashes, truth/interval membership, and aggregate counts. The generator
uses dense rate covariance and matrix-exponential sequence transitions rather
than NWKIT's simulation and pruning code. The local Python environment and
command arguments are in each `metadata.json`.
The exact tested Python sources are retained in `source.tar.gz`, including the
RADTE runner and independent generator. Recheck the source, inputs, point-fit
pairing, interval indicators, and summaries with:

```sh
PYTHONPATH=. python tools/verify_radte_profile_comparison.py \
  examples/radte/profile-comparison-20260928
```

Reproduce into fresh output directories with `tools/validate_radte_intervals.py`,
using `--generator independent --inference joint-map --likelihood exact`,
`--methods laplace studentized profile`, `--families 200`,
`--study-role validation --protocol examples/radte/profile-comparison-20260928/protocol.json`,
and the SD/seed pairs in the protocol. The validation runner now accepts
`profile` for paired comparisons; this is a research-runner addition, not a
new CLI uncertainty mode.

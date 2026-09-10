# Calibration decision — investigation closed, 2026-09-10

The final confidence-set calibration attempt did not improve detection and is
not adopted. The production finite-grid envelope remains in place. This bounded
investigation is closed: no post-result tuning, extra simulation sampling or
replacement calibration is planned within this task. Backend acceptance guards,
explicit boundary metadata, independent audits and negative evidence are
integrated in the main checkout.

## Prespecified experiment and result

The [protocol](reviews/shift-calibration-final-protocol.md) was frozen before
100 independent development blocks (seed 20260930). Each block contains four
paired observations: no shift, one shift, two distinct shifts, and two shared
shifts. Conditions share generating innovations; these are **100 independent
blocks, not 400 independent replicates**. The balanced eight-tip, no-error,
fixed-root generator uses alpha × height = 2.1, sigma² × height = 0.25 and
optimum effect 2. All searches use the unchanged 27-point alpha grid, two-shift
candidate limit and B=999 Monte Carlo draws per alpha. A separate two-block
pilot was for execution timing only.

The candidate restricts the nuisance grid using a null-only likelihood confidence
set, then adds its exclusion allowance beta to the maximum retained shift-test
P value. Beta was fixed at 0.005 and 0.01 before the experiment. The critical
screen required at least a five-point net detection gain for distinct shifts,
no greater than five-point loss for either other alternative, no failures and
null rejection no greater than 7.5%. This was a development screen; a pass would
have required fresh independent validation before adoption.

| Generating condition | Current envelope / 100 | Beta 0.005 / 100 | Beta 0.01 / 100 |
|---|---:|---:|---:|
| No shift (false selections) | 5 | 4 | 3 |
| One shift | 55 | 53 | 53 |
| Two distinct shifts | 58 | 52 | 51 |
| Two shared shifts | 47 | 45 | 41 |

All 100 blocks completed with zero failures. Neither candidate produced a newly
detected observation in any condition. Distinct-shift paired losses were six
and seven, respectively. The one-sided exact 95% upper bound on **gross paired
gain**, Bonferroni-adjusted over the two prespecified beta candidates, is 3.62%
for each in that critical condition. Since net gain cannot exceed gross gain,
this excludes a five-point improvement there under the tested design. It does
not establish impossibility for other methods or generating conditions; the
bound is not jointly adjusted across all four conditions.

Confidence sets remained broad (mean retained points for distinct shifts:
19.72 and 15.25 of 27); restriction did not yield gains sufficient to offset the
beta correction. Neither candidate passed the screen. The shared generating
condition measures any-shift detection, not a test establishing convergence.

## Mathematical and empirical limits

For a true alpha on the declared grid, let C contain points whose null-only
Monte Carlo P value exceeds beta. The restricted probability is
`min(1, beta + max(p_shift[a] for a in C))`, with empty-set maximum zero.
Rejection implies either exclusion of the true alpha (probability at most beta)
or its shift-test P value at most level minus beta. A union bound therefore
controls the finite-grid no-error null, assuming each pointwise Monte Carlo test
is valid. Reusing simulation draws does not require independence between these
P values. This adapts the confidence-set nuisance construction of
[Berger and Boos (1994)](https://doi.org/10.1080/01621459.1994.10476836);
see also [Dufour's Monte Carlo test analysis](https://jeanmariedufour.github.io/Dufour_1995_MCT_W.pdf).
The OU-specific confidence statistic is implemented and checked here.

This is no theorem for continuous alpha, measurement errors or later-stage
mean selection. Nor does the experiment validate those settings. The earlier
[response validation](SHIFT_RESPONSE_VALIDATION.md) passed its restricted
12,000-dataset null gate (largest simultaneous upper bound 6.23%, gate 7.5%),
but its five-point power-noninferiority criterion was not established. Those
limitations remain. The present negative experiment closes the improvement
attempt rather than turning either study into complete scientific acceptance.

## Audit and reproduction

[Development evidence](examples/shift/confidence-development/summary.json)
retains all observations, truths, seeds, pointwise exceedance counts, confidence
sets and source hashes. The [full replay audit](examples/shift/confidence-development/audit.json)
regenerated all inputs, independently recomputed the selected and 27 null
likelihoods, reconstructed decisions and summaries, and replayed all 100 complete
simulation banks exactly. Maximum independent likelihood error was 2.67e-14.

```sh
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
python tools/try_shift_confidence_calibration.py --phase development \
  --workers 4 --output /tmp/shift-confidence-new
python tools/verify_shift_confidence_calibration.py /tmp/shift-confidence-new --workers 4
python tools/verify_shift_null_contract.py examples/shift/null-contract-validation \
  --frozen-engine --replay-stride 1000
```

Use a new output directory. The prototype is research tooling and has no CLI
selection mode. The archived null-contract audits explicitly load their
hash-verified historical engine; they do not claim to validate the current CLI.
The final confidence experiment uses the current main-checkout fitting engine.
Existing mixed-zero-error guards and separately developed ASR work are preserved.

The integration suite passed **3752 tests**, with 27 tests skipped,
and **85% branch coverage** (required minimum 80%). Ruff, mypy, dependency
consistency, Bandit, dependency vulnerability scanning and complexity limits
passed. A subsequent change preventing archive-loader bytecode writes passed
all 49 focused tests, including original/corrected R backends and an input-file
fingerprint assertion; production package code is unchanged after the full run.
These checks used an isolated copy of the base commit plus this response,
excluding other work in progress. Other Python versions and CI operating systems
were not run.

[Integration check records](examples/shift/integration-audits/software-checks.json)
retain exact outcomes and logs. The final distribution outcome is recorded in
`reviews/shift-integration-distribution.json` in the repository, outside the
source archive. No GitHub push or release is included.

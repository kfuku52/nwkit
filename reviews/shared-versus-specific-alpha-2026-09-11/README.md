# Estimated shared versus trait-specific alpha: results

See [PROTOCOL.md](PROTOCOL.md) for the frozen model, budgets, seeds and limitations. Both fitted models estimate alpha and full evolutionary covariance. This is a Docker experiment on 100 tips, not a SIF validation.

## Fixed-layout timing

Serial warmup plus three repeated fits on the true shifted layout. Times exclude imports and simulation. A timeout is a failed attempt to obtain a complete timing series; it is not a measured median.

| Traits | Generating alpha | Shared, s | Trait-specific, s | Ratio | LL gain, specific − shared |
|---:|---|---:|---:|---:|---:|
| 2 | shared | 0.0630 | 6.2451 | 99.1× | 0.398520 |
| 2 | different | 0.0569 | 21.7070 | 381.2× | 7.947779 |
| 5 | shared | 0.0712 | >120 (timeout) | — | — |
| 5 | different | 0.0630 | >120 (timeout) | — | — |
| 10 | shared | 0.0765 | >120 (timeout) | — | — |
| 10 | different | 0.0926 | >120 (timeout) | — | — |

`shared` generating alpha means all ones; `different` means geometric spacing from 0.25 to 4. Likelihood gain is not a speed-correctness check between identical models: the models differ. A materially negative gain instead warns that numerical optimization failed to recover the nested shared submodel. Alpha estimates and optimizer diagnostics are retained in JSON.

## Approximate-search pilot

Ten paired datasets per cell; full covariance, four candidate branches and five refits. AIC/BIC selection is not a calibrated significance test. Reported detection fractions are conditional on completed runs; failures/timeouts remain visible. These small cells cannot establish general detection performance.

| True alpha | Scenario | Criterion | Model | Complete / planned | Any shift | Exact branch | False branch | Failed / timeout |
|---|---|---|---|---:|---:|---:|---:|---:|
| shared | null | AIC | shared | 10/10 | 10/10 | 0/10 | 10/10 | 0/0 |
| shared | null | AIC | trait-specific | 10/10 | 10/10 | 0/10 | 10/10 | 0/0 |
| shared | null | BIC | shared | 10/10 | 1/10 | 0/10 | 1/10 | 0/0 |
| shared | null | BIC | trait-specific | 10/10 | 1/10 | 0/10 | 1/10 | 0/0 |
| shared | aligned | AIC | shared | 10/10 | 10/10 | 1/10 | 9/10 | 0/0 |
| shared | aligned | AIC | trait-specific | 10/10 | 10/10 | 0/10 | 10/10 | 0/0 |
| shared | aligned | BIC | shared | 10/10 | 3/10 | 1/10 | 2/10 | 0/0 |
| shared | aligned | BIC | trait-specific | 10/10 | 2/10 | 0/10 | 2/10 | 0/0 |
| shared | opposed | AIC | shared | 10/10 | 10/10 | 10/10 | 0/10 | 0/0 |
| shared | opposed | AIC | trait-specific | 9/10 | 9/9 | 9/9 | 0/9 | 1/0 |
| shared | opposed | BIC | shared | 10/10 | 10/10 | 10/10 | 0/10 | 0/0 |
| shared | opposed | BIC | trait-specific | 9/10 | 9/9 | 9/9 | 0/9 | 1/0 |
| different | null | AIC | shared | 10/10 | 10/10 | 0/10 | 10/10 | 0/0 |
| different | null | AIC | trait-specific | 10/10 | 10/10 | 0/10 | 10/10 | 0/0 |
| different | null | BIC | shared | 10/10 | 4/10 | 0/10 | 4/10 | 0/0 |
| different | null | BIC | trait-specific | 10/10 | 2/10 | 0/10 | 2/10 | 0/0 |
| different | aligned | AIC | shared | 10/10 | 10/10 | 1/10 | 9/10 | 0/0 |
| different | aligned | AIC | trait-specific | 10/10 | 10/10 | 2/10 | 8/10 | 0/0 |
| different | aligned | BIC | shared | 10/10 | 2/10 | 1/10 | 1/10 | 0/0 |
| different | aligned | BIC | trait-specific | 10/10 | 2/10 | 2/10 | 0/10 | 0/0 |
| different | opposed | AIC | shared | 10/10 | 10/10 | 6/10 | 4/10 | 0/0 |
| different | opposed | AIC | trait-specific | 10/10 | 10/10 | 10/10 | 0/10 | 0/0 |
| different | opposed | BIC | shared | 10/10 | 7/10 | 6/10 | 1/10 | 0/0 |
| different | opposed | BIC | trait-specific | 10/10 | 10/10 | 10/10 | 0/10 | 0/0 |

Completed paired searches: 59/60. The per-job JSON records the actual candidate sets, nuisance estimates and numerical-mode diagnostics. Conditional success fractions must not be compared without the failure counts.

## Reproduction

Use the pinned local image recorded in the protocol, one BLAS/OpenMP thread, and run timing before parallel search jobs:

```sh
python tools/benchmark_shift_alpha_models.py --part timing --traits 2,5,10 --repeats 3 --timeout 120 --output results/timing
python tools/benchmark_shift_alpha_models.py --part pilot --traits 2 --replicates 10 --workers 4 --timeout 180 --output results/pilot
python tools/summarize_shift_alpha_models.py results
```


## Observation-error timing supplement

Two traits, 100 tips, true shifted layout. Known sampling SE 0.1; generating extra measurement variance 0.04. Both models estimate the extra diagonal measurement variances as well as alpha/full process covariance. Serial warmup plus three fits, 120-second limit per fit. Both models use general vector pruning here.

| Generating alpha | Shared, s | Trait-specific, s | Ratio | All evaluated modes succeeded |
|---|---:|---:|---:|---|
| shared | 30.4590 | 40.6511 | 1.33× | True / True |
| different | 14.7412 | 55.9770 | 3.80× | True / True |

This supplement does not test detection power with observation errors. Run `python tools/benchmark_shift_alpha_errors.py --output results/timing_errors` only after parallel jobs have finished.

# Direct comparison of estimated shared and trait-specific alpha

**In this implementation, trait-specific alpha costs more computation, but it does
not uniformly reduce shift detection.** With unequal generating rates, the small
pilot found cases where allowing separate rates improved exact-branch detection
and reduced null detections. The large computational advantage of shared alpha
on error-free data becomes much smaller when observation errors require both
models to use general numerical fitting.

Both models estimate alpha from the data and fit **full evolutionary covariance**.
Generating alpha was either `[1, 1]` or `[0.25, 4]` for two traits. The tree had
100 tips and height one. These are synthetic fixed-root OU data; fitting never
receives the generating alpha. See the [prespecified protocol](PROTOCOL.md),
[complete tables](README.md), and [validation checks](validation.json).

## Computation

Fixed true shift layout, one dataset per generating-alpha condition, one warmup
and three serial timed fits. Values are medians in seconds, including parameter
optimization and excluding imports/data generation. These are fixed-layout fit
times, not whole-search or bootstrap times.

| Generating alpha | Observation model | Shared alpha | Trait-specific alpha | Time ratio |
|---|---|---:|---:|---:|
| `[1, 1]` | No observation error | 0.0630 | 6.2451 | 99.1× |
| `[0.25, 4]` | No observation error | 0.0569 | 21.7070 | 381.2× |
| `[1, 1]` | Known + estimated observation error | 30.4590 | 40.6511 | 1.33× |
| `[0.25, 4]` | Known + estimated observation error | 14.7412 | 55.9770 | 3.80× |

The error supplement uses known sampling SE 0.1 for every tip/trait, generating
additional measurement variance 0.04 per trait, and **estimates** the additional
diagonal measurement variances in both fitted models. There are no missing
coordinates. Both models use vector-tree pruning in this supplement, and every
evaluated alpha mode passed its numerical convergence checks. Thus the 99–381×
error-free ratios should not be applied to an analysis that includes sampling or
additional measurement error.

With 5 or 10 traits, all four trait-specific-alpha timing attempts exceeded the
prespecified 120-second limit during their first warmup fit. Shared-alpha fits
completed in 0.063–0.093 seconds. No completed trait-specific median or extrapolated
speed ratio is claimed for those larger cases.

## Shift discovery pilot

Ten datasets in each generating-alpha/scenario cell, paired between the fitted
models. Both models use the same approximate search budget: four candidate
branches and five layout refits, at most one shift. Each model generates its own
covariance-aware candidates; the true branch is never supplied to the search.
The following table uses BIC. AIC results are also retained in the complete tables.

Each cell below reports **shared-alpha fit → trait-specific-alpha fit**, on the
same completed paired datasets. “Correct” requires the exact generating branch.

| Generating alpha | Null: any false shift | Aligned shift: correct branch | Opposed shift: correct branch |
|---|---:|---:|---:|
| `[1, 1]` | 1/10 → 1/10 | 1/10 → 0/10 | 9/9 → 9/9 |
| `[0.25, 4]` | 4/10 → 2/10 | 1/10 → 2/10 | 6/10 → 10/10 |

There were 119 completed model runs out of 120: shared alpha completed 60/60,
trait-specific alpha 59/60. One trait-specific run failed its numerical convergence
check in the common-alpha/opposed-shift cell. It was not redrawn or treated as a
non-detection; that cell's paired comparison uses the remaining nine datasets.
No search reached the 180-second timeout. The failed outcome and optimizer
messages are preserved in
[pilot-p2-shared-opposed-1-trait-specific.json](pilot/pilot-p2-shared-opposed-1-trait-specific.json).

These counts are **pilot observations**, not precise power estimates or proof of
5% false-positive control. BIC/AIC are research selection scores here; this study
does not run the full fitted-null bootstrap. Candidate budgets are small, only
one tree is used, and the search pilot has two traits without observation error.
Within each generating condition, the paired fits see exactly the same data.
Across generating-alpha conditions, diffusion correlation is held at 0.8 and
marginal process tip variance at one; off-diagonal tip correlation consequently
changes with alpha. Broader statistical claims would require more trees,
replicates, effect sizes, errors and a calibrated procedure.

## Verification and reproduction

- All expected outcomes are present: 12 error-free timing jobs, 120 search jobs
  and 4 observation-error timing jobs, including timeouts/failures.
- Paired observed-data hashes match. An independent formula check verifies the
  intended unit process tip variances and diffusion correlation 0.8.
- Repeated successful fits have matching likelihoods within 1e-7. In 246 matched
  dataset/layout comparisons, the more general trait-specific model
  never had a materially smaller likelihood (tolerance 1e-5). The fixed-layout
  timing comparisons with two completed fits also obey this nesting check. These checks are not a proof
  of global numerical optimization; existing mixed-boundary limitations remain.
- The preserved package source archive and benchmark scripts match every saved
  runtime/source manifest. No NWKIT fitting implementation was changed for this
  comparison. Ruff lint/format checks pass for the three new benchmark/report
  scripts. No full package regression rerun is implied by this benchmark work.
- Apple M2 Max, Linux aarch64 GeneGalleon Docker runtime, one BLAS/OpenMP thread.
  Timing runs were serial and separated from the four-process discovery pilot.
  Host/container resources are in [host-environment.json](host-environment.json).
  This does not establish SIF compatibility.

The frozen local image is
`sha256:c8c69075b5d591499bf3195ec912ef5b5e1aa1cd163f03864639ed6a669e9429`.
Raw results and manifests are in `timing/`, `pilot/`, and `timing_errors/`.
`source/` preserves the scripts and package Python source used. The initial
`preflight/` is excluded development evidence, as explained in the protocol.

Inside that runtime, set `OPENBLAS_NUM_THREADS=1`, `OMP_NUM_THREADS=1`, and
`MKL_NUM_THREADS=1`, then run these sequentially:

```sh
python tools/benchmark_shift_alpha_models.py --part timing --traits 2,5,10 --repeats 3 --timeout 120 --output results/timing
python tools/benchmark_shift_alpha_models.py --part pilot --traits 2 --replicates 10 --workers 4 --timeout 180 --output results/pilot
python tools/benchmark_shift_alpha_errors.py --output results/timing_errors
python tools/summarize_shift_alpha_models.py results --require-complete
```

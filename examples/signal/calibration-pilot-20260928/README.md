# Signal lambda calibration pilot

The [study protocol](../../../docs/validation/SIGNAL_CALIBRATION.md) defines
independent Gaussian trait generation and the prespecified 5,000-dataset
confirmation criterion. This is a **200-dataset-per-case pilot**, not that
confirmation. The no-SE source hashes, individual records, and summaries are
under `no-se/`; each case used 199 bootstrap draws and master seed 20260928.

| Case | Conventional χ² rejections / 200 | Bootstrap rejections / 200 | Profile CI covers true λ / 200 |
| --- | ---: | ---: | ---: |
| Balanced 8-tip null | 0 | 9 | 200 |
| Pectinate 8-tip null | 0 | 10 | 200 |
| Balanced 8-tip null, six observed | 0 | 16 | 200 |
| Balanced 32-tip null | 2 | 20 | 198 |
| Balanced 8-tip, true λ = 0.6 | 0 | 19 | 200 |

All 1,000 lambda fits and bootstrap P-values were available. The conventional
χ² approximation is conservative in these cases. The 32-tip bootstrap cell
has 20 rejections, twice the nominal expectation; its Wilson 95% interval
does not establish 5% calibration. This is a reason to replicate with a new
seed and inspect the computation, not to tune a cutoff from this pilot.
The alternative case shows weak detection at this effect size. The profile
confidence interval still uses the χ² cutoff, so its null coverage here does
not validate it over the full λ range.

An independent-seed [32-tip follow-up](32-tip-independent/summary.json) with
1,000 new null datasets and 199 draws per dataset rejected 43/1,000 (4.3%)
by bootstrap, with a Wilson 95% interval of 3.21–5.74%. The conventional
χ² test rejected 6/1,000 (0.6%). All fits and P-values were available.
This reduces concern that the first 20/200 bootstrap result represents a
persistent 10% false-positive rate. It still falls short of the prespecified
5,000-dataset confirmation and its 4–6% Wilson containment criterion.

The no-SE likelihood ratio is mathematically invariant to a nonzero linear
rescaling and translation of the observations; a separate numerical test
checks that property on these tree shapes. This supports the finite-simulation
null reference under the stated Gaussian model but cannot rule out numerical
errors or model misspecification. Known-SE results and larger confirmation
samples must be evaluated separately.

The separate [known-SE coarse pilot](known-se-coarse/summary.json) kept the
same heterogeneous eight-tip SE vector and master seed, with 200 datasets
and **19** bootstrap draws. All 200 fits and P-values were available;
the χ² approximation rejected 0/200 and bootstrap rejected 9/200 (4.5%).
Because the Monte Carlo P-value changes in steps of 0.05 at B=19, this can
screen for gross issues but cannot establish the B=199 calibration criterion.
It took 499.7 s with eight workers even after the zero-LR shortcut; 191/200
observed fits had λ at a boundary. A 5,000-dataset known-SE study at B=199
requires further optimization or a larger compute budget.

The production bootstrap now returns the exact P-value 1 when the observed
likelihood ratio is exactly zero, without refitting simulated traits. On this
pilot, all 1,000 individual JSONL records were byte-for-byte identical before
and after that change. A [matched single-case timing](zero-lr-benchmark.json)
with known SEs and nine bootstrap draws retained P=1 and reduced the median
of three measured runs from 23.91 s to 0.0081 s (one warmup excluded).
This is a zero-LR microbenchmark, not a speedup estimate for positive-LR
traits or the entire command.

`source.tar.gz` retains the three source files named by the protocols. Audit
any case directory from the repository root with:

```sh
PYTHONPATH=. python tools/verify_signal_calibration.py \
  examples/signal/calibration-pilot-20260928/no-se \
  --source-archive examples/signal/calibration-pilot-20260928/source.tar.gz
```

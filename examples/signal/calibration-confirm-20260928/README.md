# Independent Pagel-lambda null confirmation

This run used 5,000 **new** independently generated Gaussian datasets per
case, master seed 20261101, 199 bootstrap draws per dataset, and the frozen
source in [source.tar.gz](source.tar.gz). It covers balanced and pectinate
eight-tip trees, a fixed six-observed-tip missing pattern, and a balanced
32-tip tree, all with exact observations and true lambda zero. The source
hashes, exact settings, 20,000 individual records, and aggregate counts are
in `protocol.json`, `records.jsonl`, and [summary.json](summary.json).

| Null condition | χ² rejections / 5,000 | Bootstrap rejections / 5,000 | Bootstrap Wilson 95% interval | Prespecified 4–6% containment |
| --- | ---: | ---: | ---: | --- |
| Balanced 8 tips | 0 | 251 (5.02%) | 4.45–5.66% | Met |
| Pectinate 8 tips | 0 | 274 (5.48%) | 4.88–6.15% | Not met |
| Balanced 8 tips, six observed | 0 | 236 (4.72%) | 4.17–5.34% | Met |
| Balanced 32 tips | 38 (0.76%) | 221 (4.42%) | 3.88–5.03% | Not met |

All 20,000 observed fits and both P-values were available. The conventional
χ² tail is strongly conservative in these finite-tree cases. Bootstrap
rejection is much closer to the nominal 5%, and an independent 1,000-case
32-tip pilot also rejected 4.3%. However, **only two of four** prespecified
Wilson-containment gates passed; the pectinate upper bound and 32-tip lower
bound overlap the acceptance boundary. No cutoff was tuned, no extra families
were appended to a failing cell, and the default `--lambda-test chi2` was not
changed. These data do not validate known SEs, alternative power, profile
intervals, misspecified trees, or estimated measurement error.

With zero SEs, the LR is mathematically invariant to trait location and scale,
so the fitted-null draws give a finite-simulation null reference under the
specified Gaussian tree model. The empirical gate above is deliberately
stricter as an implementation check. The zero-LR early return was active in
this run; it returns P=1 exactly and never changes a nonzero-LR bootstrap.
The run took 1,418 s with eight workers and BLAS threads fixed to one.

Verify all records, Monte Carlo grids, Wilson intervals, and source hashes:

```sh
PYTHONPATH=. python tools/verify_signal_calibration.py \
  examples/signal/calibration-confirm-20260928 \
  --source-archive examples/signal/calibration-confirm-20260928/source.tar.gz
```

# Pagel-lambda profile interval, interior-value confirmation

This independent study used 5,000 Gaussian datasets with true lambda 0.6 on
the balanced eight-tip covariance, master seed 20261115, and the existing
95% chi-square profile interval. `--inner 0 --profile-ci` skipped bootstrap
P-values so this run measured interval coverage directly. The source archive,
exact protocol, per-dataset records, and summary are retained here.

All 5,000 fits and intervals were available. **All 5,000 intervals contained
the true lambda 0.6** (Wilson 95% Monte Carlo interval 99.92–100%). The
prespecified coverage criterion required the Wilson interval to lie within
93–97%, so it failed decisively. The maximum was at a lambda boundary in
4,429/5,000 fits; the conventional chi-square test rejected lambda=0 in
0/5,000. This simulated tree and effect have weak signal for lambda, and the
current profile limits are highly conservative in this condition. Coverage
at lambda=0, or a broad interval containing the truth, does not establish
calibration across the parameter range.

No profile cutoff or default test was changed to fit this result. Different
tree shapes, stronger covariance, known SEs, and model misspecification need
separate evaluation before any general interval claim. The run took 8.4 s
with eight workers and no bootstrap refits.

Audit the 5,000 records, summary, and archived source:

```sh
PYTHONPATH=. python tools/verify_signal_calibration.py \
  examples/signal/profile-ci-confirm-20260928 \
  --source-archive examples/signal/profile-ci-confirm-20260928/source.tar.gz
```

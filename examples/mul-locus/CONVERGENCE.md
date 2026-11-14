# Fixed-Data Monte Carlo Budget Follow-Up

Predeclared on 2026-10-03, before either additional run. This diagnostic does
not overwrite the conditional-reference probe or change the adoption decision.

Control: `/tmp/nwkit-locus-conditional-probe-v2-20261003`, 2,000 raw draws per
bank. Follow-ups use **4,000 and 8,000** draws, all three scenarios and all seven
truth cases, ten families per dataset and one dataset per case. Biological
parameters, exact hidden-tip limit four, work cap 100,000, Chernoff-KL intervals,
conditional independent reference, and 19 null replicates at each of four null
grid points are unchanged. All three methods and all failures are retained.

Data seed 20261111, bank seed base 20261112, calibration seed base 20261113
are fixed. The paired-bank per-history namespaces produce nested prefixes
when stratum budgets increase. Observed datasets must match the stored control
exactly before a fixed-data comparison is reported. Null sampling uses the
same generating grid and seed namespaces, with a complete search at each budget.

Both directories must be new. Freeze the current scientific source snapshot
before the first run. Execute the existing `integration.py` separately at each
budget, with the same arguments as CONDITIONAL.md except `--samples` and
`--output`. Preserve its nonzero exit when histogram support failures remain.

Report completed/planned analyses, failures, parent point rankings, calibrated
point and score-MC P-values, and contrast interval widths (including infinite
bounds). Verify all observed input hashes and independently recompute all
completed null-calibration endpoints. Do not interpret two additional budgets
or one dataset per truth case as proof of convergence, empirical power, a 5%
error guarantee, or validation of biological node-origin probabilities.

# Conditional Reference And KL Interval Probe

Frozen before execution on 2026-10-03. This follow-up does not overwrite the
earlier integration study or change its adoption decision.

The independent topology-only reference retains the global SSA locus process,
but samples conditional branch-end lineage counts using independently built
dense pure-death CTMC exponentials and uniform pair mergers. It imports no
production coalescent transitions/messages/sampling or detection functions.
Full hidden loci, nested daughter conditions and the observed-size selection
are retained. Work/underflow failures abort; there is no history retry or
switch to a different sampler after failure. The old msprime rejection and
sequence/NJ reference remains available unchanged. The new reference does
not provide coalescent times or calibrate inferred gene trees.

Weighted contribution intervals invert `n*kl(mean,p) <= log(2/alpha)`.
This is Theorem 3 of [Foong, Bruinsma and Burt (2022)](https://arxiv.org/pdf/2205.07880),
valid for IID observations in [0,1], not fractional Clopper-Pearson counts.
Each endpoint uses an outward numerical bracket. The existing full-universe,
bank, stratum and selection-denominator allocation is preserved.

Predeclared comparison: all three historical scenarios, all seven truth
cases, one dataset per case, ten families, 2,000 raw draws per bank, 19 null
replicates per grid point; hidden exact-integration limit four and work cap
100,000. Seeds are new and fixed. Histogram, detection RB and hybrid RB use
paired banks/data/null draws. Failures remain in all planned denominators.
Small-history full-forest oracles and rare daughter cases precede execution.
This remains a functional/support probe, not a power or type-I-error study.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
  "$PY" examples/mul-locus/integration.py \
  --output /tmp/nwkit-locus-conditional-probe-20261003 \
  --samples 2000 --families 10 --bootstrap 19 --replicates 1 \
  --seed 20261111 --bank-seed 20261112 --calibration-seed 20261113 \
  --reference conditional --interval chernoff-kl
```

The protocol JSON freezes source hashes before any simulation. Historical
results are controls, not resimulated or relabeled with the new reference.
Only opt-in research integration into GeneGalleon is authorized by this
probe; D+L remains the default and no automatic WGD call is introduced.

The first execution failed on tiny `expm` roundoff in theoretically impossible
lineage-count increases. Its complete 63 failure records and source snapshot
are preserved. A second fixed execution uses the same budgets and seeds in
`/tmp/nwkit-locus-conditional-probe-v2-20261003`, restoring only exact
pure-death structural zeros. Negative feasible transitions still fail. This
is an implementation repair, not a change to the biological model or budget.

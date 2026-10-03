# Conditional MSC fitting pilot

This is a small parent-candidate experiment, not WGD detection or a release
gate. Read [the model guide](../../docs/guides/MUL_MSC.md) before interpreting it.
No known Ne or attachment age is passed to fitting. Species topology/ages,
the polyploid clade, candidate set and broad search bounds are supplied.

## Protocol

All cases use `(((A:2,X:2):1,B:3):2,C:5);` in generations and direct
attachment age 0.7. Standard cases cross H2 A/B/C with diploid/subgenome Ne
0.5 and 5. Extra H2 B/Ne 5 cases introduce 25% random homoeolog missingness,
an erroneous A/X divergence age of 2.4, a ghost-parent proxy with donor
separation 1.3 instead of hybridization 0.7, or omit true H2 B from the search.
The ghost proxy holds constant Ne and no founder population jump: it changes
the observable donor separation, not a separately identifiable hybridization
time. It is not a full ghost-network simulator.

Each family is independently simulated with msprime's standard coalescent on
the unfolded population tree. Global ploidy 2 sets the diploid time scale;
each `SampleSet` uses ploidy 1 to draw exactly one genome per diploid/subgenome.
Homoeolog labels are independently flipped. JC69 mutations generate 600-site
alignments at rate 0.003/site/generation, including invariant sites.
Biopython NJ uses JC69-corrected distances and midpoint rooting, without
consulting the true gene root. This is an estimated-tree pilot, not an IQ-TREE
ML or sequence-marginal likelihood comparison. Saturated distances or failed
root/binary validation abort explicitly, not by dropping families.

True and estimated trees from each dataset are separately passed to both
methods. Both use the same families and candidate list. D+L ignores lengths;
fitted MSC receives fixed species ages and bounds `(0.1,1.9)` for attachment
age and `(0.1,10)` for Ne, never the true Ne/age. Missing copies are conditioned
on by MSC. Wrong ages and omitted candidates deliberately violate or restrict
its assumptions. The ghost proxy probes observable donor separation, not a
separate ghost-network model or a model-adequacy test.

Report parent recovery, numerical ties, withheld parameter estimates,
parameter errors when reported, and gene topology error, retaining each
replicate rather than just pooled accuracy. A unique maximum is not calibrated
parent support. The omitted-parent case must not be interpreted as evidence
that a wrong candidate is the true parent.

The seed is 20261018, distinct from analytic/training fixtures. Defaults are
two replicates of 100 families per case: a bounded pilot with very imprecise
recovery-rate estimates. No finite-sample accuracy threshold is tuned into CI.
Do not adjust the protocol after seeing results to obtain a desired comparison.

## Run

Use an isolated research environment with NWKIT plus `msprime` and `biopython`;
the latter two are not new NWKIT runtime dependencies. From the repository root:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
  python examples/mul-msc/conditional_pilot.py --output /tmp/new-msc-pilot
```

The output directory must not exist. Outputs include inputs, complete replicate
scores/fits, raw alignments, family audits, protocol/environment/source hashes,
and summary TSV. Failures are retained in the output with nonzero exit status;
datasets/replicates are never silently discarded. Repeated calls use fresh
directories, preserving previous evidence. Larger independent-seed studies,
ML gene trees, linkage, SSD/loss, Ne heterogeneity and external empirical data
remain later validation work.

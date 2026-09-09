# Synthetic species-age uncertainty example

This five-gene, three-species example has two paralogous speciation nodes sharing
species event AB at age 10, a gene duplication at age 15, and a species root at
age 20. It illustrates behavior; it is not an empirical confidence assessment.

The 20 supplied chronograms scale both species ages together by factors evenly
spaced from 0.9 to 1.1, preserving their dependence. The external 95% percentile
intervals are the exact empirical quantiles of these input samples. The separate
hard-bound file allows AB to vary from 8 to 12 and keeps root R fixed at 20 in
the bounded analysis. No probability prior is inferred from those bounds.

From the NWKIT repository root, in an environment with NWKIT's dependencies:

```sh
bash examples/radte/species-uncertainty/run.sh
```

An optional first argument selects a different output directory. Outputs default
to `output/pdf/species-uncertainty/`:

* `fixed.pdf`, `bounded.pdf`, `ensemble.pdf`: gene/species trees with distinct
  external and fitted intervals.
* `comparison.pdf`: three columns of gene/species trees for fixed, bounded and
  ensemble ages, using one common age scale. Page two aligns input-refit SD and
  mean within-chronogram conditional interval width with the reference tree nodes.
* Each run's numerical tables and manifest, plus comparison TSVs and manifest.

The example requests conditional profile intervals for every ensemble sample.
Inspect statuses and diagnostics: hard-bound-limited intervals and local-optimum
warnings are deliberately retained. The comparison is not a combined posterior
or a guarantee of coverage. For schemas and bootstrap options, see
[RADTE species uncertainty](../../../RADTE_SPECIES_UNCERTAINTY.md).

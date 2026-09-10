# Individual-level multivariate BM example

The synthetic `individuals.tsv` contains 28 individuals from eight species,
with two traits (`x`, `y`), unequal replication (three/four individuals), and
one missing `y` value for individual `C/2`. It was generated with NumPy's RNG
seed 197, BM Sigma `[[1.1,0.45],[0.45,0.8]]`, common individual W
`[[0.5,0.18],[0.18,0.3]]`, and a root mean `[4,7]`.

Run from the repository root:

```bash
nwkit asr --infile examples/individual_asr/tree.nwk \
  --trait examples/individual_asr/individuals.tsv --state-column x,y \
  --model MV-BM --within-species-covariance full --covariance-method REML \
  --target all --outfile /tmp/species-and-ancestors.tsv \
  --individual-out /tmp/individual-predictions.tsv --model-out /tmp/joint-model.tsv
```

The REML log likelihood is approximately `-68.3308836542`. Estimated Sigma is
`[[1.510481,0.883442],[0.883442,0.613049]]` and W is
`[[0.583255,0.226964],[0.226964,0.392897]]`. A single small dataset need not
recover the simulation parameters exactly.

The primary table reports latent species means and ancestors. The individual
table predicts `C/2/y` using the measured `C/2/x` coordinate as well as other
individuals and the tree. All intervals condition on estimated Sigma and W;
they do not include covariance-parameter uncertainty.

`Rscript examples/individual_asr/reference.R` reproduces the independent
Rphylopars comparisons when that optional R package is installed.
See [ASR_INDIVIDUALS.md](../../ASR_INDIVIDUALS.md) for assumptions and schemas.

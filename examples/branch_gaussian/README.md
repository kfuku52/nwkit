# Fixed branch-specific Gaussian example

Run from the repository root after installing NWKIT. The tree and parameters
match [mixed_process.py](mixed_process.py): background BM, OU on A, and BM plus
an independent Gaussian end jump on B. The explicit Gaussian root has mean 0.3
and variance 0.9. Parameters are supplied, not estimated.

ASR, with normalized models and complete run settings:

```sh
nwkit asr --model BRANCH-GAUSSIAN -i examples/branch_gaussian/tree.nwk --input-rooted yes \
  --branch-models examples/branch_gaussian/branch_models.tsv \
  --root-prior gaussian --root-mean 0.3 --root-variance 0.9 \
  --trait examples/branch_gaussian/traits.tsv --state-column x \
  --standard-error-column se --output summary -o branch-asr.tsv \
  --process-out branch-model.json --branch-models-out branch-normalized.tsv
```

Equivalent regime tables, evaluated by likelihood only:

```sh
nwkit asr --model BRANCH-GAUSSIAN -i examples/branch_gaussian/tree.nwk --input-rooted yes \
  --branch-regimes examples/branch_gaussian/branch_regimes.tsv \
  --regime-models examples/branch_gaussian/regime_models.tsv \
  --root-prior gaussian --root-mean 0.3 --root-variance 0.9 \
  --trait examples/branch_gaussian/traits.tsv --state-column x \
  --standard-error-column se --output likelihood -o branch-likelihood.tsv
```

The log likelihood is approximately **−5.487990618**. Direct and regime tables
produce the same process. `branch-normalized.tsv` can replace the original
`--branch-models` input without changing results.

Reproducible latent prior simulation:

```sh
nwkit asr --model BRANCH-GAUSSIAN -i examples/branch_gaussian/tree.nwk --input-rooted yes \
  --branch-models examples/branch_gaussian/branch_models.tsv \
  --root-prior gaussian --root-mean 0.3 --root-variance 0.9 \
  --output prior-samples --prior-samples 3 --seed 42 -o branch-simulations.tsv \
  --process-out branch-simulation-model.json
```

This gives 15 rows (3 draws × 5 nodes), including the root, and does not add
measurement noise. `simulation` is one-based; `branch_id` is the original
level-order tree ID. Reordering the children changes those IDs, so assignments
must be regenerated for a reordered tree. Root ID 0 never belongs in a model
or regime assignment TSV.

See [BRANCH_GAUSSIAN.md](../../BRANCH_GAUSSIAN.md) for input validation, explicit
root treatments and output schemas.

The [eight-tip plotting example](plot/README.md) includes posterior histories,
prior distributions, a missing tip and multiple OU optima. The existing files
in `output/` preserve an earlier standalone CLI demonstration; reproduce new
outputs with the integrated ASR commands above.

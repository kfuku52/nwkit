# Phylogenetic signal example

The synthetic eight-tip tree has equal root-to-tip distances of two.
`clustered` differs mainly between the two root clades; `mixed` varies within
sister pairs. SE columns are synthetic known sampling errors, not estimates
from replicates.

From the repository root:

```sh
nwkit signal -i examples/signal/tree.nwk \
  --trait examples/signal/traits.tsv --columns clustered,mixed \
  --n-sim 999 --seed 1 -o signal.tsv
```

Without SEs, K is approximately 2.164516 for `clustered` and 0.476939 for
`mixed`; their constrained lambda estimates are 1 and 0 respectively.
The boundary status is expected. See [the signal guide](../../SIGNAL.md) for
error-aware fits, confidence intervals, test interpretation and output columns.

# Branch-specific scalar Gaussian processes

The Python API in `nwkit.branch_gaussian` builds a single process with fixed
BM, OU and prescribed Gaussian jump parameters assigned to individual branches.
It returns the existing `GaussianTreeProcess`, usable directly for likelihood,
ancestral-state conditioning, covariance and simulation. Parameters and branch
assignments are supplied by the caller; fitting, automatic model search and a
CLI/table importer are outside this API.

## Model and units

Each incoming branch obeys `child = a * parent + b + independent error`, with
error variance `q`. For branch length `t`:

| Diffusion | a | b | q |
| --- | --- | --- | --- |
| `BrownianBranch(variance_rate=sigma2)` | 1 | 0 | sigma2 × t |
| `OUBranch(alpha, variance_rate=sigma2, optimum=theta)` | exp(−alpha × t) | (1−a) × theta | sigma2 × (1−a²) / (2 × alpha) |
| No diffusion (`None`) | 1 | 0 | 0 |

OU with `alpha=0` is exactly BM; small positive alpha uses a stable finite-time
formula. Alpha has inverse tree-time units, sigma2 has trait²/tree-time units,
and theta has trait units. Finite, nonnegative lengths, alpha and variance
parameters are required. Zero variance is allowed, including deterministic
branches. Unrepresentable positive diffusion variances raise an error.

`BranchGaussianModel(diffusion, jump=GaussianJump(mean=m, variance=v))`
adds an independent event **after diffusion, at the branch end**: `b += m`
and `q += v`. Jump variance is per event and is not multiplied by branch length.
A jump on a zero-length branch still occurs. `diffusion=None` gives a pure jump
with identity propagation between parent and child. A model must specify at
least a diffusion or a jump. This differs from the existing CLI `JUMP-BM`, which
integrates latent Poisson jump counts and is not marginally Gaussian.

## Branch assignments and root

Use IDs from `nwkit.util.assign_branch_ids(tree)`: level-order traversal with
root ID 0, matching the existing branch-ID convention. Supply **every non-root
ID exactly once**. Root ID 0 is excluded because the root prior is a separate,
required `GaussianRootPrior`. Integer IDs are required; strings and booleans
are rejected. Pass a mapping, not a list of pairs (which can contain duplicate
IDs). The supplied tree must be a root node; detach an attached subtree first.
Recompute assignments after topology or child-order changes.
Repeated node names do not affect ID-based assignment.

The root can be fixed (zero variance), Gaussian (positive variance), or flat
(`variance=None`). The `stationary` label is also accepted with an explicitly
supplied mean and positive variance; heterogeneous branches have no automatically
chosen stationary root. Root variance and jump variances are absolute trait²
values, independent of each branch's diffusion rate. Unlike the global
`variance_scale` in `build_evolutionary_process`, this builder does not rescale
the root. All innovations and jumps are mutually independent and independent
of the proper root draw.

## Example

Run the complete [mixed-model example](examples/branch_gaussian/mixed_process.py):

```sh
python examples/branch_gaussian/mixed_process.py
```

A minimal assignment is:

```python
from nwkit.branch_gaussian import (
    BranchGaussianModel, BrownianBranch, GaussianJump, OUBranch,
    build_branch_gaussian_process,
)
from nwkit.gaussian_tree import GaussianRootPrior
from nwkit.util import assign_branch_ids

ids = assign_branch_ids(tree)
background = BranchGaussianModel(BrownianBranch(variance_rate=0.8))
models = {branch_id: background for node, branch_id in ids.items() if not node.is_root}
selected_id = next(branch_id for node, branch_id in ids.items() if node.name == "A")
models[selected_id] = BranchGaussianModel(
    OUBranch(alpha=0.6, variance_rate=1.3, optimum=0.9),
    GaussianJump(mean=0.2, variance=0.15),
)
process = build_branch_gaussian_process(
    tree, models, root=GaussianRootPrior("gaussian", mean=0.3, variance=0.9),
)
```

Shared regime parameters can be expressed as a dictionary of named
`BranchGaussianModel` instances, then expanded using a branch-ID-to-regime map;
the runnable example shows this pattern. Model objects are immutable.

Pass `process` to `gaussian_tree_likelihood`, `condition_gaussian_tree`,
`simulate_gaussian_process`, or `sample_gaussian_posterior` in
`nwkit.gaussian_inference`. Likelihood and conditioning accept leaf-name values
and `standard_errors` (standard deviations, not variances); missing values use
the existing inference API's `None` convention. Conditional intervals describe
state uncertainty given the supplied parameters, excluding parameter-selection
uncertainty. Prior simulation returns latent node values, without observation
error; add independent measurement noise separately if needed. Flat-root prior
simulation requires explicit `root_values`, and unconditional covariance requires
a proper root. Exact deterministic constraints follow the existing inference
engine's singular-Gaussian likelihood conventions.

Numerical checks also compare 153 extreme parameter combinations with a
750-digit Decimal oracle, including weak selection whose unscaled alpha × time
product underflows.

Tests compare mixed-process covariance, likelihood and all-node conditioning
with an independently assembled structural-equation oracle, simulation moments,
uniform BM/OU reductions, the alpha→0 limit, and zero-duration end jumps.

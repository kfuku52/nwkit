"""Run from the repository root: python examples/branch_gaussian/mixed_process.py."""

from nwkit.branch_gaussian import (
    BranchGaussianModel,
    BrownianBranch,
    GaussianJump,
    OUBranch,
    build_branch_gaussian_process,
)
from nwkit.gaussian_inference import (
    condition_gaussian_tree,
    gaussian_tree_likelihood,
    simulate_gaussian_process,
)
from nwkit.gaussian_tree import GaussianRootPrior
from nwkit.util import assign_branch_ids, read_tree


def main():
    tree = read_tree(
        "((A:0.3,B:0.7)I:0.4,C:1.1)R;", "1", True, quiet=True, rooted="yes"
    )
    regimes = {
        "background": BranchGaussianModel(BrownianBranch(0.8)),
        "selected": BranchGaussianModel(
            OUBranch(alpha=0.6, variance_rate=1.3, optimum=0.9)
        ),
        "jumped": BranchGaussianModel(
            BrownianBranch(0.8), GaussianJump(mean=0.2, variance=0.15)
        ),
    }
    regime_by_name = {
        "I": "background",
        "A": "selected",
        "B": "jumped",
        "C": "background",
    }
    identifiers = assign_branch_ids(tree)
    models = {
        identifier: regimes[regime_by_name[node.name]]
        for node, identifier in identifiers.items()
        if not node.is_root
    }
    process = build_branch_gaussian_process(
        tree, models, root=GaussianRootPrior("gaussian", mean=0.3, variance=0.9)
    )
    values = {"A": 1.2, "B": -0.4, "C": 2.1}
    errors = {"A": 0.2, "B": 0.1, "C": 0.4}
    likelihood = gaussian_tree_likelihood(process, values, standard_errors=errors)
    posterior = condition_gaussian_tree(process, values, standard_errors=errors)
    print(f"log_likelihood\t{likelihood.log_likelihood:.9f}")
    print("branch_id\tname\tposterior_mean\tposterior_variance")
    for node, identifier in identifiers.items():
        marginal = posterior.marginals[node]
        print(
            f"{identifier}\t{node.name}\t{marginal.mean:.9f}\t{marginal.variance:.9f}"
        )
    samples = simulate_gaussian_process(process, num_samples=3, seed=42)
    print("simulation_shape\t" + str(samples.values.shape))


if __name__ == "__main__":
    main()

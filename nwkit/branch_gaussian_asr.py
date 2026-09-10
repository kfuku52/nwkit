"""Fixed branch processes exposed through the common ASR fit interface."""

from dataclasses import dataclass
from typing import Any

from nwkit.asr_regimes import RegimeAssignment
from nwkit.branch_gaussian import (
    BrownianBranch,
    OUBranch,
    build_branch_gaussian_process,
)
from nwkit.branch_gaussian_input import (
    BranchGaussianAssignment,
    load_branch_gaussian_assignment,
)
from nwkit.branch_gaussian_options import branch_root
from nwkit.gaussian_inference import (
    GaussianConditioningResult,
    GaussianLikelihoodResult,
    condition_gaussian_tree,
    gaussian_tree_likelihood,
)
from nwkit.gaussian_tree import GaussianTreeProcess
from nwkit.util import assign_branch_ids


@dataclass(frozen=True)
class BranchGaussianFit:
    process: GaussianTreeProcess
    branch_assignment: BranchGaussianAssignment
    display_assignment: RegimeAssignment
    theta_by_regime: dict[str, float]
    jump_nodes: tuple[Any, ...]
    log_likelihood: float | None
    num_observed: int
    num_observed_positions: int
    num_effective_observations: int
    fit_status: str = "fixed"
    optimizer_success: bool = True
    summary_kind: str = "posterior"
    prior_root_value: float | None = None

    @property
    def root_prior(self):
        return self.process.root.mode

    @property
    def root_mean(self):
        return None if self.root_prior == "flat" else self.process.root.mean

    @property
    def root_variance(self):
        return self.process.root.variance


def load_branch_assignment(tree, args):
    return load_branch_gaussian_assignment(
        tree,
        branch_models=getattr(args, "branch_models", None),
        branch_regimes=getattr(args, "branch_regimes", None),
        regime_models=getattr(args, "regime_models", None),
    )


def _display_assignment(tree, assignment):
    ids = assign_branch_ids(tree)
    groups = {}
    seen: dict[str, int] = {}
    labels = ["Root prior"]
    by_node, theta, jumps = {tree: "Root prior"}, {}, []
    for node, identifier in ids.items():
        if node.is_root:
            continue
        model = assignment.models_by_branch_id[identifier]
        kind = (
            "BM"
            if isinstance(model.diffusion, BrownianBranch)
            else "OU"
            if isinstance(model.diffusion, OUBranch)
            else "Jump"
        )
        if model.jump is not None and model.diffusion is not None:
            kind += " + jump"
        key = (
            model
            if assignment.regime_by_branch_id is None
            else assignment.regime_by_branch_id[identifier]
        )
        if key not in groups:
            seen[kind] = seen.get(kind, 0) + 1
            name = kind if assignment.regime_by_branch_id is None else f"{key} ({kind})"
            if assignment.regime_by_branch_id is None and seen[kind] > 1:
                name += f" {seen[kind]}"
            groups[key] = name
            labels.append(name)
        label = groups[key]
        by_node[node] = label
        if isinstance(model.diffusion, OUBranch) and model.diffusion.alpha > 0:
            theta[label] = model.diffusion.optimum
        if model.jump is not None:
            jumps.append(node)
    return (
        RegimeAssignment(tuple(labels), by_node, "fixed branch models"),
        theta,
        tuple(jumps),
    )


def branch_fit(
    process, assignment, result=None, *, summary_kind="posterior", prior_root_value=None
):
    display, theta, jumps = _display_assignment(process.tree, assignment)
    return BranchGaussianFit(
        process,
        assignment,
        display,
        theta,
        jumps,
        None if result is None else result.log_likelihood,
        0 if result is None else result.num_observed,
        0 if result is None else result.num_observed_positions,
        0 if result is None else result.likelihood_rank,
        summary_kind=summary_kind,
        prior_root_value=prior_root_value,
    )


def compute_branch_marginals(
    tree,
    observed,
    errors,
    *,
    args=None,
    assignment=None,
    root=None,
    compute_posterior=True,
):
    if assignment is None:
        assignment = load_branch_assignment(tree, args)
    if root is None:
        root = branch_root(args)
    process = build_branch_gaussian_process(
        tree, assignment.models_by_branch_id, root=root
    )
    result: GaussianConditioningResult | GaussianLikelihoodResult
    if compute_posterior:
        result = condition_gaussian_tree(process, observed, standard_errors=errors)
        posterior = result.marginals
    else:
        result = gaussian_tree_likelihood(process, observed, standard_errors=errors)
        posterior = {}
    return posterior, branch_fit(process, assignment, result)


def branch_model_table(fit, args, ci_level):
    import pandas as pd

    return pd.DataFrame(
        [
            {
                "trait_type": "continuous",
                "trait_type_requested": getattr(args, "trait_type", "auto"),
                "trait": args.state_column,
                "model": "BRANCH-GAUSSIAN",
                "root_prior": fit.root_prior,
                "root_mean": fit.root_mean,
                "root_variance": fit.root_variance,
                "estimation_method": "fixed",
                "num_parameters_estimated": 0,
                "log_likelihood": fit.log_likelihood,
                "likelihood_kind": "flat_root_integrated"
                if fit.root_prior == "flat"
                else "proper_root_ml",
                "num_observed": fit.num_observed,
                "num_observed_positions": fit.num_observed_positions,
                "likelihood_rank": fit.num_effective_observations,
                "fit_status": fit.fit_status,
                "ci_level": ci_level,
                "interval_kind": "prior"
                if fit.summary_kind == "prior"
                else "conditional_on_parameters",
                "parameter_uncertainty_included": False,
                "tree_uncertainty_included": False,
            }
        ]
    )

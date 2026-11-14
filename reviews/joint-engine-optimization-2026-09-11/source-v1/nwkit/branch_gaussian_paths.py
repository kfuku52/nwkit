"""Exact finite-grid histories with one explicit end jump per original branch."""

import numpy as np
from ete4 import Tree

from nwkit.asr_paths import (
    _MAX_GRID_NODES,
    _MAX_PATH_VALUES,
    BranchPath,
    SimulationPaths,
)
from nwkit.branch_gaussian import BranchGaussianModel
from nwkit.gaussian_inference import (
    sample_gaussian_posterior,
    simulate_gaussian_process,
)
from nwkit.gaussian_tree import GaussianTransition, GaussianTreeProcess
from nwkit.rooting_state import copy_rooting_info
from nwkit.util import assign_branch_ids


def refine_branch_process(tree, fit, counts):
    refined = Tree()
    refined.name = str(tree.name or "")
    copy_rooting_info(tree, refined)
    endpoints, chains, transitions = {tree: refined}, {}, {}
    ids = assign_branch_ids(tree)
    for node in tree.traverse("preorder"):
        if node.is_root:
            continue
        model = fit.branch_assignment.models_by_branch_id[ids[node]]
        cursor = endpoints[node.up]
        chain = [cursor]
        grid = np.linspace(0.0, float(node.dist), counts[node] + 1)
        lengths = np.diff(grid)
        if float(node.dist) > 0 and np.any(lengths <= 0):
            raise ValueError(
                "Path simulation grid underflowed; rescale branch lengths."
            )
        for length in lengths:
            cursor = cursor.add_child(dist=float(length))
            transitions[cursor] = (
                GaussianTransition(1.0, 0.0, 0.0)
                if model.diffusion is None
                else BranchGaussianModel(model.diffusion).transition(float(length))
            )
            chain.append(cursor)
        if model.jump is not None:
            cursor = cursor.add_child(dist=0.0)
            transitions[cursor] = GaussianTransition(
                1.0, model.jump.mean, model.jump.variance
            )
            chain.append(cursor)
            grid = np.append(grid, float(node.dist))
        cursor.name = str(node.name or "")
        endpoints[node] = cursor
        chains[node] = (chain, grid)
    return GaussianTreeProcess(
        refined, transitions, fit.process.root, "branch-gaussian-grid"
    ), chains


def simulate_branch_paths(
    tree, observed, errors, posterior, fit, counts, count, steps, mode, seed
):
    if tree is not fit.process.tree:
        raise ValueError("Branch paths require the fitted process's original tree.")
    nodes_count = 1 + sum(counts.values()) + len(fit.jump_nodes)
    if nodes_count > _MAX_GRID_NODES or nodes_count * count > _MAX_PATH_VALUES:
        raise ValueError(
            "Figure simulation grid including end jumps is too large; reduce steps or simulations."
        )
    process, chains = refine_branch_process(tree, fit, counts)
    if mode == "conditional":
        if fit.summary_kind == "prior":
            raise ValueError("Prior-only output cannot generate conditional histories.")
        samples = sample_gaussian_posterior(
            process, observed, standard_errors=errors, num_samples=count, seed=seed
        )
        description = (
            "Conditioned on observed tips, including their measurement errors."
        )
    else:
        roots = None
        if process.root.mode == "flat":
            roots = (
                fit.prior_root_value
                if fit.summary_kind == "prior"
                else posterior[tree].mean
            )
        samples = simulate_gaussian_process(
            process, num_samples=count, seed=seed, root_values=roots
        )
        description = (
            "Root fixed at the supplied starting value (flat root prior)."
            if roots is not None and fit.summary_kind == "prior"
            else "Root fixed at inferred mean (flat root prior)."
            if roots is not None
            else "Root fixed at the model's specified value."
            if process.root.mode == "fixed"
            else f"Root drawn from the {process.root.mode} prior."
        )
    values = samples.values[:, :, None]
    if not np.isfinite(values).all():
        raise ValueError(
            "A simulated path exceeds floating-point range; rescale trait units."
        )
    indices = {node: index for index, node in enumerate(samples.nodes)}
    branches = {
        node: BranchPath(grid, values[:, [indices[point] for point in chain], :])
        for node, (chain, grid) in chains.items()
    }
    return SimulationPaths(
        branches, values[:, indices[process.tree], :], mode, description, steps
    )

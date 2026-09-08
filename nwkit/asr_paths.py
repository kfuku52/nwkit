"""Finite-grid Gaussian histories on a fitted ASR tree.

Subdividing a homogeneous BM/OU branch preserves its exact endpoint transition.
Unary nodes are introduced on a private tree; fitting and the input tree are
unchanged. Paths share the same sampled value at every branching point.
"""

import math
from dataclasses import dataclass
from typing import Any

import numpy as np
from ete4 import Tree

from nwkit.asr_regimes import RegimeAssignment
from nwkit.rooting_state import copy_rooting_info

PATH_MODELS = frozenset(
    {
        "BM",
        "BM-DRIFT",
        "BMS",
        "BMS-DRIFT",
        "OU",
        "OUM",
        "OUMA",
        "OUMV",
        "OUMVA",
        "MV-BM",
        "MV-OU",
        "MV-OU-DIAG",
        "MV-OU-FULL",
    }
)
_VECTOR_MODELS = {"MV-BM", "MV-OU", "MV-OU-DIAG", "MV-OU-FULL"}
_MAX_GRID_NODES = 100_000
_MAX_PATH_VALUES = 2_000_000
_MAX_PATH_TRACES = 50_000
_MAX_VECTOR_GRID_BYTES = 512 * 1024**2


@dataclass(frozen=True)
class BranchPath:
    elapsed: np.ndarray
    values: np.ndarray  # history, grid position (including parent), trait


@dataclass(frozen=True)
class SimulationPaths:
    branches: dict[Any, BranchPath]
    root_values: np.ndarray  # history, trait
    mode: str
    root_description: str
    steps: int

    @property
    def count(self):
        return self.root_values.shape[0]


def _integer(value, name, minimum):
    if (
        isinstance(value, (bool, np.bool_))
        or not isinstance(value, (int, np.integer))
        or value < minimum
    ):
        raise ValueError(f"{name} must be an integer >= {minimum}.")
    return int(value)


def simulation_options(args):
    count = _integer(getattr(args, "figure_simulations", 0), "--figure-simulations", 0)
    steps_arg = getattr(args, "figure_simulation_steps", None)
    steps = _integer(
        200 if steps_arg is None else steps_arg, "--figure-simulation-steps", 1
    )
    mode_arg = getattr(args, "figure_simulation_mode", None)
    mode = "unconditional" if mode_arg is None else mode_arg
    if mode not in {"unconditional", "conditional"}:
        raise ValueError(
            "--figure-simulation-mode must be unconditional or conditional."
        )
    if count and getattr(args, "figure_out", None) in (None, ""):
        raise ValueError("--figure-simulations requires --figure-out.")
    if not count and (steps_arg is not None or mode_arg is not None):
        raise ValueError("Figure simulation settings require --figure-simulations > 0.")
    return count, steps, mode


def validate_path_model(args, model):
    if getattr(args, "figure_simulations", 0) and model not in PATH_MODELS:
        raise ValueError(
            f"--figure-simulations does not support --model {model}; "
            "use a BM/OU model without a branch-length transformation. "
            "ASR node-summary figures remain available for this model."
        )


def _grid_counts(tree, steps, count, dimension, *, vector=False, mode="unconditional"):
    lengths = {node: float(node.dist) for node in tree.traverse() if not node.is_root}
    if any(not math.isfinite(length) or length < 0 for length in lengths.values()):
        raise ValueError("Path simulation requires finite nonnegative branch lengths.")
    longest = max(lengths.values(), default=0.0) or 1.0
    if steps > _MAX_GRID_NODES and any(lengths.values()):
        raise ValueError(
            "Figure simulation grid is too large; reduce --figure-simulation-steps."
        )
    if len(lengths) * count > _MAX_PATH_TRACES:
        raise ValueError(
            "Too many simulated branch traces; reduce --figure-simulations (maximum 50,000 traces)."
        )
    counts = {
        node: max(1, math.ceil((length / longest) * steps))
        for node, length in lengths.items()
    }
    nodes = 1 + sum(counts.values())
    if nodes > _MAX_GRID_NODES or nodes * count * dimension > _MAX_PATH_VALUES:
        raise ValueError(
            "Figure simulation grid is too large; reduce --figure-simulation-steps "
            "or --figure-simulations (at most 100,000 grid nodes and 2,000,000 trait values)."
        )
    if vector:
        # Transitions retain two matrices per node; conditioning retains four
        # more. Include working matrices and both sampled/output path arrays.
        matrices = 8 if mode == "conditional" else 3
        estimated = 8 * nodes * (matrices * dimension**2 + 2 * count * dimension)
        if estimated > _MAX_VECTOR_GRID_BYTES:
            raise ValueError(
                "Figure simulation matrix memory exceeds 512 MiB; reduce "
                "--figure-simulation-steps, --figure-simulations, or trait count."
            )
    return counts


def _refined_tree(tree, counts, assignment):
    refined = Tree()
    refined.name = str(tree.name or "")
    copy_rooting_info(tree, refined)
    endpoints = {tree: refined}
    paths = {}
    regimes = {} if assignment is None else {refined: assignment.by_node[tree]}
    for node in tree.traverse("preorder"):
        if node.is_root:
            continue
        cursor = endpoints[node.up]
        chain = [cursor]
        grid = np.linspace(0.0, float(node.dist), counts[node] + 1)
        lengths = np.diff(grid)
        if float(node.dist) > 0 and np.any(lengths <= 0):
            raise ValueError(
                "Path simulation grid underflowed; rescale the branch lengths."
            )
        for length in lengths:
            cursor = cursor.add_child(dist=float(length))
            chain.append(cursor)
            if assignment is not None:
                regimes[cursor] = assignment.by_node[node]
        cursor.name = str(node.name or "")
        endpoints[node] = cursor
        paths[node] = (tuple(chain), grid)
    mapped = (
        None
        if assignment is None
        else RegimeAssignment(assignment.regimes, regimes, assignment.source)
    )
    return refined, paths, mapped


def _scalar_samples(
    tree,
    observed,
    errors,
    fit,
    model,
    root_prior,
    assignment,
    root_mean,
    count,
    mode,
    seed,
):
    from nwkit.continuous_asr_process import fitted_scalar_process
    from nwkit.gaussian_inference import (
        sample_gaussian_posterior,
        simulate_gaussian_process,
    )

    process = fitted_scalar_process(
        tree, model, fit, root_prior=root_prior, regime_assignment=assignment
    )
    if mode == "conditional":
        samples = sample_gaussian_posterior(
            process, observed, standard_errors=errors, num_samples=count, seed=seed
        )
        description = (
            "Conditioned on observed tips, including their measurement errors."
        )
    else:
        roots = float(root_mean) if process.root.mode == "flat" else None
        samples = simulate_gaussian_process(
            process, num_samples=count, seed=seed, root_values=roots
        )
        description = (
            "Root fixed at inferred mean (flat root prior)."
            if roots is not None
            else "Root fixed at the model's specified value."
            if process.root.mode == "fixed"
            else f"Root drawn from the {process.root.mode} prior."
        )
    return samples.nodes, samples.values[:, :, None], description


def _vector_samples(tree, observed, errors, fit, model, root_mean, count, mode, seed):
    from nwkit.asr_multivariate_diagnostics import (
        fitted_vector_process,
        vector_error_covariances,
    )
    from nwkit.vector_gaussian import condition_vector_tree
    from nwkit.vector_simulation import simulate_vector_process

    process = fitted_vector_process(tree, model, fit)
    if mode == "conditional":
        posterior = condition_vector_tree(
            process,
            observed,
            error_covariances=vector_error_covariances(observed, errors, fit),
        )
        nodes, values = posterior.nodes, posterior.sample(count, seed=seed)
        description = (
            "Conditioned on observed tips, including their measurement errors."
        )
    else:
        roots = root_mean if process.root_mean is None else None
        nodes, values = simulate_vector_process(
            process, count, seed=seed, root_values=roots
        )
        description = (
            "Root fixed at inferred mean (flat root prior)."
            if roots is not None
            else "Root drawn from the stationary prior."
        )
    return nodes, values, description


def simulate_fitted_paths(
    tree,
    observed,
    errors,
    posterior,
    *,
    fit,
    model,
    root_prior=None,
    regime_assignment=None,
    count=1,
    steps=200,
    mode="unconditional",
    seed=None,
):
    """Sample exact Gaussian transitions at a finite grid along every branch.

    ``unconditional`` generates new latent histories with the fitted parameters;
    improper-root models start at the inferred root mean. ``conditional`` draws
    from the joint posterior, retaining noisy and missing tip semantics.
    """
    count = _integer(count, "Path count", 1)
    steps = _integer(steps, "Path steps", 1)
    if mode not in {"unconditional", "conditional"}:
        raise ValueError("Path mode must be unconditional or conditional.")
    if model not in PATH_MODELS:
        raise ValueError(f"Path simulation does not support model {model}.")
    dimension = len(fit.trait_names) if model in _VECTOR_MODELS else 1
    counts = _grid_counts(
        tree, steps, count, dimension, vector=model in _VECTOR_MODELS, mode=mode
    )
    refined, chains, assignment = _refined_tree(tree, counts, regime_assignment)
    if model in _VECTOR_MODELS:
        nodes, values, description = _vector_samples(
            refined,
            observed,
            errors,
            fit,
            model,
            posterior[tree].mean,
            count,
            mode,
            seed,
        )
    else:
        nodes, values, description = _scalar_samples(
            refined,
            observed,
            errors,
            fit,
            model,
            root_prior,
            assignment,
            posterior[tree].mean,
            count,
            mode,
            seed,
        )
    if not np.isfinite(values).all():
        raise ValueError(
            "A simulated path exceeds floating-point range; rescale trait units."
        )
    indices = {node: index for index, node in enumerate(nodes)}
    branches = {
        node: BranchPath(grid, values[:, [indices[point] for point in chain], :])
        for node, (chain, grid) in chains.items()
    }
    return SimulationPaths(
        branches, values[:, indices[refined], :], mode, description, steps
    )

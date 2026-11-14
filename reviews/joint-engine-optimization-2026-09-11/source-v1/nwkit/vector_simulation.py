"""Joint vector process simulation with explicit treatment of an improper root."""

import numpy as np

from nwkit.compiled_tree import CompiledTree


def simulate_vector_process(process, num_samples=1, *, seed=None, root_values=None):
    """Simulate all nodes; a flat root requires explicitly supplied root values."""
    if (
        isinstance(num_samples, bool)
        or not isinstance(num_samples, (int, np.integer))
        or num_samples < 1
    ):
        raise ValueError("num_samples must be a positive integer.")
    compiled = CompiledTree.from_tree(process.tree)
    rng = np.random.default_rng(seed)
    dimension = process.dimension
    result = np.empty((num_samples, len(compiled.nodes), dimension))
    if root_values is None:
        if process.root_mean is None or process.root_covariance is None:
            raise ValueError(
                "A flat root requires explicit root_values for simulation."
            )
        result[:, 0] = rng.multivariate_normal(
            process.root_mean,
            process.root_covariance,
            size=num_samples,
            check_valid="raise",
        )
    else:
        roots = np.asarray(root_values, dtype=float)
        if (
            roots.shape not in {(dimension,), (num_samples, dimension)}
            or not np.isfinite(roots).all()
        ):
            raise ValueError(
                "root_values must be a finite vector or one vector per sample."
            )
        result[:, 0] = roots
    indices = {node: index for index, node in enumerate(compiled.nodes)}
    for index, node in enumerate(compiled.nodes[1:], 1):
        transition = process.transitions[node]
        result[:, index] = (
            result[:, indices[node.up]] @ transition.slope.T
            + transition.intercept
            + rng.multivariate_normal(
                np.zeros(dimension),
                transition.covariance,
                size=num_samples,
                check_valid="raise",
            )
        )
    return compiled.nodes, result


def simulated_vector_observations(
    process, observed, errors, num_samples, seed, *, root_values=None
):
    """Simulate observation noise jointly across traits and retain the missing mask."""
    process_seed, error_seed = np.random.SeedSequence(seed).spawn(2)
    nodes, latent = simulate_vector_process(
        process,
        num_samples,
        seed=int(process_seed.generate_state(1)[0]),
        root_values=root_values,
    )
    rng = np.random.default_rng(error_seed)
    result = [dict.fromkeys(observed) for _ in range(num_samples)]
    for index, node in enumerate(nodes):
        name = str(node.name)
        if not node.is_leaf or observed.get(name) is None:
            continue
        values = latent[:, index].copy()
        if errors is not None:
            values += rng.multivariate_normal(
                np.zeros(process.dimension),
                errors[name],
                size=num_samples,
                check_valid="raise",
            )
        for replicate in range(num_samples):
            result[replicate][name] = [
                None if original is None else float(value)
                for original, value in zip(
                    observed[name], values[replicate], strict=True
                )
            ]
    return result

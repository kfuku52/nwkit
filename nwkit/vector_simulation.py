"""Joint vector process simulation with explicit treatment of an improper root."""

import numpy as np

from nwkit.compiled_tree import CompiledTree


def _zero_mean_draws(rng, dimension, count):
    """Reuse the SVD for consecutive draws with identical covariance."""
    previous = None
    factor = None
    zero = np.zeros(dimension)

    def draw(covariance):
        nonlocal previous, factor
        covariance = np.asarray(covariance)
        if covariance.shape != (dimension, dimension) or covariance.dtype != np.float64:
            return rng.multivariate_normal(
                zero, covariance, size=count, check_valid="raise"
            )
        bits = covariance.view(np.uint64)
        if previous is None or not np.array_equal(bits, previous):
            # Match Generator.multivariate_normal's default SVD and validation
            # so the same standard-normal draws use the same factor.
            vectors, singular, right = np.linalg.svd(covariance)
            tolerance = 1e-8
            if not np.allclose(
                (right.T * singular) @ right, covariance, rtol=tolerance, atol=tolerance
            ):
                raise ValueError("covariance is not symmetric positive-semidefinite.")
            factor = vectors * np.sqrt(singular)
            previous = bits.copy()
        draws = rng.standard_normal((count, dimension)) @ factor.T
        draws += zero
        return draws

    return draw


def simulate_vector_process(
    process, num_samples=1, *, seed=None, root_values=None, compiled_tree=None
):
    """Simulate all nodes; a flat root requires explicitly supplied root values."""
    if (
        isinstance(num_samples, bool)
        or not isinstance(num_samples, (int, np.integer))
        or num_samples < 1
    ):
        raise ValueError("num_samples must be a positive integer.")
    compiled = (
        CompiledTree.from_tree(process.tree) if compiled_tree is None else compiled_tree
    )
    compiled.require_tree(process.tree)
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
    draw = _zero_mean_draws(rng, dimension, num_samples)
    for index, node in enumerate(compiled.nodes[1:], 1):
        transition = process.transitions[node]
        result[:, index] = (
            result[:, compiled.parents[index]] @ transition.slope.T
            + transition.intercept
            + draw(transition.covariance)
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
    draw = _zero_mean_draws(rng, process.dimension, num_samples)
    result = [dict.fromkeys(observed) for _ in range(num_samples)]
    for index, node in enumerate(nodes):
        name = str(node.name)
        if not node.is_leaf or observed.get(name) is None:
            continue
        values = latent[:, index]
        if errors is not None:
            values = values + draw(errors[name])
        for replicate in range(num_samples):
            result[replicate][name] = [
                None if original is None else float(value)
                for original, value in zip(
                    observed[name], values[replicate], strict=True
                )
            ]
    return result

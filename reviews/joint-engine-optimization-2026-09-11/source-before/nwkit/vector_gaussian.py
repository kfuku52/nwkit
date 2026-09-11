"""Vector affine-Gaussian tree pruning, smoothing, and joint posterior draws.

Storage is O(nodes * traits**2), independent of observed-pair count. Branch
innovation matrices must be positive definite. Exact coordinate observations
are eliminated algebraically, without adding a numerical noise floor.
"""

import math
from dataclasses import dataclass

import numpy as np
from scipy.linalg import cho_factor, cho_solve

from nwkit.compiled_tree import CompiledTree


@dataclass(frozen=True)
class VectorTransition:
    slope: np.ndarray
    intercept: np.ndarray
    covariance: np.ndarray


@dataclass(frozen=True)
class VectorProcess:
    tree: object
    dimension: int
    transitions: dict
    root_mean: np.ndarray | None = None
    root_covariance: np.ndarray | None = None


@dataclass(frozen=True)
class VectorConditioning:
    nodes: tuple
    means: np.ndarray
    covariances: np.ndarray
    log_likelihood: float
    conditional_slopes: np.ndarray
    conditional_intercepts: np.ndarray
    conditional_covariances: np.ndarray
    parents: tuple

    def sample(self, num_samples=1000, seed=None):
        """Draw jointly across nodes, preserving parent-child dependence."""
        if (
            isinstance(num_samples, bool)
            or not isinstance(num_samples, (int, np.integer))
            or num_samples < 1
        ):
            raise ValueError("num_samples must be a positive integer.")
        rng = np.random.default_rng(seed)
        result = np.empty((num_samples, len(self.nodes), self.means.shape[1]))
        for index in range(len(self.nodes)):
            covariance = self.conditional_covariances[index]
            eigenvalues, eigenvectors = np.linalg.eigh((covariance + covariance.T) / 2)
            tolerance = (
                np.finfo(float).eps * max(1.0, float(np.max(np.abs(eigenvalues)))) * 100
            )
            if np.min(eigenvalues) < -tolerance:
                raise ValueError(
                    "Posterior conditional covariance is not positive semidefinite."
                )
            factor = eigenvectors * np.sqrt(np.maximum(eigenvalues, 0))
            location = self.conditional_intercepts[index]
            if index:
                location = (
                    location
                    + result[:, self.parents[index]] @ self.conditional_slopes[index].T
                )
            result[:, index] = (
                location
                + rng.normal(size=(num_samples, self.means.shape[1])) @ factor.T
            )
        return result


def _array(value, shape, label):
    result = np.asarray(value, dtype=float)
    if result.shape != shape or not np.isfinite(result).all():
        raise ValueError(f"{label} must have shape {shape} and finite values.")
    return result


def _inverse(matrix, label):
    if not np.allclose(matrix, matrix.T, rtol=1e-12, atol=1e-14):
        raise ValueError(f"{label} must be symmetric.")
    try:
        factor = cho_factor(matrix, lower=True, check_finite=False)
    except np.linalg.LinAlgError as exc:
        raise ValueError(f"{label} must be positive definite.") from exc
    return cho_solve(factor, np.eye(len(matrix))), 2 * float(
        np.log(np.diag(factor[0])).sum()
    )


def _observations(
    process, compiled, observed, error_covariances
) -> list[tuple[np.ndarray, np.ndarray, float, dict[int, float]]]:
    d = process.dimension
    factors: list[tuple[np.ndarray, np.ndarray, float, dict[int, float]]] = [
        (np.zeros((d, d)), np.zeros(d), 0.0, {}) for _ in compiled.nodes
    ]
    unknown = set(observed) - set(compiled.leaf_index_by_name)
    if unknown:
        raise ValueError(
            "Unknown vector observation tips: " + ", ".join(sorted(unknown))
        )
    for name, values in observed.items():
        if values is None:
            continue
        if len(values) != d:
            raise ValueError(f"Trait dimension mismatch for tip '{name}'.")
        indices = np.array(
            [i for i, value in enumerate(values) if value is not None], dtype=int
        )
        if not len(indices):
            continue
        y = _array([values[i] for i in indices], (len(indices),), "Observations")
        if error_covariances is not None and name not in error_covariances:
            raise ValueError(f"Observation covariance is required for tip '{name}'.")
        covariance = (
            np.zeros((d, d))
            if error_covariances is None
            else _array(error_covariances[name], (d, d), "Observation covariance")
        )
        covariance = covariance[np.ix_(indices, indices)]
        exact = np.flatnonzero(np.diag(covariance) == 0)
        if len(exact) and np.any(covariance[exact] != 0):
            raise ValueError(
                "An exact observed coordinate cannot have nonzero error covariance."
            )
        noisy = np.flatnonzero(np.diag(covariance) != 0)
        precision, linear, constant, fixed = factors[compiled.leaf_index_by_name[name]]
        fixed.update({int(indices[i]): float(y[i]) for i in exact})
        if len(noisy):
            inverse, logdet = _inverse(
                covariance[np.ix_(noisy, noisy)], "Observation covariance"
            )
            selected = indices[noisy]
            precision[np.ix_(selected, selected)] += inverse
            linear[selected] += inverse @ y[noisy]
            constant = -0.5 * (
                len(noisy) * math.log(2 * math.pi)
                + logdet
                + y[noisy] @ inverse @ y[noisy]
            )
        factors[compiled.leaf_index_by_name[name]] = precision, linear, constant, fixed
    return factors


def _restrict(precision, linear, constant, fixed):
    dimension = len(linear)
    coordinates = np.array([fixed.get(i, 0.0) for i in range(dimension)])
    free = np.array([i for i in range(dimension) if i not in fixed], dtype=int)
    reduced_linear = (linear - precision @ coordinates)[free]
    constant += linear @ coordinates - 0.5 * coordinates @ precision @ coordinates
    return free, coordinates, reduced_linear, float(constant)


def _propagate(precision, linear, constant, fixed, transition, dimension):
    slope = _array(transition.slope, (dimension, dimension), "Branch slope")
    intercept = _array(transition.intercept, (dimension,), "Branch intercept")
    covariance = _array(
        transition.covariance, (dimension, dimension), "Branch covariance"
    )
    _inverse(covariance, "Branch covariance")
    free, coordinates, reduced_linear, constant = _restrict(
        precision, linear, constant, fixed
    )
    message_precision = np.zeros((dimension, dimension))
    message_linear = np.zeros(dimension)
    constrained = np.array(sorted(fixed), dtype=int)
    conditional_prior_slope = slope[free].copy()
    conditional_prior_mean = intercept[free].copy()
    conditional_prior_covariance = covariance[np.ix_(free, free)].copy()
    if len(constrained):
        exact_inverse, exact_logdet = _inverse(
            covariance[np.ix_(constrained, constrained)], "Exact-coordinate innovation"
        )
        delta = coordinates[constrained] - intercept[constrained]
        exact_slope = slope[constrained]
        message_precision = exact_slope.T @ exact_inverse @ exact_slope
        message_linear = exact_slope.T @ exact_inverse @ delta
        constant -= 0.5 * (
            len(constrained) * math.log(2 * math.pi)
            + exact_logdet
            + delta @ exact_inverse @ delta
        )
        gain = covariance[np.ix_(free, constrained)] @ exact_inverse
        conditional_prior_slope -= gain @ exact_slope
        conditional_prior_mean += gain @ delta
        conditional_prior_covariance -= gain @ covariance[np.ix_(constrained, free)]
    conditional_covariance = np.zeros((dimension, dimension))
    conditional_slope = np.zeros((dimension, dimension))
    location = coordinates.copy()
    if len(free):
        local_precision = precision[np.ix_(free, free)]
        update = np.eye(len(free)) + conditional_prior_covariance @ local_precision
        # Solve the covariance-form update directly: subtracting two nearly
        # equal innovation precisions loses weak/missing-subtree information.
        gain = np.linalg.solve(update, np.eye(len(free)))
        posterior_covariance = gain @ conditional_prior_covariance
        posterior_covariance = (posterior_covariance + posterior_covariance.T) / 2
        effective_precision = local_precision @ gain
        effective_linear = gain.T @ reduced_linear
        sign, logdet = np.linalg.slogdet(update)
        if sign <= 0:
            raise ValueError("Invalid vector Gaussian covariance update.")
        mean = conditional_prior_mean
        matrix = conditional_prior_slope
        message_precision += matrix.T @ effective_precision @ matrix
        message_linear += matrix.T @ (effective_linear - effective_precision @ mean)
        constant += (
            0.5
            * (
                reduced_linear @ posterior_covariance @ reduced_linear
                - logdet
                - mean @ effective_precision @ mean
            )
            + effective_linear @ mean
        )
        conditional_covariance[np.ix_(free, free)] = posterior_covariance
        conditional_slope[free] = gain @ matrix
        location[free] = gain @ mean + posterior_covariance @ reduced_linear
    message_precision = (message_precision + message_precision.T) / 2
    return (message_precision, message_linear, constant), (
        conditional_slope,
        location,
        conditional_covariance,
    )


def _root(precision, linear, constant, fixed, process):
    d = process.dimension
    if process.root_mean is not None:
        mean = _array(process.root_mean, (d,), "Root mean")
        covariance = _array(process.root_covariance, (d, d), "Root covariance")
        # Integrate in covariance form. Adding a nearly singular root's huge
        # precision and then subtracting its quadratic form can create a
        # spurious positive likelihood during covariance optimization.
        message, conditional = _propagate(
            precision,
            linear,
            constant,
            fixed,
            VectorTransition(np.zeros((d, d)), mean, covariance),
            d,
        )
        return conditional[1], conditional[2], float(message[2])
    elif process.root_covariance is not None:
        raise ValueError("Root covariance requires a root mean.")
    free, mean, reduced_linear, constant = _restrict(precision, linear, constant, fixed)
    covariance = np.zeros((d, d))
    if len(free):
        inverse, logdet = _inverse(
            precision[np.ix_(free, free)], "Root posterior precision"
        )
        mean[free] = inverse @ reduced_linear
        covariance[np.ix_(free, free)] = inverse
        constant += 0.5 * (
            len(free) * math.log(2 * math.pi) - logdet + reduced_linear @ mean[free]
        )
    return mean, covariance, float(constant)


def condition_vector_tree(process, observed, *, error_covariances=None):
    """Condition vector tips, allowing missing and exact/noisy coordinates."""
    d = process.dimension
    if isinstance(d, bool) or not isinstance(d, int) or d < 1:
        raise ValueError("Vector dimension must be a positive integer.")
    compiled = CompiledTree.from_tree(process.tree)
    local = _observations(process, compiled, observed, error_covariances)
    n = len(compiled.nodes)
    slopes = np.zeros((n, d, d))
    intercepts = np.zeros((n, d))
    covariances = np.zeros((n, d, d))
    for index in compiled.postorder:
        if index == 0:
            break
        precision, linear, constant, fixed = local[index]
        parent = compiled.parents[index]
        parent_precision, parent_linear, parent_constant, parent_fixed = local[parent]
        transition = process.transitions[compiled.nodes[index]]
        if np.array_equal(transition.covariance, np.zeros((d, d))):
            if not np.array_equal(transition.slope, np.eye(d)) or not np.array_equal(
                transition.intercept, np.zeros(d)
            ):
                raise ValueError(
                    "Zero-variance vector branches must be identity transitions."
                )
            for trait, value in fixed.items():
                if trait in parent_fixed and parent_fixed[trait] != value:
                    raise ValueError(
                        "Conflicting exact observations across zero-length branches."
                    )
                parent_fixed[trait] = value
            message = precision, linear, constant
            slopes[index] = np.eye(d)
        else:
            message, conditional = _propagate(
                precision, linear, constant, fixed, transition, d
            )
            slopes[index], intercepts[index], covariances[index] = conditional
        local[parent] = (
            parent_precision + message[0],
            parent_linear + message[1],
            parent_constant + message[2],
            parent_fixed,
        )
    precision, linear, constant, fixed = local[0]
    intercepts[0], covariances[0], likelihood = _root(
        precision, linear, constant, fixed, process
    )
    means = intercepts.copy()
    marginal_covariances = covariances.copy()
    for index in range(1, n):
        parent = compiled.parents[index]
        means[index] += slopes[index] @ means[parent]
        marginal_covariances[index] += (
            slopes[index] @ marginal_covariances[parent] @ slopes[index].T
        )
    return VectorConditioning(
        compiled.nodes,
        means,
        marginal_covariances,
        likelihood,
        slopes,
        intercepts,
        covariances,
        compiled.parents,
    )

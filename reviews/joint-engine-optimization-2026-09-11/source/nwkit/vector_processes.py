"""Brownian and stable full-attraction Ornstein-Uhlenbeck vector processes."""

import numpy as np
from scipy.linalg import expm, solve_continuous_lyapunov

from nwkit.compiled_tree import CompiledTree
from nwkit.vector_gaussian import VectorProcess, VectorTransition, _array, _inverse


def vector_brownian_process(tree, sigma):
    """Build flat-root correlated Brownian motion, including zero-length edges."""
    sigma = np.asarray(sigma, dtype=float)
    if sigma.ndim != 2 or sigma.shape[0] != sigma.shape[1] or not len(sigma):
        raise ValueError("Diffusion covariance must be a nonempty square matrix.")
    d = len(sigma)
    sigma = _array(sigma, (d, d), "Diffusion covariance")
    _inverse(sigma, "Diffusion covariance")
    transitions = {}
    for node in CompiledTree.from_tree(tree).nodes[1:]:
        time = float(node.dist)
        if not np.isfinite(time) or time < 0:
            raise ValueError("Vector branch lengths must be finite and nonnegative.")
        transitions[node] = VectorTransition(np.eye(d), np.zeros(d), sigma * time)
    return VectorProcess(tree, d, transitions)


def vector_ou_process(tree, attraction, diffusion, theta):
    """Build stationary OU with a general (possibly nonsymmetric) stable drift.

    SDE: dX = -A(X-theta)dt + LdW, LL' = diffusion. Positive real
    eigenvalues of A are required. No diagonal/eigenvector approximation is
    made. A block exponential avoids subtractive cancellation on short edges.
    """
    theta = np.asarray(theta, dtype=float)
    if theta.ndim != 1 or not len(theta):
        raise ValueError("OU theta must be a nonempty vector.")
    d = len(theta)
    theta = _array(theta, (d,), "OU theta")
    attraction = _array(attraction, (d, d), "OU attraction")
    diffusion = _array(diffusion, (d, d), "OU diffusion")
    _inverse(diffusion, "OU diffusion")
    # Solve in diffusion-standardized trait units: a valid similarity-scaled
    # attraction can otherwise defeat the Lyapunov solver's absolute tolerance.
    scales = np.sqrt(np.diag(diffusion))
    attraction = attraction * scales[None, :] / scales[:, None]
    diffusion = diffusion / np.outer(scales, scales)
    if np.min(np.linalg.eigvals(attraction).real) <= 0:
        raise ValueError("OU attraction must have strictly positive real eigenvalues.")
    stationary = solve_continuous_lyapunov(attraction, diffusion)
    stationary = (stationary + stationary.T) / 2
    _inverse(stationary, "OU stationary covariance")
    transitions = {}
    block = np.block([[-attraction, diffusion], [np.zeros((d, d)), attraction.T]])
    for node in CompiledTree.from_tree(tree).nodes[1:]:
        time = float(node.dist)
        if not np.isfinite(time) or time < 0:
            raise ValueError("Vector branch lengths must be finite and nonnegative.")
        if time * np.linalg.norm(attraction, ord=1) < 0.5:
            exponential = expm(block * time)
            slope = exponential[:d, :d]
            covariance = exponential[:d, d:] @ slope.T
        else:
            slope = expm(-attraction * time)
            covariance = stationary - slope @ stationary @ slope.T
        slope = slope * scales[:, None] / scales[None, :]
        covariance = covariance * np.outer(scales, scales)
        transitions[node] = VectorTransition(
            slope, theta - slope @ theta, (covariance + covariance.T) / 2
        )
    return VectorProcess(
        tree, d, transitions, theta, stationary * np.outer(scales, scales)
    )

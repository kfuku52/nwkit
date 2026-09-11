"""Importance integration of latent histories with conditional Gaussian pruning.

Parameters are fixed, not estimated. Independent histories are proposed from the
process prior (or the CTMC posterior given discrete tips), then reweighted by the
continuous likelihood. The finite mixture is an approximation, not an MCMC fit.
"""

import math
from dataclasses import dataclass
from typing import Any

import numpy as np
from scipy.special import logsumexp

from nwkit.asr_averaging import gaussian_mixture_summary
from nwkit.compiled_tree import CompiledTree
from nwkit.gaussian_inference import (
    GaussianSamples,
    condition_gaussian_tree,
    sample_gaussian_posterior,
)


@dataclass(frozen=True)
class MixtureMarginal:
    means: np.ndarray
    variances: np.ndarray
    weights: np.ndarray

    @property
    def mean(self):
        return float(self.weights @ self.means)

    @property
    def variance(self):
        return float(self.weights @ (self.variances + (self.means - self.mean) ** 2))

    def interval(self, level):
        summary = gaussian_mixture_summary(
            self.means, self.variances, self.weights, level=level
        )
        return summary.lower, summary.upper


@dataclass
class LatentGaussianFit:
    model: str
    nodes: tuple
    processes: list
    latent: list
    weights: np.ndarray
    marginals: dict
    log_likelihood: float
    continuous_log_likelihood: float
    discrete_log_likelihood: float
    effective_sample_size: float
    relative_mcse: float
    num_observed: int
    parameters: dict

    @property
    def fit_status(self):
        return (
            "importance_degenerate"
            if self.effective_sample_size < max(20, 0.1 * len(self.weights))
            else "monte_carlo"
        )

    def sample(self, observed, errors, count, seed=None):
        """Select a history once per joint draw, preserving node dependence."""
        count = positive_integer(count, "posterior sample count")
        rng = np.random.default_rng(seed)
        assignments = rng.choice(len(self.weights), count, p=self.weights)
        values = np.empty((count, len(self.nodes)))
        for index in np.unique(assignments):
            selected = np.flatnonzero(assignments == index)
            draws = sample_gaussian_posterior(
                self.processes[index],
                observed,
                standard_errors=errors,
                num_samples=len(selected),
                seed=int(rng.integers(2**32)),
            )
            values[selected] = draws.values
        return GaussianSamples(self.nodes, values)


def positive_integer(value, name):
    if isinstance(value, (bool, np.bool_)):
        raise ValueError(f"{name} must be a positive integer.")
    try:
        integer = int(value)
        valid = integer == float(value) and integer > 0
    except (ValueError, TypeError, OverflowError):
        valid = False
    if not valid:
        raise ValueError(f"{name} must be a positive integer.")
    return integer


def integrate_histories(
    tree,
    observed,
    errors,
    propose,
    *,
    model,
    samples,
    seed,
    parameters,
    discrete_log_likelihood=0.0,
):
    """Self-normalized importance posterior and unbiased likelihood estimator.

    The log of the likelihood estimator is biased. ``relative_mcse`` is the
    delta-method standard error of its logarithm, reliable only with adequate
    importance overlap; ESS alone does not establish convergence.
    """
    samples = positive_integer(samples, "history samples")
    if samples < 2:
        raise ValueError("At least two history samples are required for MC error.")
    compiled = CompiledTree.from_tree(tree)
    if samples * len(compiled.nodes) > 1_000_000:
        raise ValueError("Latent integration exceeds 1,000,000 history-node pairs.")
    rng = np.random.default_rng(seed)
    means = np.empty((samples, len(compiled.nodes)))
    variances = np.empty_like(means)
    log_weights = np.empty(samples)
    processes: list[Any] = []
    latent: list[Any] = []
    for index in range(samples):
        process, history = propose(rng)
        conditioned = condition_gaussian_tree(
            process, observed, standard_errors=errors, compiled_tree=compiled
        )
        log_weights[index] = conditioned.log_likelihood
        if not math.isfinite(log_weights[index]):
            raise ValueError("A latent history has a nonfinite Gaussian likelihood.")
        for node_index, node in enumerate(compiled.nodes):
            means[index, node_index] = conditioned.marginals[node].mean
            variances[index, node_index] = conditioned.marginals[node].variance
        processes.append(process)
        latent.append(history)
    normalizer = float(logsumexp(log_weights))
    weights = np.exp(log_weights - normalizer)
    weights /= weights.sum()
    ess = float(1 / (weights @ weights))
    relative_mcse = math.sqrt(max(0.0, (samples / ess - 1) / (samples - 1)))
    marginals = {
        node: MixtureMarginal(means[:, index], variances[:, index], weights)
        for index, node in enumerate(compiled.nodes)
    }
    continuous_ll = normalizer - math.log(samples)
    return LatentGaussianFit(
        model,
        compiled.nodes,
        processes,
        latent,
        weights,
        marginals,
        continuous_ll + discrete_log_likelihood,
        continuous_ll,
        discrete_log_likelihood,
        ess,
        relative_mcse,
        sum(value is not None for value in observed.values()),
        parameters,
    )

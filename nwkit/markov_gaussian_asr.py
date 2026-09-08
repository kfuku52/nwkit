"""CTMC-modulated scalar BM/OU, integrating uncertain regime histories."""

import math

import numpy as np

from nwkit.gaussian_tree import (
    GaussianRootPrior,
    GaussianTransition,
    GaussianTreeProcess,
)
from nwkit.jump_asr import finite_parameter
from nwkit.latent_gaussian import integrate_histories


def validate_regime_config(config, model):
    """Strict, explicitly fixed CTMC and per-state diffusion/OU parameters."""
    required = {"states", "q", "sigma2"}
    if model == "MM-OU":
        required |= {"alpha", "theta"}
    if not isinstance(config, dict) or set(config) != required:
        raise ValueError(f"{model} config must contain exactly {sorted(required)}.")
    states = config["states"]
    if (
        not isinstance(states, list)
        or len(states) < 2
        or any(not isinstance(state, str) or not state for state in states)
        or len(set(states)) != len(states)
    ):
        raise ValueError("Regime states must be at least two unique nonempty strings.")
    q = np.asarray(config["q"], dtype=float)
    if q.shape != (len(states), len(states)) or not np.isfinite(q).all():
        raise ValueError("Regime q must be a finite square matrix matching states.")
    off_diagonal = q.copy()
    np.fill_diagonal(off_diagonal, 0)
    scale = max(float(np.abs(q).max()), np.finfo(float).tiny)
    if np.any(off_diagonal < 0) or not np.allclose(
        q.sum(axis=1), 0, atol=1e-12 * scale, rtol=0
    ):
        raise ValueError("Regime q needs nonnegative off-diagonals and zero row sums.")
    vectors = {}
    for key in required - {"states", "q"}:
        vector = np.asarray(config[key], dtype=float)
        if vector.shape != (len(states),) or not np.isfinite(vector).all():
            raise ValueError(f"Regime {key} must be finite and match states.")
        if key != "theta" and np.any(vector <= 0):
            raise ValueError(f"Regime {key} must be strictly positive.")
        vectors[key] = vector
    return tuple(states), q, vectors


def sample_bridge_segments(start, end, length, context, rng):
    """Uniformization bridge including event times and self-transition durations."""
    from nwkit.asr import _bridge_probabilities, _draw_bridge_event_count

    count, backward = _draw_bridge_event_count(start, end, context, rng)
    times = np.concatenate(([0.0], np.sort(rng.uniform(0, length, count)), [length]))
    current = start
    segments = []
    for index in range(count):
        segments.append((current, float(times[index + 1] - times[index])))
        probabilities = _bridge_probabilities(
            context["r_matrix"][current] * backward[count - index - 1]
        )
        current = int(rng.choice(len(probabilities), p=probabilities))
    segments.append((current, float(times[-1] - times[-2])))
    if current != end:
        raise ValueError("Regime bridge failed to reach its conditioned endpoint.")
    return tuple(segments)


def history_transition(segments, parameters):
    slope, intercept, variance = 1.0, 0.0, 0.0
    for state, duration in segments:
        sigma2 = parameters["sigma2"][state]
        if "alpha" in parameters:
            alpha = parameters["alpha"][state]
            exponent = alpha * duration
            attenuation = math.exp(-exponent)
            shift = -math.expm1(-exponent) * parameters["theta"][state]
            innovation = sigma2 * -math.expm1(-2 * exponent) / (2 * alpha)
        else:
            attenuation, shift, innovation = 1.0, 0.0, sigma2 * duration
        slope *= attenuation
        intercept = attenuation * intercept + shift
        variance = attenuation**2 * variance + innovation
    if slope == 0:
        raise ValueError(
            "MM-OU attenuation underflow; rescale rates or branch lengths."
        )
    return GaussianTransition(slope, intercept, variance)


def fit_markov_gaussian(
    tree,
    observed,
    regime_observed,
    config,
    *,
    model="MM-BM",
    standard_errors=None,
    history_samples=1000,
    seed=None,
):
    """Propose histories conditional on discrete tips; reweight by continuous data.

    The regime root is uniform; the continuous root is flat, independent of the
    regime before conditioning. This is not a stationary-root switching OU.
    """
    from nwkit.asr import (
        _build_uniformization_context,
        _sample_node_states,
        compute_mk_marginals,
    )

    if model not in {"MM-BM", "MM-OU"}:
        raise ValueError("Markov Gaussian model must be MM-BM or MM-OU.")
    states, q, parameters = validate_regime_config(config, model)
    names = set(tree.leaf_names())
    if set(regime_observed) - names:
        raise ValueError("Regime observations include tips absent from the tree.")
    if any(
        value is not None and value not in states for value in regime_observed.values()
    ):
        raise ValueError("Regime observations contain a state absent from config.")
    likelihoods = {
        name: np.ones(len(states))
        if regime_observed.get(name) is None
        else np.asarray(
            [state == regime_observed[name] for state in states], dtype=float
        )
        for name in names
    }
    _, discrete_fit = compute_mk_marginals(
        tree, states, regime_observed, likelihoods, model="CUSTOM", fixed_rate_matrix=q
    )
    branches = tuple(node for node in tree.traverse() if not node.is_root)
    lengths = {node: finite_parameter(node.dist, "branch length") for node in branches}
    contexts = {
        length: _build_uniformization_context(q, length)
        for length in set(lengths.values())
    }

    def propose(rng):
        node_states = _sample_node_states(tree, states, discrete_fit, rng)
        segments = {
            node: sample_bridge_segments(
                node_states[node.up],
                node_states[node],
                lengths[node],
                contexts[lengths[node]],
                rng,
            )
            for node in branches
        }
        process = GaussianTreeProcess(
            tree,
            {
                node: history_transition(value, parameters)
                for node, value in segments.items()
            },
            GaussianRootPrior("flat", variance=None),
            model,
        )
        return process, {"node_states": node_states, "segments": segments}

    return integrate_histories(
        tree,
        observed,
        standard_errors,
        propose,
        model=model,
        samples=history_samples,
        seed=seed,
        parameters=config,
        discrete_log_likelihood=discrete_fit["log_likelihood"],
    )

"""State-space construction for discrete hidden-rates ASR models."""

import numpy as np

from nwkit.discrete_asr_models import complete_transition_graph


def expanded_state_labels(states, hidden_categories):
    return [
        f"hidden_{hidden_index + 1}__state_{state.encode('utf-8').hex()}"
        for hidden_index in range(hidden_categories)
        for state in states
    ]


def expanded_transition_graph(observed_graph, hidden_categories):
    """Allow observed changes within class and hidden changes within state."""

    observed_graph = np.asarray(observed_graph, dtype=bool)
    num_states = observed_graph.shape[0]
    if observed_graph.shape != (num_states, num_states):
        raise ValueError("HRM observed transition graph must be square.")
    expanded_size = num_states * hidden_categories
    graph = np.zeros((expanded_size, expanded_size), dtype=bool)
    hidden_graph = complete_transition_graph(hidden_categories)
    for hidden_index in range(hidden_categories):
        offset = hidden_index * num_states
        graph[offset : offset + num_states, offset : offset + num_states] = (
            observed_graph
        )
    for state_index in range(num_states):
        for from_hidden in range(hidden_categories):
            for to_hidden in range(hidden_categories):
                if hidden_graph[from_hidden, to_hidden]:
                    source = from_hidden * num_states + state_index
                    target = to_hidden * num_states + state_index
                    graph[source, target] = True
    return graph


def expand_tip_likelihoods(likelihood_by_leaf, hidden_categories):
    return {
        leaf: np.tile(np.asarray(likelihood, dtype=float), hidden_categories)
        for leaf, likelihood in likelihood_by_leaf.items()
    }


def state_projection(num_states, hidden_categories):
    return np.tile(np.arange(num_states, dtype=int), hidden_categories)


def aggregate_probabilities(probabilities, num_states, hidden_categories):
    values = np.asarray(probabilities, dtype=float).reshape(
        hidden_categories, num_states
    )
    return values.sum(axis=0)


def covarion_rate_matrix(
    observed_graph,
    hidden_categories,
    base_rate,
    log_rate_spread,
    switching_rate,
    *,
    effective_rate_bounds=None,
):
    """Build an identifiable ordered hidden-rate covarion generator.

    Hidden-class multipliers are log-spaced around a geometric mean of one.
    This removes the arbitrary class labels and the scale confounding between
    the base observed rate and class multipliers.
    """

    observed_graph = np.asarray(observed_graph, dtype=bool)
    num_states = observed_graph.shape[0]
    if observed_graph.shape != (num_states, num_states):
        raise ValueError("A covarion observed transition graph must be square.")
    if np.any(np.diag(observed_graph)):
        raise ValueError("A covarion transition graph cannot contain self edges.")
    if hidden_categories < 2:
        raise ValueError("A covarion model requires at least two hidden classes.")
    values = (base_rate, log_rate_spread, switching_rate)
    if not all(np.isfinite(value) for value in values):
        raise ValueError("Covarion parameters must be finite.")
    if base_rate < 0.0 or log_rate_spread < 0.0 or switching_rate < 0.0:
        raise ValueError("Covarion parameters must be non-negative.")
    multipliers = np.exp(
        np.linspace(-log_rate_spread, log_rate_spread, hidden_categories)
    )
    effective_rates = base_rate * multipliers
    if effective_rate_bounds is not None:
        lower, upper = (float(value) for value in effective_rate_bounds)
        if (
            not np.isfinite(lower)
            or not np.isfinite(upper)
            or lower <= 0.0
            or lower >= upper
        ):
            raise ValueError(
                "Covarion effective-rate bounds must be positive and increasing."
            )
        lower_tolerance = np.finfo(float).eps * max(1.0, lower) * 64.0
        upper_tolerance = np.finfo(float).eps * max(1.0, upper) * 64.0
        if (
            float(effective_rates[0]) < lower - lower_tolerance
            or float(effective_rates[-1]) > upper + upper_tolerance
        ):
            raise ValueError(
                "Covarion hidden-class rates must remain within the fitted rate bounds."
            )
        # The constrained parameterization can land a few ulps beyond an
        # endpoint.  Keep both the reported multipliers and the generator
        # strictly inside the public fitted-rate contract.
        effective_rates = np.clip(effective_rates, lower, upper)
        if base_rate > 0.0:
            multipliers = effective_rates / base_rate
    expanded_size = num_states * hidden_categories
    matrix = np.zeros((expanded_size, expanded_size), dtype=float)
    for hidden_index, effective_rate in enumerate(effective_rates):
        start = hidden_index * num_states
        block = matrix[start : start + num_states, start : start + num_states]
        block[observed_graph] = effective_rate
    hidden_rate = switching_rate / float(hidden_categories - 1)
    for state_index in range(num_states):
        for source_hidden in range(hidden_categories):
            source = source_hidden * num_states + state_index
            for target_hidden in range(hidden_categories):
                if target_hidden != source_hidden:
                    target = target_hidden * num_states + state_index
                    matrix[source, target] = hidden_rate
    np.fill_diagonal(matrix, -matrix.sum(axis=1))
    return matrix, multipliers

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

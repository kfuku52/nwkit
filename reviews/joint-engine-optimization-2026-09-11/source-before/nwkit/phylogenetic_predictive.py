"""Topology-aware predictive discrepancies shared by scalar and vector ASR."""

import itertools

import numpy as np


def sister_clade_contrast(tree, values_by_leaf):
    """Mean squared differences of observed sister-clade means.

    Every child-clade pair at every branching node contributes once. This is
    a topology discrepancy, not an independent-contrasts estimator: branch
    lengths and measurement errors are reproduced by the simulation model.
    Missing tips are omitted, identically in observed and replicated datasets.
    """
    summaries = {}
    contrasts: list[float] = []
    for node in tree.traverse("postorder"):
        if node.is_leaf:
            value = values_by_leaf.get(str(node.name))
            summaries[node] = (0.0, 0) if value is None else (float(value), 1)
            continue
        children = [summaries[child] for child in node.children]
        means = [total / count for total, count in children if count]
        contrasts.extend((a - b) ** 2 for a, b in itertools.combinations(means, 2))
        summaries[node] = (
            sum(total for total, _ in children),
            sum(count for _, count in children),
        )
    if not contrasts:
        raise ValueError("Phylogenetic predictive checks require two observed tips.")
    return float(np.mean(contrasts))


def predictive_summary(statistic, observed, replicated):
    """Summarize a scalar discrepancy without interpreting tail areas as tests."""
    values = np.asarray(replicated, dtype=float)
    if values.ndim != 1 or not len(values) or not np.isfinite(values).all():
        raise ValueError("Predictive discrepancies must be finite and nonempty.")
    lower = (1 + np.sum(values <= observed)) / (len(values) + 1)
    upper = (1 + np.sum(values >= observed)) / (len(values) + 1)
    return {
        "statistic": statistic,
        "observed": float(observed),
        "replicate_mean": float(np.mean(values)),
        "replicate_sd": float(np.std(values, ddof=1)) if len(values) > 1 else 0.0,
        "replicate_q025": float(np.quantile(values, 0.025)),
        "replicate_q975": float(np.quantile(values, 0.975)),
        "p_lower": float(lower),
        "p_upper": float(upper),
        "p_two_sided": min(1.0, 2 * min(float(lower), float(upper))),
        "num_simulations": len(values),
    }

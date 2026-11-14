"""Likelihood-compatible averaging of discrete and Gaussian reconstructions."""

import json
import math
from dataclasses import dataclass

import numpy as np
from scipy.special import ndtr


def comparison_average_table(context, table, criterion):
    """Use only ranked, non-duplicate models in one existing comparison group."""
    from nwkit.asr_compare import _custom_discrete_data, _single_discrete_data
    from nwkit.asr_tree_ensemble import (
        summarize_tree_ensemble,
        summarize_vector_tree_ensemble,
    )

    eligible = table.loc[
        table["rankable"].astype(str).str.lower().isin({"yes", "true", "1"})
        & table[criterion].notna()
        & table["equivalent_to"].fillna("").eq("")
    ]
    if eligible.empty:
        raise ValueError(
            "Model averaging requires at least one rankable model with a finite criterion."
        )
    weights = compatible_model_weights(
        eligible[criterion], eligible["comparison_group"]
    )
    states = None
    if context.trait_type == "discrete":
        data = (
            _custom_discrete_data(context)
            if set(eligible.model) == {"CUSTOM"}
            else _single_discrete_data(context)
        )
        states = data[0]
    posteriors = context.cache.get("averaging_posteriors", {})
    components = []
    for row in eligible.itertuples():
        posterior = posteriors[row.model_id]
        if row.model in {"HRM", "COVARION"}:
            if states is None:
                raise ValueError(
                    "Hidden-state averaging requires discrete state labels."
                )
            posterior = {
                node: value.reshape(-1, len(states)).sum(axis=0)
                for node, value in posterior.items()
            }
        components.append((context.tree, posterior))
    if context.trait_type == "continuous" and len(context.trait_columns) > 1:
        result = summarize_vector_tree_ensemble(
            context.tree, components, weights=weights, trait_names=context.trait_columns
        )
    else:
        result = summarize_tree_ensemble(
            context.tree, components, weights=weights, states=states
        )
    result = result.drop(
        columns=["mapping", "num_matched_trees", "matched_tree_weight"]
    )
    result = result.rename(
        columns={
            "num_trees": "num_models",
            "within_tree_variance": "within_model_variance",
            "between_tree_variance": "between_model_variance",
        }
    )
    result["comparison_group"] = eligible.comparison_group.iloc[0]
    result["criterion"] = criterion
    result["model_weights"] = json.dumps(
        dict(zip(eligible.model_id, weights.tolist(), strict=True)), sort_keys=True
    )
    return result


def normalized_weights(weights):
    result = np.asarray(weights, dtype=float)
    if (
        result.ndim != 1
        or not len(result)
        or not np.isfinite(result).all()
        or np.any(result < 0)
    ):
        raise ValueError(
            "Mixture weights must be a finite nonempty nonnegative vector."
        )
    maximum = float(np.max(result))
    if maximum == 0:
        raise ValueError("At least one mixture weight must be positive.")
    result = result / maximum
    return result / result.sum()


def compatible_model_weights(scores, groups):
    """IC weights for a single validated likelihood/root/data comparison group."""
    scores = np.asarray(scores, dtype=float)
    groups = tuple(groups)
    if (
        scores.ndim != 1
        or len(scores) != len(groups)
        or not len(scores)
        or not np.isfinite(scores).all()
    ):
        raise ValueError("Model averaging requires finite scores for every component.")
    if any(group != groups[0] for group in groups[1:]) or groups[0] in (None, ""):
        raise ValueError(
            "Model averaging cannot mix incompatible likelihood/root/data groups."
        )
    return normalized_weights(np.exp(-0.5 * (scores - scores.min())))


@dataclass(frozen=True)
class GaussianMixtureSummary:
    mean: float
    variance: float
    within_variance: float
    between_variance: float
    lower: float
    upper: float
    level: float


def gaussian_mixture_summary(means, variances, weights, *, level=0.95):
    """Exact mixture moments and numerical mixture quantiles, including atoms."""
    means = np.asarray(means, dtype=float)
    variances = np.asarray(variances, dtype=float)
    weights = normalized_weights(weights)
    if (
        means.shape != weights.shape
        or variances.shape != weights.shape
        or not np.isfinite(means).all()
        or not np.isfinite(variances).all()
        or np.any(variances < 0)
    ):
        raise ValueError(
            "Gaussian components require matching finite means and nonnegative variances."
        )
    if not math.isfinite(level) or not 0 < level < 1:
        raise ValueError("Mixture interval level must be between zero and one.")
    positive = weights > 0
    means, variances, weights = means[positive], variances[positive], weights[positive]
    mean = float(weights @ means)
    within = float(weights @ variances)
    between = float(weights @ ((means - mean) ** 2))
    standard_deviations = np.sqrt(variances)
    continuous = standard_deviations > 0

    def tail_probability(value, upper_tail):
        probabilities = np.asarray(
            value < means if upper_tail else value >= means, dtype=float
        )
        sign = -1 if upper_tail else 1
        probabilities[continuous] = ndtr(
            sign * (value - means[continuous]) / standard_deviations[continuous]
        )
        return float(weights @ probabilities)

    def quantile(probability, upper_tail=False):
        lower = float(np.min(means - 40 * standard_deviations))
        upper = float(np.max(means + 40 * standard_deviations))
        if lower == upper:
            return lower
        for _ in range(200):
            middle = lower / 2 + upper / 2
            if middle in (lower, upper):
                break
            tail = tail_probability(middle, upper_tail)
            if (tail <= probability) if upper_tail else (tail >= probability):
                upper = middle
            else:
                lower = middle
        return upper

    return GaussianMixtureSummary(
        mean,
        within + between,
        within,
        between,
        quantile((1 - level) / 2),
        quantile((1 - level) / 2, upper_tail=True),
        level,
    )


def average_state_probabilities(probabilities, weights):
    weights = normalized_weights(weights)
    values = np.asarray(probabilities, dtype=float)
    if (
        values.ndim != 2
        or len(values) != len(weights)
        or not values.shape[1]
        or not np.isfinite(values).all()
        or np.any(values < 0)
        or not np.allclose(values.sum(axis=1), 1, rtol=1e-10, atol=1e-12)
    ):
        raise ValueError(
            "Discrete mixture components must be normalized probabilities in a common state order."
        )
    return weights @ values

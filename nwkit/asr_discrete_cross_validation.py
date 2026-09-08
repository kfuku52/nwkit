"""Refitted tip/clade predictive probabilities for single-character CTMCs."""

import json
import math
from copy import copy

import numpy as np
import pandas as pd

from nwkit.asr_cross_validation import holdout_groups


def validate_discrete_cv(args, trait_type, model):
    if trait_type != "discrete":
        return
    output = getattr(args, "cross_validation_out", None)
    if getattr(args, "cross_validation_unit", None) and not output:
        raise ValueError("--cross-validation-unit requires --cross-validation-out.")
    if output and model in {
        "THRESHOLD",
        "MK-MIXTURE",
        "PAGEL-INDEPENDENT",
        "PAGEL-DEPENDENT",
    }:
        raise ValueError(
            "Discrete cross-validation requires a single-character CTMC model."
        )


def discrete_cross_validation(
    tree, states, observed, likelihoods, refit, *, mode="tip"
):
    """Remove held-out likelihoods before fitting, then score their observations.

    Log scores use P(observation|training), including known error/ambiguity.
    Brier scores are supplied only for one-hot likelihoods (known latent states).
    Clade output contains marginal tip scores, not a joint clade score.
    """
    from nwkit.asr import _is_informative_tip_likelihood

    informative = {
        name: True if _is_informative_tip_likelihood(vector) else None
        for name, vector in likelihoods.items()
    }
    leaves = {str(node.name): node for node in tree.leaves()}
    rows = []
    for fold, held_out in enumerate(holdout_groups(tree, informative, mode)):
        training_observed = dict(observed)
        training_likelihoods = {
            name: np.asarray(vector).copy() for name, vector in likelihoods.items()
        }
        for name in held_out:
            training_observed[name] = None
            training_likelihoods[name] = np.ones(len(states))
        fit = refit(training_observed, training_likelihoods)
        for name in held_out:
            posterior = np.asarray(fit["posterior_by_node"][leaves[name]])
            # Hidden classes are nuisance states, summed before observable scoring.
            probabilities = posterior.reshape(-1, len(states)).sum(axis=0)
            probabilities /= probabilities.sum()
            likelihood = np.asarray(likelihoods[name])
            probability = float(probabilities @ likelihood)
            exact = np.count_nonzero(likelihood) == 1 and float(likelihood.max()) == 1
            rows.append(
                {
                    "fold": fold,
                    "holdout": mode,
                    "name": name,
                    "num_training": sum(
                        _is_informative_tip_likelihood(value)
                        for value in training_likelihoods.values()
                    ),
                    "num_held_out": len(held_out),
                    "observed": observed.get(name) or "",
                    "observation_probability": probability,
                    "log_score": -math.inf
                    if probability == 0
                    else math.log(probability),
                    "brier_score": float(np.sum((probabilities - likelihood) ** 2))
                    if exact
                    else None,
                    "predicted_state": states[int(np.argmax(probabilities))],
                    "state_probabilities": json.dumps(
                        dict(zip(states, probabilities.tolist(), strict=True))
                    ),
                    "fit_status": fit.get("fit_status", ""),
                    "optimizer_success": fit.get("optimizer_success", ""),
                }
            )
    return pd.DataFrame(rows)


def write_discrete_cross_validation(
    tree, frame, states, observed, likelihoods, args, settings, fixed_matrix
):
    if not getattr(args, "cross_validation_out", None):
        return
    from nwkit.asr import _write_table
    from nwkit.asr_compare import (
        ComparisonCandidate,
        ComparisonContext,
        _fit_single_discrete,
    )

    training_args = copy(args)
    # Observations have already been converted to likelihoods. Reapplying a
    # file-backed error model here would reintroduce held-out data.
    training_args.tip_likelihoods = None
    training_args.misclassification_matrix = None
    candidate = ComparisonCandidate(
        settings.model, settings.root_prior, "cross_validation"
    )

    def refit(training_observed, training_likelihoods):
        context = ComparisonContext(
            tree, frame, "discrete", (args.state_column,), None, training_args
        )
        data = (states, training_observed, training_likelihoods)
        if settings.model == "CUSTOM":
            context.cache["custom_discrete_data"] = (*data, fixed_matrix)
        else:
            context.cache["single_discrete_data"] = data
        return _fit_single_discrete(context, candidate)

    table = discrete_cross_validation(
        tree,
        states,
        observed,
        likelihoods,
        refit,
        mode=getattr(args, "cross_validation_unit", None) or "tip",
    )
    _write_table(table, args.cross_validation_out)

"""Parametric-bootstrap prediction-error intervals for latent node values."""

import math

import numpy as np
import pandas as pd

from nwkit.gaussian_inference import condition_gaussian_tree
from nwkit.util import assign_branch_ids, get_node_class


def bootstrap_prediction_intervals(
    process,
    observed,
    simulate,
    refit_process,
    *,
    errors=None,
    num_simulations=100,
    seed=None,
    level=0.95,
):
    """Calibrate plug-in latent predictions using simulated reconstruction errors.

    ``simulate(seed)`` returns an observed mapping and joint latent GaussianSamples;
    ``refit_process(mapping)`` refits the original fixed/free parameter contract.
    This is a frequentist bootstrap, not parameter-posterior sampling. Failed
    refits are fatal so the interval is never based on selected successes.
    """
    if (
        isinstance(num_simulations, bool)
        or not isinstance(num_simulations, (int, np.integer))
        or num_simulations < 2
    ):
        raise ValueError(
            "Bootstrap prediction intervals require at least two simulations."
        )
    if not math.isfinite(level) or not 0 < level < 1:
        raise ValueError("Bootstrap interval level must be between zero and one.")
    original = condition_gaussian_tree(process, observed, standard_errors=errors)
    nodes = original.compiled_tree.nodes
    differences = np.empty((num_simulations, len(nodes)))
    for index, child_seed in enumerate(
        np.random.SeedSequence(seed).spawn(num_simulations)
    ):
        data, latent = simulate(int(child_seed.generate_state(1)[0]))
        if tuple(latent.nodes) != nodes or latent.values.shape != (1, len(nodes)):
            raise ValueError(
                "Bootstrap latent samples must follow the input tree order."
            )
        fitted = condition_gaussian_tree(
            refit_process(data), data, standard_errors=errors
        )
        differences[index] = latent.values[0] - np.array(
            [fitted.marginals[node].mean for node in nodes]
        )
    if not np.isfinite(differences).all():
        raise ValueError("Bootstrap reconstruction errors must be finite.")
    tail = (1 - level) / 2
    lower, upper = np.quantile(differences, [tail, 1 - tail], axis=0)
    ids = assign_branch_ids(process.tree)
    return pd.DataFrame(
        [
            {
                "branch_id": ids[node],
                "node_class": get_node_class(node),
                "name": "" if node.name in (None, "") else str(node.name),
                "mean": original.marginals[node].mean,
                "lower": original.marginals[node].mean + lower[index],
                "upper": original.marginals[node].mean + upper[index],
                "interval_level": level,
                "num_simulations": num_simulations,
                "method": "parametric_bootstrap_prediction_error",
                "error_bias": float(np.mean(differences[:, index])),
                "error_sd": float(np.std(differences[:, index], ddof=1)),
            }
            for index, node in enumerate(nodes)
        ]
    )

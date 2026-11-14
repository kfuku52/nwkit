"""Refitted, held-out predictive checks for scalar Gaussian ASR.

Scores are marginal observation densities, not a joint clade likelihood.
Evolutionary parameters are refitted without held-out trait values.
"""

import math

import pandas as pd
from scipy.stats import norm

from nwkit.compiled_tree import CompiledTree
from nwkit.gaussian_inference import condition_gaussian_tree


def holdout_groups(tree, observed, mode="tip"):
    """Partition observed tips into singleton or root-child clade folds."""
    compiled = CompiledTree.from_tree(tree)
    names = {
        name for name in compiled.leaf_index_by_name if observed.get(name) is not None
    }
    if mode == "tip":
        groups: list[tuple[str, ...]] = [(name,) for name in sorted(names)]
    elif mode == "clade":
        groups = [
            tuple(
                sorted(
                    str(leaf.name) for leaf in child.leaves() if str(leaf.name) in names
                )
            )
            for child in tree.children
        ]
        groups = [group for group in groups if group]
    else:
        raise ValueError("Cross-validation mode must be tip or clade.")
    if len(groups) < 2:
        raise ValueError("Cross-validation requires at least two nonempty folds.")
    return groups


def gaussian_cross_validation(
    tree, observed, refit_process, *, errors=None, mode="tip", level=0.95
):
    """Refit each fold and report log scores, PIT and predictive coverage.

    ``refit_process`` receives the training mapping with held-out entries set
    to None and must return a GaussianTreeProcess on the same input tree.
    Failed folds raise instead of producing selectively filtered scores.
    """
    if not 0 < level < 1:
        raise ValueError("Predictive interval level must be between zero and one.")
    compiled = CompiledTree.from_tree(tree)
    quantile = float(norm.ppf((1 + level) / 2))
    rows = []
    for fold, held_out in enumerate(holdout_groups(tree, observed, mode)):
        training = dict(observed)
        training.update(dict.fromkeys(held_out))
        process = refit_process(training)
        posterior = condition_gaussian_tree(process, training, standard_errors=errors)
        for name in held_out:
            node = compiled.nodes[compiled.leaf_index_by_name[name]]
            marginal = posterior.marginals[node]
            error = 0.0 if errors is None else float(errors[name])
            variance = marginal.variance + error * error
            if not math.isfinite(variance) or variance <= 0:
                raise ValueError(
                    f"Predictive variance for held-out tip '{name}' must be positive."
                )
            value = float(observed[name])
            sd = math.sqrt(variance)
            z = (value - marginal.mean) / sd
            rows.append(
                {
                    "fold": fold,
                    "holdout": mode,
                    "name": name,
                    "num_training": sum(
                        value is not None for value in training.values()
                    ),
                    "num_held_out": len(held_out),
                    "observed": value,
                    "predicted_mean": marginal.mean,
                    "predicted_sd": sd,
                    "log_score": float(norm.logpdf(z) - math.log(sd)),
                    "pit": float(norm.cdf(z)),
                    "interval_level": level,
                    "lower": marginal.mean - quantile * sd,
                    "upper": marginal.mean + quantile * sd,
                    "covered": abs(z) <= quantile,
                }
            )
    return pd.DataFrame(rows)

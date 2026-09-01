"""Continuous ASR summaries and NHX serialization."""

import math
from copy import copy

import pandas as pd
from scipy.special import erfinv

from nwkit.util import assign_branch_ids, get_node_class, write_tree


def _summary(marginal, ci_level):
    sd = math.sqrt(marginal.variance)
    # This centered form retains both tiny levels, where 1-level rounds to one,
    # and the largest representable level below one.
    width = math.sqrt(2.0) * float(erfinv(ci_level)) * sd
    lower, upper = marginal.mean - width, marginal.mean + width
    if not math.isfinite(lower) or not math.isfinite(upper):
        raise ValueError(
            "An ancestral interval exceeds floating-point range; rescale the trait units."
        )
    return {
        "mean": marginal.mean,
        "variance": marginal.variance,
        "sd": sd,
        "ci_lower": lower,
        "ci_upper": upper,
        "ci_level": ci_level,
    }


def continuous_output_table(
    tree, selected_nodes, observed, errors, posterior, *, trait, ci_level
):
    ids = assign_branch_ids(tree)
    rows = []
    for node in selected_nodes:
        value = observed.get(str(node.name)) if node.is_leaf else None
        error = 0.0 if errors is None else errors.get(str(node.name))
        rows.append(
            {
                "branch_id": ids[node],
                "parent": -1 if node.is_root else ids[node.up],
                "node_class": get_node_class(node),
                "name": "" if node.name in (None, "") else str(node.name),
                "trait": trait,
                "observed_value": "" if value is None else value,
                "observed_se": "" if value is None else error,
                "is_imputed": bool(node.is_leaf and value is None),
                **_summary(posterior[node], ci_level),
            }
        )
    return pd.DataFrame(
        rows,
        columns=[
            "branch_id",
            "parent",
            "node_class",
            "name",
            "trait",
            "observed_value",
            "observed_se",
            "is_imputed",
            "mean",
            "variance",
            "sd",
            "ci_lower",
            "ci_upper",
            "ci_level",
        ],
    )


def continuous_model_table(fit, args, ci_level):
    return pd.DataFrame(
        [
            {
                "trait_type": "continuous",
                "trait_type_requested": getattr(args, "trait_type", "auto"),
                "trait": args.state_column,
                "model": "BM",
                "root_prior": "flat",
                "sigma2": fit.sigma2,
                "sigma2_estimated": fit.sigma2_estimated,
                "estimation_method": "REML" if fit.sigma2_estimated else "fixed",
                "restricted_log_likelihood": fit.restricted_log_likelihood,
                "likelihood_kind": "flat_root_integrated",
                "num_observed": fit.num_observed,
                "num_effective_observations": fit.num_effective_observations,
                "residual_df": fit.residual_df,
                "fit_status": fit.fit_status,
                "standard_error_column": getattr(args, "standard_error_column", None)
                or "",
                "ci_level": ci_level,
                "interval_kind": "conditional_on_sigma2",
                "parameter_uncertainty_included": False,
                "tree_uncertainty_included": False,
            }
        ]
    )


def write_continuous_tree(tree, observed, errors, posterior, args, ci_level):
    if getattr(args, "tree_out", None) in (None, ""):
        return
    for node in tree.traverse():
        summary = _summary(posterior[node], ci_level)
        node.add_props(**{"asr_" + key: value for key, value in summary.items()})
        node.add_props(
            asr_trait_type="continuous", asr_interval_kind="conditional_on_sigma2"
        )
        if node.is_leaf:
            value = observed.get(str(node.name))
            error = 0.0 if errors is None else errors.get(str(node.name))
            node.add_props(
                asr_observed_value="" if value is None else value,
                asr_observed_se="" if value is None else error,
                asr_is_imputed="yes" if value is None else "no",
            )
    props = ["asr_mean", "asr_trait_type"]
    if args.tree_annotation != "mean":
        props += [
            "asr_variance",
            "asr_sd",
            "asr_ci_lower",
            "asr_ci_upper",
            "asr_ci_level",
            "asr_interval_kind",
        ]
    if args.tree_annotation == "all":
        props += ["asr_observed_value", "asr_observed_se", "asr_is_imputed"]
    output_args = copy(args)
    output_args.outfile = args.tree_out
    write_tree(
        tree,
        output_args,
        format=getattr(args, "tree_outformat", "auto"),
        quiet=True,
        props=props,
    )

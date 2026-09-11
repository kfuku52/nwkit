"""Tree uncertainty summaries aligned by descendant clades, never branch IDs."""

import json
from copy import copy
from typing import Any

import pandas as pd

from nwkit.asr_averaging import (
    average_state_probabilities,
    gaussian_mixture_summary,
    normalized_weights,
)
from nwkit.util import assign_branch_ids, get_node_class


def validate_tree_ensemble_options(args, settings=None):
    source = getattr(args, "tree_ensemble", None)
    output = getattr(args, "tree_ensemble_out", None)
    if bool(source) != bool(output):
        raise ValueError("--tree-ensemble and --tree-ensemble-out require each other.")
    if not source:
        if getattr(args, "tree_ensemble_weights", None) or getattr(
            args, "tree_ensemble_mapping", None
        ):
            raise ValueError("Tree ensemble options require --tree-ensemble.")
        return
    if settings is None:
        return
    if settings.model in {
        "THRESHOLD",
        "MK-MIXTURE",
        "PAGEL-INDEPENDENT",
        "PAGEL-DEPENDENT",
    }:
        raise ValueError(
            "Tree ensemble summaries require a single Gaussian or CTMC trait."
        )
    if getattr(args, "regime_map", None):
        raise ValueError(
            "Tree ensembles cannot reuse a branch-ID --regime-map across trees."
        )


def _fit_ensemble_tree(tree, trait_df, args, settings):
    from nwkit.asr import (
        _continuous_observations,
        _fit_continuous_model,
    )
    from nwkit.asr_compare import (
        ComparisonCandidate,
        ComparisonContext,
        _custom_discrete_data,
        _fit_single_discrete,
        _single_discrete_data,
    )

    if settings.trait_type == "continuous":
        columns, observed, errors = _continuous_observations(
            tree, trait_df, args, settings
        )
        posterior, _ = _fit_continuous_model(
            tree, observed, errors, columns, args, settings, None
        )
        return posterior, None
    context = ComparisonContext(
        tree, trait_df, "discrete", (args.state_column,), None, args
    )
    data = (
        _custom_discrete_data(context)
        if settings.model == "CUSTOM"
        else _single_discrete_data(context)
    )
    states = data[0]
    fit = _fit_single_discrete(
        context, ComparisonCandidate(settings.model, settings.root_prior, "ensemble")
    )
    posterior = fit["posterior_by_node"]
    if settings.model in {"HRM", "COVARION"}:
        posterior = {
            node: value.reshape(-1, len(states)).sum(axis=0)
            for node, value in posterior.items()
        }
    return posterior, states


def write_tree_ensemble(reference, trait_df, args, settings):
    """Refit each sampled tree independently; preserve reference-only outputs."""
    from nwkit.asr import _validate_tree_for_asr, _write_table
    from nwkit.util import read_trees

    if not getattr(args, "tree_ensemble", None):
        return
    trees = read_trees(
        args.tree_ensemble,
        args.format,
        args.quoted_node_names,
        rooted=getattr(args, "input_rooted", "auto"),
    )
    mapping = getattr(args, "tree_ensemble_mapping", None) or "clade"
    for tree in trees:
        _validate_tree_for_asr(tree)
        align_ensemble_nodes(reference, tree, mapping=mapping)
    raw_weights = getattr(args, "tree_ensemble_weights", None)
    weights = normalized_weights(
        [float(value) for value in raw_weights.split(",")]
        if raw_weights
        else [1] * len(trees)
    )
    if len(weights) != len(trees):
        raise ValueError("Tree ensemble weights must match the tree count.")
    fit_args = copy(args)
    fit_args.profile_ci_level = None
    fitted = [_fit_ensemble_tree(tree, trait_df, fit_args, settings) for tree in trees]
    states = fitted[0][1]
    if any(labels != states for _, labels in fitted):
        raise ValueError("Ensemble fits have inconsistent state ordering.")
    summarize = (
        summarize_vector_tree_ensemble
        if settings.model.startswith("MV-")
        else summarize_tree_ensemble
    )
    options = (
        {"trait_names": tuple(args.state_column.split(","))}
        if settings.model.startswith("MV-")
        else {"states": states}
    )
    table = summarize(
        reference,
        [(tree, result[0]) for tree, result in zip(trees, fitted, strict=True)],
        weights=weights,
        mapping=mapping,
        **options,
        level=settings.ci_level or 0.95,
    )
    _write_table(table, args.tree_ensemble_out)


def summarize_vector_tree_ensemble(
    reference, tree_posteriors, *, trait_names, **options
):
    """Marginal mixtures for each trait; retain both variance sources."""
    from nwkit.continuous_asr import GaussianMarginal

    tables = []
    for index, name in enumerate(trait_names):
        scalar = [
            (
                tree,
                {
                    node: GaussianMarginal(
                        float(value.mean[index]), float(value.covariance[index, index])
                    )
                    for node, value in posterior.items()
                },
            )
            for tree, posterior in tree_posteriors
        ]
        table = summarize_tree_ensemble(reference, scalar, **options)
        table["trait"] = name.strip()
        tables.append(table)
    return pd.concat(tables, ignore_index=True)


def descendant_clades(tree):
    """Map nodes to immutable descendant-tip sets; reject duplicate tip labels."""
    names = list(tree.leaf_names())
    if len(names) != len(set(names)) or any(name in (None, "") for name in names):
        raise ValueError("Tree ensembles require unique nonempty tip names.")
    clades: dict[Any, frozenset[str]] = {}
    for node in tree.traverse("postorder"):
        clades[node] = (
            frozenset([str(node.name)])
            if node.is_leaf
            else frozenset().union(*(clades[child] for child in node.children))
        )
    return clades


def align_ensemble_nodes(reference, tree, *, mapping="clade"):
    """Match exact clades or the MRCA of each reference node's descendant tips.

    Unary duplicate clades map to their deepest node, i.e. their MRCA. The
    reference root maps to the sampled root, retaining the root estimand.
    """
    if mapping not in {"clade", "mrca"}:
        raise ValueError("Tree ensemble mapping must be clade or mrca.")
    reference_clades = descendant_clades(reference)
    sampled_clades = descendant_clades(tree)
    if reference_clades[reference] != sampled_clades[tree]:
        raise ValueError(
            "Every ensemble tree must contain exactly the reference tip set."
        )
    exact: dict[frozenset[str], Any] = {}
    for node in tree.traverse("postorder"):
        exact.setdefault(sampled_clades[node], node)
    result = {}
    for node, clade in reference_clades.items():
        if node is reference:
            result[node] = tree
        elif clade in exact:
            result[node] = exact[clade]
        elif mapping == "mrca":
            result[node] = tree.common_ancestor([str(name) for name in sorted(clade)])
    return result


def summarize_tree_ensemble(
    reference,
    tree_posteriors,
    *,
    weights=None,
    mapping="clade",
    states=None,
    level=0.95,
):
    """Average scalar Gaussian or discrete posteriors over an explicit tree sample.

    Missing clades are excluded, with their original weight reported as missing
    support. A clade-conditioned summary is not an unconditional tree posterior.
    """
    tree_posteriors = tuple(tree_posteriors)
    weights = normalized_weights(
        [1] * len(tree_posteriors) if weights is None else weights
    )
    if len(weights) != len(tree_posteriors):
        raise ValueError("Tree ensemble weights must match the tree count.")
    mappings = [
        align_ensemble_nodes(reference, tree, mapping=mapping)
        for tree, _ in tree_posteriors
    ]
    clades = descendant_clades(reference)
    ids = assign_branch_ids(reference)
    rows = []
    for node in reference.traverse():
        included = [
            i
            for i, aligned in enumerate(mappings)
            if node in aligned and weights[i] > 0
        ]
        support = float(sum(weights[i] for i in included))
        row = {
            "branch_id": ids[node],
            "parent": -1 if node.is_root else ids[node.up],
            "node_class": get_node_class(node),
            "name": str(node.name or ""),
            "descendant_tips": json.dumps(sorted(clades[node]), ensure_ascii=False),
            "mapping": mapping,
            "num_trees": len(tree_posteriors),
            "num_matched_trees": len(included),
            "matched_tree_weight": support,
        }
        if included:
            selected = [tree_posteriors[i][1][mappings[i][node]] for i in included]
            selected_weights = [weights[i] for i in included]
            if states is None:
                summary = gaussian_mixture_summary(
                    [value.mean for value in selected],
                    [value.variance for value in selected],
                    selected_weights,
                    level=level,
                )
                row.update(
                    mean=summary.mean,
                    variance=summary.variance,
                    within_tree_variance=summary.within_variance,
                    between_tree_variance=summary.between_variance,
                    lower=summary.lower,
                    upper=summary.upper,
                    level=level,
                )
            else:
                from nwkit.asr import _safe_column_state

                probabilities = average_state_probabilities(selected, selected_weights)
                if len(probabilities) != len(states):
                    raise ValueError(
                        "Ensemble state labels must match posterior dimension."
                    )
                row.update(
                    {
                        "p_" + _safe_column_state(state): float(value)
                        for state, value in zip(states, probabilities, strict=True)
                    }
                )
        rows.append(row)
    return pd.DataFrame(rows)

"""Diagnostics for delayed trait origins in reconciled speciation contrasts."""

from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd

from nwkit.asr import _simulate_stochastic_maps, compute_mk_marginals
from nwkit.clade_index import CladeIndex
from nwkit.model_matrix import CategoricalObservation
from nwkit.util import assign_branch_ids

ORIGIN_DIAGNOSTIC_COLUMNS = [
    "trait",
    "origin_id",
    "branch_id",
    "parent",
    "branch_clade_id",
    "descendant_taxa",
    "descendant_species_event_ids",
    "from_state",
    "to_state",
    "total_count",
    "mean_count",
    "posterior_frequency",
    "num_simulations",
    "mk_model",
    "mk_rates",
    "log_likelihood",
    "credible",
]


def _categorical_states(values_by_leaf):
    states: set[str] = set()
    for value in values_by_leaf.values():
        if isinstance(value, CategoricalObservation):
            states.update(str(state) for state in value.probabilities)
        else:
            states.add(str(value))
    if len(states) < 2:
        raise ValueError(
            "Categorical origin diagnostics require at least two observed states."
        )
    return sorted(states)


def _tip_likelihoods(tree, states, values_by_leaf):
    observed = {}
    likelihoods = {}
    state_index = {state: index for index, state in enumerate(states)}
    for leaf in tree.leaves():
        name = str(leaf.name)
        value = values_by_leaf[name]
        probabilities = np.zeros(len(states), dtype=float)
        if isinstance(value, CategoricalObservation):
            for state, probability in value.probabilities.items():
                probabilities[state_index[str(state)]] = float(probability)
            observed[name] = states[int(np.argmax(probabilities))]
        else:
            state = str(value)
            probabilities[state_index[state]] = 1.0
            observed[name] = state
        total = float(probabilities.sum())
        if total <= 0.0 or not np.isfinite(probabilities).all():
            raise ValueError(
                "Categorical origin likelihoods must be finite with positive mass."
            )
        likelihoods[name] = probabilities / total
    return observed, likelihoods


def _descendant_event_ids(node, clades, available_event_ids):
    return sorted(
        clades.clade_id_for_node(descendant)
        for descendant in node.traverse()
        if clades.clade_id_for_node(descendant) in available_event_ids
    )


def _origin_identifier(trait, record):
    return "{}:b{}:{}->{}".format(
        trait,
        record["branch_id"],
        record["from_state"],
        record["to_state"],
    )


def _origin_label(trait, record):
    return "{} {}->{} branch {}".format(
        trait,
        record["from_state"],
        record["to_state"],
        record["branch_id"],
    )


def _format_mk_rates(rates):
    return ",".join("{:.12g}".format(float(rate)) for rate in rates)


def _is_credible_origin(record, minimum_posterior):
    return float(record["posterior_frequency"]) >= minimum_posterior


def build_categorical_origin_diagnostics(
    species_tree,
    values_by_trait: dict[str, dict[str, Any]],
    categorical_traits,
    species_event_ids,
    *,
    num_simulations=200,
    minimum_posterior=0.5,
    seed=1,
    threads=1,
):
    """Map categorical transitions and define descendant-event sensitivity sets."""
    if not 0.0 <= minimum_posterior <= 1.0:
        raise ValueError("Origin minimum posterior must lie in [0, 1].")
    if not isinstance(num_simulations, int) or num_simulations <= 0:
        raise ValueError("Origin map replicates must be a positive integer.")
    clades = CladeIndex(species_tree)
    branch_ids = assign_branch_ids(species_tree)
    node_by_branch = {branch_id: node for node, branch_id in branch_ids.items()}
    available_event_ids = set(str(value) for value in species_event_ids)
    rows = []
    omissions = []
    for trait_index, trait in enumerate(categorical_traits):
        values_by_leaf = values_by_trait[trait]
        states = _categorical_states(values_by_leaf)
        observed, likelihoods = _tip_likelihoods(species_tree, states, values_by_leaf)
        _, fit = compute_mk_marginals(
            species_tree,
            states,
            observed,
            likelihoods,
            model="ER",
            root_prior_mode="equal",
        )
        mapped = _simulate_stochastic_maps(
            species_tree,
            states,
            fit,
            num_simulations,
            seed=seed + trait_index,
            threads=threads,
        )
        rates = _format_mk_rates(fit["rates"])
        for record in mapped.to_dict("records"):
            node = node_by_branch[int(record["branch_id"])]
            branch_clade_id = clades.clade_id_for_node(node)
            event_ids = _descendant_event_ids(node, clades, available_event_ids)
            origin_id = _origin_identifier(trait, record)
            credible = _is_credible_origin(record, minimum_posterior)
            rows.append(
                {
                    **record,
                    "trait": trait,
                    "origin_id": origin_id,
                    "branch_clade_id": branch_clade_id,
                    "descendant_taxa": ",".join(
                        sorted(str(leaf.name) for leaf in node.leaves())
                    ),
                    "descendant_species_event_ids": ",".join(event_ids),
                    "mk_model": "ER",
                    "mk_rates": rates,
                    "log_likelihood": float(fit["log_likelihood"]),
                    "credible": "yes" if credible else "no",
                }
            )
            if credible and event_ids:
                omissions.append(
                    {
                        "analysis_type": "trait-origin-leave-one-out",
                        "group_id": origin_id,
                        "group_label": _origin_label(trait, record),
                        "event_ids": set(event_ids),
                    }
                )
    return (
        pd.DataFrame(rows, columns=ORIGIN_DIAGNOSTIC_COLUMNS),
        omissions,
    )

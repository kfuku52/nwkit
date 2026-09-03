"""Pagel independent/dependent evolution for two binary traits."""

import json
from dataclasses import dataclass

import numpy as np

from nwkit.discrete_asr_models import RateDesign, rate_design_from_edges

PAGEL_MODELS = frozenset({"PAGEL-INDEPENDENT", "PAGEL-DEPENDENT"})


class PartialPagelObservation(str):
    """Displayable joint state whose likelihood has one missing component."""

    is_partial_missing = True


@dataclass(frozen=True, slots=True)
class PagelData:
    trait_columns: tuple[str, str]
    trait_states: tuple[tuple[str, str], tuple[str, str]]
    joint_states: tuple[str, ...]
    observed_state_by_leaf: dict[str, str | None]
    likelihood_by_leaf: dict[str, np.ndarray]


def _explicit_trait_states(value):
    if value in (None, ""):
        return None
    groups = str(value).split(";")
    if len(groups) != 2:
        raise ValueError(
            "Pagel models require --states 'A0,A1;B0,B1' when binary state "
            "spaces are supplied explicitly."
        )
    result = []
    for trait_index, group in enumerate(groups, start=1):
        states = tuple(item.strip() for item in group.split(","))
        if len(states) != 2 or any(state == "" for state in states):
            raise ValueError(
                f"Pagel trait {trait_index} requires exactly two non-empty states."
            )
        if len(set(states)) != 2:
            raise ValueError(f"Pagel trait {trait_index} contains duplicated states.")
        result.append(states)
    return tuple(result)


def _joint_state_label(first, second):
    return json.dumps(
        [str(first), str(second)], ensure_ascii=False, separators=(",", ":")
    )


def prepare_pagel_data(
    trait_path,
    trait_columns,
    tree_leaf_names,
    *,
    states_arg=None,
    missing_values_arg=None,
    ambiguous_separator="|",
    unmatched="warn",
    trait_df=None,
):
    """Read two binary characters and form their four-state tip likelihoods."""

    from nwkit.asr import _read_tip_states

    columns = tuple(str(column) for column in trait_columns)
    if len(columns) != 2 or columns[0] == columns[1]:
        raise ValueError(
            "Pagel models require exactly two distinct comma-separated "
            "--state-column values."
        )
    explicit = _explicit_trait_states(states_arg)
    trait_states = []
    observed_by_trait = []
    likelihoods_by_trait = []
    for index, column in enumerate(columns):
        states, observed, likelihoods = _read_tip_states(
            trait_path,
            column,
            tree_leaf_names,
            states_arg=None if explicit is None else explicit[index],
            missing_values_arg=missing_values_arg,
            ambiguous_separator=ambiguous_separator,
            unmatched=unmatched,
            trait_df=trait_df,
            state_source=(
                "--states" if explicit is not None else f"trait column '{column}'"
            ),
        )
        if len(states) != 2:
            raise ValueError(
                f"Pagel trait column '{column}' must have exactly two model states; "
                "use --states 'A0,A1;B0,B1' when an unobserved state is required."
            )
        trait_states.append((states[0], states[1]))
        observed_by_trait.append(observed)
        likelihoods_by_trait.append(likelihoods)

    joint_states = tuple(
        _joint_state_label(first, second)
        for first in trait_states[0]
        for second in trait_states[1]
    )
    observed_joint: dict[str, str | None] = {}
    likelihood_joint: dict[str, np.ndarray] = {}
    for leaf_name in tree_leaf_names:
        first_observed = observed_by_trait[0].get(leaf_name)
        second_observed = observed_by_trait[1].get(leaf_name)
        if first_observed is None and second_observed is None:
            observed_joint[leaf_name] = None
        else:
            label = json.dumps(
                [
                    None if first_observed is None else str(first_observed),
                    None if second_observed is None else str(second_observed),
                ],
                ensure_ascii=False,
                separators=(",", ":"),
            )
            observed_joint[leaf_name] = (
                PartialPagelObservation(label)
                if first_observed is None or second_observed is None
                else label
            )
        likelihood_joint[leaf_name] = np.kron(
            likelihoods_by_trait[0][leaf_name],
            likelihoods_by_trait[1][leaf_name],
        )
    return PagelData(
        (columns[0], columns[1]),
        (trait_states[0], trait_states[1]),
        joint_states,
        observed_joint,
        likelihood_joint,
    )


def pagel_rate_design(data, model):
    """Return Pagel's four- or eight-rate direct-edge partition."""

    if model not in PAGEL_MODELS:
        raise ValueError(f"Unsupported Pagel model: {model}.")
    edges = []
    for first in range(2):
        for second in range(2):
            source = data.joint_states[2 * first + second]
            target_first = data.joint_states[2 * (1 - first) + second]
            target_second = data.joint_states[2 * first + (1 - second)]
            if model == "PAGEL-INDEPENDENT":
                first_class = f"trait1_{first}_to_{1 - first}"
                second_class = f"trait2_{second}_to_{1 - second}"
            else:
                first_class = f"trait1_{first}_to_{1 - first}_when_trait2_{second}"
                second_class = f"trait2_{second}_to_{1 - second}_when_trait1_{first}"
            edges.append((source, target_first, first_class))
            edges.append((source, target_second, second_class))
    return rate_design_from_edges(
        data.joint_states,
        edges,
        source=model,
    )


def compute_pagel_marginals(
    tree,
    data,
    *,
    model,
    rate=None,
    root_prior_mode="equal",
    rate_bounds=None,
):
    """Fit one Pagel binary-pair model and return joint ancestral marginals."""

    from nwkit.asr import compute_mk_marginals

    design: RateDesign = pagel_rate_design(data, model)
    posterior, fit = compute_mk_marginals(
        tree,
        data.joint_states,
        data.observed_state_by_leaf,
        data.likelihood_by_leaf,
        model="MK-DESIGN",
        rate=rate,
        root_prior_mode=root_prior_mode,
        rate_bounds=rate_bounds,
        rate_design=design,
    )
    fit.update(
        {
            "model": model,
            "likelihood_kind": "pagel_joint_ml",
            "pagel_trait_columns": data.trait_columns,
            "pagel_trait_states": data.trait_states,
            "pagel_joint_states_json": json.dumps(
                data.joint_states, ensure_ascii=False, separators=(",", ":")
            ),
            "rate_class_names": design.class_names,
            "rate_design_source": design.source,
        }
    )
    return posterior, fit

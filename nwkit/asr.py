import math
import multiprocessing
import sys
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from copy import copy
from typing import Any

import numpy as np
import pandas as pd
from scipy.linalg import expm
from scipy.optimize import minimize
from scipy.special import logsumexp
from scipy.stats import poisson

from nwkit.asr_input import (
    AsrSettings,
    asr_trait_columns,
    continuous_tip_values,
    continuous_tip_vectors,
    effective_asr_args,
    read_asr_table,
    resolve_trait_type,
)
from nwkit.asr_models import model_names
from nwkit.discrete_asr_models import (
    FREQUENCY_RATIO_BOUNDS,
    read_rate_matrix,
    read_transition_graph,
    stationary_distribution,
    validate_rate_matrix,
)
from nwkit.discrete_asr_models import build_rate_matrix as _build_structured_rate_matrix
from nwkit.discrete_asr_models import initial_parameters as _initial_model_parameters
from nwkit.discrete_asr_models import parameter_kinds as _structured_parameter_kinds
from nwkit.discrete_asr_models import parameter_labels as _structured_parameter_labels
from nwkit.rooting_state import require_rooted
from nwkit.util import (
    assign_branch_ids,
    get_node_class,
    is_missing_table_value,
    parse_table_missing_values,
    read_tree,
    validate_distinct_output_paths,
    validate_unique_named_leaves,
    write_tree,
)

DEFAULT_RATE_BOUNDS = (10**-9, 10**3)
DEFAULT_AMBIGUOUS_SEPARATOR = "|"
DEFAULT_TARGET = "all"
SUPPORTED_MODELS = model_names("discrete")
_MAX_UNIFORMIZATION_TERMS = 2_000_000
_MAX_CACHED_BACKWARD_BYTES = 256 * 1024


def _parse_comma_list(value, option_name):
    if value in ["", None]:
        return list()
    items = [item.strip() for item in str(value).split(",")]
    if any(item == "" for item in items):
        raise ValueError("'{}' contains an empty item.".format(option_name))
    return items


def _parse_missing_values(value):
    return parse_table_missing_values(value)


def _validate_states(states):
    if len(states) != len(set(states)):
        raise ValueError("'--states' contains duplicated states.")
    return states


def _parse_states(value):
    return _validate_states(_parse_comma_list(value, "--states"))


def _parse_state_argument(value):
    if isinstance(value, (list, tuple)):
        return _validate_states([str(state) for state in value])
    return _parse_states(value)


def _parse_rate_bounds(value):
    if value in ["", None]:
        return DEFAULT_RATE_BOUNDS
    items = _parse_comma_list(value, "--rate-bounds")
    if len(items) != 2:
        raise ValueError(
            "'--rate-bounds' must contain exactly two comma-separated values."
        )
    try:
        lower = float(items[0])
        upper = float(items[1])
    except ValueError as exc:
        raise ValueError("'--rate-bounds' values must be numeric.") from exc
    if (
        (not math.isfinite(lower))
        or (not math.isfinite(upper))
        or lower <= 0.0
        or upper <= 0.0
    ):
        raise ValueError("'--rate-bounds' values must be positive finite numbers.")
    if lower >= upper:
        raise ValueError(
            "'--rate-bounds' lower bound must be smaller than the upper bound."
        )
    return lower, upper


def _parse_targets(value):
    if value in ["", None]:
        value = DEFAULT_TARGET
    targets = set()
    aliases = {
        "tip": "leaf",
        "tips": "leaf",
        "leaves": "leaf",
        "missing_tip": "missing-leaf",
        "missing-tip": "missing-leaf",
        "missing_leaf": "missing-leaf",
    }
    valid_targets = {"all", "intnode", "leaf", "missing-leaf"}
    for raw_target in _parse_comma_list(value, "--target"):
        target = aliases.get(raw_target, raw_target)
        if target not in valid_targets:
            raise ValueError("Unknown '--target': {}".format(raw_target))
        targets.add(target)
    if "all" in targets:
        return {"all"}
    return targets


def _is_missing_trait_value(value, missing_values):
    return is_missing_table_value(value, missing_values)


def _split_state_value(value, ambiguous_separator):
    state_text = str(value).strip()
    if ambiguous_separator in ["", None]:
        return [state_text]
    states = [state.strip() for state in state_text.split(ambiguous_separator)]
    if any(state == "" for state in states):
        raise ValueError(
            "Ambiguous state value contains an empty state: {}".format(state_text)
        )
    if len(states) != len(set(states)):
        raise ValueError(
            "Ambiguous state value contains duplicated states: {}".format(state_text)
        )
    return states


def _read_tip_states(
    trait_path,
    state_column,
    tree_leaf_names,
    states_arg=None,
    missing_values_arg=None,
    ambiguous_separator=DEFAULT_AMBIGUOUS_SEPARATOR,
    unmatched="warn",
    *,
    trait_df=None,
    state_source="--states",
):
    if trait_df is None:
        trait_df = read_asr_table(
            trait_path,
            state_column,
            tree_leaf_names,
            unmatched=unmatched,
            missing_values=missing_values_arg,
        )

    missing_values = _parse_missing_values(missing_values_arg)
    observed_state_by_leaf_input: dict[str, str | None] = {}
    state_set_by_leaf: dict[str, tuple[str, ...] | None] = {}
    observed_state_sets: list[list[str]] = []
    for _, row in trait_df.iterrows():
        leaf_name = str(row["leaf_name"])
        raw_state = row[state_column]
        if _is_missing_trait_value(raw_state, missing_values):
            observed_state_by_leaf_input[leaf_name] = None
            state_set_by_leaf[leaf_name] = None
            continue
        state_set = _split_state_value(raw_state, ambiguous_separator)
        state_set_by_leaf[leaf_name] = tuple(state_set)
        separator = (
            "" if ambiguous_separator in ["", None] else str(ambiguous_separator)
        )
        observed_state_by_leaf_input[leaf_name] = separator.join(state_set)
        observed_state_sets.append(state_set)

    states = _parse_state_argument(states_arg)
    observed_states = [
        state for state_set in observed_state_sets for state in state_set
    ]
    if states:
        unknown_states = sorted(set(observed_states) - set(states))
        if unknown_states:
            raise ValueError(
                "State(s) in '--trait' are absent from the state space defined by "
                "'{}': {}".format(state_source, ", ".join(unknown_states))
            )
    else:
        seen_states = set()
        for state in observed_states:
            if state not in seen_states:
                states.append(state)
                seen_states.add(state)
    if len(states) == 0:
        raise ValueError(
            "At least one observed or explicitly listed state is required."
        )

    state_to_index = {state: index for index, state in enumerate(states)}
    observed_state_by_leaf = {
        leaf_name: observed_state_by_leaf_input.get(leaf_name)
        for leaf_name in tree_leaf_names
    }
    likelihood_by_leaf = dict()
    for leaf_name in tree_leaf_names:
        state_set = state_set_by_leaf.get(leaf_name)
        likelihood = np.ones(len(states), dtype=float)
        if state_set is not None:
            likelihood = np.zeros(len(states), dtype=float)
            for state in state_set:
                likelihood[state_to_index[state]] = 1.0
        likelihood_by_leaf[leaf_name] = likelihood
    return states, observed_state_by_leaf, likelihood_by_leaf


def _validate_tree_for_asr(tree):
    validate_unique_named_leaves(tree, option_name="--infile", context=" for 'asr'")
    require_rooted(tree, "ASR requires a rooted tree.")
    for node in tree.traverse():
        if node.is_root:
            continue
        if node.dist is None:
            raise ValueError("ASR requires branch lengths for all non-root nodes.")
        try:
            dist = float(node.dist)
        except (TypeError, ValueError) as exc:
            raise ValueError("Branch lengths must be numeric.") from exc
        if (not math.isfinite(dist)) or dist < 0.0:
            raise ValueError("Branch lengths must be non-negative finite numbers.")


def _er_transition_matrix(branch_length, rate, num_states):
    if num_states == 1:
        return np.ones((1, 1), dtype=float)
    if rate < 0.0:
        raise ValueError("Mk transition rate must be non-negative.")
    decay = math.exp(-float(num_states) * float(rate) * float(branch_length))
    off_diagonal = (1.0 - decay) / float(num_states)
    matrix = np.full((num_states, num_states), off_diagonal, dtype=float)
    diagonal = off_diagonal + decay
    np.fill_diagonal(matrix, diagonal)
    return matrix


def _get_root_prior(root_prior, states, observed_state_by_leaf, likelihood_by_leaf):
    num_states = len(states)
    if root_prior == "equal":
        return np.full(num_states, 1.0 / float(num_states), dtype=float)
    if root_prior == "empirical":
        counts = np.zeros(num_states, dtype=float)
        for leaf_name, observed_state in observed_state_by_leaf.items():
            if observed_state is None:
                continue
            likelihood = likelihood_by_leaf[leaf_name]
            counts += likelihood / likelihood.sum()
        if counts.sum() == 0.0:
            return np.full(num_states, 1.0 / float(num_states), dtype=float)
        return counts / counts.sum()
    raise ValueError("Unsupported '--root-prior': {}".format(root_prior))


def _root_prior_for_matrix(
    root_prior, states, observed_state_by_leaf, likelihood_by_leaf, rate_matrix
):
    if root_prior == "stationary":
        return stationary_distribution(rate_matrix)
    return _get_root_prior(
        root_prior, states, observed_state_by_leaf, likelihood_by_leaf
    )


def _num_rate_parameters(model, num_states, transition_graph=None):
    states = [str(index) for index in range(num_states)]
    return len(_structured_parameter_labels(model, states, transition_graph))


def _rate_parameter_labels(model, states, transition_graph=None):
    return _structured_parameter_labels(model, states, transition_graph)


def _build_rate_matrix(model, states, rates, transition_graph=None):
    return _build_structured_rate_matrix(model, states, rates, transition_graph)


def _transition_matrix(rate_matrix, branch_length):
    branch_length = float(branch_length)
    if not math.isfinite(branch_length) or branch_length < 0.0:
        raise ValueError("Mk branch lengths must be non-negative and finite.")
    if branch_length == 0.0:
        return np.eye(rate_matrix.shape[0], dtype=float)
    er_rate = _get_er_rate_from_matrix(rate_matrix)
    if er_rate is not None:
        return _er_transition_matrix(branch_length, er_rate, rate_matrix.shape[0])
    matrix = expm(rate_matrix * branch_length)
    exponent_norm = float(np.linalg.norm(rate_matrix, ord=np.inf)) * branch_length
    return _validated_transition_matrix(matrix, exponent_norm)


def _validated_transition_matrix(matrix, exponent_norm):
    if np.any(~np.isfinite(matrix)):
        raise ValueError("Failed to calculate a finite Mk transition matrix.")
    tolerance = min(
        1e-8,
        max(
            1e-12,
            512.0 * np.finfo(float).eps * max(1.0, exponent_norm),
        ),
    )
    if float(np.min(matrix)) < -tolerance:
        raise ValueError(
            "Mk matrix exponentiation produced materially negative probabilities."
        )
    row_sums = matrix.sum(axis=1)
    if np.any(row_sums <= 0.0) or np.any(np.abs(row_sums - 1.0) > tolerance):
        raise ValueError(
            "Mk matrix exponentiation produced invalid transition-probability rows."
        )
    matrix = np.maximum(matrix, 0.0)
    matrix /= matrix.sum(axis=1, keepdims=True)
    return matrix


def _get_er_rate_from_matrix(rate_matrix):
    num_states = rate_matrix.shape[0]
    if num_states <= 1:
        return 0.0
    off_diagonal_mask = ~np.eye(num_states, dtype=bool)
    off_diagonal_rates = rate_matrix[off_diagonal_mask]
    if len(off_diagonal_rates) == 0:
        return 0.0
    rate = float(off_diagonal_rates[0])
    if rate < 0.0:
        return None
    if not np.allclose(off_diagonal_rates, rate, rtol=10**-12, atol=0.0):
        return None
    expected_diagonal = np.full(num_states, -rate * float(num_states - 1), dtype=float)
    if not np.allclose(np.diag(rate_matrix), expected_diagonal, rtol=10**-12, atol=0.0):
        return None
    return rate


def _log_probabilities(values):
    with np.errstate(divide="ignore"):
        return np.log(values)


def _probabilities_from_logs(values):
    normalizer = float(logsumexp(values))
    if not math.isfinite(normalizer):
        raise ValueError("Failed to calculate posterior state probabilities.")
    return np.exp(values - normalizer)


def _compute_inside_likelihoods(
    tree, likelihood_by_leaf, rate_matrix, *, log_space=False
):
    inside, log_scales, child_terms, matrices = _compute_log_inside_likelihoods(
        tree, likelihood_by_leaf, rate_matrix
    )
    if not log_space:
        inside = {node: np.exp(values) for node, values in inside.items()}
        child_terms = {
            node: {child: np.exp(values) for child, values in terms.items()}
            for node, terms in child_terms.items()
        }
    return inside, log_scales, child_terms, matrices


def _compute_log_inside_likelihoods(tree, likelihood_by_leaf, rate_matrix):
    if isinstance(rate_matrix, dict):
        sample_matrix = next(iter(rate_matrix.values()))
    else:
        sample_matrix = rate_matrix
    num_states = sample_matrix.shape[0]
    inside = dict()
    log_scales = dict()
    child_terms: dict[Any, Any] = {}
    transition_matrices = dict()
    transition_matrix_cache = dict()
    for node in tree.traverse(strategy="postorder"):
        if node.is_leaf:
            inside[node] = _log_probabilities(
                likelihood_by_leaf[node.name].astype(float)
            )
            log_scales[node] = 0.0
            continue

        likelihood = np.zeros(num_states, dtype=float)
        log_scale = 0.0
        child_terms[node] = dict()
        for child in node.get_children():
            branch_length = float(child.dist)
            child_rate_matrix = (
                rate_matrix[child] if isinstance(rate_matrix, dict) else rate_matrix
            )
            cache_key = (id(child_rate_matrix), branch_length)
            if cache_key not in transition_matrix_cache:
                matrix = _transition_matrix(child_rate_matrix, branch_length)
                transition_matrix_cache[cache_key] = (
                    matrix,
                    _log_probabilities(matrix),
                )
            matrix, log_matrix = transition_matrix_cache[cache_key]
            transition_matrices[child] = matrix
            term = logsumexp(log_matrix + inside[child][None, :], axis=1)
            child_terms[node][child] = term
            likelihood += term
            log_scale += log_scales[child]
        scale = float(np.max(likelihood))
        if not math.isfinite(scale):
            inside[node] = likelihood
            log_scales[node] = -math.inf
            continue
        inside[node] = likelihood - scale
        log_scales[node] = log_scale + scale
    return inside, log_scales, child_terms, transition_matrices


def _log_likelihood(tree, likelihood_by_leaf, root_prior, rate_matrix):
    inside, log_scales, _, _ = _compute_inside_likelihoods(
        tree=tree,
        likelihood_by_leaf=likelihood_by_leaf,
        rate_matrix=rate_matrix,
        log_space=True,
    )
    root_term = float(logsumexp(_log_probabilities(root_prior) + inside[tree]))
    return log_scales[tree] + root_term


def _initial_rate_value(tree, rate, rate_bounds):
    if rate is not None:
        value = float(rate)
    else:
        branch_lengths = [
            float(node.dist)
            for node in tree.traverse()
            if (not node.is_root) and node.dist is not None and float(node.dist) > 0.0
        ]
        if branch_lengths:
            scale = max(branch_lengths)
            mean_branch_length = scale * (
                math.fsum(branch_length / scale for branch_length in branch_lengths)
                / len(branch_lengths)
            )
        else:
            mean_branch_length = 1.0
        value = 1.0 / max(mean_branch_length, 1.0)
    lower, upper = rate_bounds
    return min(max(value, lower), upper)


def _fit_custom_rate_matrix(
    tree,
    states,
    likelihood_by_leaf,
    rate_bounds,
    fixed_rate_matrix,
    prior_for,
):
    if fixed_rate_matrix is None:
        raise ValueError("--model CUSTOM requires --rate-matrix.")
    rate_matrix = validate_rate_matrix(fixed_rate_matrix, num_states=len(states))
    final_prior = prior_for(rate_matrix)
    return {
        "model": "CUSTOM",
        "rates": np.array([], dtype=float),
        "rate_matrix": rate_matrix,
        "log_likelihood": _log_likelihood(
            tree, likelihood_by_leaf, final_prior, rate_matrix
        ),
        "rate_estimated": False,
        "rate_bounds": rate_bounds,
        "root_prior": final_prior,
        "fit_status": "fixed",
    }


def _fit_parametric_rate_matrix(
    tree,
    model,
    states,
    likelihood_by_leaf,
    root_prior,
    rate=None,
    rate_bounds=None,
    transition_graph=None,
    root_prior_factory=None,
):
    def prior_for(matrix):
        if root_prior_factory is None:
            return root_prior
        return root_prior_factory(matrix)

    num_params = _num_rate_parameters(model, len(states), transition_graph)
    if num_params == 0:
        rate_matrix = _build_rate_matrix(model, states, [], transition_graph)
        final_prior = prior_for(rate_matrix)
        log_likelihood = _log_likelihood(
            tree, likelihood_by_leaf, final_prior, rate_matrix
        )
        return {
            "model": model,
            "rates": np.array([], dtype=float),
            "rate_matrix": rate_matrix,
            "log_likelihood": log_likelihood,
            "rate_estimated": False,
            "rate_bounds": rate_bounds,
            "root_prior": final_prior,
            "fit_status": "fixed",
        }

    if model == "ER" and rate is not None:
        fixed_rate = float(rate)
        if (not math.isfinite(fixed_rate)) or fixed_rate < 0.0:
            raise ValueError("'--rate' must be a non-negative finite number.")
        rate_matrix = _build_rate_matrix(model, states, [fixed_rate], transition_graph)
        final_prior = prior_for(rate_matrix)
        log_likelihood = _log_likelihood(
            tree, likelihood_by_leaf, final_prior, rate_matrix
        )
        return {
            "model": model,
            "rates": np.array([fixed_rate], dtype=float),
            "rate_matrix": rate_matrix,
            "log_likelihood": log_likelihood,
            "rate_estimated": False,
            "rate_bounds": rate_bounds,
            "root_prior": final_prior,
            "fit_status": "fixed",
        }

    initial_rate = _initial_rate_value(tree, rate, rate_bounds)
    parameter_kinds = _structured_parameter_kinds(
        model, states, transition_graph
    )
    parameter_bounds = [
        FREQUENCY_RATIO_BOUNDS if kind == "frequency_ratio" else rate_bounds
        for kind in parameter_kinds
    ]
    lower_logs = np.log([bounds[0] for bounds in parameter_bounds])
    upper_logs = np.log([bounds[1] for bounds in parameter_bounds])
    initial_values = _initial_model_parameters(
        model, states, initial_rate, transition_graph
    )
    initial_logs = np.log(initial_values)

    initial_matrix = _build_rate_matrix(
        model, states, initial_values, transition_graph
    )
    try:
        initial_prior = prior_for(initial_matrix)
        initial_log_likelihood = _log_likelihood(
            tree, likelihood_by_leaf, initial_prior, initial_matrix
        )
    except (ValueError, ArithmeticError):
        initial_log_likelihood = None
    if initial_log_likelihood is not None and not math.isfinite(initial_log_likelihood):
        raise ValueError(
            "Observed states have zero likelihood under the selected transition "
            "graph and root prior."
        )

    def objective(log_rates):
        try:
            rates = np.exp(log_rates)
            rate_matrix = _build_rate_matrix(model, states, rates, transition_graph)
            current_prior = prior_for(rate_matrix)
            log_likelihood = _log_likelihood(
                tree, likelihood_by_leaf, current_prior, rate_matrix
            )
        except (ValueError, ArithmeticError):
            return 1e100
        if not math.isfinite(log_likelihood):
            return 1e100
        return -log_likelihood

    global_fractions = np.linspace(0.0, 1.0, 7)
    starts = [
        lower_logs + fraction * (upper_logs - lower_logs)
        for fraction in global_fractions
    ]
    starts.append(initial_logs)
    if num_params > 1:
        lower_quartiles = lower_logs + 0.25 * (upper_logs - lower_logs)
        upper_quartiles = lower_logs + 0.75 * (upper_logs - lower_logs)
        starts.extend(
            [
                np.asarray(
                    [
                        lower_quartiles[index]
                        if index % 2 == 0
                        else upper_quartiles[index]
                        for index in range(num_params)
                    ],
                    dtype=float,
                ),
                np.asarray(
                    [
                        upper_quartiles[index]
                        if index % 2 == 0
                        else lower_quartiles[index]
                        for index in range(num_params)
                    ],
                    dtype=float,
                ),
                lower_quartiles
                + np.linspace(0.0, 1.0, num_params)
                * (upper_quartiles - lower_quartiles),
                upper_quartiles
                - np.linspace(0.0, 1.0, num_params)
                * (upper_quartiles - lower_quartiles),
            ]
        )
        if num_params <= 4:
            for index in range(num_params):
                for value in (lower_quartiles[index], upper_quartiles[index]):
                    start = initial_logs.copy()
                    start[index] = value
                    starts.append(start)
    unique_starts = []
    seen_starts = set()
    for start in starts:
        key = tuple(float(value) for value in start)
        if key not in seen_starts:
            unique_starts.append(start)
            seen_starts.add(key)

    candidates = []
    failed_messages = []
    for start in unique_starts:
        result = minimize(
            objective,
            start,
            method="L-BFGS-B",
            bounds=list(zip(lower_logs, upper_logs, strict=True)),
        )
        messages = [str(result.message)]
        success = (
            bool(result.success)
            and math.isfinite(float(result.fun))
            and np.all(np.isfinite(result.x))
            and float(result.fun) < 1e99
        )
        if not success:
            fallback = minimize(
                objective,
                start,
                method="Powell",
                bounds=list(zip(lower_logs, upper_logs, strict=True)),
                options={"xtol": 1e-8, "ftol": 1e-10, "maxiter": 500},
            )
            messages.append("Powell fallback: {}".format(fallback.message))
            if (
                bool(fallback.success)
                and math.isfinite(float(fallback.fun))
                and np.all(np.isfinite(fallback.x))
                and float(fallback.fun) < 1e99
            ):
                result = fallback
                success = True
        if success:
            candidates.append((float(result.fun), result, "; ".join(messages)))
        else:
            failed_messages.append("; ".join(messages))
    if not candidates:
        failure_summary = "; ".join(dict.fromkeys(failed_messages))
        raise ValueError(
            "Failed to estimate Mk model parameters: {}".format(
                failure_summary or "no finite fit"
            )
        )
    _, result, selected_message = min(candidates, key=lambda item: item[0])
    if (not math.isfinite(float(result.fun))) or np.any(~np.isfinite(result.x)):
        raise ValueError("Failed to estimate finite Mk model parameters.")
    rates = np.exp(result.x)
    rate_matrix = _build_rate_matrix(model, states, rates, transition_graph)
    final_prior = prior_for(rate_matrix)
    log_likelihood = _log_likelihood(tree, likelihood_by_leaf, final_prior, rate_matrix)
    tolerance = 1e-5
    rate_mask = np.asarray([kind == "rate" for kind in parameter_kinds], dtype=bool)
    frequency_mask = ~rate_mask
    num_rates_at_lower_bound = int(
        np.sum(rate_mask & (rates <= rate_bounds[0] * (1.0 + tolerance)))
    )
    num_rates_at_upper_bound = int(
        np.sum(rate_mask & (rates >= rate_bounds[1] * (1.0 - tolerance)))
    )
    num_frequencies_at_lower_bound = int(
        np.sum(
            frequency_mask
            & (rates <= FREQUENCY_RATIO_BOUNDS[0] * (1.0 + tolerance))
        )
    )
    num_frequencies_at_upper_bound = int(
        np.sum(
            frequency_mask
            & (rates >= FREQUENCY_RATIO_BOUNDS[1] * (1.0 - tolerance))
        )
    )
    statuses = []
    if num_rates_at_lower_bound:
        statuses.append("rate_lower_boundary")
    if num_rates_at_upper_bound:
        statuses.append("rate_upper_boundary")
    if num_frequencies_at_lower_bound:
        statuses.append("frequency_lower_boundary")
    if num_frequencies_at_upper_bound:
        statuses.append("frequency_upper_boundary")
    return {
        "model": model,
        "rates": rates,
        "rate_matrix": rate_matrix,
        "log_likelihood": log_likelihood,
        "rate_estimated": True,
        "rate_bounds": rate_bounds,
        "optimizer_success": bool(result.success),
        "optimizer_message": (
            f"deterministic multistart: {len(candidates)}/{len(unique_starts)} "
            f"starts converged; {selected_message}"
        ),
        "optimizer_starts": len(unique_starts),
        "optimizer_converged_starts": len(candidates),
        "optimizer_failed_starts": len(failed_messages),
        "num_rates_at_lower_bound": num_rates_at_lower_bound,
        "num_rates_at_upper_bound": num_rates_at_upper_bound,
        "num_frequencies_at_lower_bound": num_frequencies_at_lower_bound,
        "num_frequencies_at_upper_bound": num_frequencies_at_upper_bound,
        "fit_status": "+".join(statuses) if statuses else "ok",
        "root_prior": final_prior,
    }


def _fit_rate_matrix(
    tree,
    model,
    states,
    likelihood_by_leaf,
    root_prior,
    rate=None,
    rate_bounds=None,
    transition_graph=None,
    fixed_rate_matrix=None,
    root_prior_factory=None,
    regime_assignment=None,
    regime_model="ER",
):
    if model not in SUPPORTED_MODELS:
        raise ValueError("Unsupported '--model': {}".format(model))
    rate_bounds = DEFAULT_RATE_BOUNDS if rate_bounds is None else rate_bounds

    def prior_for(matrix):
        if root_prior_factory is None:
            return root_prior
        return root_prior_factory(matrix)

    if model == "CUSTOM":
        return _fit_custom_rate_matrix(
            tree,
            states,
            likelihood_by_leaf,
            rate_bounds,
            fixed_rate_matrix,
            prior_for,
        )
    if model == "MK-REGIME":
        return _fit_regime_rate_matrices(
            tree=tree,
            states=states,
            likelihood_by_leaf=likelihood_by_leaf,
            root_prior=root_prior,
            rate=rate,
            rate_bounds=rate_bounds,
            transition_graph=transition_graph,
            root_prior_factory=root_prior_factory,
            regime_assignment=regime_assignment,
            regime_model=regime_model,
        )
    return _fit_parametric_rate_matrix(
        tree=tree,
        model=model,
        states=states,
        likelihood_by_leaf=likelihood_by_leaf,
        root_prior=root_prior,
        rate=rate,
        rate_bounds=rate_bounds,
        transition_graph=transition_graph,
        root_prior_factory=root_prior_factory,
    )


def _fit_regime_rate_matrices(
    tree,
    states,
    likelihood_by_leaf,
    root_prior,
    rate,
    rate_bounds,
    transition_graph,
    root_prior_factory,
    regime_assignment,
    regime_model,
):
    if regime_assignment is None:
        raise ValueError("--model MK-REGIME requires --regime-map.")
    if regime_model not in {"ER", "SYM", "ARD", "F81", "GTR"}:
        raise ValueError(
            "--regime-model must be one of ER, SYM, ARD, F81, or GTR."
        )
    regimes = regime_assignment.regimes
    represented_edges = {
        regime_assignment.by_node[node]
        for node in tree.traverse()
        if not node.is_root
    }
    root_only = sorted(set(regimes) - represented_edges)
    if root_only:
        raise ValueError(
            "Every estimated regime must label at least one non-root branch; "
            "root-only regime(s): " + ", ".join(root_only)
        )
    labels = _rate_parameter_labels(regime_model, states, transition_graph)
    if not labels:
        raise ValueError("MK-REGIME requires at least one rate per regime.")
    kinds = _structured_parameter_kinds(regime_model, states, transition_graph)
    per_regime = len(labels)
    initial_rate = _initial_rate_value(tree, rate, rate_bounds)
    base_initial = _initial_model_parameters(
        regime_model, states, initial_rate, transition_graph
    )
    initial_values = np.tile(base_initial, len(regimes))
    all_kinds = kinds * len(regimes)
    bounds = [
        FREQUENCY_RATIO_BOUNDS if kind == "frequency_ratio" else rate_bounds
        for kind in all_kinds
    ]
    lower_logs = np.log([item[0] for item in bounds])
    upper_logs = np.log([item[1] for item in bounds])

    def matrices_for(parameters):
        by_regime = {}
        for regime_index, regime in enumerate(regimes):
            start = regime_index * per_regime
            by_regime[regime] = _build_rate_matrix(
                regime_model,
                states,
                parameters[start : start + per_regime],
                transition_graph,
            )
        by_node = {
            node: by_regime[regime_assignment.by_node[node]]
            for node in tree.traverse()
            if not node.is_root
        }
        return by_regime, by_node

    def prior_for(by_regime):
        root_matrix = by_regime[regime_assignment.root_regime]
        return (
            root_prior
            if root_prior_factory is None
            else root_prior_factory(root_matrix)
        )

    def objective(log_parameters):
        try:
            parameters = np.exp(log_parameters)
            by_regime, by_node = matrices_for(parameters)
            current_prior = prior_for(by_regime)
            likelihood = _log_likelihood(
                tree, likelihood_by_leaf, current_prior, by_node
            )
        except (ValueError, ArithmeticError):
            return 1e100
        return -likelihood if math.isfinite(likelihood) else 1e100

    initial_logs = np.log(initial_values)
    starts = [initial_logs]
    for fraction in np.linspace(0.0, 1.0, 7):
        starts.append(lower_logs + fraction * (upper_logs - lower_logs))
    if len(regimes) > 1:
        for fraction in (0.25, 0.75):
            start = initial_logs.copy()
            for regime_index in range(len(regimes)):
                if regime_index % 2:
                    left = regime_index * per_regime
                    right = left + per_regime
                    start[left:right] = lower_logs[left:right] + fraction * (
                        upper_logs[left:right] - lower_logs[left:right]
                    )
            starts.append(start)
    unique_starts = []
    seen = set()
    for start in starts:
        key = tuple(float(value) for value in start)
        if key not in seen:
            seen.add(key)
            unique_starts.append(start)
    candidates = []
    failures = []
    optimizer_bounds = list(zip(lower_logs, upper_logs, strict=True))
    for start in unique_starts:
        result = minimize(
            objective, start, method="L-BFGS-B", bounds=optimizer_bounds
        )
        messages = [str(result.message)]
        success = (
            bool(result.success)
            and math.isfinite(float(result.fun))
            and float(result.fun) < 1e99
            and np.all(np.isfinite(result.x))
        )
        if not success:
            fallback = minimize(
                objective,
                start,
                method="Powell",
                bounds=optimizer_bounds,
                options={"xtol": 1e-8, "ftol": 1e-10, "maxiter": 800},
            )
            messages.append(f"Powell fallback: {fallback.message}")
            if (
                bool(fallback.success)
                and math.isfinite(float(fallback.fun))
                and float(fallback.fun) < 1e99
                and np.all(np.isfinite(fallback.x))
            ):
                result, success = fallback, True
        if success:
            candidates.append((float(result.fun), result, "; ".join(messages)))
        else:
            failures.append("; ".join(messages))
    if not candidates:
        raise ValueError(
            "Failed to estimate MK-REGIME parameters: "
            + ("; ".join(dict.fromkeys(failures)) or "no finite fit")
        )
    _, selected, selected_message = min(candidates, key=lambda item: item[0])
    parameters = np.exp(selected.x)
    by_regime, by_node = matrices_for(parameters)
    final_prior = prior_for(by_regime)
    log_likelihood = _log_likelihood(
        tree, likelihood_by_leaf, final_prior, by_node
    )
    tolerance = 1e-5
    lower = np.asarray([item[0] for item in bounds])
    upper = np.asarray([item[1] for item in bounds])
    rate_mask = np.asarray([kind == "rate" for kind in all_kinds], dtype=bool)
    frequency_mask = ~rate_mask
    lower_hits = parameters <= lower * (1.0 + tolerance)
    upper_hits = parameters >= upper * (1.0 - tolerance)
    num_rates_at_lower_bound = int(np.sum(rate_mask & lower_hits))
    num_rates_at_upper_bound = int(np.sum(rate_mask & upper_hits))
    num_frequencies_at_lower_bound = int(np.sum(frequency_mask & lower_hits))
    num_frequencies_at_upper_bound = int(np.sum(frequency_mask & upper_hits))
    statuses = []
    if num_rates_at_lower_bound:
        statuses.append("rate_lower_boundary")
    if num_rates_at_upper_bound:
        statuses.append("rate_upper_boundary")
    if num_frequencies_at_lower_bound:
        statuses.append("frequency_lower_boundary")
    if num_frequencies_at_upper_bound:
        statuses.append("frequency_upper_boundary")
    rates_by_regime = {
        regime: parameters[index * per_regime : (index + 1) * per_regime]
        for index, regime in enumerate(regimes)
    }
    return {
        "model": "MK-REGIME",
        "regime_model": regime_model,
        "regimes": regimes,
        "rates": parameters,
        "rates_by_regime": rates_by_regime,
        "rate_matrix": by_regime[regime_assignment.root_regime],
        "rate_matrices_by_regime": by_regime,
        "rate_matrix_by_node": by_node,
        "regime_by_node": regime_assignment.by_node,
        "regime_map_source": regime_assignment.source,
        "root_regime": regime_assignment.root_regime,
        "log_likelihood": log_likelihood,
        "rate_estimated": True,
        "rate_bounds": rate_bounds,
        "optimizer_success": bool(selected.success),
        "optimizer_message": (
            f"deterministic multistart: {len(candidates)}/{len(unique_starts)} "
            f"starts converged; {selected_message}"
        ),
        "optimizer_starts": len(unique_starts),
        "optimizer_converged_starts": len(candidates),
        "optimizer_failed_starts": len(failures),
        "num_rates_at_lower_bound": num_rates_at_lower_bound,
        "num_rates_at_upper_bound": num_rates_at_upper_bound,
        "num_frequencies_at_lower_bound": num_frequencies_at_lower_bound,
        "num_frequencies_at_upper_bound": num_frequencies_at_upper_bound,
        "fit_status": "+".join(statuses) if statuses else "ok",
        "root_prior": final_prior,
    }


def _compute_outside_likelihoods(
    tree, child_terms, transition_matrices, root_prior, *, log_space=False
):
    if not log_space:
        child_terms = {
            node: {child: _log_probabilities(term) for child, term in terms.items()}
            for node, terms in child_terms.items()
        }
    outside = _compute_log_outside_likelihoods(
        tree, child_terms, transition_matrices, root_prior
    )
    return (
        outside
        if log_space
        else {node: np.exp(values) for node, values in outside.items()}
    )


def _compute_log_outside_likelihoods(
    tree, child_terms, transition_matrices, root_prior
):
    outside = {tree: _log_probabilities(root_prior)}
    for node in tree.traverse(strategy="preorder"):
        children = list(node.get_children())
        if len(children) == 0:
            continue
        terms = [child_terms[node][child] for child in children]
        prefix_log_sums = [np.zeros_like(outside[node], dtype=float)]
        for term in terms:
            prefix_log_sums.append(prefix_log_sums[-1] + term)
        suffix_log_sum = np.zeros_like(outside[node], dtype=float)
        for child_index in range(len(children) - 1, -1, -1):
            child = children[child_index]
            sibling_log_sum = prefix_log_sums[child_index] + suffix_log_sum
            parent_log_weight = outside[node] + sibling_log_sum
            matrix = _log_probabilities(transition_matrices[child])
            outside[child] = logsumexp(matrix.T + parent_log_weight[None, :], axis=1)
            total = float(logsumexp(outside[child]))
            if math.isfinite(total):
                outside[child] -= total
            suffix_log_sum = suffix_log_sum + terms[child_index]
    return outside


def _root_prior_configuration(
    root_prior_mode, states, observed_state_by_leaf, likelihood_by_leaf
):
    if root_prior_mode != "stationary":
        return (
            _get_root_prior(
                root_prior_mode,
                states,
                observed_state_by_leaf,
                likelihood_by_leaf,
            ),
            None,
        )
    root_prior = np.full(len(states), 1.0 / float(len(states)), dtype=float)

    def root_prior_factory(matrix):
        return _root_prior_for_matrix(
            root_prior_mode,
            states,
            observed_state_by_leaf,
            likelihood_by_leaf,
            matrix,
        )

    return root_prior, root_prior_factory


def compute_mk_marginals(
    tree,
    states,
    observed_state_by_leaf,
    likelihood_by_leaf,
    model="ER",
    rate=None,
    root_prior_mode="equal",
    rate_bounds=None,
    transition_graph=None,
    fixed_rate_matrix=None,
    regime_assignment=None,
    regime_model="ER",
):
    root_prior, root_prior_factory = _root_prior_configuration(
        root_prior_mode, states, observed_state_by_leaf, likelihood_by_leaf
    )
    fit = _fit_rate_matrix(
        tree=tree,
        model=model,
        states=states,
        likelihood_by_leaf=likelihood_by_leaf,
        root_prior=root_prior,
        rate=rate,
        rate_bounds=DEFAULT_RATE_BOUNDS if rate_bounds is None else rate_bounds,
        transition_graph=transition_graph,
        fixed_rate_matrix=fixed_rate_matrix,
        root_prior_factory=root_prior_factory,
        regime_assignment=regime_assignment,
        regime_model=regime_model,
    )
    root_prior = fit.get("root_prior", root_prior)
    if not math.isfinite(fit["log_likelihood"]):
        raise ValueError(
            "The observed tip states have zero likelihood under the Mk model."
        )

    inside, _, child_terms, transition_matrices = _compute_inside_likelihoods(
        tree=tree,
        likelihood_by_leaf=likelihood_by_leaf,
        rate_matrix=fit.get("rate_matrix_by_node", fit["rate_matrix"]),
        log_space=True,
    )
    outside = _compute_outside_likelihoods(
        tree=tree,
        child_terms=child_terms,
        transition_matrices=transition_matrices,
        root_prior=root_prior,
        log_space=True,
    )
    posterior_by_node = dict()
    for node in tree.traverse():
        posterior_by_node[node] = _probabilities_from_logs(inside[node] + outside[node])
    fit.update(
        {
            "root_prior": root_prior,
            "inside": {node: np.exp(values) for node, values in inside.items()},
            "log_inside": inside,
            "transition_matrices": transition_matrices,
            "posterior_by_node": posterior_by_node,
        }
    )
    return posterior_by_node, fit


def compute_hrm_marginals(
    tree,
    states,
    observed_state_by_leaf,
    likelihood_by_leaf,
    *,
    hidden_categories=2,
    rate=None,
    root_prior_mode="equal",
    rate_bounds=None,
    transition_graph=None,
):
    """Fit a fully parameterized hidden-rates CTMC and marginalize classes."""

    from nwkit.hidden_rate_asr import (
        aggregate_probabilities,
        expand_tip_likelihoods,
        expanded_state_labels,
        expanded_transition_graph,
        state_projection,
    )

    try:
        hidden_categories = int(hidden_categories)
    except (TypeError, ValueError) as exc:
        raise ValueError("--hidden-categories must be an integer.") from exc
    if hidden_categories < 2:
        raise ValueError("--hidden-categories must be at least 2 for HRM.")
    expanded_states = expanded_state_labels(states, hidden_categories)
    observed_graph = (
        np.ones((len(states), len(states)), dtype=bool)
        if transition_graph is None
        else np.asarray(transition_graph, dtype=bool)
    )
    np.fill_diagonal(observed_graph, False)
    expanded_graph = expanded_transition_graph(observed_graph, hidden_categories)
    num_parameters = int(np.sum(expanded_graph))
    if num_parameters > 256:
        raise ValueError(
            "HRM would require more than 256 free transition rates; reduce the "
            "observed state space or --hidden-categories."
        )
    expanded_likelihoods = expand_tip_likelihoods(
        likelihood_by_leaf, hidden_categories
    )
    expanded_posterior, fit = compute_mk_marginals(
        tree,
        expanded_states,
        observed_state_by_leaf,
        expanded_likelihoods,
        model="ARD",
        rate=rate,
        root_prior_mode=root_prior_mode,
        rate_bounds=rate_bounds,
        transition_graph=expanded_graph,
    )
    projection = state_projection(len(states), hidden_categories)
    posterior = {
        node: aggregate_probabilities(
            probabilities, len(states), hidden_categories
        )
        for node, probabilities in expanded_posterior.items()
    }
    expanded_root_prior = fit["root_prior"]
    fit.update(
        {
            "model": "HRM",
            "hidden_categories": hidden_categories,
            "expanded_states": expanded_states,
            "state_projection": projection,
            "expanded_root_prior": expanded_root_prior,
            "root_prior": aggregate_probabilities(
                expanded_root_prior, len(states), hidden_categories
            ),
            "expanded_posterior_by_node": expanded_posterior,
            # Sampling must retain the hidden-class posterior.
            "posterior_by_node": expanded_posterior,
            "observed_transition_graph": observed_graph,
        }
    )
    return posterior, fit


def _is_missing_tip(node, observed_state_by_leaf):
    return node.is_leaf and observed_state_by_leaf.get(node.name) is None


def _should_output_node(node, observed_state_by_leaf, targets):
    if "all" in targets:
        return True
    if (not node.is_leaf) and ("intnode" in targets):
        return True
    if node.is_leaf and ("leaf" in targets):
        return True
    if _is_missing_tip(node, observed_state_by_leaf) and ("missing-leaf" in targets):
        return True
    return False


def _build_output_table(
    tree, states, observed_state_by_leaf, posterior_by_node, targets, output_mode
):
    node_to_branch_id = assign_branch_ids(tree)
    state_ids = [_safe_column_state(state) for state in states]
    rows = list()
    for node in tree.traverse():
        if not _should_output_node(node, observed_state_by_leaf, targets):
            continue
        posterior = posterior_by_node[node]
        map_index = int(np.argmax(posterior))
        observed_state = observed_state_by_leaf.get(node.name) if node.is_leaf else None
        row = {
            "branch_id": node_to_branch_id[node],
            "parent": -1 if node.is_root else node_to_branch_id[node.up],
            "node_class": get_node_class(node),
            "name": "" if node.name in [None, ""] else str(node.name),
            "observed_state": "" if observed_state is None else observed_state,
            "is_imputed": bool(_is_missing_tip(node, observed_state_by_leaf)),
        }
        if output_mode == "probabilities":
            row["map_state"] = states[map_index]
            row["map_probability"] = float(posterior[map_index])
            for state_id, probability in zip(state_ids, posterior, strict=True):
                row["p_{}".format(state_id)] = float(probability)
        elif output_mode == "map":
            row["state"] = states[map_index]
            row["probability"] = float(posterior[map_index])
        else:
            raise ValueError("Unsupported '--output': {}".format(output_mode))
        rows.append(row)
    if output_mode == "probabilities":
        columns = [
            "branch_id",
            "parent",
            "node_class",
            "name",
            "observed_state",
            "is_imputed",
            "map_state",
            "map_probability",
        ] + ["p_{}".format(state_id) for state_id in state_ids]
    else:
        columns = [
            "branch_id",
            "parent",
            "node_class",
            "name",
            "observed_state",
            "is_imputed",
            "state",
            "probability",
        ]
    return pd.DataFrame(rows, columns=columns)


def _safe_column_state(state):
    state_text = str(state)
    escape_prefix = "state_"
    is_legacy_safe = bool(state_text) and all(
        char.isalnum() or char == "_" for char in state_text
    )
    if is_legacy_safe and not state_text.startswith(escape_prefix):
        return state_text
    return "{}{}".format(escape_prefix, state_text.encode("utf-8").hex())


def _model_rate_bounds_text(fit):
    if fit["model"] == "CUSTOM":
        return ""
    return "{},{}".format(fit["rate_bounds"][0], fit["rate_bounds"][1])


def _build_model_table(states, root_prior_mode, fit):
    state_ids = [_safe_column_state(state) for state in states]
    row = {
        "trait_type": "discrete",
        "trait_type_requested": fit.get("trait_type_requested", "discrete"),
        "model": fit["model"],
        "num_states": len(states),
        "states": ",".join(states),
        "root_prior": root_prior_mode,
        "root_prior_values": ",".join(str(float(value)) for value in fit["root_prior"]),
        "rate_estimated": bool(fit["rate_estimated"]),
        "log_likelihood": float(fit["log_likelihood"]),
        "rate_bounds": _model_rate_bounds_text(fit),
        "transition_graph": fit.get("transition_graph", "complete"),
        "regime_model": fit.get("regime_model", ""),
        "regimes": ",".join(fit.get("regimes", ())),
        "root_regime": fit.get("root_regime", ""),
        "regime_map": fit.get("regime_map_source", ""),
        "hidden_categories": fit.get("hidden_categories", ""),
        "q_source": fit.get("rate_matrix_source", "estimated"),
        "fit_status": fit.get("fit_status", ""),
        "optimizer_success": fit.get("optimizer_success", ""),
        "optimizer_message": fit.get("optimizer_message", ""),
        "optimizer_starts": fit.get("optimizer_starts", ""),
        "optimizer_converged_starts": fit.get("optimizer_converged_starts", ""),
        "optimizer_failed_starts": fit.get("optimizer_failed_starts", ""),
        "num_rates_at_lower_bound": fit.get("num_rates_at_lower_bound", 0),
        "num_rates_at_upper_bound": fit.get("num_rates_at_upper_bound", 0),
        "frequency_ratio_bounds": (
            f"{FREQUENCY_RATIO_BOUNDS[0]},{FREQUENCY_RATIO_BOUNDS[1]}"
            if fit["model"] == "GTR"
            or (
                fit["model"] == "MK-REGIME"
                and fit.get("regime_model") == "GTR"
            )
            else ""
        ),
        "num_frequencies_at_lower_bound": fit.get(
            "num_frequencies_at_lower_bound", 0
        ),
        "num_frequencies_at_upper_bound": fit.get(
            "num_frequencies_at_upper_bound", 0
        ),
    }
    if fit["model"] == "ER" and len(fit["rates"]) == 1:
        row["rate"] = float(fit["rates"][0])
    if fit["model"] in {"F81", "GTR"}:
        equilibrium = stationary_distribution(fit["rate_matrix"])
        for state_id, probability in zip(state_ids, equilibrium, strict=True):
            row[f"equilibrium_{state_id}"] = float(probability)
    if fit["model"] == "MK-REGIME":
        for regime, matrix in fit["rate_matrices_by_regime"].items():
            regime_id = _safe_column_state(regime)
            if fit.get("regime_model") in {"F81", "GTR"}:
                equilibrium = stationary_distribution(matrix)
                for state_id, probability in zip(
                    state_ids, equilibrium, strict=True
                ):
                    row[f"equilibrium_{regime_id}_{state_id}"] = float(probability)
            for from_index, from_state_id in enumerate(state_ids):
                for to_index, to_state_id in enumerate(state_ids):
                    if from_index != to_index:
                        row[
                            f"rate_{regime_id}_{from_state_id}_to_{to_state_id}"
                        ] = float(matrix[from_index, to_index])
    if fit["model"] == "HRM":
        expanded_ids = [
            _safe_column_state(state) for state in fit["expanded_states"]
        ]
        row["num_expanded_states"] = len(fit["expanded_states"])
        row["expanded_states"] = ",".join(fit["expanded_states"])
        for from_index, from_state_id in enumerate(expanded_ids):
            for to_index, to_state_id in enumerate(expanded_ids):
                if from_index != to_index:
                    row[f"rate_{from_state_id}_to_{to_state_id}"] = float(
                        fit["rate_matrix"][from_index, to_index]
                    )
    elif fit["model"] != "MK-REGIME":
        for from_index, from_state_id in enumerate(state_ids):
            for to_index, to_state_id in enumerate(state_ids):
                if from_index == to_index:
                    continue
                row["rate_{}_to_{}".format(from_state_id, to_state_id)] = float(
                    fit["rate_matrix"][from_index, to_index]
                )
    return pd.DataFrame([row])


def _write_table(table, outfile):
    if outfile == "-":
        print(table.to_csv(sep="\t", index=False), end="")
    else:
        table.to_csv(outfile, sep="\t", index=False)


def _write_model_table(states, root_prior_mode, fit, outfile):
    if outfile in ["", None]:
        return
    table = _build_model_table(states, root_prior_mode, fit)
    _write_table(table, outfile)


def _annotate_tree(tree, states, posterior_by_node, observed_state_by_leaf):
    state_ids = [_safe_column_state(state) for state in states]
    for node in tree.traverse():
        posterior = posterior_by_node[node]
        map_index = int(np.argmax(posterior))
        node.add_props(
            asr_state=states[map_index],
            asr_probability=float(posterior[map_index]),
        )
        if node.is_leaf:
            observed_state = observed_state_by_leaf.get(node.name)
            node.add_props(
                asr_observed_state="" if observed_state is None else observed_state
            )
        for state_id, probability in zip(state_ids, posterior, strict=True):
            node.add_props(**{"asr_p_{}".format(state_id): float(probability)})


def _write_annotated_tree(
    tree, states, posterior_by_node, observed_state_by_leaf, args
):
    tree_out = getattr(args, "tree_out", None)
    if tree_out in ["", None]:
        return
    _annotate_tree(tree, states, posterior_by_node, observed_state_by_leaf)
    tree_annotation = getattr(args, "tree_annotation", "map")
    if tree_annotation == "state":
        props = ["asr_state"]
    elif tree_annotation == "probability":
        props = ["asr_probability"]
    elif tree_annotation == "map":
        props = ["asr_state", "asr_probability"]
    elif tree_annotation == "all":
        props = ["asr_state", "asr_probability"]
        props += ["asr_observed_state"] + [
            "asr_p_{}".format(_safe_column_state(state)) for state in states
        ]
    else:
        raise ValueError("Unsupported '--tree-annotation': {}".format(tree_annotation))
    output_args = copy(args)
    output_args.outfile = tree_out
    write_tree(
        tree,
        output_args,
        format=getattr(args, "tree_outformat", "auto"),
        quiet=True,
        props=props,
    )


def _normalize_probability_vector(weights):
    weights = np.asarray(weights, dtype=float)
    weights = np.maximum(weights, 0.0)
    total = float(weights.sum())
    if total <= 0.0:
        return np.full(len(weights), 1.0 / float(len(weights)), dtype=float)
    return weights / total


def _build_uniformization_context(rate_matrix, branch_length):
    branch_length = float(branch_length)
    num_states = rate_matrix.shape[0]
    if branch_length == 0.0 or not np.any(rate_matrix):
        return {"no_events": True}
    omega = float(np.max(-np.diag(rate_matrix)))
    if omega <= 0.0:
        return {"no_events": True}
    lam = omega * branch_length
    if not 0.0 <= lam < math.inf:
        raise ValueError(
            "Stochastic mapping requires a finite non-negative rate × time."
        )
    poisson_limit = float(poisson.ppf(1.0 - 10**-12, lam))
    if not math.isfinite(poisson_limit):
        raise ValueError(
            "Stochastic mapping requires too many uniformization terms for "
            "floating-point indexing."
        )
    max_n = int(max(10, num_states - 1, poisson_limit))
    if max_n > _MAX_UNIFORMIZATION_TERMS:
        raise ValueError(
            "Stochastic mapping would require more than "
            f"{_MAX_UNIFORMIZATION_TERMS:,} potential events on one branch; "
            "reduce the rate/time scale or fitted rate bounds."
        )
    r_matrix = np.eye(num_states, dtype=float) + (rate_matrix / omega)
    r_matrix = np.maximum(r_matrix, 0.0)
    r_matrix = r_matrix / r_matrix.sum(axis=1, keepdims=True)
    event_counts = np.arange(max_n + 1, dtype=np.int64)
    log_poisson = poisson.logpmf(event_counts, lam)
    full_cache_bytes = 2 * (max_n + 1) * num_states * num_states * 8
    cache_small_context = full_cache_bytes <= _MAX_CACHED_BACKWARD_BYTES
    return {
        "no_events": False,
        "r_matrix": r_matrix,
        "event_counts": event_counts,
        "log_poisson": log_poisson,
        "backward_by_end": {} if cache_small_context else None,
        "bridge_probability_cache": {} if cache_small_context else None,
    }


def _build_uniformization_contexts(rate_matrix, branch_lengths):
    contexts = dict()
    for branch_length in sorted(set(float(length) for length in branch_lengths)):
        contexts[branch_length] = _build_uniformization_context(
            rate_matrix, branch_length
        )
    return contexts


def _sample_bridge_event_count(
    start_state,
    end_state,
    branch_length,
    rate_matrix,
    rng,
    uniformization_contexts=None,
):
    branch_length = float(branch_length)
    if uniformization_contexts is None:
        context = _build_uniformization_context(rate_matrix, branch_length)
    else:
        context = uniformization_contexts[branch_length]
    event_count, _ = _draw_bridge_event_count(start_state, end_state, context, rng)
    return event_count


def _constant_bridge_event_count(start_state, end_state):
    if start_state != end_state:
        raise ValueError("A zero-rate or zero-length branch cannot change state.")
    return 0


def _bridge_probabilities(weights):
    if not np.any(weights > 0.0):
        raise ValueError("Stochastic bridge has no feasible next state.")
    return _normalize_probability_vector(weights)


def _backward_bridge_messages(context, end_state):
    cache = context.get("backward_by_end")
    if cache is not None and end_state in cache:
        return cache[end_state]
    r_matrix = context["r_matrix"]
    num_states = r_matrix.shape[0]
    max_n = len(context["event_counts"]) - 1
    backward = np.empty((max_n + 1, num_states), dtype=float)
    backward[0, :] = 0.0
    backward[0, end_state] = 1.0
    for event_count in range(max_n):
        backward[event_count + 1, :] = r_matrix.dot(backward[event_count, :])
    if cache is not None:
        cache[end_state] = backward
    return backward


def _draw_bridge_event_count(start_state, end_state, context, rng):
    if context.get("no_events"):
        return _constant_bridge_event_count(start_state, end_state), None
    cache = context.get("bridge_probability_cache")
    cache_key = (start_state, end_state)
    probabilities = None if cache is None else cache.get(cache_key)
    backward = _backward_bridge_messages(context, end_state)
    if probabilities is None:
        bridge_probabilities = backward[:, start_state]
        valid = bridge_probabilities > 0.0
        if not np.any(valid):
            raise ValueError(
                "No feasible stochastic bridge connects the requested states."
            )
        log_weights = np.full(len(bridge_probabilities), -math.inf, dtype=float)
        log_weights[valid] = context["log_poisson"][valid] + np.log(
            bridge_probabilities[valid]
        )
        normalizer = float(logsumexp(log_weights))
        if not math.isfinite(normalizer):
            raise ValueError(
                "No feasible stochastic bridge connects the requested states."
            )
        probabilities = _normalize_probability_vector(np.exp(log_weights - normalizer))
        if cache is not None:
            cache[cache_key] = probabilities
    event_count = int(rng.choice(context["event_counts"], p=probabilities))
    return event_count, backward


def _sample_bridge_transition_counts(
    start_state,
    end_state,
    branch_length,
    rate_matrix,
    rng,
    uniformization_contexts=None,
    uniformization_context=None,
):
    branch_length = float(branch_length)
    if uniformization_context is not None:
        context = uniformization_context
    elif uniformization_contexts is None:
        context = _build_uniformization_context(rate_matrix, branch_length)
    else:
        context = uniformization_contexts[branch_length]
    event_count, backward = _draw_bridge_event_count(
        start_state, end_state, context, rng
    )
    num_states = rate_matrix.shape[0]
    counts: defaultdict[Any, int] = defaultdict(int)
    if event_count == 0:
        return counts
    r_matrix = context["r_matrix"]
    current_state = start_state
    for step_index in range(event_count):
        remaining = event_count - step_index - 1
        weights = r_matrix[current_state, :] * backward[remaining, :]
        probabilities = _bridge_probabilities(weights)
        next_state = int(rng.choice(np.arange(num_states), p=probabilities))
        if next_state != current_state:
            counts[(current_state, next_state)] += 1
        current_state = next_state
    return counts


def _sample_node_states(tree, states, fit, rng):
    sampled_states = dict()
    num_states = len(fit["posterior_by_node"][tree])
    sampled_states[tree] = int(
        rng.choice(np.arange(num_states), p=fit["posterior_by_node"][tree])
    )
    for node in tree.traverse(strategy="preorder"):
        for child in node.get_children():
            parent_state = sampled_states[node]
            transition_matrix = fit["transition_matrices"][child]
            weights = (
                _log_probabilities(transition_matrix[parent_state, :])
                + fit["log_inside"][child]
            )
            probabilities = _probabilities_from_logs(weights)
            sampled_states[child] = int(
                rng.choice(np.arange(num_states), p=probabilities)
            )
    return sampled_states


def _simulation_seed_sequence(seed, num_simulations):
    if seed is None:
        seed_sequence = np.random.SeedSequence()
    else:
        seed_sequence = np.random.SeedSequence(int(seed))
    return seed_sequence.spawn(num_simulations)


def _build_stochastic_map_spec(tree, fit, node_to_branch_id, uniformization_contexts):
    nodes = list(tree.traverse(strategy="preorder"))
    node_to_index = {node: index for index, node in enumerate(nodes)}
    children_by_index: list[list[int]] = [[] for _ in nodes]
    parent_indices = [-1] * len(nodes)
    branch_ids = [-1] * len(nodes)
    branch_lengths = [0.0] * len(nodes)
    transition_matrices = np.zeros(
        (len(nodes), fit["rate_matrix"].shape[0], fit["rate_matrix"].shape[1]),
        dtype=float,
    )
    rate_matrices = np.zeros_like(transition_matrices)
    contexts = [None] * len(nodes)
    inside = np.zeros((len(nodes), fit["rate_matrix"].shape[0]), dtype=float)
    for node, node_index in node_to_index.items():
        inside[node_index, :] = fit["log_inside"][node]
        if node.is_root:
            continue
        parent_index = node_to_index[node.up]
        parent_indices[node_index] = parent_index
        children_by_index[parent_index].append(node_index)
        branch_ids[node_index] = node_to_branch_id[node]
        branch_lengths[node_index] = float(node.dist)
        transition_matrices[node_index, :, :] = fit["transition_matrices"][node]
        rate_matrices[node_index, :, :] = fit.get("rate_matrix_by_node", {}).get(
            node, fit["rate_matrix"]
        )
        contexts[node_index] = (
            uniformization_contexts[node]
            if node in uniformization_contexts
            else uniformization_contexts[float(node.dist)]
        )
    return {
        "children_by_index": children_by_index,
        "parent_indices": parent_indices,
        "branch_ids": branch_ids,
        "branch_lengths": branch_lengths,
        "log_inside": inside,
        "log_transition_matrices": _log_probabilities(transition_matrices),
        "rate_matrix": fit["rate_matrix"],
        "rate_matrices": rate_matrices,
        "root_posterior": fit["posterior_by_node"][tree],
        "uniformization_contexts": contexts,
        "num_states": fit["rate_matrix"].shape[0],
        "state_projection": fit.get("state_projection"),
    }


def _sample_node_states_from_spec(spec, rng):
    sampled_states = np.zeros(len(spec["parent_indices"]), dtype=int)
    state_indices = np.arange(spec["num_states"])
    sampled_states[0] = int(rng.choice(state_indices, p=spec["root_posterior"]))
    for parent_index, child_indices in enumerate(spec["children_by_index"]):
        parent_state = sampled_states[parent_index]
        for child_index in child_indices:
            transition_matrix = spec["log_transition_matrices"][child_index]
            weights = (
                transition_matrix[parent_state, :] + spec["log_inside"][child_index]
            )
            probabilities = _probabilities_from_logs(weights)
            sampled_states[child_index] = int(
                rng.choice(state_indices, p=probabilities)
            )
    return sampled_states


def _simulate_stochastic_map_once(spec, seed_sequence):
    rng = np.random.default_rng(seed_sequence)
    sampled_states = _sample_node_states_from_spec(spec, rng)
    simulation_counts: defaultdict[Any, int] = defaultdict(int)
    for node_index, parent_index in enumerate(spec["parent_indices"]):
        if parent_index < 0:
            continue
        branch_counts = _sample_bridge_transition_counts(
            int(sampled_states[parent_index]),
            int(sampled_states[node_index]),
            spec["branch_lengths"][node_index],
            spec["rate_matrices"][node_index],
            rng,
            uniformization_context=spec["uniformization_contexts"][node_index],
        )
        for pair, count in branch_counts.items():
            projection = spec.get("state_projection")
            from_state, to_state = pair
            if projection is not None:
                from_state = int(projection[from_state])
                to_state = int(projection[to_state])
            if from_state != to_state:
                simulation_counts[
                    (spec["branch_ids"][node_index], from_state, to_state)
                ] += count
    return simulation_counts


def _simulate_stochastic_map_chunk(spec, seed_sequences):
    return [
        _simulate_stochastic_map_once(spec, seed_sequence)
        for seed_sequence in seed_sequences
    ]


def _simulate_stochastic_map_chunk_worker(payload):
    spec, seed_sequences = payload
    return _merge_simulation_counts(
        _simulate_stochastic_map_once(spec, seed_sequence)
        for seed_sequence in seed_sequences
    )


def _merge_simulation_counts(simulation_counts_list):
    total_counts: defaultdict[Any, int] = defaultdict(int)
    any_counts: defaultdict[Any, int] = defaultdict(int)
    for simulation_counts in simulation_counts_list:
        for key, count in simulation_counts.items():
            total_counts[key] += count
            any_counts[key] += 1
    return total_counts, any_counts


def _merge_count_summaries(count_summaries):
    total_counts: defaultdict[Any, int] = defaultdict(int)
    any_counts: defaultdict[Any, int] = defaultdict(int)
    for chunk_total_counts, chunk_any_counts in count_summaries:
        for key, count in chunk_total_counts.items():
            total_counts[key] += count
        for key, count in chunk_any_counts.items():
            any_counts[key] += count
    return total_counts, any_counts


def _get_process_pool_context():
    try:
        return multiprocessing.get_context("forkserver")
    except ValueError:
        return None


def _simulate_stochastic_maps(tree, states, fit, num_simulations, seed=None, threads=1):
    if num_simulations <= 0:
        raise ValueError(
            "'--n-sim' must be positive when '--stochastic-map-out' is specified."
        )
    try:
        threads = int(threads)
    except (TypeError, ValueError) as exc:
        raise ValueError("'--threads' must be an integer.") from exc
    if threads <= 0:
        raise ValueError("'--threads' must be positive.")
    node_to_branch_id = assign_branch_ids(tree)
    branches = [node for node in tree.traverse() if not node.is_root]
    branch_lengths = [float(node.dist) for node in branches]
    if "rate_matrix_by_node" in fit:
        uniformization_contexts = {}
        cache = {}
        for node in branches:
            matrix = fit["rate_matrix_by_node"][node]
            key = (id(matrix), float(node.dist))
            if key not in cache:
                cache[key] = _build_uniformization_context(matrix, node.dist)
            uniformization_contexts[node] = cache[key]
    else:
        uniformization_contexts = _build_uniformization_contexts(
            fit["rate_matrix"], branch_lengths
        )
    spec = _build_stochastic_map_spec(
        tree, fit, node_to_branch_id, uniformization_contexts
    )
    seed_sequences = _simulation_seed_sequence(seed, num_simulations)
    if threads == 1 or num_simulations == 1:
        total_counts, any_counts = _merge_simulation_counts(
            _simulate_stochastic_map_once(spec, seed_sequence)
            for seed_sequence in seed_sequences
        )
    else:
        max_workers = min(threads, num_simulations)
        seed_sequence_chunks = [
            seed_sequences[worker_index::max_workers]
            for worker_index in range(max_workers)
        ]
        process_pool_context = _get_process_pool_context()
        executor_kwargs = {"max_workers": max_workers}
        if process_pool_context is not None:
            executor_kwargs["mp_context"] = process_pool_context
        with ProcessPoolExecutor(**executor_kwargs) as executor:
            total_counts, any_counts = _merge_count_summaries(
                executor.map(
                    _simulate_stochastic_map_chunk_worker,
                    (
                        (spec, seed_sequence_chunk)
                        for seed_sequence_chunk in seed_sequence_chunks
                    ),
                )
            )
    rows = list()
    for node in tree.traverse():
        if node.is_root:
            continue
        branch_id = node_to_branch_id[node]
        for from_index, from_state in enumerate(states):
            for to_index, to_state in enumerate(states):
                if from_index == to_index:
                    continue
                count_key = (branch_id, from_index, to_index)
                total_count = total_counts.get(count_key, 0)
                any_count = any_counts.get(count_key, 0)
                rows.append(
                    {
                        "branch_id": branch_id,
                        "parent": node_to_branch_id[node.up],
                        "node_class": get_node_class(node),
                        "name": "" if node.name in [None, ""] else str(node.name),
                        "from_state": from_state,
                        "to_state": to_state,
                        "total_count": total_count,
                        "mean_count": total_count / float(num_simulations),
                        "posterior_frequency": any_count / float(num_simulations),
                        "num_simulations": num_simulations,
                    }
                )
    columns = (
        "branch_id",
        "parent",
        "node_class",
        "name",
        "from_state",
        "to_state",
        "total_count",
        "mean_count",
        "posterior_frequency",
        "num_simulations",
    )
    return pd.DataFrame(rows, columns=columns)


def _write_stochastic_map(tree, states, fit, args):
    stochastic_map_out = getattr(args, "stochastic_map_out", None)
    if stochastic_map_out in ["", None]:
        return
    table = _simulate_stochastic_maps(
        tree=tree,
        states=states,
        fit=fit,
        num_simulations=getattr(args, "n_sim", 100),
        seed=getattr(args, "seed", None),
        threads=getattr(args, "threads", 1),
    )
    _write_table(table, stochastic_map_out)


def _validate_asr_output_paths(args):
    auxiliary_outputs = {
        "--model-out": getattr(args, "model_out", None),
        "--tree-out": getattr(args, "tree_out", None),
        "--stochastic-map-out": getattr(args, "stochastic_map_out", None),
        "--covariance-out": getattr(args, "covariance_out", None),
    }
    stdout_auxiliary_outputs = [
        option_name for option_name, path in auxiliary_outputs.items() if path == "-"
    ]
    if stdout_auxiliary_outputs:
        raise ValueError(
            "Auxiliary outputs require file paths, not STDOUT: {}".format(
                ", ".join(stdout_auxiliary_outputs)
            )
        )
    validate_distinct_output_paths(
        [
            ("--outfile", getattr(args, "outfile", None)),
            ("--model-out", getattr(args, "model_out", None)),
            ("--tree-out", getattr(args, "tree_out", None)),
            ("--stochastic-map-out", getattr(args, "stochastic_map_out", None)),
            ("--covariance-out", getattr(args, "covariance_out", None)),
        ]
    )


def asr_main(args):
    _validate_asr_output_paths(args)
    tree = read_tree(
        args.infile,
        args.format,
        args.quoted_node_names,
        rooted=getattr(args, "input_rooted", "auto"),
    )
    _validate_tree_for_asr(tree)
    trait_columns = asr_trait_columns(
        args.state_column, getattr(args, "model", None)
    )
    state_column_input = (
        trait_columns if getattr(args, "model", None) == "MV-BM" else args.state_column
    )
    trait_df = read_asr_table(
        args.trait,
        state_column_input,
        list(tree.leaf_names()),
        missing_values=getattr(args, "missing_values", None),
        unmatched=getattr(args, "unmatched", "warn"),
        standard_error_column=getattr(args, "standard_error_column", None),
    )
    requested_type = getattr(args, "trait_type", "auto")
    trait_type = resolve_trait_type(requested_type, trait_df, state_column_input)
    settings = AsrSettings.from_args(args, trait_type)
    effective = effective_asr_args(args, settings)
    sys.stderr.write(f"ASR trait type: {trait_type} ({requested_type}).\n")
    if requested_type == "auto" and trait_type == "continuous":
        sys.stderr.write(
            "ASR auto-detection treats numeric columns as continuous; use "
            "--trait-type discrete for numeric category codes.\n"
        )
    targets = _parse_targets(getattr(args, "target", DEFAULT_TARGET))
    handlers = {"discrete": _run_discrete_asr, "continuous": _run_continuous_asr}
    handlers[trait_type](tree, trait_df, effective, settings, targets)


def _run_discrete_asr(tree, trait_df, args, settings, targets):
    from nwkit.asr_regimes import read_regime_map

    leaf_names = list(tree.leaf_names())
    fixed_rate_matrix = None
    states_arg = getattr(args, "states", None)
    rate_matrix_source = (
        "fixed:--rate"
        if settings.model == "ER" and getattr(args, "rate", None) is not None
        else "estimated"
    )
    if settings.model == "CUSTOM":
        matrix_states, fixed_rate_matrix = read_rate_matrix(args.rate_matrix)
        if (
            states_arg not in (None, "")
            and _parse_state_argument(states_arg) != matrix_states
        ):
            raise ValueError(
                "--states must exactly match the state order in --rate-matrix."
            )
        states_arg = matrix_states
        rate_matrix_source = f"fixed:{args.rate_matrix}"
    states, observed_state_by_leaf, likelihood_by_leaf = _read_tip_states(
        trait_path=args.trait,
        state_column=args.state_column,
        tree_leaf_names=leaf_names,
        states_arg=states_arg,
        missing_values_arg=getattr(args, "missing_values", None),
        ambiguous_separator=getattr(
            args, "ambiguous_separator", DEFAULT_AMBIGUOUS_SEPARATOR
        ),
        unmatched=getattr(args, "unmatched", "warn"),
        trait_df=trait_df,
        state_source="--rate-matrix" if settings.model == "CUSTOM" else "--states",
    )
    if settings.model == "CUSTOM":
        transition_graph = None
        transition_graph_source = "fixed-Q"
    else:
        transition_graph, transition_graph_source = read_transition_graph(
            getattr(args, "transition_graph", None),
            states,
            state_source="--states" if states_arg not in (None, "") else "--trait",
        )
    regime_assignment = (
        read_regime_map(getattr(args, "regime_map", None), tree)
        if settings.model == "MK-REGIME"
        else None
    )
    rate_bounds = _parse_rate_bounds(getattr(args, "rate_bounds", None))
    common_options = {
        "tree": tree,
        "states": states,
        "observed_state_by_leaf": observed_state_by_leaf,
        "likelihood_by_leaf": likelihood_by_leaf,
        "rate": getattr(args, "rate", None),
        "root_prior_mode": settings.root_prior,
        "rate_bounds": rate_bounds,
        "transition_graph": transition_graph,
    }
    if settings.model == "HRM":
        posterior_by_node, fit = compute_hrm_marginals(
            **common_options,
            hidden_categories=getattr(args, "hidden_categories", 2),
        )
    else:
        posterior_by_node, fit = compute_mk_marginals(
            **common_options,
            model=settings.model,
            fixed_rate_matrix=fixed_rate_matrix,
            regime_assignment=regime_assignment,
            regime_model=getattr(args, "regime_model", None) or "ER",
        )
    fit["transition_graph"] = transition_graph_source
    fit["rate_matrix_source"] = rate_matrix_source
    fit["trait_type_requested"] = getattr(args, "trait_type", "auto")
    if fit["rate_estimated"] and (
        fit.get("fit_status") != "ok" or fit.get("optimizer_failed_starts", 0)
    ):
        sys.stderr.write(
            "Discrete ASR: Mk fit status={}; {}/{} optimizer starts converged.\n".format(
                fit.get("fit_status", "unknown"),
                fit.get("optimizer_converged_starts", 0),
                fit.get("optimizer_starts", 0),
            )
        )
    table = _build_output_table(
        tree=tree,
        states=states,
        observed_state_by_leaf=observed_state_by_leaf,
        posterior_by_node=posterior_by_node,
        targets=targets,
        output_mode=settings.output,
    )
    _write_table(table, args.outfile)
    _write_model_table(
        states,
        settings.root_prior,
        fit,
        getattr(args, "model_out", None),
    )
    _write_annotated_tree(tree, states, posterior_by_node, observed_state_by_leaf, args)
    _write_stochastic_map(tree, states, fit, args)


def _run_continuous_asr(tree, trait_df, args, settings, targets):
    from nwkit.asr_regimes import read_regime_map, read_regime_parameters
    from nwkit.continuous_asr import compute_bm_marginals
    from nwkit.continuous_asr_io import (
        continuous_model_table,
        continuous_output_table,
        write_continuous_tree,
    )

    if settings.model == "MV-BM":
        trait_columns = asr_trait_columns(args.state_column, settings.model)
        observed = continuous_tip_vectors(
            trait_df, trait_columns, list(tree.leaf_names())
        )
        errors = None
    else:
        trait_columns = (args.state_column,)
        observed, errors = continuous_tip_values(
            trait_df,
            args.state_column,
            list(tree.leaf_names()),
            getattr(args, "standard_error_column", None),
        )
    regime_assignment = (
        read_regime_map(getattr(args, "regime_map", None), tree)
        if settings.model in {"BMS", "OUM"}
        else None
    )
    if settings.model == "BM":
        posterior, fit = compute_bm_marginals(
            tree,
            observed,
            sigma2=getattr(args, "sigma2", None),
            standard_errors=errors,
            _tree_validated=True,
        )
    elif settings.model == "BMS":
        from nwkit.regime_continuous_asr import compute_bms_marginals

        assert regime_assignment is not None
        fixed_parameters = read_regime_parameters(
            getattr(args, "regime_parameters", None),
            regime_assignment.regimes,
            ("sigma2",),
        )
        parameter_values = None if fixed_parameters is None else fixed_parameters[0]
        posterior, fit = compute_bms_marginals(
            tree,
            observed,
            regime_assignment,
            sigma2=getattr(args, "sigma2", None),
            sigma2_by_regime=None
            if parameter_values is None
            else {
                regime: values["sigma2"]
                for regime, values in parameter_values.items()
            },
            regime_parameters_source=None
            if fixed_parameters is None
            else fixed_parameters[1],
            standard_errors=errors,
            _tree_validated=True,
        )
    elif settings.model == "OU":
        from nwkit.ou_asr import compute_ou_marginals, parse_alpha_bounds

        posterior, fit = compute_ou_marginals(
            tree,
            observed,
            alpha=getattr(args, "alpha", None),
            sigma2=getattr(args, "sigma2", None),
            theta=getattr(args, "theta", None),
            alpha_bounds=parse_alpha_bounds(getattr(args, "alpha_bounds", None), tree),
            standard_errors=errors,
            _tree_validated=True,
        )
    elif settings.model == "OUM":
        from nwkit.ou_asr import parse_alpha_bounds
        from nwkit.regime_continuous_asr import compute_oum_marginals

        assert regime_assignment is not None
        fixed_parameters = read_regime_parameters(
            getattr(args, "regime_parameters", None),
            regime_assignment.regimes,
            ("theta",),
        )
        parameter_values = None if fixed_parameters is None else fixed_parameters[0]
        posterior, fit = compute_oum_marginals(
            tree,
            observed,
            regime_assignment,
            alpha=getattr(args, "alpha", None),
            sigma2=getattr(args, "sigma2", None),
            theta=getattr(args, "theta", None),
            theta_by_regime=None
            if parameter_values is None
            else {
                regime: values["theta"]
                for regime, values in parameter_values.items()
            },
            regime_parameters_source=None
            if fixed_parameters is None
            else fixed_parameters[1],
            alpha_bounds=parse_alpha_bounds(
                getattr(args, "alpha_bounds", None), tree
            ),
            standard_errors=errors,
            _tree_validated=True,
        )
    elif settings.model == "EB":
        from nwkit.nonstationary_continuous_asr import (
            compute_eb_marginals,
            parse_eb_rate_bounds,
        )

        posterior, fit = compute_eb_marginals(
            tree,
            observed,
            sigma2=getattr(args, "sigma2", None),
            eb_rate=getattr(args, "eb_rate", None),
            eb_rate_bounds=parse_eb_rate_bounds(
                getattr(args, "eb_rate_bounds", None), tree
            ),
            standard_errors=errors,
            _tree_validated=True,
        )
    elif settings.model == "BM-DRIFT":
        from nwkit.nonstationary_continuous_asr import compute_bm_drift_marginals

        posterior, fit = compute_bm_drift_marginals(
            tree,
            observed,
            sigma2=getattr(args, "sigma2", None),
            drift=getattr(args, "drift", None),
            standard_errors=errors,
            _tree_validated=True,
        )
    elif settings.model == "MV-BM":
        from nwkit.multivariate_asr import compute_mvbm_marginals

        posterior, fit = compute_mvbm_marginals(
            tree,
            observed,
            trait_columns,
            _tree_validated=True,
        )
    else:
        raise ValueError(f"Unsupported continuous ASR model: {settings.model}.")
    selected = [
        node for node in tree.traverse() if _should_output_node(node, observed, targets)
    ]
    if settings.model == "MV-BM":
        from nwkit.multivariate_asr import (
            multivariate_covariance_table,
            multivariate_model_table,
            multivariate_output_table,
            write_multivariate_tree,
        )

        table = multivariate_output_table(
            tree,
            selected,
            observed,
            posterior,
            trait_columns,
            settings.ci_level,
        )
    else:
        table = continuous_output_table(
            tree,
            selected,
            observed,
            errors,
            posterior,
            trait=args.state_column,
            ci_level=settings.ci_level,
        )
    if settings.model == "BM" and fit.sigma2_estimated and fit.fit_status != "ok":
        sys.stderr.write(
            f"Continuous ASR: sigma2=0 ({fit.fit_status}); intervals condition on this rate "
            "and exclude rate-estimation uncertainty.\n"
        )
    elif settings.model == "MV-BM" and fit.fit_status != "ok":
        sys.stderr.write(
            f"Continuous ASR: MV-BM fit status={fit.fit_status}; marginal "
            "reconstruction is available, but the ordinary multivariate likelihood is not.\n"
        )
    elif settings.model in {"BMS", "OU", "OUM", "EB", "BM-DRIFT"} and (
        fit.fit_status != "ok" or not fit.optimizer_success
    ):
        sys.stderr.write(
            f"Continuous ASR: {settings.model} fit status={fit.fit_status}; intervals condition "
            "on fitted parameters and exclude parameter-estimation uncertainty.\n"
        )
    _write_table(table, args.outfile)
    if settings.model == "MV-BM" and getattr(args, "covariance_out", None) not in (
        None,
        "",
    ):
        _write_table(
            multivariate_covariance_table(
                tree, selected, posterior, trait_columns
            ),
            args.covariance_out,
        )
    if getattr(args, "model_out", None) not in (None, ""):
        model_table = (
            multivariate_model_table(fit, args, settings.ci_level)
            if settings.model == "MV-BM"
            else continuous_model_table(fit, args, settings.ci_level)
        )
        _write_table(model_table, args.model_out)
    if settings.model == "MV-BM":
        write_multivariate_tree(
            tree, observed, posterior, trait_columns, args, settings.ci_level
        )
    else:
        write_continuous_tree(
            tree, observed, errors, posterior, args, settings.ci_level
        )

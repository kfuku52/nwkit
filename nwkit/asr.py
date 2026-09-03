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
from scipy.stats import gamma as gamma_distribution
from scipy.stats import poisson

from nwkit.asr_input import (
    AsrSettings,
    asr_standard_error_columns,
    asr_trait_columns,
    continuous_tip_values,
    continuous_tip_vector_errors,
    continuous_tip_vectors,
    effective_asr_args,
    read_asr_table,
    resolve_trait_type,
)
from nwkit.asr_models import model_names
from nwkit.discrete_asr_models import (
    FREQUENCY_RATIO_BOUNDS,
    model_equivalence_family,
    read_rate_matrix,
    read_transition_graph,
    stationary_distribution,
    validate_rate_matrix,
)
from nwkit.discrete_asr_models import build_rate_matrix as _build_structured_rate_matrix
from nwkit.discrete_asr_models import initial_parameters as _initial_model_parameters
from nwkit.discrete_asr_models import parameter_kinds as _structured_parameter_kinds
from nwkit.discrete_asr_models import parameter_labels as _structured_parameter_labels
from nwkit.optimization import deterministic_multistart
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
_MAX_STOCHASTIC_MAP_WORK = 2_000_000
_MAX_CACHED_BACKWARD_BYTES = 256 * 1024
_MAX_BACKWARD_BYTES = 256 * 1024 * 1024
_MAX_DIRECT_EXPONENT_NORM = 32.0
_MAX_HIDDEN_EXPANDED_STATES = 64
_MAX_COVARION_LOG_SPREAD = math.log(1e4)
_MAX_FREE_TRANSITION_PARAMETERS = 256


def _parse_comma_list(value, option_name):
    if value in ["", None]:
        return list()
    items = [item.strip() for item in str(value).split(",")]
    if any(item == "" for item in items):
        raise ValueError("'{}' contains an empty item.".format(option_name))
    return items


def _integer_option(value, option_name):
    if isinstance(value, (bool, np.bool_)):
        raise ValueError(f"{option_name} must be an integer.")
    try:
        number = float(value)
        integer = int(number)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError(f"{option_name} must be an integer.") from exc
    if not math.isfinite(number) or integer != number:
        raise ValueError(f"{option_name} must be an integer.")
    return integer


def _parse_missing_values(value):
    return parse_table_missing_values(value)


def _validate_states(states):
    if len(states) != len(set(states)):
        raise ValueError("'--states' contains duplicated states.")
    return states


def _require_multiple_discrete_states(states):
    if len(states) < 2:
        raise ValueError(
            "Discrete ASR requires at least two states; provide an explicit "
            "multi-state space with --states when only one state is observed."
        )


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
    _require_multiple_discrete_states(states)

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
    exponent = -float(num_states) * float(rate) * float(branch_length)
    decay = math.exp(exponent)
    # exp() rounds to one for tiny negative exponents.  expm1() preserves the
    # correspondingly tiny but non-zero transition probabilities.
    off_diagonal = -math.expm1(exponent) / float(num_states)
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
    branch_length = _validated_branch_length(branch_length)
    if branch_length == 0.0:
        return np.eye(rate_matrix.shape[0], dtype=float)
    er_rate = _get_er_rate_from_matrix(rate_matrix)
    if er_rate is not None:
        return _er_transition_matrix(branch_length, er_rate, rate_matrix.shape[0])
    return _general_transition_matrix(rate_matrix, branch_length)


def _validated_branch_length(branch_length):
    try:
        value = float(branch_length)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError("Mk branch lengths must be non-negative and finite.") from exc
    if not math.isfinite(value) or value < 0.0:
        raise ValueError("Mk branch lengths must be non-negative and finite.")
    return value


def _general_transition_matrix(rate_matrix, branch_length):
    absolute_rates = np.abs(rate_matrix)
    rate_scale = float(np.max(absolute_rates))
    if rate_scale == 0.0:
        return np.eye(rate_matrix.shape[0], dtype=float)
    scaled_rate_norm = float(np.max(np.sum(absolute_rates / rate_scale, axis=1)))
    log2_rate_norm = math.log2(rate_scale) + math.log2(scaled_rate_norm)
    log2_exponent_norm = log2_rate_norm + math.log2(branch_length)
    if log2_exponent_norm <= math.log2(_MAX_DIRECT_EXPONENT_NORM):
        exponent_norm = 2.0**log2_exponent_norm
        return _validated_transition_matrix(
            expm(rate_matrix * branch_length), exponent_norm
        )

    # scipy.linalg.expm is backward stable, but on very long CTMC branches its
    # generic scaling-and-squaring path can accumulate row-sum error large
    # enough to destroy stochasticity.  Exponentiate a moderate time interval,
    # then use the Markov semigroup while projecting only roundoff-sized row
    # errors after every square.  Computing the scale in log space also avoids
    # overflowing rate_norm * branch_length for finite inputs.
    squarings = max(
        1,
        int(math.ceil(log2_exponent_norm - math.log2(_MAX_DIRECT_EXPONENT_NORM))),
    )
    scaled_length = math.ldexp(branch_length, -squarings)
    scaled_norm = rate_scale * scaled_length * scaled_rate_norm
    matrix = _validated_transition_matrix(
        expm(rate_matrix * scaled_length), scaled_norm
    )
    exponent_norm = 2.0**log2_exponent_norm if log2_exponent_norm < 1023.0 else math.inf
    for _ in range(squarings):
        previous = matrix
        matrix = _validated_transition_matrix(matrix @ matrix, exponent_norm)
        # Once a finite-state CTMC has reached its idempotent long-time limit,
        # further squaring cannot change the result beyond floating-point
        # roundoff.  This also bounds work for extreme finite rate-time scales.
        if np.array_equal(matrix, previous):
            break
    return matrix


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
    tree,
    likelihood_by_leaf,
    rate_matrix,
    *,
    log_space=False,
    fixed_transition_matrices=None,
):
    inside, log_scales, child_terms, matrices = _compute_log_inside_likelihoods(
        tree,
        likelihood_by_leaf,
        rate_matrix,
        fixed_transition_matrices=fixed_transition_matrices,
    )
    if not log_space:
        inside = {node: np.exp(values) for node, values in inside.items()}
        child_terms = {
            node: {child: np.exp(values) for child, values in terms.items()}
            for node, terms in child_terms.items()
        }
    return inside, log_scales, child_terms, matrices


def _compute_log_inside_likelihoods(
    tree, likelihood_by_leaf, rate_matrix, *, fixed_transition_matrices=None
):
    if isinstance(rate_matrix, dict):
        sample_matrix = next(iter(rate_matrix.values()))
    else:
        sample_matrix = rate_matrix
        if fixed_transition_matrices is None:
            fixed_transition_matrices = _transition_matrices_for_tree(tree, rate_matrix)
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
            if fixed_transition_matrices is not None:
                matrix = fixed_transition_matrices[child]
                log_matrix = _log_probabilities(matrix)
            else:
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


def _log_likelihood(
    tree,
    likelihood_by_leaf,
    root_prior,
    rate_matrix,
    *,
    fixed_transition_matrices=None,
):
    inside, log_scales, _, _ = _compute_inside_likelihoods(
        tree=tree,
        likelihood_by_leaf=likelihood_by_leaf,
        rate_matrix=rate_matrix,
        log_space=True,
        fixed_transition_matrices=fixed_transition_matrices,
    )
    root_term = float(logsumexp(_log_probabilities(root_prior) + inside[tree]))
    return log_scales[tree] + root_term


def _is_informative_tip_likelihood(likelihood):
    """Return whether a tip likelihood distinguishes at least two states."""

    values = np.asarray(likelihood, dtype=float)
    return values.ndim == 1 and len(values) > 1 and bool(np.any(values != values[0]))


def _has_informative_tip_likelihood(likelihood_by_leaf):
    return any(
        _is_informative_tip_likelihood(likelihood)
        for likelihood in likelihood_by_leaf.values()
    )


def _informative_tip_likelihood_count(likelihood_by_leaf):
    return sum(
        _is_informative_tip_likelihood(likelihood)
        for likelihood in likelihood_by_leaf.values()
    )


def _informative_descendants(tree, likelihood_by_leaf):
    informative = {}
    for node in tree.traverse(strategy="postorder"):
        if node.is_leaf:
            informative[node] = _is_informative_tip_likelihood(
                likelihood_by_leaf[node.name]
            )
        else:
            informative[node] = any(informative[child] for child in node.children)
    return informative


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


def _rate_root_prior(root_prior, root_prior_factory, matrix):
    return root_prior if root_prior_factory is None else root_prior_factory(matrix)


def _fixed_parametric_rate_fit(
    tree,
    model,
    states,
    likelihood_by_leaf,
    rates,
    rate_bounds,
    transition_graph,
    prior_for,
):
    rate_matrix = _build_rate_matrix(model, states, rates, transition_graph)
    final_prior = prior_for(rate_matrix)
    return {
        "model": model,
        "rates": np.asarray(rates, dtype=float),
        "rate_matrix": rate_matrix,
        "log_likelihood": _log_likelihood(
            tree, likelihood_by_leaf, final_prior, rate_matrix
        ),
        "rate_estimated": False,
        "rate_bounds": rate_bounds,
        "root_prior": final_prior,
        "fit_status": "fixed",
    }


def _validate_initial_rate_likelihood(
    tree, states, likelihood_by_leaf, model, values, transition_graph, prior_for
):
    matrix = _build_rate_matrix(model, states, values, transition_graph)
    try:
        value = _log_likelihood(tree, likelihood_by_leaf, prior_for(matrix), matrix)
    except (ValueError, ArithmeticError):
        return
    if not math.isfinite(value):
        raise ValueError(
            "Observed states have zero likelihood under the selected transition "
            "graph and root prior."
        )


def _parametric_rate_starts(lower_logs, upper_logs, initial_logs):
    num_params = len(initial_logs)
    starts = [
        lower_logs + fraction * (upper_logs - lower_logs)
        for fraction in np.linspace(0.0, 1.0, 7)
    ]
    if num_params <= 1:
        return starts
    lower_quartiles = lower_logs + 0.25 * (upper_logs - lower_logs)
    upper_quartiles = lower_logs + 0.75 * (upper_logs - lower_logs)
    starts.extend(
        [
            np.where(np.arange(num_params) % 2 == 0, lower_quartiles, upper_quartiles),
            np.where(np.arange(num_params) % 2 == 0, upper_quartiles, lower_quartiles),
            lower_quartiles
            + np.linspace(0.0, 1.0, num_params) * (upper_quartiles - lower_quartiles),
            upper_quartiles
            - np.linspace(0.0, 1.0, num_params) * (upper_quartiles - lower_quartiles),
        ]
    )
    if num_params <= 4:
        for index in range(num_params):
            for value in (lower_quartiles[index], upper_quartiles[index]):
                start = initial_logs.copy()
                start[index] = value
                starts.append(start)
    return starts


def _parameter_boundary_summary(parameters, parameter_kinds, parameter_bounds):
    tolerance = 1e-5
    values = np.asarray(parameters, dtype=float)
    lower = np.asarray([item[0] for item in parameter_bounds], dtype=float)
    upper = np.asarray([item[1] for item in parameter_bounds], dtype=float)
    rate_mask = np.asarray([kind == "rate" for kind in parameter_kinds], dtype=bool)
    frequency_mask = ~rate_mask
    lower_hits = values <= lower * (1.0 + tolerance)
    upper_hits = values >= upper * (1.0 - tolerance)
    counts: dict[str, int | str] = {
        "num_rates_at_lower_bound": int(np.sum(rate_mask & lower_hits)),
        "num_rates_at_upper_bound": int(np.sum(rate_mask & upper_hits)),
        "num_frequencies_at_lower_bound": int(np.sum(frequency_mask & lower_hits)),
        "num_frequencies_at_upper_bound": int(np.sum(frequency_mask & upper_hits)),
    }
    status_labels = (
        ("num_rates_at_lower_bound", "rate_lower_boundary"),
        ("num_rates_at_upper_bound", "rate_upper_boundary"),
        ("num_frequencies_at_lower_bound", "frequency_lower_boundary"),
        ("num_frequencies_at_upper_bound", "frequency_upper_boundary"),
    )
    statuses = [label for name, label in status_labels if counts[name]]
    counts["fit_status"] = "+".join(statuses) if statuses else "ok"
    return counts


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
        return _rate_root_prior(root_prior, root_prior_factory, matrix)

    num_params = _num_rate_parameters(model, len(states), transition_graph)
    if num_params > _MAX_FREE_TRANSITION_PARAMETERS:
        raise ValueError(
            f"{model} would require more than {_MAX_FREE_TRANSITION_PARAMETERS} "
            "free transition parameters; reduce the state space or constrain "
            "the transition graph."
        )
    if num_params == 0:
        return _fixed_parametric_rate_fit(
            tree,
            model,
            states,
            likelihood_by_leaf,
            [],
            rate_bounds,
            transition_graph,
            prior_for,
        )

    if model == "ER" and rate is not None:
        fixed_rate = float(rate)
        if (not math.isfinite(fixed_rate)) or fixed_rate < 0.0:
            raise ValueError("'--rate' must be a non-negative finite number.")
        return _fixed_parametric_rate_fit(
            tree,
            model,
            states,
            likelihood_by_leaf,
            [fixed_rate],
            rate_bounds,
            transition_graph,
            prior_for,
        )

    initial_rate = _initial_rate_value(tree, rate, rate_bounds)
    parameter_kinds = _structured_parameter_kinds(model, states, transition_graph)
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

    _validate_initial_rate_likelihood(
        tree,
        states,
        likelihood_by_leaf,
        model,
        initial_values,
        transition_graph,
        prior_for,
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

    try:
        optimized = deterministic_multistart(
            objective,
            initial_logs,
            list(zip(lower_logs, upper_logs, strict=True)),
            maxiter=500,
            minimizer=minimize,
            ftol=None,
            additional_starts=_parametric_rate_starts(
                lower_logs, upper_logs, initial_logs
            ),
        )
    except ValueError as exc:
        raise ValueError(f"Failed to estimate Mk model parameters: {exc}") from exc
    rates = np.exp(optimized.x)
    rate_matrix = _build_rate_matrix(model, states, rates, transition_graph)
    final_prior = prior_for(rate_matrix)
    log_likelihood = _log_likelihood(tree, likelihood_by_leaf, final_prior, rate_matrix)
    result = {
        "model": model,
        "rates": rates,
        "rate_matrix": rate_matrix,
        "log_likelihood": log_likelihood,
        "rate_estimated": True,
        "rate_bounds": rate_bounds,
        "optimizer_success": optimized.success,
        "optimizer_message": optimized.message,
        "optimizer_starts": optimized.starts,
        "optimizer_converged_starts": optimized.converged_starts,
        "optimizer_failed_starts": optimized.failed_starts,
        "root_prior": final_prior,
    }
    result.update(_parameter_boundary_summary(rates, parameter_kinds, parameter_bounds))
    return result


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
    num_parameters = (
        1
        if model == "MK-REGIME"
        else _num_rate_parameters(model, len(states), transition_graph)
    )
    fully_fixed = num_parameters == 0 or (model == "ER" and rate is not None)
    if not _has_informative_tip_likelihood(likelihood_by_leaf) and not fully_fixed:
        raise ValueError(
            "Prior-only discrete ASR requires a fully fixed transition process; "
            "use ER with --rate or CUSTOM with --rate-matrix."
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


def _validate_informative_regimes(tree, likelihood_by_leaf, regime_assignment):
    informative_descendant = _informative_descendants(tree, likelihood_by_leaf)
    represented_edges = {
        regime_assignment.by_node[node]
        for node in tree.traverse()
        if not node.is_root and float(node.dist) > 0.0 and informative_descendant[node]
    }
    uninformative = sorted(set(regime_assignment.regimes) - represented_edges)
    if uninformative:
        raise ValueError(
            "Every estimated regime must label at least one positive-length "
            "branch ancestral to an informative observation; uninformative "
            "regime(s): " + ", ".join(uninformative)
        )


def _regime_rate_starts(initial_logs, lower_logs, upper_logs, num_regimes, per_regime):
    starts = [
        lower_logs + fraction * (upper_logs - lower_logs)
        for fraction in np.linspace(0.0, 1.0, 7)
    ]
    if num_regimes <= 1:
        return starts
    for fraction in (0.25, 0.75):
        start = initial_logs.copy()
        for regime_index in range(num_regimes):
            if regime_index % 2:
                left = regime_index * per_regime
                right = left + per_regime
                start[left:right] = lower_logs[left:right] + fraction * (
                    upper_logs[left:right] - lower_logs[left:right]
                )
        starts.append(start)
    return starts


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
        raise ValueError("--regime-model must be one of ER, SYM, ARD, F81, or GTR.")
    regimes = regime_assignment.regimes
    _validate_informative_regimes(tree, likelihood_by_leaf, regime_assignment)
    labels = _rate_parameter_labels(regime_model, states, transition_graph)
    if not labels:
        raise ValueError("MK-REGIME requires at least one rate per regime.")
    kinds = _structured_parameter_kinds(regime_model, states, transition_graph)
    per_regime = len(labels)
    total_parameters = per_regime * len(regimes)
    if total_parameters > _MAX_FREE_TRANSITION_PARAMETERS:
        raise ValueError(
            "MK-REGIME would require more than "
            f"{_MAX_FREE_TRANSITION_PARAMETERS} free transition parameters; reduce "
            "the state space, number of regimes, or transition graph."
        )
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
    optimizer_bounds = list(zip(lower_logs, upper_logs, strict=True))
    try:
        optimized = deterministic_multistart(
            objective,
            initial_logs,
            optimizer_bounds,
            maxiter=800,
            minimizer=minimize,
            ftol=None,
            additional_starts=_regime_rate_starts(
                initial_logs,
                lower_logs,
                upper_logs,
                len(regimes),
                per_regime,
            ),
        )
    except ValueError as exc:
        raise ValueError(f"Failed to estimate MK-REGIME parameters: {exc}") from exc
    parameters = np.exp(optimized.x)
    by_regime, by_node = matrices_for(parameters)
    final_prior = prior_for(by_regime)
    log_likelihood = _log_likelihood(tree, likelihood_by_leaf, final_prior, by_node)
    rates_by_regime = {
        regime: parameters[index * per_regime : (index + 1) * per_regime]
        for index, regime in enumerate(regimes)
    }
    result = {
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
        "optimizer_success": optimized.success,
        "optimizer_message": optimized.message,
        "optimizer_starts": optimized.starts,
        "optimizer_converged_starts": optimized.converged_starts,
        "optimizer_failed_starts": optimized.failed_starts,
        "root_prior": final_prior,
    }
    result.update(_parameter_boundary_summary(parameters, all_kinds, bounds))
    return result


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
    _require_multiple_discrete_states(states)
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
            "sample_size": _informative_tip_likelihood_count(likelihood_by_leaf),
            "inside": {node: np.exp(values) for node, values in inside.items()},
            "log_inside": inside,
            "transition_matrices": transition_matrices,
            "posterior_by_node": posterior_by_node,
        }
    )
    return posterior_by_node, fit


def _posterior_for_rate_matrix(
    tree,
    likelihood_by_leaf,
    root_prior,
    rate_matrix,
    *,
    fixed_transition_matrices=None,
):
    inside, _, child_terms, transition_matrices = _compute_inside_likelihoods(
        tree=tree,
        likelihood_by_leaf=likelihood_by_leaf,
        rate_matrix=rate_matrix,
        log_space=True,
        fixed_transition_matrices=fixed_transition_matrices,
    )
    outside = _compute_outside_likelihoods(
        tree=tree,
        child_terms=child_terms,
        transition_matrices=transition_matrices,
        root_prior=root_prior,
        log_space=True,
    )
    return {
        node: _probabilities_from_logs(inside[node] + outside[node])
        for node in tree.traverse()
    }


def _transition_matrices_for_tree(tree, rate_matrix):
    """Build each unique branch-length transition once for one generator."""

    lengths = sorted({float(node.dist) for node in tree.traverse() if not node.is_root})
    # Preserve the closed-form ER path: it is both cheaper and more accurate
    # than a generic eigendecomposition for equal-rates generators.
    by_length = (
        None
        if _get_er_rate_from_matrix(rate_matrix) is not None
        else _symmetric_transition_matrices(rate_matrix, lengths)
    )
    if by_length is None:
        by_length = {
            length: _transition_matrix(rate_matrix, length) for length in lengths
        }
    by_node = {}
    for node in tree.traverse(strategy="preorder"):
        if node.is_root:
            continue
        length = float(node.dist)
        by_node[node] = by_length[length]
    return by_node


def _symmetric_transition_matrices(rate_matrix, lengths):
    """Exponentiate an exactly reversible symmetric generator once per Q."""

    matrix = np.asarray(rate_matrix, dtype=float)
    scale = float(np.max(np.abs(matrix), initial=0.0))
    tolerance = np.finfo(float).eps * max(1.0, scale) * max(100, len(matrix))
    if not np.allclose(matrix, matrix.T, rtol=0.0, atol=tolerance):
        return None
    eigenvalues, eigenvectors = np.linalg.eigh((matrix + matrix.T) / 2.0)
    if float(np.max(eigenvalues)) > tolerance:
        return None
    eigenvalues = np.minimum(eigenvalues, 0.0)
    row_norm = float(np.max(np.sum(np.abs(matrix), axis=1), initial=0.0))
    result = {}
    for raw_length in lengths:
        length = _validated_branch_length(raw_length)
        if length == 0.0:
            result[length] = np.eye(len(matrix), dtype=float)
            continue
        with np.errstate(over="ignore", under="ignore"):
            transition = (eigenvectors * np.exp(eigenvalues * length)) @ eigenvectors.T
        if row_norm == 0.0:
            exponent_norm = 0.0
        else:
            log_exponent_norm = math.log(row_norm) + math.log(length)
            exponent_norm = (
                math.exp(log_exponent_norm)
                if log_exponent_norm <= math.log(np.finfo(float).max)
                else math.inf
            )
        result[length] = _validated_transition_matrix(transition, exponent_norm)
    return result


def _discrete_gamma_category_rates(shape, categories):
    """Return Yang's equal-probability, conditional-mean gamma categories."""

    boundaries = np.concatenate(
        (
            np.zeros(1, dtype=float),
            gamma_distribution.ppf(
                np.arange(1, categories, dtype=float) / categories,
                a=shape,
                scale=1.0 / shape,
            ),
            np.full(1, np.inf, dtype=float),
        )
    )
    # For X ~ Gamma(shape, scale=1/shape), the partial first moment through x
    # is the CDF of Gamma(shape + 1, scale=1/shape), because E[X] = 1.
    first_moments = gamma_distribution.cdf(boundaries, a=shape + 1.0, scale=1.0 / shape)
    rates = categories * np.diff(first_moments)
    rates /= float(np.mean(rates))
    return rates


def compute_mk_mixture_marginals(
    tree,
    states,
    observed_by_character,
    likelihoods_by_character,
    *,
    base_model="ER",
    mixture="gamma",
    categories=4,
    rate=None,
    root_prior_mode="equal",
    rate_bounds=None,
    transition_graph=None,
):
    """Jointly fit several Mk characters with gamma or ordered free rates."""

    from nwkit.optimization import deterministic_multistart

    _require_multiple_discrete_states(states)
    if base_model not in {"ER", "SYM", "ARD", "F81", "GTR"}:
        raise ValueError("MK-MIXTURE base model must be ER, SYM, ARD, F81, or GTR.")
    if mixture not in {"gamma", "free"}:
        raise ValueError("MK-MIXTURE rate mixture must be gamma or free.")
    categories = _integer_option(categories, "--rate-categories")
    if not 2 <= categories <= 8:
        raise ValueError("--rate-categories must be between 2 and 8.")
    if len(likelihoods_by_character) < 2:
        raise ValueError("MK-MIXTURE requires at least two character columns.")
    rate_bounds = DEFAULT_RATE_BOUNDS if rate_bounds is None else rate_bounds
    labels = _rate_parameter_labels(base_model, states, transition_graph)
    if len(labels) > _MAX_FREE_TRANSITION_PARAMETERS:
        raise ValueError(
            f"MK-MIXTURE {base_model} would require more than "
            f"{_MAX_FREE_TRANSITION_PARAMETERS} free transition parameters; reduce "
            "the state space or constrain the transition graph."
        )
    kinds = _structured_parameter_kinds(base_model, states, transition_graph)
    if not labels:
        raise ValueError("MK-MIXTURE requires at least one base-rate parameter.")
    fixed_base = None
    if base_model == "ER" and rate is not None:
        fixed_rate = float(rate)
        if not math.isfinite(fixed_rate) or fixed_rate <= 0.0:
            raise ValueError(
                "'--rate' must be positive for MK-MIXTURE because a zero base "
                "rate makes the mixture distribution unidentifiable."
            )
        fixed_base = np.asarray([fixed_rate], dtype=float)

    initial_rate = _initial_rate_value(tree, rate, rate_bounds)
    initial_base = _initial_model_parameters(
        base_model, states, initial_rate, transition_graph
    )
    initial: list[float] = []
    bounds: list[tuple[float | None, float | None]] = []
    if fixed_base is None:
        initial.extend(np.log(initial_base))
        bounds.extend(
            (
                math.log(
                    FREQUENCY_RATIO_BOUNDS[0]
                    if kind == "frequency_ratio"
                    else rate_bounds[0]
                ),
                math.log(
                    FREQUENCY_RATIO_BOUNDS[1]
                    if kind == "frequency_ratio"
                    else rate_bounds[1]
                ),
            )
            for kind in kinds
        )
    if mixture == "gamma":
        mixture_offset = len(initial)
        initial.append(0.0)
        bounds.append((math.log(0.05), math.log(50.0)))
    else:
        mixture_offset = len(initial)
        initial_gap = math.log(2.0 / float(categories - 1))
        initial.extend([initial_gap] * (categories - 1))
        bounds.extend([(-5.0, 2.0)] * (categories - 1))
        initial.extend([0.0] * (categories - 1))
        bounds.extend([(-8.0, 8.0)] * (categories - 1))

    def mixture_parameters(parameters):
        if mixture == "gamma":
            shape = math.exp(float(parameters[mixture_offset]))
            category_rates = _discrete_gamma_category_rates(shape, categories)
            weights = np.full(categories, 1.0 / categories, dtype=float)
            return category_rates, weights, shape
        gap_parameters = np.asarray(
            parameters[mixture_offset : mixture_offset + categories - 1],
            dtype=float,
        )
        weight_logits = np.concatenate(
            (
                np.asarray(
                    parameters[
                        mixture_offset + categories - 1 : mixture_offset
                        + 2 * (categories - 1)
                    ],
                    dtype=float,
                ),
                np.zeros(1),
            )
        )
        log_rates = np.concatenate((np.zeros(1), np.cumsum(np.exp(gap_parameters))))
        weights = np.exp(weight_logits - float(np.max(weight_logits)))
        weights /= float(weights.sum())
        category_rates = np.exp(log_rates - float(np.max(log_rates)))
        category_rates /= float(weights @ category_rates)
        return category_rates, weights, None

    def unpack(parameters):
        base = (
            fixed_base
            if fixed_base is not None
            else np.exp(np.asarray(parameters[: len(labels)], dtype=float))
        )
        matrix = _build_rate_matrix(base_model, states, base, transition_graph)
        category_rates, weights, shape = mixture_parameters(parameters)
        return base, matrix, category_rates, weights, shape

    def prior_for(matrix):
        if root_prior_mode == "stationary":
            return stationary_distribution(matrix)
        if root_prior_mode == "equal":
            return np.full(len(states), 1.0 / len(states), dtype=float)
        if root_prior_mode != "empirical":
            raise ValueError(f"Unsupported MK-MIXTURE root prior: {root_prior_mode}.")
        counts = np.zeros(len(states), dtype=float)
        for character, likelihood_by_leaf in likelihoods_by_character.items():
            observed = observed_by_character[character]
            for leaf_name, value in observed.items():
                if value is None:
                    continue
                likelihood = likelihood_by_leaf[leaf_name]
                counts += likelihood / likelihood.sum()
        return (
            counts / counts.sum()
            if counts.sum() > 0.0
            else np.full(len(states), 1.0 / len(states), dtype=float)
        )

    def character_scores(
        matrix,
        category_rates,
        weights,
        root_prior,
        transition_matrices_by_category=None,
    ):
        if transition_matrices_by_category is None:
            transition_matrices_by_category = [
                _transition_matrices_for_tree(tree, matrix * category_rate)
                for category_rate in category_rates
            ]
        scores = {}
        total = 0.0
        for character, likelihood_by_leaf in likelihoods_by_character.items():
            values = np.asarray(
                [
                    math.log(weight)
                    + _log_likelihood(
                        tree,
                        likelihood_by_leaf,
                        root_prior,
                        matrix * category_rate,
                        fixed_transition_matrices=transition_matrices,
                    )
                    for category_rate, weight, transition_matrices in zip(
                        category_rates,
                        weights,
                        transition_matrices_by_category,
                        strict=True,
                    )
                ]
            )
            score = float(logsumexp(values))
            scores[character] = values
            total += score
        return total, scores, transition_matrices_by_category

    def objective(parameters):
        try:
            _, matrix, category_rates, weights, _ = unpack(parameters)
            likelihood, _, _ = character_scores(
                matrix, category_rates, weights, prior_for(matrix)
            )
        except (ValueError, ArithmeticError, OverflowError):
            return 1e100
        return -likelihood if math.isfinite(likelihood) else 1e100

    optimized = deterministic_multistart(objective, initial, bounds, maxiter=1500)
    base, matrix, category_rates, weights, shape = unpack(optimized.x)
    root_prior = prior_for(matrix)
    log_likelihood, scores, transition_matrices_by_category = character_scores(
        matrix, category_rates, weights, root_prior
    )
    posterior_by_character = {}
    category_posterior_by_character = {}
    for character, likelihood_by_leaf in likelihoods_by_character.items():
        category_probabilities = np.exp(
            scores[character] - logsumexp(scores[character])
        )
        category_posterior_by_character[character] = category_probabilities
        category_posteriors = [
            _posterior_for_rate_matrix(
                tree,
                likelihood_by_leaf,
                root_prior,
                matrix * category_rate,
                fixed_transition_matrices=transition_matrices,
            )
            for category_rate, transition_matrices in zip(
                category_rates, transition_matrices_by_category, strict=True
            )
        ]
        posterior_by_character[character] = {
            node: sum(
                probability * posterior[node]
                for probability, posterior in zip(
                    category_probabilities, category_posteriors, strict=True
                )
            )
            for node in tree.traverse()
        }

    tolerance = 1e-5
    statuses = []
    base_boundary = {}
    if fixed_base is None:
        base_bounds = [
            FREQUENCY_RATIO_BOUNDS if kind == "frequency_ratio" else rate_bounds
            for kind in kinds
        ]
        base_boundary = _parameter_boundary_summary(base, kinds, base_bounds)
        base_status = str(base_boundary.pop("fit_status"))
        if base_status != "ok":
            statuses.extend(base_status.split("+"))
    if mixture == "gamma" and shape is not None:
        if shape <= 0.05 * (1.0 + tolerance):
            statuses.append("gamma_shape_lower_boundary")
        elif shape >= 50.0 * (1.0 - tolerance):
            statuses.append("gamma_shape_upper_boundary")
    elif mixture == "free":
        gap_values = np.asarray(
            optimized.x[mixture_offset : mixture_offset + categories - 1]
        )
        weight_values = np.asarray(
            optimized.x[
                mixture_offset + categories - 1 : mixture_offset + 2 * (categories - 1)
            ]
        )
        if np.any(gap_values <= -5.0 + 7.0 * tolerance) or np.any(
            gap_values >= 2.0 - 7.0 * tolerance
        ):
            statuses.append("free_rate_gap_boundary")
        if np.any(weight_values <= -8.0 + 16.0 * tolerance) or np.any(
            weight_values >= 8.0 - 16.0 * tolerance
        ):
            statuses.append("mixture_weight_boundary")
    fit = {
        "model": "MK-MIXTURE",
        "mixture_model": base_model,
        "rate_mixture": mixture,
        "rate_categories": categories,
        "rates": base,
        "rate_matrix": matrix,
        "rate_estimated": True,
        "base_rate_estimated": fixed_base is None,
        "rate_bounds": rate_bounds,
        "mixture_category_rates": category_rates,
        "mixture_weights": weights,
        "gamma_shape": shape,
        "log_likelihood": log_likelihood,
        "sample_size": sum(
            _informative_tip_likelihood_count(likelihood_by_leaf)
            for likelihood_by_leaf in likelihoods_by_character.values()
        ),
        "root_prior": root_prior,
        "posterior_by_character": posterior_by_character,
        "category_posterior_by_character": category_posterior_by_character,
        "optimizer_success": optimized.success,
        "optimizer_message": optimized.message,
        "optimizer_starts": optimized.starts,
        "optimizer_converged_starts": optimized.converged_starts,
        "optimizer_failed_starts": optimized.failed_starts,
        "fit_status": "+".join(statuses) if statuses else "ok",
    }
    fit.update(base_boundary)
    return posterior_by_character, fit


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

    _require_multiple_discrete_states(states)
    hidden_categories = _integer_option(hidden_categories, "--hidden-categories")
    if hidden_categories < 2:
        raise ValueError("--hidden-categories must be at least 2 for HRM.")
    observed_graph = (
        np.ones((len(states), len(states)), dtype=bool)
        if transition_graph is None
        else np.array(transition_graph, dtype=bool, copy=True)
    )
    np.fill_diagonal(observed_graph, False)
    num_parameters = hidden_categories * int(np.sum(observed_graph)) + len(
        states
    ) * hidden_categories * (hidden_categories - 1)
    if num_parameters > _MAX_FREE_TRANSITION_PARAMETERS:
        raise ValueError(
            "HRM would require more than "
            f"{_MAX_FREE_TRANSITION_PARAMETERS} free transition rates; reduce the "
            "observed state space or --hidden-categories."
        )
    expanded_size = len(states) * hidden_categories
    if expanded_size > _MAX_HIDDEN_EXPANDED_STATES:
        raise ValueError(
            "HRM would require more than "
            f"{_MAX_HIDDEN_EXPANDED_STATES} expanded states; reduce the observed "
            "state space or --hidden-categories."
        )
    expanded_states = expanded_state_labels(states, hidden_categories)
    expanded_graph = expanded_transition_graph(observed_graph, hidden_categories)
    expanded_likelihoods = expand_tip_likelihoods(likelihood_by_leaf, hidden_categories)
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
        node: aggregate_probabilities(probabilities, len(states), hidden_categories)
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


def _covarion_parameter_setup(rate_bounds, initial_rate, fixed_rate):
    lower, upper = (float(value) for value in rate_bounds)
    if (
        not math.isfinite(lower)
        or not math.isfinite(upper)
        or lower <= 0.0
        or lower >= upper
    ):
        raise ValueError("COVARION rate bounds must be positive and increasing.")
    if fixed_rate is not None and not lower < fixed_rate < upper:
        raise ValueError(
            "A fixed COVARION --rate must lie strictly inside --rate-bounds so "
            "hidden-class rate spread remains identifiable."
        )
    log_bounds = math.log(lower), math.log(upper)
    if not log_bounds[0] < log_bounds[1]:
        raise ValueError(
            "COVARION rate bounds must be distinguishable on the fitted log scale."
        )
    center = math.log(fixed_rate if fixed_rate is not None else initial_rate)
    spread_limit = min(
        _MAX_COVARION_LOG_SPREAD,
        center - log_bounds[0],
        log_bounds[1] - center,
    )
    if spread_limit <= 0.0:
        center = (log_bounds[0] + log_bounds[1]) / 2.0
        spread_limit = min(
            _MAX_COVARION_LOG_SPREAD,
            center - log_bounds[0],
            log_bounds[1] - center,
        )
    spread_fraction = min(0.9, math.log(2.0) / spread_limit)
    names = []
    initial = []
    bounds = []
    if fixed_rate is None:
        names.append("log_base_rate")
        initial.append(center)
        bounds.append(log_bounds)
    names.extend(("spread_fraction", "log_switching_rate"))
    initial.extend((spread_fraction, math.log(initial_rate)))
    bounds.extend(((0.0, 1.0), log_bounds))
    return names, initial, bounds, log_bounds


def _unpack_covarion_parameters(parameters, names, fixed_rate, log_bounds):
    supplied = dict(zip(names, parameters, strict=True))
    base_log = (
        math.log(fixed_rate)
        if fixed_rate is not None
        else float(supplied["log_base_rate"])
    )
    spread_limit = max(
        0.0,
        min(
            _MAX_COVARION_LOG_SPREAD,
            base_log - log_bounds[0],
            log_bounds[1] - base_log,
        ),
    )
    spread_fraction = float(supplied["spread_fraction"])
    return (
        math.exp(base_log),
        spread_fraction * spread_limit,
        math.exp(float(supplied["log_switching_rate"])),
        spread_fraction,
    )


def compute_covarion_marginals(
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
    """Fit a parsimonious ordered hidden-rate covarion CTMC."""

    from nwkit.hidden_rate_asr import (
        aggregate_probabilities,
        covarion_rate_matrix,
        expand_tip_likelihoods,
        expanded_state_labels,
        state_projection,
    )
    from nwkit.optimization import deterministic_multistart

    _require_multiple_discrete_states(states)
    hidden_categories = _integer_option(hidden_categories, "--hidden-categories")
    if hidden_categories < 2:
        raise ValueError("--hidden-categories must be at least 2 for COVARION.")
    expanded_size = len(states) * hidden_categories
    if expanded_size > _MAX_HIDDEN_EXPANDED_STATES:
        raise ValueError(
            "COVARION would require more than "
            f"{_MAX_HIDDEN_EXPANDED_STATES} expanded states; reduce the observed "
            "state space or --hidden-categories."
        )
    observed_graph = (
        np.ones((len(states), len(states)), dtype=bool)
        if transition_graph is None
        else np.array(transition_graph, dtype=bool, copy=True)
    )
    np.fill_diagonal(observed_graph, False)
    if not np.any(observed_graph):
        raise ValueError("COVARION requires at least one observed-state transition.")
    rate_bounds = DEFAULT_RATE_BOUNDS if rate_bounds is None else rate_bounds
    initial_rate = _initial_rate_value(tree, rate, rate_bounds)
    fixed_rate = None
    if rate is not None:
        try:
            fixed_rate = float(rate)
        except (TypeError, ValueError, OverflowError) as exc:
            raise ValueError("'--rate' must be a non-negative finite number.") from exc
        if not math.isfinite(fixed_rate) or fixed_rate <= 0.0:
            raise ValueError(
                "'--rate' must be positive for COVARION because a zero base "
                "rate makes hidden-rate spread unidentifiable."
            )

    expanded_states = expanded_state_labels(states, hidden_categories)
    expanded_likelihoods = expand_tip_likelihoods(likelihood_by_leaf, hidden_categories)
    root_prior, root_prior_factory = _root_prior_configuration(
        root_prior_mode,
        expanded_states,
        observed_state_by_leaf,
        expanded_likelihoods,
    )
    names, initial, bounds, log_bounds = _covarion_parameter_setup(
        rate_bounds, initial_rate, fixed_rate
    )

    def unpack(parameters):
        return _unpack_covarion_parameters(parameters, names, fixed_rate, log_bounds)

    def matrix_for(parameters):
        return covarion_rate_matrix(
            observed_graph,
            hidden_categories,
            *unpack(parameters)[:3],
            effective_rate_bounds=rate_bounds,
        )

    def objective(parameters):
        try:
            matrix, _ = matrix_for(parameters)
            prior = (
                root_prior if root_prior_factory is None else root_prior_factory(matrix)
            )
            likelihood = _log_likelihood(tree, expanded_likelihoods, prior, matrix)
        except (ValueError, ArithmeticError, OverflowError):
            return 1e100
        return -likelihood if math.isfinite(likelihood) else 1e100

    optimized = deterministic_multistart(objective, initial, bounds, maxiter=1200)
    matrix, multipliers = matrix_for(optimized.x)
    base_rate, spread, switching_rate, spread_fraction = unpack(optimized.x)
    effective_rates = base_rate * multipliers
    expanded_posterior, fit = compute_mk_marginals(
        tree,
        expanded_states,
        observed_state_by_leaf,
        expanded_likelihoods,
        model="CUSTOM",
        root_prior_mode=root_prior_mode,
        fixed_rate_matrix=matrix,
    )
    posterior = {
        node: aggregate_probabilities(probabilities, len(states), hidden_categories)
        for node, probabilities in expanded_posterior.items()
    }
    tolerance = 1e-5
    statuses = []
    if spread <= tolerance:
        statuses.append("hidden_rate_homogeneous_boundary")
    elif spread_fraction >= 1.0 - tolerance:
        statuses.append("hidden_rate_spread_upper_boundary")
    if fixed_rate is None:
        if base_rate <= rate_bounds[0] * (1.0 + tolerance):
            statuses.append("base_rate_lower_boundary")
        elif base_rate >= rate_bounds[1] * (1.0 - tolerance):
            statuses.append("base_rate_upper_boundary")
    if switching_rate <= rate_bounds[0] * (1.0 + tolerance):
        statuses.append("switching_rate_lower_boundary")
    elif switching_rate >= rate_bounds[1] * (1.0 - tolerance):
        statuses.append("switching_rate_upper_boundary")
    if effective_rates[0] <= rate_bounds[0] * (1.0 + tolerance):
        statuses.append("effective_rate_lower_boundary")
    if effective_rates[-1] >= rate_bounds[1] * (1.0 - tolerance):
        statuses.append("effective_rate_upper_boundary")
    projection = state_projection(len(states), hidden_categories)
    expanded_root_prior = fit["root_prior"]
    fit.update(
        {
            "model": "COVARION",
            "rates": np.asarray([base_rate, spread, switching_rate]),
            "rate_estimated": True,
            "rate_bounds": rate_bounds,
            "base_rate": base_rate,
            "base_rate_estimated": fixed_rate is None,
            "log_rate_spread": spread,
            "hidden_rate_multipliers": multipliers,
            "hidden_rate_effective_rates": effective_rates,
            "switching_rate": switching_rate,
            "hidden_categories": hidden_categories,
            "expanded_states": expanded_states,
            "state_projection": projection,
            "expanded_root_prior": expanded_root_prior,
            "root_prior": aggregate_probabilities(
                expanded_root_prior, len(states), hidden_categories
            ),
            "expanded_posterior_by_node": expanded_posterior,
            "posterior_by_node": expanded_posterior,
            "observed_transition_graph": observed_graph,
            "optimizer_success": optimized.success,
            "optimizer_message": optimized.message,
            "optimizer_starts": optimized.starts,
            "optimizer_converged_starts": optimized.converged_starts,
            "optimizer_failed_starts": optimized.failed_starts,
            "fit_status": "+".join(statuses) if statuses else "ok",
        }
    )
    return posterior, fit


def _is_missing_tip(node, observed_state_by_leaf):
    return node.is_leaf and _is_missing_observation(
        observed_state_by_leaf.get(node.name)
    )


def _is_missing_observation(value):
    if value is None:
        return True
    return isinstance(value, (tuple, list)) and any(item is None for item in value)


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


def _threshold_model_row(states, root_prior_mode, fit):
    threshold_fit = fit["threshold_fit"]
    return {
        "trait_type": "discrete",
        "trait_type_requested": fit.get("trait_type_requested", "discrete"),
        "model": "THRESHOLD",
        "num_states": len(states),
        "states": ",".join(states),
        "root_prior": root_prior_mode,
        "liability_process": "standardized_brownian",
        "liability_root_mean": 0.0,
        "liability_root_variance": 1.0,
        "liability_sigma2": 1.0,
        "thresholds": ",".join(str(value) for value in threshold_fit.thresholds),
        "thresholds_estimated": threshold_fit.thresholds_estimated,
        "estimation_method": "data_augmentation_mcmc",
        "likelihood_kind": "posterior_sampling_no_marginal_likelihood",
        "mcmc_samples_per_chain": threshold_fit.num_samples,
        "mcmc_burnin": threshold_fit.burnin,
        "mcmc_thin": threshold_fit.thin,
        "mcmc_chains": threshold_fit.chains,
        "mcmc_seed": "" if threshold_fit.seed is None else threshold_fit.seed,
        "mcmc_rhat_max": threshold_fit.rhat_max,
        "mcmc_ess_min": threshold_fit.ess_min,
        "fit_status": threshold_fit.fit_status,
    }


def _base_discrete_model_row(states, root_prior_mode, fit):
    frequency_model = (
        fit["model"] == "GTR"
        or (fit["model"] == "MK-REGIME" and fit.get("regime_model") == "GTR")
        or (fit["model"] == "MK-MIXTURE" and fit.get("mixture_model") == "GTR")
    )
    return {
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
        "frequency_ratio_bounds": f"{FREQUENCY_RATIO_BOUNDS[0]},{FREQUENCY_RATIO_BOUNDS[1]}"
        if frequency_model
        else "",
        "num_frequencies_at_lower_bound": fit.get("num_frequencies_at_lower_bound", 0),
        "num_frequencies_at_upper_bound": fit.get("num_frequencies_at_upper_bound", 0),
    }


def _add_matrix_rates(row, matrix, state_ids, prefix=""):
    for from_index, from_state_id in enumerate(state_ids):
        for to_index, to_state_id in enumerate(state_ids):
            if from_index != to_index:
                row[f"rate_{prefix}{from_state_id}_to_{to_state_id}"] = float(
                    matrix[from_index, to_index]
                )


def _add_mixture_metadata(row, fit):
    row["mixture_model"] = fit["mixture_model"]
    row["rate_mixture"] = fit["rate_mixture"]
    row["rate_categories"] = fit["rate_categories"]
    row["mixture_category_rates"] = ",".join(
        str(float(value)) for value in fit["mixture_category_rates"]
    )
    row["mixture_weights"] = ",".join(
        str(float(value)) for value in fit["mixture_weights"]
    )
    row["gamma_shape"] = "" if fit["gamma_shape"] is None else fit["gamma_shape"]
    row["base_rate_estimated"] = fit["base_rate_estimated"]


def _add_regime_matrix_metadata(row, fit, state_ids):
    for regime, matrix in fit["rate_matrices_by_regime"].items():
        regime_id = _safe_column_state(regime)
        if fit.get("regime_model") in {"F81", "GTR"}:
            equilibrium = stationary_distribution(matrix)
            for state_id, probability in zip(state_ids, equilibrium, strict=True):
                row[f"equilibrium_{regime_id}_{state_id}"] = float(probability)
        _add_matrix_rates(row, matrix, state_ids, prefix=f"{regime_id}_")


def _add_hidden_model_metadata(row, fit):
    expanded_ids = [_safe_column_state(state) for state in fit["expanded_states"]]
    row["num_expanded_states"] = len(fit["expanded_states"])
    row["expanded_states"] = ",".join(fit["expanded_states"])
    if fit["model"] == "COVARION":
        row["base_rate"] = fit["base_rate"]
        row["base_rate_estimated"] = fit["base_rate_estimated"]
        row["hidden_rate_multipliers"] = ",".join(
            str(float(value)) for value in fit["hidden_rate_multipliers"]
        )
        row["hidden_rate_effective_rates"] = ",".join(
            str(float(value)) for value in fit["hidden_rate_effective_rates"]
        )
        row["switching_rate"] = fit["switching_rate"]
    _add_matrix_rates(row, fit["rate_matrix"], expanded_ids)


def _augment_discrete_model_row(row, fit, state_ids):
    if fit["model"] == "ER" and len(fit["rates"]) == 1:
        row["rate"] = float(fit["rates"][0])
    if fit["model"] in {"F81", "GTR"} or (
        fit["model"] == "MK-MIXTURE" and fit.get("mixture_model") in {"F81", "GTR"}
    ):
        equilibrium = stationary_distribution(fit["rate_matrix"])
        for state_id, probability in zip(state_ids, equilibrium, strict=True):
            row[f"equilibrium_{state_id}"] = float(probability)
    if fit["model"] == "MK-MIXTURE":
        _add_mixture_metadata(row, fit)
    if fit["model"] == "MK-REGIME":
        _add_regime_matrix_metadata(row, fit, state_ids)
    if fit["model"] in {"HRM", "COVARION"}:
        _add_hidden_model_metadata(row, fit)
    elif fit["model"] != "MK-REGIME":
        _add_matrix_rates(row, fit["rate_matrix"], state_ids)


def _build_model_table(states, root_prior_mode, fit):
    if fit["model"] == "THRESHOLD":
        return pd.DataFrame([_threshold_model_row(states, root_prior_mode, fit)])
    state_ids = [_safe_column_state(state) for state in states]
    row = _base_discrete_model_row(states, root_prior_mode, fit)
    _augment_discrete_model_row(row, fit, state_ids)
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


def _uniformization_parameters(rate_matrix, branch_length):
    branch_length = float(branch_length)
    num_states = rate_matrix.shape[0]
    if branch_length == 0.0 or not np.any(rate_matrix):
        return 0.0, 0.0, 0
    omega = float(np.max(-np.diag(rate_matrix)))
    if omega <= 0.0:
        return 0.0, 0.0, 0
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
    return omega, lam, max_n


def _build_uniformization_context(rate_matrix, branch_length):
    num_states = rate_matrix.shape[0]
    omega, lam, max_n = _uniformization_parameters(rate_matrix, branch_length)
    if max_n == 0:
        return {"no_events": True}
    backward_bytes = (max_n + 1) * num_states * np.dtype(float).itemsize
    if backward_bytes > _MAX_BACKWARD_BYTES:
        raise ValueError(
            "Stochastic mapping would require more than 256 MiB for one "
            "endpoint-specific bridge history; reduce the state space or "
            "rate/time scale."
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
        "max_n": max_n,
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
    context = _resolve_uniformization_context(
        rate_matrix,
        float(branch_length),
        uniformization_contexts,
        uniformization_context,
    )
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


def _resolve_uniformization_context(
    rate_matrix, branch_length, contexts, explicit_context
):
    if explicit_context is not None:
        return explicit_context
    if contexts is not None:
        return contexts[branch_length]
    return _build_uniformization_context(rate_matrix, branch_length)


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
        contexts[node_index] = _stochastic_map_context(uniformization_contexts, node)
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


def _stochastic_map_context(uniformization_contexts, node):
    if node in uniformization_contexts:
        return uniformization_contexts[node]
    return uniformization_contexts[float(node.dist)]


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
        _merge_projected_branch_counts(
            simulation_counts,
            branch_counts,
            spec["branch_ids"][node_index],
            spec.get("state_projection"),
        )
    return simulation_counts


def _merge_projected_branch_counts(
    simulation_counts, branch_counts, branch_id, projection
):
    for (from_state, to_state), count in branch_counts.items():
        if projection is not None:
            from_state = int(projection[from_state])
            to_state = int(projection[to_state])
        if from_state != to_state:
            simulation_counts[(branch_id, from_state, to_state)] += count


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


def _validated_simulation_threads(num_simulations, threads):
    num_simulations = _integer_option(num_simulations, "'--n-sim'")
    if num_simulations <= 0:
        raise ValueError(
            "'--n-sim' must be positive when '--stochastic-map-out' is specified."
        )
    threads = _integer_option(threads, "'--threads'")
    if threads <= 0:
        raise ValueError("'--threads' must be positive.")
    return num_simulations, threads


def _stochastic_uniformization_contexts(fit, branches):
    if "rate_matrix_by_node" in fit:
        uniformization_contexts = {}
        cache = {}
        for node in branches:
            matrix = fit["rate_matrix_by_node"][node]
            key = (id(matrix), float(node.dist))
            if key not in cache:
                cache[key] = _build_uniformization_context(matrix, node.dist)
            uniformization_contexts[node] = cache[key]
        return uniformization_contexts
    return _build_uniformization_contexts(
        fit["rate_matrix"], [float(node.dist) for node in branches]
    )


def _validate_stochastic_map_work(fit, branches, num_simulations):
    """Reject uniformization jobs whose bounded loop count is impractical."""

    by_node = fit.get("rate_matrix_by_node", {})
    cached_terms = {}
    per_simulation = 0
    preparation = 0
    for node in branches:
        matrix = by_node.get(node, fit["rate_matrix"])
        key = (id(matrix), float(node.dist))
        if key not in cached_terms:
            max_n = _uniformization_parameters(matrix, node.dist)[2]
            cached_terms[key] = max_n
            preparation += max_n * min(matrix.shape[0], num_simulations)
        per_simulation += cached_terms[key]
    work = preparation + num_simulations * per_simulation
    if work > _MAX_STOCHASTIC_MAP_WORK:
        raise ValueError(
            "Stochastic mapping would require more than "
            f"{_MAX_STOCHASTIC_MAP_WORK:,} bounded uniformization steps; reduce "
            "--n-sim, the rate/time scale, or fitted rate bounds."
        )


def _run_stochastic_simulations(spec, seed_sequences, threads):
    num_simulations = len(seed_sequences)
    if threads == 1 or num_simulations == 1:
        return _merge_simulation_counts(
            _simulate_stochastic_map_once(spec, seed_sequence)
            for seed_sequence in seed_sequences
        )
    max_workers = min(threads, num_simulations)
    chunks = [
        seed_sequences[worker_index::max_workers] for worker_index in range(max_workers)
    ]
    context = _get_process_pool_context()
    executor_kwargs = {"max_workers": max_workers}
    if context is not None:
        executor_kwargs["mp_context"] = context
    with ProcessPoolExecutor(**executor_kwargs) as executor:
        return _merge_count_summaries(
            executor.map(
                _simulate_stochastic_map_chunk_worker,
                ((spec, chunk) for chunk in chunks),
            )
        )


def _stochastic_map_rows(
    tree, states, node_to_branch_id, total_counts, any_counts, num_simulations
):
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
    return rows


def _simulate_stochastic_maps(tree, states, fit, num_simulations, seed=None, threads=1):
    num_simulations, threads = _validated_simulation_threads(num_simulations, threads)
    node_to_branch_id = assign_branch_ids(tree)
    branches = [node for node in tree.traverse() if not node.is_root]
    _validate_stochastic_map_work(fit, branches, num_simulations)
    contexts = _stochastic_uniformization_contexts(fit, branches)
    spec = _build_stochastic_map_spec(tree, fit, node_to_branch_id, contexts)
    seed_sequences = _simulation_seed_sequence(seed, num_simulations)
    total_counts, any_counts = _run_stochastic_simulations(
        spec, seed_sequences, threads
    )
    rows = _stochastic_map_rows(
        tree,
        states,
        node_to_branch_id,
        total_counts,
        any_counts,
        num_simulations,
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
        "--liability-out": getattr(args, "liability_out", None),
        "--posterior-samples-out": getattr(args, "posterior_samples_out", None),
        "--posterior-predictive-out": getattr(args, "posterior_predictive_out", None),
        "--bootstrap-out": getattr(args, "bootstrap_out", None),
        "--model-comparison-out": getattr(args, "model_comparison_out", None),
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
            ("--liability-out", getattr(args, "liability_out", None)),
            (
                "--posterior-samples-out",
                getattr(args, "posterior_samples_out", None),
            ),
            (
                "--posterior-predictive-out",
                getattr(args, "posterior_predictive_out", None),
            ),
            ("--bootstrap-out", getattr(args, "bootstrap_out", None)),
            (
                "--model-comparison-out",
                getattr(args, "model_comparison_out", None),
            ),
        ]
    )
    compare_models = getattr(args, "compare_models", None)
    comparison_out = getattr(args, "model_comparison_out", None)
    if compare_models not in (None, "") and comparison_out in (None, ""):
        raise ValueError("--compare-models requires --model-comparison-out.")
    if comparison_out not in (None, "") and compare_models in (None, ""):
        raise ValueError("--model-comparison-out requires --compare-models.")


def asr_main(args):
    _validate_asr_output_paths(args)
    tree = read_tree(
        args.infile,
        args.format,
        args.quoted_node_names,
        rooted=getattr(args, "input_rooted", "auto"),
    )
    _validate_tree_for_asr(tree)
    model = getattr(args, "model", None)
    trait_columns = asr_trait_columns(args.state_column, model)
    standard_error_columns = asr_standard_error_columns(
        getattr(args, "standard_error_column", None), model, trait_columns
    )
    state_column_input = (
        trait_columns
        if model in {"MV-BM", "MV-OU", "MK-MIXTURE"}
        else args.state_column
    )
    trait_df = read_asr_table(
        args.trait,
        state_column_input,
        list(tree.leaf_names()),
        missing_values=getattr(args, "missing_values", None),
        unmatched=getattr(args, "unmatched", "warn"),
        standard_error_column=standard_error_columns,
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

    if settings.model == "MK-MIXTURE":
        return _run_discrete_mixture_asr(tree, trait_df, args, settings, targets)
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
    if settings.model == "THRESHOLD":
        from nwkit.threshold_asr import (
            compute_threshold_marginals,
            threshold_liability_table,
        )

        posterior_by_node, threshold_fit = compute_threshold_marginals(
            tree,
            states,
            observed_state_by_leaf,
            likelihood_by_leaf,
            thresholds=getattr(args, "thresholds", None),
            num_samples=getattr(args, "liability_samples", 1000),
            burnin=getattr(args, "liability_burnin", 500),
            thin=getattr(args, "liability_thin", 1),
            chains=getattr(args, "liability_chains", 4),
            seed=getattr(args, "seed", None),
            _tree_validated=True,
        )
        fit = {
            "model": "THRESHOLD",
            "trait_type_requested": getattr(args, "trait_type", "auto"),
            "threshold_fit": threshold_fit,
        }
        table = _build_output_table(
            tree,
            states,
            observed_state_by_leaf,
            posterior_by_node,
            targets,
            settings.output,
        )
        _write_table(table, args.outfile)
        _write_model_table(states, settings.root_prior, fit, args.model_out)
        _write_annotated_tree(
            tree, states, posterior_by_node, observed_state_by_leaf, args
        )
        if getattr(args, "liability_out", None) not in (None, ""):
            _write_table(
                threshold_liability_table(tree, threshold_fit), args.liability_out
            )
        if threshold_fit.fit_status != "ok":
            sys.stderr.write(
                "Threshold ASR MCMC diagnostics: "
                f"status={threshold_fit.fit_status}, "
                f"max R-hat={threshold_fit.rhat_max:.4g}, "
                f"min ESS={threshold_fit.ess_min:.4g}.\n"
            )
        return
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
    elif settings.model == "COVARION":
        posterior_by_node, fit = compute_covarion_marginals(
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
    _write_discrete_model_comparison(
        tree,
        states,
        observed_state_by_leaf,
        likelihood_by_leaf,
        transition_graph,
        args,
        settings,
        fit,
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


def _comparison_model_names(value):
    names = tuple(item.strip().upper() for item in str(value).split(","))
    if len(names) < 2 or any(not item for item in names):
        raise ValueError("--compare-models requires at least two model names.")
    if len(names) != len(set(names)):
        raise ValueError("--compare-models contains duplicated model names.")
    return names


def _write_discrete_model_comparison(
    tree,
    states,
    observed_state_by_leaf,
    likelihood_by_leaf,
    transition_graph,
    args,
    settings,
    primary_fit,
):
    value = getattr(args, "compare_models", None)
    if value in (None, ""):
        return
    from nwkit.asr_comparison import model_comparison_table, summarize_fit

    supported = {"ER", "SYM", "ARD", "F81", "GTR", "COVARION"}
    models = _comparison_model_names(value)
    unsupported = [model for model in models if model not in supported]
    if unsupported:
        raise ValueError(
            "Discrete ASR comparison currently supports ER, SYM, ARD, F81, "
            "GTR, and COVARION; unsupported: " + ", ".join(unsupported)
        )
    summaries = []
    representatives: dict[tuple[Any, ...], tuple[str, dict[str, Any]]] = {}
    has_equivalent_alias = False
    for model in models:
        equivalence_key = None
        if model != "COVARION" and not (
            model == "ER" and getattr(args, "rate", None) is not None
        ):
            family = model_equivalence_family(model, states, transition_graph)
            equivalence_key = (
                family,
                settings.root_prior,
                transition_graph.shape,
                transition_graph.tobytes(),
            )
        if equivalence_key is not None and equivalence_key in representatives:
            representative_model, fit = representatives[equivalence_key]
            summary = summarize_fit(model, fit, trait_type="discrete")
            summary.update(
                {
                    "status": "equivalent",
                    "rankable": "no",
                    "equivalent_to": representative_model,
                    "message": (
                        f"Statistically equivalent to {representative_model} for "
                        "the same binary transition/root contract; excluded from "
                        "duplicate IC weighting."
                    ),
                }
            )
            summaries.append(summary)
            has_equivalent_alias = True
            continue
        if model == settings.model:
            fit = primary_fit
        elif model == "COVARION":
            hidden_categories = getattr(args, "hidden_categories", None)
            _, fit = compute_covarion_marginals(
                tree,
                states,
                observed_state_by_leaf,
                likelihood_by_leaf,
                hidden_categories=(
                    2 if hidden_categories is None else hidden_categories
                ),
                rate=getattr(args, "rate", None),
                root_prior_mode=settings.root_prior,
                rate_bounds=_parse_rate_bounds(getattr(args, "rate_bounds", None)),
                transition_graph=transition_graph,
            )
        else:
            _, fit = compute_mk_marginals(
                tree,
                states,
                observed_state_by_leaf,
                likelihood_by_leaf,
                model=model,
                rate=getattr(args, "rate", None),
                root_prior_mode=settings.root_prior,
                rate_bounds=_parse_rate_bounds(getattr(args, "rate_bounds", None)),
                transition_graph=transition_graph,
            )
        if equivalence_key is not None:
            representatives[equivalence_key] = (model, fit)
        summary = summarize_fit(model, fit, trait_type="discrete")
        summary.update(
            {"status": "ok", "rankable": "yes", "equivalent_to": "", "message": ""}
        )
        summaries.append(summary)
    if not has_equivalent_alias:
        for summary in summaries:
            for column in ("status", "rankable", "equivalent_to", "message"):
                summary.pop(column, None)
    _write_table(
        model_comparison_table(summaries),
        args.model_comparison_out,
    )


def _write_mixture_tree(
    tree,
    states,
    observed_by_character,
    posterior_by_character,
    args,
):
    tree_out = getattr(args, "tree_out", None)
    if tree_out in (None, ""):
        return
    props = ["asr_model"]
    for character, posterior_by_node in posterior_by_character.items():
        character_id = str(character).encode("utf-8").hex()
        state_property = f"asr_state_{character_id}"
        probability_property = f"asr_probability_{character_id}"
        props.extend((state_property, probability_property))
        for node in tree.traverse():
            posterior = posterior_by_node[node]
            map_index = int(np.argmax(posterior))
            node.add_props(
                asr_model="MK-MIXTURE",
                **{
                    state_property: states[map_index],
                    probability_property: float(posterior[map_index]),
                },
            )
            if args.tree_annotation == "all":
                observed_property = f"asr_observed_state_{character_id}"
                node.add_props(
                    **{
                        observed_property: ""
                        if not node.is_leaf
                        or observed_by_character[character].get(str(node.name)) is None
                        else observed_by_character[character][str(node.name)]
                    }
                )
                if observed_property not in props:
                    props.append(observed_property)
                for state_index, state in enumerate(states):
                    state_id = _safe_column_state(state)
                    probability_prop = f"asr_p_{character_id}_{state_id}"
                    node.add_props(**{probability_prop: float(posterior[state_index])})
                    if probability_prop not in props:
                        props.append(probability_prop)
    if args.tree_annotation == "state":
        props = [prop for prop in props if "state_" in prop or prop == "asr_model"]
    elif args.tree_annotation == "probability":
        props = [
            prop for prop in props if "probability_" in prop or prop == "asr_model"
        ]
    output_args = copy(args)
    output_args.outfile = tree_out
    write_tree(
        tree,
        output_args,
        format=getattr(args, "tree_outformat", "auto"),
        quiet=True,
        props=props,
    )


def _run_discrete_mixture_asr(tree, trait_df, args, settings, targets):
    leaf_names = list(tree.leaf_names())
    characters = asr_trait_columns(args.state_column, settings.model)
    states_arg = getattr(args, "states", None)
    if states_arg in (None, ""):
        states = []
        seen = set()
        for character in characters:
            character_states, _, _ = _read_tip_states(
                args.trait,
                character,
                leaf_names,
                missing_values_arg=getattr(args, "missing_values", None),
                ambiguous_separator=getattr(
                    args, "ambiguous_separator", DEFAULT_AMBIGUOUS_SEPARATOR
                ),
                unmatched=getattr(args, "unmatched", "warn"),
                trait_df=trait_df,
            )
            for state in character_states:
                if state not in seen:
                    states.append(state)
                    seen.add(state)
    else:
        states = _parse_state_argument(states_arg)
    observed_by_character = {}
    likelihoods_by_character = {}
    for character in characters:
        _, observed, likelihoods = _read_tip_states(
            args.trait,
            character,
            leaf_names,
            states_arg=states,
            missing_values_arg=getattr(args, "missing_values", None),
            ambiguous_separator=getattr(
                args, "ambiguous_separator", DEFAULT_AMBIGUOUS_SEPARATOR
            ),
            unmatched=getattr(args, "unmatched", "warn"),
            trait_df=trait_df,
            state_source="--states"
            if states_arg not in (None, "")
            else "joint characters",
        )
        observed_by_character[character] = observed
        likelihoods_by_character[character] = likelihoods
    transition_graph, transition_graph_source = read_transition_graph(
        getattr(args, "transition_graph", None),
        states,
        state_source="--states" if states_arg not in (None, "") else "joint characters",
    )
    posterior_by_character, fit = compute_mk_mixture_marginals(
        tree,
        states,
        observed_by_character,
        likelihoods_by_character,
        base_model=getattr(args, "mixture_model", None) or "ER",
        mixture=getattr(args, "rate_mixture", None) or "gamma",
        categories=(
            4
            if getattr(args, "rate_categories", None) is None
            else args.rate_categories
        ),
        rate=getattr(args, "rate", None),
        root_prior_mode=settings.root_prior,
        rate_bounds=_parse_rate_bounds(getattr(args, "rate_bounds", None)),
        transition_graph=transition_graph,
    )
    fit["transition_graph"] = transition_graph_source
    fit["rate_matrix_source"] = "estimated"
    fit["trait_type_requested"] = getattr(args, "trait_type", "auto")
    tables = []
    for character in characters:
        table = _build_output_table(
            tree,
            states,
            observed_by_character[character],
            posterior_by_character[character],
            targets,
            settings.output,
        )
        table.insert(4, "trait", character)
        tables.append(table)
    _write_table(pd.concat(tables, ignore_index=True), args.outfile)
    _write_model_table(
        states, settings.root_prior, fit, getattr(args, "model_out", None)
    )
    _write_mixture_tree(
        tree,
        states,
        observed_by_character,
        posterior_by_character,
        args,
    )


def _continuous_observations(tree, trait_df, args, settings):
    if settings.model in {"MV-BM", "MV-OU"}:
        trait_columns = asr_trait_columns(args.state_column, settings.model)
        observed = continuous_tip_vectors(
            trait_df, trait_columns, list(tree.leaf_names())
        )
        error_columns = asr_standard_error_columns(
            getattr(args, "standard_error_column", None),
            settings.model,
            trait_columns,
        )
        errors = continuous_tip_vector_errors(
            trait_df,
            error_columns,
            trait_columns,
            list(tree.leaf_names()),
        )
        return trait_columns, observed, errors
    observed, errors = continuous_tip_values(
        trait_df,
        args.state_column,
        list(tree.leaf_names()),
        getattr(args, "standard_error_column", None),
    )
    return (args.state_column,), observed, errors


def _continuous_regime_assignment(tree, args, model):
    if model not in {"BMS", "BMS-DRIFT", "OUM", "OUMA", "OUMV", "OUMVA"}:
        return None
    from nwkit.asr_regimes import read_regime_map

    return read_regime_map(getattr(args, "regime_map", None), tree)


def _fit_continuous_model(
    tree, observed, errors, trait_columns, args, settings, regime_assignment
):
    from nwkit.asr_regimes import read_regime_parameters
    from nwkit.continuous_asr import compute_bm_marginals

    model = settings.model
    if model == "BM":
        return compute_bm_marginals(
            tree,
            observed,
            sigma2=getattr(args, "sigma2", None),
            standard_errors=errors,
            _tree_validated=True,
        )
    if model == "BMS":
        from nwkit.regime_continuous_asr import compute_bms_marginals

        assert regime_assignment is not None
        fixed = read_regime_parameters(
            getattr(args, "regime_parameters", None),
            regime_assignment.regimes,
            ("sigma2",),
        )
        values = None if fixed is None else fixed[0]
        return compute_bms_marginals(
            tree,
            observed,
            regime_assignment,
            sigma2=getattr(args, "sigma2", None),
            sigma2_by_regime=None
            if values is None
            else {regime: item["sigma2"] for regime, item in values.items()},
            regime_parameters_source=None if fixed is None else fixed[1],
            standard_errors=errors,
            _tree_validated=True,
        )
    if model == "BMS-DRIFT":
        from nwkit.regime_drift_asr import compute_regime_drift_marginals

        assert regime_assignment is not None
        fixed = read_regime_parameters(
            getattr(args, "regime_parameters", None),
            regime_assignment.regimes,
            ("sigma2", "drift"),
        )
        return compute_regime_drift_marginals(
            tree,
            observed,
            regime_assignment,
            sigma2=getattr(args, "sigma2", None),
            drift=getattr(args, "drift", None),
            regime_parameters=None if fixed is None else fixed[0],
            regime_parameters_source=None if fixed is None else fixed[1],
            standard_errors=errors,
            _tree_validated=True,
        )
    if model == "OU":
        return _fit_ou_model(tree, observed, errors, args, settings)
    if model in {"OUM", "OUMA", "OUMV", "OUMVA"}:
        return _fit_regime_ou_model(
            tree, observed, errors, args, model, regime_assignment
        )
    if model in {"LAMBDA", "KAPPA", "DELTA", "EB", "ACDC"}:
        return _fit_transformed_continuous_model(tree, observed, errors, args, model)
    if model == "BM-DRIFT":
        from nwkit.nonstationary_continuous_asr import compute_bm_drift_marginals

        return compute_bm_drift_marginals(
            tree,
            observed,
            sigma2=getattr(args, "sigma2", None),
            drift=getattr(args, "drift", None),
            standard_errors=errors,
            _tree_validated=True,
        )
    if model == "MV-BM":
        from nwkit.multivariate_asr import compute_mvbm_marginals

        return compute_mvbm_marginals(
            tree,
            observed,
            trait_columns,
            standard_errors=errors,
            _tree_validated=True,
        )
    if model == "MV-OU":
        from nwkit.multivariate_gaussian_asr import fit_dense_mvou
        from nwkit.ou_asr import parse_alpha_bounds

        return fit_dense_mvou(
            tree,
            observed,
            trait_columns,
            alpha=getattr(args, "alpha", None),
            alpha_bounds=parse_alpha_bounds(getattr(args, "alpha_bounds", None), tree),
            standard_errors=errors,
        )
    raise ValueError(f"Unsupported continuous ASR model: {model}.")


def _fit_ou_model(tree, observed, errors, args, settings):
    from nwkit.ou_asr import compute_ou_marginals, parse_alpha_bounds

    alpha_bounds = parse_alpha_bounds(getattr(args, "alpha_bounds", None), tree)
    if settings.root_prior == "stationary":
        return compute_ou_marginals(
            tree,
            observed,
            alpha=getattr(args, "alpha", None),
            sigma2=getattr(args, "sigma2", None),
            theta=getattr(args, "theta", None),
            alpha_bounds=alpha_bounds,
            standard_errors=errors,
            _tree_validated=True,
        )
    from nwkit.general_ou_asr import compute_general_ou_marginals

    return compute_general_ou_marginals(
        tree,
        observed,
        root_prior=settings.root_prior,
        root_mean=getattr(args, "root_mean", None),
        root_variance=getattr(args, "root_variance", None),
        alpha=getattr(args, "alpha", None),
        sigma2=getattr(args, "sigma2", None),
        theta=getattr(args, "theta", None),
        alpha_bounds=alpha_bounds,
        standard_errors=errors,
    )


def _fit_regime_ou_model(tree, observed, errors, args, model, regime_assignment):
    from nwkit.asr_regimes import read_regime_parameters
    from nwkit.ou_asr import parse_alpha_bounds
    from nwkit.regime_gaussian_asr import (
        compute_regime_ou_marginals,
        regime_parameter_columns,
    )

    assert regime_assignment is not None
    fixed = read_regime_parameters(
        getattr(args, "regime_parameters", None),
        regime_assignment.regimes,
        regime_parameter_columns(model),
    )
    return compute_regime_ou_marginals(
        tree,
        observed,
        regime_assignment,
        model=model,
        alpha=getattr(args, "alpha", None),
        sigma2=getattr(args, "sigma2", None),
        theta=getattr(args, "theta", None),
        regime_parameters=None if fixed is None else fixed[0],
        regime_parameters_source=None if fixed is None else fixed[1],
        alpha_bounds=parse_alpha_bounds(getattr(args, "alpha_bounds", None), tree),
        standard_errors=errors,
        _tree_validated=True,
    )


def _fit_transformed_continuous_model(tree, observed, errors, args, model):
    from nwkit.transformed_continuous_asr import (
        compute_transformed_bm_marginals,
        parse_parameter_bounds,
    )

    parameter = getattr(args, "evolution_parameter", None)
    parameter_bounds = getattr(args, "evolution_parameter_bounds", None)
    if model in {"EB", "ACDC"}:
        parameter = (
            parameter if parameter is not None else getattr(args, "eb_rate", None)
        )
        parameter_bounds = (
            parameter_bounds
            if parameter_bounds is not None
            else getattr(args, "eb_rate_bounds", None)
        )
    process_model = model.lower()
    return compute_transformed_bm_marginals(
        tree,
        observed,
        model=process_model,
        sigma2=getattr(args, "sigma2", None),
        evolution_parameter=parameter,
        evolution_parameter_bounds=parse_parameter_bounds(
            parameter_bounds, tree, process_model
        ),
        profile_ci_level=getattr(args, "profile_ci_level", None),
        standard_errors=errors,
        _tree_validated=True,
    )


def _multivariate_fit_status_message(model, fit):
    likelihood_available = any(
        value is not None and math.isfinite(float(value))
        for value in (
            getattr(fit, "log_likelihood", None),
            getattr(fit, "restricted_log_likelihood", None),
        )
    )
    if likelihood_available:
        return (
            f"Continuous ASR: {model} fit status={fit.fit_status}; intervals "
            "condition on fitted parameters and exclude parameter-estimation "
            "uncertainty.\n"
        )
    return (
        f"Continuous ASR: {model} fit status={fit.fit_status}; marginal "
        "reconstruction is available, but the ordinary multivariate likelihood "
        "is not.\n"
    )


def _run_continuous_asr(tree, trait_df, args, settings, targets):
    from nwkit.continuous_asr_io import (
        continuous_model_table,
        continuous_output_table,
        write_continuous_tree,
    )

    trait_columns, observed, errors = _continuous_observations(
        tree, trait_df, args, settings
    )
    regime_assignment = _continuous_regime_assignment(tree, args, settings.model)
    posterior, fit = _fit_continuous_model(
        tree,
        observed,
        errors,
        trait_columns,
        args,
        settings,
        regime_assignment,
    )
    _write_continuous_model_comparison(tree, observed, errors, args, settings, fit)
    from nwkit.asr_continuous_diagnostics import write_continuous_diagnostics

    write_continuous_diagnostics(
        tree,
        observed,
        errors,
        args,
        settings,
        fit,
        regime_assignment,
    )
    selected = [
        node for node in tree.traverse() if _should_output_node(node, observed, targets)
    ]
    if settings.model in {"MV-BM", "MV-OU"}:
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
            errors=errors,
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
    elif settings.model in {"MV-BM", "MV-OU"} and fit.fit_status != "ok":
        sys.stderr.write(_multivariate_fit_status_message(settings.model, fit))
    elif settings.model in {
        "BMS",
        "BMS-DRIFT",
        "LAMBDA",
        "KAPPA",
        "DELTA",
        "OU",
        "OUM",
        "OUMA",
        "OUMV",
        "OUMVA",
        "EB",
        "ACDC",
        "BM-DRIFT",
    } and (fit.fit_status != "ok" or not fit.optimizer_success):
        sys.stderr.write(
            f"Continuous ASR: {settings.model} fit status={fit.fit_status}; intervals condition "
            "on fitted parameters and exclude parameter-estimation uncertainty.\n"
        )
    _write_table(table, args.outfile)
    if settings.model in {"MV-BM", "MV-OU"} and getattr(
        args, "covariance_out", None
    ) not in (None, ""):
        _write_table(
            multivariate_covariance_table(tree, selected, posterior, trait_columns),
            args.covariance_out,
        )
    if getattr(args, "model_out", None) not in (None, ""):
        model_table = (
            multivariate_model_table(fit, args, settings.ci_level)
            if settings.model in {"MV-BM", "MV-OU"}
            else continuous_model_table(fit, args, settings.ci_level)
        )
        _write_table(model_table, args.model_out)
    if settings.model in {"MV-BM", "MV-OU"}:
        write_multivariate_tree(
            tree,
            observed,
            posterior,
            trait_columns,
            args,
            settings.ci_level,
            errors=errors,
        )
    else:
        write_continuous_tree(
            tree, observed, errors, posterior, args, settings.ci_level
        )


def _write_continuous_model_comparison(
    tree, observed, errors, args, settings, primary_fit
):
    value = getattr(args, "compare_models", None)
    if value in (None, ""):
        return
    from nwkit.asr_comparison import model_comparison_table, summarize_fit
    from nwkit.continuous_asr import compute_bm_marginals
    from nwkit.nonstationary_continuous_asr import compute_bm_drift_marginals
    from nwkit.transformed_continuous_asr import compute_transformed_bm_marginals

    supported = {"BM", "LAMBDA", "KAPPA", "DELTA", "EB", "ACDC", "BM-DRIFT"}
    models = _comparison_model_names(value)
    unsupported = [model for model in models if model not in supported]
    if unsupported:
        raise ValueError(
            "Continuous flat-root ASR comparison currently supports BM, LAMBDA, "
            "KAPPA, DELTA, EB, ACDC, and BM-DRIFT; unsupported: "
            + ", ".join(unsupported)
        )
    if settings.root_prior != "flat":
        raise ValueError(
            "The requested ASR comparison uses flat-root integrated likelihoods; "
            "select a flat-root primary model."
        )
    summaries = []
    for model in models:
        if model == settings.model:
            candidate_fit = primary_fit
        elif model == "BM":
            _, candidate_fit = compute_bm_marginals(
                tree,
                observed,
                sigma2=getattr(args, "sigma2", None),
                standard_errors=errors,
                _tree_validated=True,
            )
        elif model == "BM-DRIFT":
            _, candidate_fit = compute_bm_drift_marginals(
                tree,
                observed,
                sigma2=getattr(args, "sigma2", None),
                standard_errors=errors,
                _tree_validated=True,
            )
        else:
            _, candidate_fit = compute_transformed_bm_marginals(
                tree,
                observed,
                model=model.lower(),
                sigma2=getattr(args, "sigma2", None),
                standard_errors=errors,
                _tree_validated=True,
            )
        summaries.append(summarize_fit(model, candidate_fit, trait_type="continuous"))
    _write_table(
        model_comparison_table(summaries),
        args.model_comparison_out,
    )

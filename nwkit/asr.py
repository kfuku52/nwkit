import math
from collections import defaultdict

import numpy as np
import pandas as pd
from scipy.linalg import expm
from scipy.optimize import minimize
from scipy.special import logsumexp
from scipy.stats import poisson

from nwkit.util import (
    TREE_FORMAT_PROP,
    assign_branch_ids,
    is_rooted,
    read_tree,
    read_tsv_preserving_leaf_name,
    validate_unique_named_leaves,
)


DEFAULT_MISSING_VALUES = ('', 'NA', 'NaN', 'nan', '?', 'missing', 'unknown')
DEFAULT_RATE_BOUNDS = (10 ** -9, 10 ** 3)
DEFAULT_AMBIGUOUS_SEPARATOR = '|'
SUPPORTED_MODELS = ('ER', 'SYM', 'ARD')


def _parse_comma_list(value, option_name):
    if value in ['', None]:
        return list()
    items = [item.strip() for item in str(value).split(',')]
    if any(item == '' for item in items):
        raise ValueError("'{}' contains an empty item.".format(option_name))
    return items


def _parse_missing_values(value):
    if value in ['', None]:
        return set(DEFAULT_MISSING_VALUES)
    return set(_parse_comma_list(value, '--missing_values'))


def _parse_states(value):
    states = _parse_comma_list(value, '--states')
    if len(states) != len(set(states)):
        raise ValueError("'--states' contains duplicated states.")
    return states


def _parse_rate_bounds(value):
    if value in ['', None]:
        return DEFAULT_RATE_BOUNDS
    items = _parse_comma_list(value, '--rate_bounds')
    if len(items) != 2:
        raise ValueError("'--rate_bounds' must contain exactly two comma-separated values.")
    try:
        lower = float(items[0])
        upper = float(items[1])
    except ValueError as exc:
        raise ValueError("'--rate_bounds' values must be numeric.") from exc
    if (not math.isfinite(lower)) or (not math.isfinite(upper)) or lower <= 0.0 or upper <= 0.0:
        raise ValueError("'--rate_bounds' values must be positive finite numbers.")
    if lower >= upper:
        raise ValueError("'--rate_bounds' lower bound must be smaller than the upper bound.")
    return lower, upper


def _parse_targets(value):
    if value in ['', None]:
        value = 'intnode,missing_tip'
    targets = set()
    aliases = {'leaf': 'tip', 'leaves': 'tip', 'tips': 'tip'}
    valid_targets = {'all', 'intnode', 'tip', 'missing_tip'}
    for raw_target in _parse_comma_list(value, '--target'):
        target = aliases.get(raw_target, raw_target)
        if target not in valid_targets:
            raise ValueError("Unknown '--target': {}".format(raw_target))
        targets.add(target)
    if 'all' in targets:
        return {'all'}
    return targets


def _is_missing_trait_value(value, missing_values):
    if pd.isna(value):
        return True
    text = str(value).strip()
    return text in missing_values


def _split_state_value(value, ambiguous_separator):
    state_text = str(value).strip()
    if ambiguous_separator in ['', None]:
        return [state_text]
    states = [state.strip() for state in state_text.split(ambiguous_separator)]
    if any(state == '' for state in states):
        raise ValueError("Ambiguous state value contains an empty state: {}".format(state_text))
    if len(states) != len(set(states)):
        raise ValueError("Ambiguous state value contains duplicated states: {}".format(state_text))
    return states


def _read_tip_states(
    trait_path,
    state_column,
    tree_leaf_names,
    states_arg=None,
    missing_values_arg=None,
    ambiguous_separator=DEFAULT_AMBIGUOUS_SEPARATOR,
):
    if trait_path in ['', None]:
        raise ValueError("'--trait' is required.")
    if state_column in ['', None]:
        raise ValueError("'--state_column' is required.")
    trait_df = read_tsv_preserving_leaf_name(trait_path)
    if 'leaf_name' not in trait_df.columns:
        raise ValueError("Column 'leaf_name' is required in '--trait'.")
    if state_column not in trait_df.columns:
        raise ValueError("Column '{}' specified by '--state_column' was not found in '--trait'.".format(state_column))

    trait_df = trait_df.copy()
    trait_df['leaf_name'] = [str(leaf_name) for leaf_name in trait_df['leaf_name'].tolist()]
    if any(leaf_name.strip() == '' for leaf_name in trait_df['leaf_name'].tolist()):
        raise ValueError("Column 'leaf_name' in '--trait' must not contain empty values.")
    duplicated_leaf_names = trait_df.loc[
        trait_df['leaf_name'].duplicated(keep=False), 'leaf_name'
    ].unique().tolist()
    if duplicated_leaf_names:
        duplicated_leaf_names = sorted(str(name) for name in duplicated_leaf_names)
        raise ValueError("Duplicated 'leaf_name' entries in '--trait': {}".format(', '.join(duplicated_leaf_names)))

    tree_leaf_name_set = set(tree_leaf_names)
    missing_leaf_names = sorted(set(trait_df['leaf_name']) - tree_leaf_name_set)
    if missing_leaf_names:
        raise ValueError(
            "The following 'leaf_name' values in '--trait' were not found in the input tree: {}".format(
                ', '.join(missing_leaf_names)
            )
        )

    missing_values = _parse_missing_values(missing_values_arg)
    raw_state_by_leaf = dict()
    observed_state_sets = list()
    for _, row in trait_df.iterrows():
        leaf_name = str(row['leaf_name'])
        raw_state = row[state_column]
        if _is_missing_trait_value(raw_state, missing_values):
            raw_state_by_leaf[leaf_name] = None
            continue
        state_set = _split_state_value(raw_state, ambiguous_separator)
        raw_state_by_leaf[leaf_name] = '|'.join(state_set)
        observed_state_sets.append(state_set)

    states = _parse_states(states_arg)
    observed_states = [state for state_set in observed_state_sets for state in state_set]
    if states:
        unknown_states = sorted(set(observed_states) - set(states))
        if unknown_states:
            raise ValueError(
                "State(s) in '--trait' were not listed in '--states': {}".format(
                    ', '.join(unknown_states)
                )
            )
    else:
        seen_states = set()
        for state in observed_states:
            if state not in seen_states:
                states.append(state)
                seen_states.add(state)
    if len(states) == 0:
        raise ValueError("At least one observed or explicitly listed state is required.")

    state_to_index = {state: index for index, state in enumerate(states)}
    observed_state_by_leaf = {
        leaf_name: raw_state_by_leaf.get(leaf_name)
        for leaf_name in tree_leaf_names
    }
    likelihood_by_leaf = dict()
    for leaf_name in tree_leaf_names:
        raw_state = observed_state_by_leaf[leaf_name]
        likelihood = np.ones(len(states), dtype=float)
        if raw_state is not None:
            likelihood = np.zeros(len(states), dtype=float)
            for state in _split_state_value(raw_state, DEFAULT_AMBIGUOUS_SEPARATOR):
                likelihood[state_to_index[state]] = 1.0
        likelihood_by_leaf[leaf_name] = likelihood
    return states, observed_state_by_leaf, likelihood_by_leaf


def _validate_tree_for_asr(tree):
    validate_unique_named_leaves(tree, option_name='--infile', context=" for 'asr'")
    if not is_rooted(tree):
        raise ValueError("ASR requires a rooted tree.")
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
    if root_prior == 'equal':
        return np.full(num_states, 1.0 / float(num_states), dtype=float)
    if root_prior == 'empirical':
        counts = np.zeros(num_states, dtype=float)
        for leaf_name, observed_state in observed_state_by_leaf.items():
            if observed_state is None:
                continue
            likelihood = likelihood_by_leaf[leaf_name]
            counts += likelihood / likelihood.sum()
        if counts.sum() == 0.0:
            return np.full(num_states, 1.0 / float(num_states), dtype=float)
        return counts / counts.sum()
    raise ValueError("Unsupported '--root_prior': {}".format(root_prior))


def _num_rate_parameters(model, num_states):
    if num_states <= 1:
        return 0
    if model == 'ER':
        return 1
    if model == 'SYM':
        return num_states * (num_states - 1) // 2
    if model == 'ARD':
        return num_states * (num_states - 1)
    raise ValueError("Unsupported '--model': {}".format(model))


def _rate_parameter_labels(model, states):
    if len(states) <= 1:
        return list()
    labels = list()
    if model == 'ER':
        return [('all', 'all')]
    if model == 'SYM':
        for i, from_state in enumerate(states):
            for to_state in states[i + 1:]:
                labels.append((from_state, to_state))
        return labels
    if model == 'ARD':
        for from_state in states:
            for to_state in states:
                if from_state != to_state:
                    labels.append((from_state, to_state))
        return labels
    raise ValueError("Unsupported '--model': {}".format(model))


def _build_rate_matrix(model, states, rates):
    num_states = len(states)
    matrix = np.zeros((num_states, num_states), dtype=float)
    if num_states <= 1:
        return matrix
    rates = np.asarray(rates, dtype=float)
    if len(rates) != _num_rate_parameters(model, num_states):
        raise ValueError("Unexpected number of rate parameters for model '{}'.".format(model))
    if np.any(~np.isfinite(rates)) or np.any(rates < 0.0):
        raise ValueError("Mk rates must be non-negative finite numbers.")
    rate_index = 0
    if model == 'ER':
        matrix[:] = rates[0]
        np.fill_diagonal(matrix, 0.0)
    elif model == 'SYM':
        for i in range(num_states):
            for j in range(i + 1, num_states):
                matrix[i, j] = rates[rate_index]
                matrix[j, i] = rates[rate_index]
                rate_index += 1
    elif model == 'ARD':
        for i in range(num_states):
            for j in range(num_states):
                if i == j:
                    continue
                matrix[i, j] = rates[rate_index]
                rate_index += 1
    else:
        raise ValueError("Unsupported '--model': {}".format(model))
    np.fill_diagonal(matrix, -matrix.sum(axis=1))
    return matrix


def _transition_matrix(rate_matrix, branch_length):
    if float(branch_length) == 0.0:
        return np.eye(rate_matrix.shape[0], dtype=float)
    matrix = expm(rate_matrix * float(branch_length))
    matrix[matrix < 0.0] = np.maximum(matrix[matrix < 0.0], -10 ** -12)
    matrix = np.maximum(matrix, 0.0)
    row_sums = matrix.sum(axis=1)
    for row_index, row_sum in enumerate(row_sums):
        if row_sum > 0.0:
            matrix[row_index, :] /= row_sum
    return matrix


def _compute_inside_likelihoods(tree, likelihood_by_leaf, rate_matrix):
    num_states = rate_matrix.shape[0]
    inside = dict()
    log_scales = dict()
    child_terms = dict()
    transition_matrices = dict()
    for node in tree.traverse(strategy='postorder'):
        if node.is_leaf:
            inside[node] = likelihood_by_leaf[node.name].astype(float)
            log_scales[node] = 0.0
            continue

        likelihood = np.ones(num_states, dtype=float)
        log_scale = 0.0
        child_terms[node] = dict()
        for child in node.get_children():
            matrix = _transition_matrix(rate_matrix, child.dist)
            transition_matrices[child] = matrix
            term = matrix.dot(inside[child])
            child_terms[node][child] = term
            likelihood *= term
            log_scale += log_scales[child]
        scale = float(np.max(likelihood))
        if scale <= 0.0:
            inside[node] = likelihood
            log_scales[node] = -math.inf
            continue
        inside[node] = likelihood / scale
        log_scales[node] = log_scale + math.log(scale)
    return inside, log_scales, child_terms, transition_matrices


def _log_likelihood(tree, likelihood_by_leaf, root_prior, rate_matrix):
    inside, log_scales, _, _ = _compute_inside_likelihoods(
        tree=tree,
        likelihood_by_leaf=likelihood_by_leaf,
        rate_matrix=rate_matrix,
    )
    root_term = float(np.dot(root_prior, inside[tree]))
    if root_term <= 0.0:
        return -math.inf
    return log_scales[tree] + math.log(root_term)


def _initial_rate_value(tree, rate, rate_bounds):
    if rate is not None:
        value = float(rate)
    else:
        branch_lengths = [
            float(node.dist)
            for node in tree.traverse()
            if (not node.is_root) and node.dist is not None and float(node.dist) > 0.0
        ]
        mean_branch_length = sum(branch_lengths) / len(branch_lengths) if branch_lengths else 1.0
        value = 1.0 / max(mean_branch_length, 1.0)
    lower, upper = rate_bounds
    return min(max(value, lower), upper)


def _fit_rate_matrix(tree, model, states, likelihood_by_leaf, root_prior, rate=None, rate_bounds=None):
    if model not in SUPPORTED_MODELS:
        raise ValueError("Unsupported '--model': {}".format(model))
    rate_bounds = DEFAULT_RATE_BOUNDS if rate_bounds is None else rate_bounds
    num_params = _num_rate_parameters(model, len(states))
    if num_params == 0:
        rate_matrix = _build_rate_matrix(model, states, [])
        log_likelihood = _log_likelihood(tree, likelihood_by_leaf, root_prior, rate_matrix)
        return {
            'model': model,
            'rates': np.array([], dtype=float),
            'rate_matrix': rate_matrix,
            'log_likelihood': log_likelihood,
            'rate_estimated': False,
            'rate_bounds': rate_bounds,
        }

    if model == 'ER' and rate is not None:
        fixed_rate = float(rate)
        if (not math.isfinite(fixed_rate)) or fixed_rate < 0.0:
            raise ValueError("'--rate' must be a non-negative finite number.")
        rate_matrix = _build_rate_matrix(model, states, [fixed_rate])
        log_likelihood = _log_likelihood(tree, likelihood_by_leaf, root_prior, rate_matrix)
        return {
            'model': model,
            'rates': np.array([fixed_rate], dtype=float),
            'rate_matrix': rate_matrix,
            'log_likelihood': log_likelihood,
            'rate_estimated': False,
            'rate_bounds': rate_bounds,
        }

    initial_rate = _initial_rate_value(tree, rate, rate_bounds)
    lower_log = math.log(rate_bounds[0])
    upper_log = math.log(rate_bounds[1])
    initial_log_rates = np.full(num_params, math.log(initial_rate), dtype=float)

    def objective(log_rates):
        rates = np.exp(log_rates)
        rate_matrix = _build_rate_matrix(model, states, rates)
        log_likelihood = _log_likelihood(tree, likelihood_by_leaf, root_prior, rate_matrix)
        if not math.isfinite(log_likelihood):
            return 10 ** 100
        return -log_likelihood

    result = minimize(
        objective,
        initial_log_rates,
        method='L-BFGS-B',
        bounds=[(lower_log, upper_log)] * num_params,
    )
    if (not result.success) and (not math.isfinite(result.fun)):
        raise ValueError("Failed to estimate Mk model parameters: {}".format(result.message))
    rates = np.exp(result.x)
    rate_matrix = _build_rate_matrix(model, states, rates)
    log_likelihood = _log_likelihood(tree, likelihood_by_leaf, root_prior, rate_matrix)
    return {
        'model': model,
        'rates': rates,
        'rate_matrix': rate_matrix,
        'log_likelihood': log_likelihood,
        'rate_estimated': True,
        'rate_bounds': rate_bounds,
        'optimizer_success': bool(result.success),
        'optimizer_message': str(result.message),
    }


def _compute_outside_likelihoods(tree, child_terms, transition_matrices, root_prior):
    outside = {tree: root_prior.copy()}
    for node in tree.traverse(strategy='preorder'):
        children = list(node.get_children())
        if len(children) == 0:
            continue
        sibling_product_by_child = dict()
        for child in children:
            sibling_product = np.ones_like(outside[node], dtype=float)
            for sibling in children:
                if sibling is child:
                    continue
                sibling_product *= child_terms[node][sibling]
            sibling_product_by_child[child] = sibling_product
        for child in children:
            parent_weight = outside[node] * sibling_product_by_child[child]
            matrix = transition_matrices[child]
            outside[child] = matrix.T.dot(parent_weight)
            total = float(outside[child].sum())
            if total > 0.0:
                outside[child] = outside[child] / total
    return outside


def compute_mk_marginals(
    tree,
    states,
    observed_state_by_leaf,
    likelihood_by_leaf,
    model='ER',
    rate=None,
    root_prior_mode='equal',
    rate_bounds=None,
):
    root_prior = _get_root_prior(root_prior_mode, states, observed_state_by_leaf, likelihood_by_leaf)
    fit = _fit_rate_matrix(
        tree=tree,
        model=model,
        states=states,
        likelihood_by_leaf=likelihood_by_leaf,
        root_prior=root_prior,
        rate=rate,
        rate_bounds=DEFAULT_RATE_BOUNDS if rate_bounds is None else rate_bounds,
    )
    if not math.isfinite(fit['log_likelihood']):
        raise ValueError("The observed tip states have zero likelihood under the Mk model.")

    inside, _, child_terms, transition_matrices = _compute_inside_likelihoods(
        tree=tree,
        likelihood_by_leaf=likelihood_by_leaf,
        rate_matrix=fit['rate_matrix'],
    )
    outside = _compute_outside_likelihoods(
        tree=tree,
        child_terms=child_terms,
        transition_matrices=transition_matrices,
        root_prior=root_prior,
    )
    posterior_by_node = dict()
    for node in tree.traverse():
        posterior = inside[node] * outside[node]
        total = float(posterior.sum())
        if total <= 0.0:
            raise ValueError("Failed to calculate posterior state probabilities.")
        posterior_by_node[node] = posterior / total
    fit.update({
        'root_prior': root_prior,
        'inside': inside,
        'transition_matrices': transition_matrices,
        'posterior_by_node': posterior_by_node,
    })
    return posterior_by_node, fit


def _node_type(node):
    if node.is_root:
        return 'root'
    if node.is_leaf:
        return 'tip'
    return 'intnode'


def _is_missing_tip(node, observed_state_by_leaf):
    return node.is_leaf and observed_state_by_leaf.get(node.name) is None


def _should_output_node(node, observed_state_by_leaf, targets):
    if 'all' in targets:
        return True
    if (not node.is_leaf) and ('intnode' in targets):
        return True
    if node.is_leaf and ('tip' in targets):
        return True
    if _is_missing_tip(node, observed_state_by_leaf) and ('missing_tip' in targets):
        return True
    return False


def _build_output_table(tree, states, observed_state_by_leaf, posterior_by_node, targets, output_mode):
    node_to_branch_id = assign_branch_ids(tree)
    rows = list()
    for node in tree.traverse():
        if not _should_output_node(node, observed_state_by_leaf, targets):
            continue
        posterior = posterior_by_node[node]
        map_index = int(np.argmax(posterior))
        observed_state = observed_state_by_leaf.get(node.name) if node.is_leaf else None
        row = {
            'branch_id': node_to_branch_id[node],
            'parent': -1 if node.is_root else node_to_branch_id[node.up],
            'node_type': _node_type(node),
            'name': '' if node.name in [None, ''] else str(node.name),
            'observed_state': '' if observed_state is None else observed_state,
            'is_imputed': bool(_is_missing_tip(node, observed_state_by_leaf)),
        }
        if output_mode == 'probabilities':
            row['map_state'] = states[map_index]
            row['map_probability'] = float(posterior[map_index])
            for state, probability in zip(states, posterior):
                row['p_{}'.format(state)] = float(probability)
        elif output_mode == 'map':
            row['state'] = states[map_index]
            row['probability'] = float(posterior[map_index])
        else:
            raise ValueError("Unsupported '--output': {}".format(output_mode))
        rows.append(row)
    if output_mode == 'probabilities':
        columns = [
            'branch_id',
            'parent',
            'node_type',
            'name',
            'observed_state',
            'is_imputed',
            'map_state',
            'map_probability',
        ] + ['p_{}'.format(state) for state in states]
    else:
        columns = [
            'branch_id',
            'parent',
            'node_type',
            'name',
            'observed_state',
            'is_imputed',
            'state',
            'probability',
        ]
    return pd.DataFrame(rows, columns=columns)


def _safe_column_state(state):
    return ''.join(char if char.isalnum() or char == '_' else '_' for char in str(state))


def _build_model_table(states, root_prior_mode, fit):
    row = {
        'model': fit['model'],
        'num_states': len(states),
        'states': ','.join(states),
        'root_prior': root_prior_mode,
        'root_prior_values': ','.join(str(float(value)) for value in fit['root_prior']),
        'rate_estimated': bool(fit['rate_estimated']),
        'log_likelihood': float(fit['log_likelihood']),
        'rate_bounds': '{},{}'.format(fit['rate_bounds'][0], fit['rate_bounds'][1]),
    }
    if fit['model'] == 'ER' and len(fit['rates']) == 1:
        row['rate'] = float(fit['rates'][0])
    for from_index, from_state in enumerate(states):
        for to_index, to_state in enumerate(states):
            if from_index == to_index:
                continue
            row['rate_{}_to_{}'.format(_safe_column_state(from_state), _safe_column_state(to_state))] = float(
                fit['rate_matrix'][from_index, to_index]
            )
    return pd.DataFrame([row])


def _write_table(table, outfile):
    if outfile == '-':
        print(table.to_csv(sep='\t', index=False), end='')
    else:
        table.to_csv(outfile, sep='\t', index=False)


def _write_model_table(states, root_prior_mode, fit, outfile):
    if outfile in ['', None]:
        return
    table = _build_model_table(states, root_prior_mode, fit)
    _write_table(table, outfile)


def _annotate_tree(tree, states, posterior_by_node, observed_state_by_leaf):
    for node in tree.traverse():
        posterior = posterior_by_node[node]
        map_index = int(np.argmax(posterior))
        node.add_props(
            asr_state=states[map_index],
            asr_probability=float(posterior[map_index]),
        )
        if node.is_leaf:
            observed_state = observed_state_by_leaf.get(node.name)
            node.add_props(asr_observed_state='' if observed_state is None else observed_state)
        for state, probability in zip(states, posterior):
            node.add_props(**{'asr_p_{}'.format(_safe_column_state(state)): float(probability)})


def _write_annotated_tree(tree, states, posterior_by_node, observed_state_by_leaf, args):
    tree_out = getattr(args, 'tree_out', None)
    if tree_out in ['', None]:
        return
    _annotate_tree(tree, states, posterior_by_node, observed_state_by_leaf)
    tree_annotation = getattr(args, 'tree_annotation', 'map')
    if tree_annotation == 'state':
        props = ['asr_state']
    elif tree_annotation == 'probability':
        props = ['asr_probability']
    elif tree_annotation == 'map':
        props = ['asr_state', 'asr_probability']
    elif tree_annotation == 'all':
        props = ['asr_state', 'asr_probability']
        props += ['asr_observed_state'] + ['asr_p_{}'.format(_safe_column_state(state)) for state in states]
    else:
        raise ValueError("Unsupported '--tree_annotation': {}".format(tree_annotation))
    parser_format = tree.props.get(TREE_FORMAT_PROP, 0)
    if getattr(args, 'tree_outformat', 'auto') != 'auto':
        parser_format = int(getattr(args, 'tree_outformat'))
    tree_string = tree.write(props=props, parser=parser_format, format_root_node=True)
    if tree_out == '-':
        print(tree_string)
    else:
        with open(tree_out, 'w') as handle:
            handle.write(tree_string)


def _normalize_probability_vector(weights):
    weights = np.asarray(weights, dtype=float)
    weights = np.maximum(weights, 0.0)
    total = float(weights.sum())
    if total <= 0.0:
        return np.full(len(weights), 1.0 / float(len(weights)), dtype=float)
    return weights / total


def _sample_bridge_event_count(start_state, end_state, branch_length, rate_matrix, rng):
    if float(branch_length) == 0.0 or np.allclose(rate_matrix, 0.0):
        return 0
    omega = float(np.max(-np.diag(rate_matrix)))
    if omega <= 0.0:
        return 0
    transition_matrix = _transition_matrix(rate_matrix, branch_length)
    endpoint_probability = transition_matrix[start_state, end_state]
    if endpoint_probability <= 0.0:
        return 0
    lam = omega * float(branch_length)
    max_n = int(max(10, poisson.ppf(1.0 - 10 ** -12, lam)))
    if start_state != end_state:
        max_n = max(max_n, 1)
    r_matrix = np.eye(rate_matrix.shape[0], dtype=float) + (rate_matrix / omega)
    r_matrix = np.maximum(r_matrix, 0.0)
    r_matrix = r_matrix / r_matrix.sum(axis=1, keepdims=True)
    powers = [np.eye(rate_matrix.shape[0], dtype=float)]
    for _ in range(max_n):
        powers.append(powers[-1].dot(r_matrix))
    log_weights = list()
    for event_count in range(max_n + 1):
        bridge_probability = powers[event_count][start_state, end_state]
        if bridge_probability <= 0.0:
            log_weights.append(-math.inf)
            continue
        log_weights.append(poisson.logpmf(event_count, lam) + math.log(bridge_probability))
    normalizer = logsumexp(log_weights)
    if not math.isfinite(normalizer):
        return 0
    probabilities = np.exp(np.asarray(log_weights) - normalizer)
    probabilities = _normalize_probability_vector(probabilities)
    return int(rng.choice(np.arange(max_n + 1), p=probabilities))


def _sample_bridge_transition_counts(start_state, end_state, branch_length, rate_matrix, rng):
    event_count = _sample_bridge_event_count(start_state, end_state, branch_length, rate_matrix, rng)
    num_states = rate_matrix.shape[0]
    counts = defaultdict(int)
    if event_count == 0:
        return counts
    omega = float(np.max(-np.diag(rate_matrix)))
    r_matrix = np.eye(num_states, dtype=float) + (rate_matrix / omega)
    r_matrix = np.maximum(r_matrix, 0.0)
    r_matrix = r_matrix / r_matrix.sum(axis=1, keepdims=True)
    powers = [np.eye(num_states, dtype=float)]
    for _ in range(event_count):
        powers.append(powers[-1].dot(r_matrix))
    current_state = start_state
    for step_index in range(event_count):
        remaining = event_count - step_index - 1
        weights = r_matrix[current_state, :] * powers[remaining][:, end_state]
        probabilities = _normalize_probability_vector(weights)
        next_state = int(rng.choice(np.arange(num_states), p=probabilities))
        if next_state != current_state:
            counts[(current_state, next_state)] += 1
        current_state = next_state
    return counts


def _sample_node_states(tree, states, fit, rng):
    sampled_states = dict()
    sampled_states[tree] = int(rng.choice(np.arange(len(states)), p=fit['posterior_by_node'][tree]))
    for node in tree.traverse(strategy='preorder'):
        for child in node.get_children():
            parent_state = sampled_states[node]
            transition_matrix = fit['transition_matrices'][child]
            weights = transition_matrix[parent_state, :] * fit['inside'][child]
            probabilities = _normalize_probability_vector(weights)
            sampled_states[child] = int(rng.choice(np.arange(len(states)), p=probabilities))
    return sampled_states


def _simulate_stochastic_maps(tree, states, fit, num_simulations, seed=None):
    if num_simulations <= 0:
        raise ValueError("'--n_sim' must be positive when '--stochastic_map_out' is specified.")
    rng = np.random.default_rng(seed)
    node_to_branch_id = assign_branch_ids(tree)
    total_counts = defaultdict(int)
    any_counts = defaultdict(int)
    for _ in range(num_simulations):
        sampled_states = _sample_node_states(tree, states, fit, rng)
        simulation_counts = defaultdict(int)
        for node in tree.traverse():
            if node.is_root:
                continue
            branch_id = node_to_branch_id[node]
            branch_counts = _sample_bridge_transition_counts(
                sampled_states[node.up],
                sampled_states[node],
                node.dist,
                fit['rate_matrix'],
                rng,
            )
            for pair, count in branch_counts.items():
                key = (branch_id, pair[0], pair[1])
                total_counts[key] += count
                simulation_counts[key] += count
        for key in simulation_counts:
            any_counts[key] += 1
    rows = list()
    for node in tree.traverse():
        if node.is_root:
            continue
        branch_id = node_to_branch_id[node]
        for from_index, from_state in enumerate(states):
            for to_index, to_state in enumerate(states):
                if from_index == to_index:
                    continue
                key = (branch_id, from_index, to_index)
                total_count = total_counts.get(key, 0)
                any_count = any_counts.get(key, 0)
                rows.append({
                    'branch_id': branch_id,
                    'parent': node_to_branch_id[node.up],
                    'node_type': _node_type(node),
                    'name': '' if node.name in [None, ''] else str(node.name),
                    'from_state': from_state,
                    'to_state': to_state,
                    'total_count': total_count,
                    'mean_count': total_count / float(num_simulations),
                    'posterior_frequency': any_count / float(num_simulations),
                    'num_simulations': num_simulations,
                })
    return pd.DataFrame(rows)


def _write_stochastic_map(tree, states, fit, args):
    stochastic_map_out = getattr(args, 'stochastic_map_out', None)
    if stochastic_map_out in ['', None]:
        return
    table = _simulate_stochastic_maps(
        tree=tree,
        states=states,
        fit=fit,
        num_simulations=getattr(args, 'n_sim', 100),
        seed=getattr(args, 'seed', None),
    )
    _write_table(table, stochastic_map_out)


def asr_main(args):
    model = getattr(args, 'model', 'ER')
    if model not in SUPPORTED_MODELS:
        raise ValueError("Unsupported '--model': {}".format(model))
    output_mode = getattr(args, 'output', 'probabilities')
    if output_mode not in ['probabilities', 'map']:
        raise ValueError("Unsupported '--output': {}".format(output_mode))
    tree = read_tree(args.infile, args.format, args.quoted_node_names)
    _validate_tree_for_asr(tree)
    leaf_names = list(tree.leaf_names())
    states, observed_state_by_leaf, likelihood_by_leaf = _read_tip_states(
        trait_path=args.trait,
        state_column=args.state_column,
        tree_leaf_names=leaf_names,
        states_arg=getattr(args, 'states', None),
        missing_values_arg=getattr(args, 'missing_values', None),
        ambiguous_separator=getattr(args, 'ambiguous_separator', DEFAULT_AMBIGUOUS_SEPARATOR),
    )
    rate_bounds = _parse_rate_bounds(getattr(args, 'rate_bounds', None))
    posterior_by_node, fit = compute_mk_marginals(
        tree=tree,
        states=states,
        observed_state_by_leaf=observed_state_by_leaf,
        likelihood_by_leaf=likelihood_by_leaf,
        model=model,
        rate=getattr(args, 'rate', None),
        root_prior_mode=getattr(args, 'root_prior', 'equal'),
        rate_bounds=rate_bounds,
    )
    targets = _parse_targets(getattr(args, 'target', 'intnode,missing_tip'))
    table = _build_output_table(
        tree=tree,
        states=states,
        observed_state_by_leaf=observed_state_by_leaf,
        posterior_by_node=posterior_by_node,
        targets=targets,
        output_mode=output_mode,
    )
    _write_table(table, args.outfile)
    _write_model_table(states, getattr(args, 'root_prior', 'equal'), fit, getattr(args, 'model_out', None))
    _write_annotated_tree(tree, states, posterior_by_node, observed_state_by_leaf, args)
    _write_stochastic_map(tree, states, fit, args)

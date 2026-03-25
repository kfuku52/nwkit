import pandas as pd

from nwkit.util import (
    compute_node_ages,
    get_species_group_records,
    inspect_tree_text,
    is_rooted,
    read_tree,
    read_input_text,
    split_newick_stream,
    support_is_missing,
)


def _get_duplicate_leaf_names(tree):
    leaf_names = list(tree.leaf_names())
    seen = set()
    duplicates = list()
    for leaf_name in leaf_names:
        if leaf_name in seen:
            duplicates.append(str(leaf_name))
        else:
            seen.add(leaf_name)
    return sorted(set(duplicates))


def _get_empty_leaf_names(tree):
    return [
        str(leaf.name)
        for leaf in tree.leaves()
        if (leaf.name is None) or (str(leaf.name).strip() == '')
    ]


def _collect_tree_metrics(tree):
    num_nodes = 0
    num_singleton_nodes = 0
    num_multifurcation_nodes = 0
    num_zero_branch_nodes = 0
    num_negative_branch_nodes = 0
    num_missing_support_internal_nodes = 0
    for node in tree.traverse():
        num_nodes += 1
        if not node.is_leaf:
            num_children = len(node.get_children())
            if num_children == 1:
                num_singleton_nodes += 1
            elif num_children > 2:
                num_multifurcation_nodes += 1
            if (not node.is_root) and support_is_missing(node.support):
                num_missing_support_internal_nodes += 1
        if (not node.is_root) and (node.dist is not None):
            if float(node.dist) == 0.0:
                num_zero_branch_nodes += 1
            if float(node.dist) < 0.0:
                num_negative_branch_nodes += 1
    try:
        compute_node_ages(tree)
        is_ultrametric = True
    except ValueError:
        is_ultrametric = False
    return {
        'num_leaves': len(list(tree.leaves())),
        'num_nodes': num_nodes,
        'num_singleton_nodes': num_singleton_nodes,
        'num_multifurcation_nodes': num_multifurcation_nodes,
        'num_zero_branch_nodes': num_zero_branch_nodes,
        'num_negative_branch_nodes': num_negative_branch_nodes,
        'num_missing_support_internal_nodes': num_missing_support_internal_nodes,
        'is_rooted': is_rooted(tree),
        'is_ultrametric': is_ultrametric,
    }


def _validate_species_labels(tree, args):
    try:
        get_species_group_records(
            tree,
            option_name='--infile',
            context=" for 'validate'",
            args=args,
        )
    except ValueError as exc:
        message = str(exc)
        if 'could not be parsed' in message:
            return False, False
        if 'not monophyletic' in message:
            return True, False
        raise
    return True, True


def _build_parse_error_row(tree_id, inspection, issues):
    return {
        'tree_id': tree_id,
        'status': 'invalid',
        'parse_ok': False,
        'parse_error': inspection['parse_error'],
        'input_format': inspection['input_format'],
        'format_ambiguous': inspection['format_ambiguous'],
        'has_quoted_node_names': inspection['has_quoted_node_names'],
        'has_quoted_internal_node_names': inspection['has_quoted_internal_node_names'],
        'num_leaves': '',
        'num_nodes': '',
        'num_singleton_nodes': '',
        'num_multifurcation_nodes': '',
        'num_zero_branch_nodes': '',
        'num_negative_branch_nodes': '',
        'num_missing_support_internal_nodes': '',
        'is_rooted': '',
        'is_ultrametric': '',
        'leaf_names_unique': '',
        'num_duplicate_leaf_names': '',
        'num_empty_leaf_names': '',
        'leaf_set_matches_first': '',
        'rooting_matches_first': '',
        'species_parseable': '',
        'species_groups_monophyletic': '',
        'issues': ','.join(issues),
    }


def validate_main(args):
    require_rooted = getattr(args, 'require_rooted', False)
    require_ultrametric = getattr(args, 'require_ultrametric', False)
    require_same_leaf_set = getattr(args, 'require_same_leaf_set', True)
    require_same_rooting = getattr(args, 'require_same_rooting', False)
    require_binary = getattr(args, 'require_binary', False)
    require_all_support = getattr(args, 'require_all_support', False)
    require_unambiguous_format = getattr(args, 'require_unambiguous_format', False)
    require_unquoted_names = getattr(args, 'require_unquoted_names', False)
    check_species = getattr(args, 'check_species', False)
    fail_on_issue = getattr(args, 'fail_on_issue', False)
    raw_text = read_input_text(args.infile)
    tree_strings = split_newick_stream(raw_text)
    if len(tree_strings) == 0:
        raise Exception('Failed to parse the input trees.')
    rows = list()
    invalid_tree_ids = list()
    reference_leaf_set = None
    reference_rooted_state = None
    first_tree_parsed = None
    for tree_id, tree_string in enumerate(tree_strings, start=1):
        inspection = inspect_tree_text(
            newick_text=tree_string,
            format=args.format,
            quoted_node_names=args.quoted_node_names,
        )
        issues = list()
        if inspection['format_ambiguous'] and require_unambiguous_format:
            issues.append('format_ambiguous')
        if inspection['has_quoted_node_names'] and require_unquoted_names:
            issues.append('quoted_node_names')
        if not inspection['parse_ok']:
            issues.append('parse_error')
            if tree_id == 1:
                first_tree_parsed = False
            row = _build_parse_error_row(tree_id, inspection, issues)
            rows.append(row)
            invalid_tree_ids.append(tree_id)
            continue
        tree = read_tree(tree_string, args.format, args.quoted_node_names, quiet=True)
        metrics = _collect_tree_metrics(tree)
        duplicate_leaf_names = _get_duplicate_leaf_names(tree)
        empty_leaf_names = _get_empty_leaf_names(tree)
        if duplicate_leaf_names:
            issues.append('duplicate_leaf_names')
        if empty_leaf_names:
            issues.append('empty_leaf_names')
        if metrics['num_negative_branch_nodes'] > 0:
            issues.append('negative_branch_length')
        if require_rooted and (not metrics['is_rooted']):
            issues.append('not_rooted')
        if require_ultrametric and (not metrics['is_ultrametric']):
            issues.append('not_ultrametric')
        if require_binary and (
            (metrics['num_singleton_nodes'] > 0) or (metrics['num_multifurcation_nodes'] > 0)
        ):
            issues.append('not_binary')
        if require_all_support and (metrics['num_missing_support_internal_nodes'] > 0):
            issues.append('missing_support')
        leaf_name_set = set(tree.leaf_names())
        if tree_id == 1:
            leaf_set_matches_first = True
            reference_leaf_set = leaf_name_set
            first_tree_parsed = True
        elif first_tree_parsed:
            leaf_set_matches_first = (leaf_name_set == reference_leaf_set)
            if require_same_leaf_set and (not leaf_set_matches_first):
                issues.append('leaf_set_mismatch')
        else:
            leaf_set_matches_first = ''
        if tree_id == 1:
            rooting_matches_first = True
            reference_rooted_state = metrics['is_rooted']
        elif first_tree_parsed:
            rooting_matches_first = (metrics['is_rooted'] == reference_rooted_state)
            if require_same_rooting and (not rooting_matches_first):
                issues.append('rooting_mismatch')
        else:
            rooting_matches_first = ''
        species_parseable = ''
        species_groups_monophyletic = ''
        if check_species:
            species_parseable, species_groups_monophyletic = _validate_species_labels(tree, args)
            if not species_parseable:
                issues.append('species_parse_failed')
            elif not species_groups_monophyletic:
                issues.append('species_not_monophyletic')
        status = 'ok' if len(issues) == 0 else 'invalid'
        if status != 'ok':
            invalid_tree_ids.append(tree_id)
        row = {
            'tree_id': tree_id,
            'status': status,
            'parse_ok': True,
            'parse_error': '',
            'input_format': inspection['input_format'],
            'format_ambiguous': inspection['format_ambiguous'],
            'has_quoted_node_names': inspection['has_quoted_node_names'],
            'has_quoted_internal_node_names': inspection['has_quoted_internal_node_names'],
            **metrics,
            'leaf_names_unique': len(duplicate_leaf_names) == 0,
            'num_duplicate_leaf_names': len(duplicate_leaf_names),
            'num_empty_leaf_names': len(empty_leaf_names),
            'leaf_set_matches_first': leaf_set_matches_first,
            'rooting_matches_first': rooting_matches_first,
            'species_parseable': species_parseable,
            'species_groups_monophyletic': species_groups_monophyletic,
            'issues': ','.join(issues),
        }
        rows.append(row)
    out = pd.DataFrame(rows)
    if args.outfile == '-':
        print(out.to_csv(sep='\t', index=False), end='')
    else:
        out.to_csv(args.outfile, sep='\t', index=False)
    if fail_on_issue and invalid_tree_ids:
        raise ValueError(
            'Validation failed for tree_id(s): {}'.format(
                ', '.join(str(tree_id) for tree_id in invalid_tree_ids)
            )
        )

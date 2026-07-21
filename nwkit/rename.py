import csv
import re
import sys

from nwkit.util import get_target_nodes, read_tree, support_is_missing, validate_unique_named_leaves, write_tree


def read_name_tsv(path):
    handle = sys.stdin if path == '-' else open(path, newline='')
    try:
        reader = csv.DictReader(handle, delimiter='\t')
        fieldnames = reader.fieldnames or list()
        required = {'old_name', 'new_name'}
        if not required.issubset(fieldnames):
            raise ValueError("--name-tsv must contain 'old_name' and 'new_name' columns.")
        mapping = dict()
        for row in reader:
            raw_old_name = row.get('old_name')
            raw_new_name = row.get('new_name')
            old_name = '' if raw_old_name is None else str(raw_old_name)
            new_name = '' if raw_new_name is None else str(raw_new_name)
            if old_name.strip() == '':
                raise ValueError("--name-tsv contains an empty 'old_name' value.")
            if new_name.strip() == '':
                raise ValueError("--name-tsv contains an empty 'new_name' value.")
            if old_name in mapping:
                raise ValueError("Duplicated 'old_name' entries are not supported in --name-tsv: {}".format(old_name))
            mapping[old_name] = new_name
    finally:
        if handle is not sys.stdin:
            handle.close()
    if len(mapping) == 0:
        raise ValueError('--name-tsv is empty.')
    return mapping


def _rename_by_mapping(tree, args):
    mapping = read_name_tsv(args.name_tsv)
    target_nodes = get_target_nodes(tree=tree, target=args.target)
    target_name_to_nodes = dict()
    for node in target_nodes:
        node_name = str(node.name or '')
        if node_name.strip() == '':
            continue
        if node_name not in target_name_to_nodes:
            target_name_to_nodes[node_name] = list()
        target_name_to_nodes[node_name].append(node)
    for old_name, matching_nodes in target_name_to_nodes.items():
        if (old_name in mapping) and (len(matching_nodes) > 1):
            raise ValueError(
                "Matched multiple target nodes for old_name '{}' in '--target {}'.".format(
                    old_name,
                    args.target,
                )
            )
    renamed_count = 0
    used_old_names = set()
    for old_name, new_name in mapping.items():
        matching_nodes = target_name_to_nodes.get(old_name, list())
        if len(matching_nodes) == 0:
            continue
        matching_nodes[0].name = new_name
        renamed_count += 1
        used_old_names.add(old_name)
    unused_old_names = sorted(set(mapping.keys()) - used_old_names)
    if args.require_all_old_names and unused_old_names:
        raise ValueError(
            "The following old_name values in '--name-tsv' were not found among target nodes: {}".format(
                ', '.join(unused_old_names)
            )
        )
    return tree, renamed_count, len(mapping)


def _rename_by_regex(tree, args):
    target_nodes = get_target_nodes(tree=tree, target=args.target)
    renamed_count = 0
    matched_count = 0
    pattern = re.compile(args.pattern)
    for node in target_nodes:
        node_name = str(node.name or '')
        if node_name == '':
            continue
        new_name, num_subs = pattern.subn(args.replacement, node_name)
        if num_subs == 0:
            continue
        matched_count += 1
        if new_name != node_name:
            node.name = new_name
            renamed_count += 1
    if args.require_match and matched_count == 0:
        raise ValueError("No target node names matched '--pattern'.")
    return tree, renamed_count, matched_count


def _has_named_internal_nodes(tree):
    for node in tree.traverse():
        if node.is_leaf:
            continue
        if str(node.name or '').strip() != '':
            return True
    return False


def _has_meaningful_internal_support(tree):
    for node in tree.traverse():
        if node.is_root or node.is_leaf:
            continue
        if not support_is_missing(node.support):
            return True
    return False


def _resolve_output_format(tree, args):
    if args.outformat != 'auto':
        return args.outformat
    if args.target == 'leaf':
        return 'auto'
    has_internal_names = _has_named_internal_nodes(tree)
    has_internal_support = _has_meaningful_internal_support(tree)
    if has_internal_names and has_internal_support:
        sys.stderr.write(
            "Warning: Preserving the input/output format because standard Newick cannot retain "
            "both internal node names and support values in auto mode.\n"
        )
        return 'auto'
    if has_internal_names:
        return 1
    return 'auto'


def rename_main(args):
    name_tsv = getattr(args, 'name_tsv', None)
    pattern = getattr(args, 'pattern', None)
    replacement = getattr(args, 'replacement', None)
    require_match = getattr(args, 'require_match', True)
    has_name_tsv = name_tsv not in ['', None]
    has_pattern = pattern not in ['', None]
    if has_name_tsv == has_pattern:
        raise ValueError("Specify exactly one of '--name-tsv' or '--pattern/--replacement'.")
    if has_pattern and (replacement is None):
        raise ValueError("'--replacement' is required when '--pattern' is specified.")
    tree = read_tree(args.infile, args.format, args.quoted_node_names)
    if has_name_tsv:
        tree, renamed_count, matched_count = _rename_by_mapping(tree=tree, args=args)
        mapping_label = '{} mapping row(s)'.format(matched_count)
    else:
        args.pattern = pattern
        args.replacement = replacement
        args.require_match = require_match
        tree, renamed_count, matched_count = _rename_by_regex(tree=tree, args=args)
        mapping_label = '{} matched target node(s)'.format(matched_count)
    if args.check_leaf_uniqueness:
        validate_unique_named_leaves(tree, option_name='--infile', context=" after 'rename'")
    sys.stderr.write('Renamed {} target node(s) using {}.\n'.format(renamed_count, mapping_label))
    outformat = _resolve_output_format(tree=tree, args=args)
    write_tree(tree, args, format=outformat)

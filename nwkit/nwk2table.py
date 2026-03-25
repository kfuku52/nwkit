import pandas as pd

from nwkit.util import (
    assign_branch_ids,
    compute_node_ages,
    read_tree,
    support_is_missing,
)


def _support_value_or_empty(node):
    support = getattr(node, 'support', None)
    if support_is_missing(support):
        return ''
    return float(support)


def _sister_branch_id(node, node_to_branch_id):
    if node.is_root:
        return -1
    siblings = node.up.get_children()
    if len(siblings) == 2:
        sister = siblings[1] if siblings[0] is node else siblings[0]
        return node_to_branch_id[sister]
    sisters = node.get_sisters()
    if len(sisters) > 0:
        return node_to_branch_id[sisters[0]]
    return -1


def nwk2table_main(args):
    tree = read_tree(args.infile, args.format, args.quoted_node_names)
    node_to_branch_id = assign_branch_ids(tree)
    age_by_node = None
    if args.age:
        try:
            age_by_node = compute_node_ages(tree)
        except ValueError as exc:
            raise ValueError("Tree must be ultrametric when '--age yes'.") from exc
    rows = list()
    for node in tree.traverse():
        row = {
            'branch_id': node_to_branch_id[node],
            'parent': -1 if node.is_root else node_to_branch_id[node.up],
            'name': '' if node.name in [None, ''] else str(node.name),
            'dist': '' if node.dist is None else float(node.dist),
            'support': _support_value_or_empty(node),
        }
        if args.sister:
            row['sister'] = _sister_branch_id(node, node_to_branch_id)
        if args.age:
            row['age'] = age_by_node[node]
        rows.append(row)
    table = pd.DataFrame(rows).sort_values('branch_id').reset_index(drop=True)
    if args.outfile == '-':
        print(table.to_csv(sep='\t', index=False), end='')
    else:
        table.to_csv(args.outfile, sep='\t', index=False)

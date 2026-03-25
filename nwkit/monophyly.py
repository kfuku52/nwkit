import pandas as pd

from nwkit.util import get_species_group_records, read_tree


def _read_trait_groups(path, tree_leaf_name_set, group_by):
    trait_df = pd.read_csv(path, sep='\t')
    if 'leaf_name' not in trait_df.columns:
        raise ValueError("Column 'leaf_name' is required in '--trait'.")
    if group_by not in trait_df.columns:
        raise ValueError("Column '{}' specified by '--group_by' was not found in '--trait'.".format(group_by))
    trait_df = trait_df[trait_df['leaf_name'].isin(tree_leaf_name_set)].copy()
    duplicated_leaf_names = trait_df.loc[
        trait_df['leaf_name'].duplicated(keep=False), 'leaf_name'
    ].unique().tolist()
    if duplicated_leaf_names:
        duplicated_leaf_names = sorted(str(name) for name in duplicated_leaf_names)
        raise ValueError("Duplicated 'leaf_name' entries in '--trait': {}".format(', '.join(duplicated_leaf_names)))
    trait_df = trait_df[~trait_df[group_by].isna()].copy()
    trait_df[group_by] = trait_df[group_by].astype(str)
    group_to_leaf_names = dict()
    for group_name, group_df in trait_df.groupby(group_by, sort=True):
        leaf_names = sorted(group_df['leaf_name'].astype(str).tolist())
        if len(leaf_names) > 0:
            group_to_leaf_names[str(group_name)] = leaf_names
    return group_to_leaf_names


def _get_groups(tree, args):
    if args.trait not in ['', None]:
        if args.group_by in ['', None]:
            raise ValueError("'--group_by' is required when '--trait' is specified.")
        return _read_trait_groups(
            path=args.trait,
            tree_leaf_name_set=set(tree.leaf_names()),
            group_by=args.group_by,
        )
    _, species_to_leaf_names, _ = get_species_group_records(
        tree,
        option_name='--infile',
        context=" for 'monophyly'",
        args=args,
    )
    return species_to_leaf_names


def monophyly_main(args):
    tree = read_tree(args.infile, args.format, args.quoted_node_names)
    group_to_leaf_names = _get_groups(tree, args)
    rows = list()
    non_monophyletic_groups = list()
    for group_name in sorted(group_to_leaf_names.keys()):
        leaf_names = sorted(str(name) for name in group_to_leaf_names[group_name])
        if len(leaf_names) == 0:
            continue
        is_monophyletic, status, _ = tree.check_monophyly(
            values=leaf_names,
            prop='name',
            unrooted=args.unrooted,
        )
        if len(leaf_names) == 1:
            mrca = tree.search_nodes(name=leaf_names[0]).__next__()
        else:
            mrca = tree.common_ancestor(leaf_names)
        mrca_leaf_names = sorted(str(name) for name in mrca.leaf_names())
        intruder_leaf_names = sorted(set(mrca_leaf_names) - set(leaf_names))
        row = {
            'group': group_name,
            'status': status,
            'is_monophyletic': bool(is_monophyletic),
            'num_target_leaves': len(leaf_names),
            'num_mrca_leaves': len(mrca_leaf_names),
            'num_intruder_leaves': len(intruder_leaf_names),
            'target_leaves': ','.join(leaf_names),
            'intruder_leaves': ','.join(intruder_leaf_names),
        }
        rows.append(row)
        if not is_monophyletic:
            non_monophyletic_groups.append(group_name)
    out = pd.DataFrame(rows)
    if args.outfile == '-':
        print(out.to_csv(sep='\t', index=False), end='')
    else:
        out.to_csv(args.outfile, sep='\t', index=False)
    if args.fail_on_non_monophyly and non_monophyletic_groups:
        raise ValueError(
            'Non-monophyletic group(s): {}'.format(', '.join(sorted(non_monophyletic_groups)))
        )

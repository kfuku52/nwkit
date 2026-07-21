import pandas as pd
import sys

from nwkit.util import (
    get_subtree_leaf_name_sets,
    read_tip_table,
    read_tree,
    validate_unique_named_leaves,
    write_tree,
)

def read_trait(args, tree):
    if args.trait is None:
        sys.stderr.write("'--trait' not specified. Sampling leaves at random.\n")
        trait_df = pd.DataFrame({'leaf_name': list(tree.leaf_names())})
        return trait_df
    leaf_names_list = list(tree.leaf_names())
    leaf_name_set = set(leaf_names_list)
    trait_df, _, missing_leaf_names = read_tip_table(
        args.trait,
        option_name='--trait',
        tree_leaf_names=leaf_names_list,
        unmatched=getattr(args, 'unmatched', 'warn'),
        missing_values=getattr(args, 'missing_values', None),
    )
    trait_df = trait_df[trait_df['leaf_name'].isin(leaf_name_set)].copy()
    if len(missing_leaf_names) > 0:
        trait_df = pd.concat(
            [trait_df, pd.DataFrame({'leaf_name': missing_leaf_names})],
            ignore_index=True,
        )
    return trait_df

def mark_traits_to_nodes(tree, trait_df, args):
    if (args.group_by is not None) and (args.group_by not in trait_df.columns):
        raise ValueError("Column '{}' specified by '--group-by' was not found in '--trait'.".format(args.group_by))
    if args.group_by is None:
        leafname2trait = {leaf_name: None for leaf_name in trait_df['leaf_name']}
    else:
        leafname2trait = {leaf_name: trait for leaf_name, trait in zip(trait_df['leaf_name'], trait_df[args.group_by])}
    for node in tree.traverse('postorder'):
        if node.is_leaf:
            # Treat unmatched/missing leaf labels as missing trait information.
            node.add_props(trait=leafname2trait.get(node.name))
        else:
            trait = None
            for child in node.children:
                child_trait = child.props.get('trait')
                if pd.isna(child_trait):
                    continue
                elif child_trait == '_MIXED_':
                    trait = '_MIXED_'
                    break
                elif pd.isna(trait):
                    trait = child_trait
                elif trait == child_trait:
                    continue
                else:
                    trait = '_MIXED_'
                    break
            node.add_props(trait=trait)
    return tree

def iter_frontier_nodes(root, stop_condition):
    stack = [root]
    while stack:
        node = stack.pop()
        if stop_condition(node):
            yield node
            continue
        stack.extend(reversed(node.children))

def add_group_ids(trait_df, marked_tree):
    group_id = 1
    leafname2group = {}
    subtree_leaf_name_sets = get_subtree_leaf_name_sets(marked_tree)
    for node in iter_frontier_nodes(root=marked_tree.root, stop_condition=lambda node: pd.isna(node.props.get('trait')) or node.props.get('trait') != '_MIXED_'):
        if node.is_leaf:
            leafname2group[node.name] = group_id
            group_id += 1
        else:
            for leaf_name in subtree_leaf_name_sets[node]:
                leafname2group[leaf_name] = group_id
            group_id += 1
    trait_df['group'] = trait_df['leaf_name'].map(leafname2group)
    return trait_df

def sample_from_groups(trait_df, args):
    if 'group' not in trait_df.columns:
        raise ValueError("Column 'group' not found. Run group assignment before sampling.")
    if (args.filter_by is not None) and (args.filter_by not in trait_df.columns):
        raise ValueError("Column '{}' specified by '--filter-by' was not found in '--trait'.".format(args.filter_by))
    shuffled_df = trait_df.sample(frac=1, random_state=getattr(args, 'seed', None))
    if args.prioritize_non_missing and args.group_by is not None:
        if args.filter_by is None:
            sorted_df = shuffled_df.sort_values(by=[args.group_by], ascending=True, kind='mergesort')
        else:
            if args.filter_mode == 'ascending':
                sorted_df = shuffled_df.sort_values(
                    by=[args.group_by, args.filter_by],
                    ascending=[True, True],
                    kind='mergesort',
                )
            elif args.filter_mode == 'descending':
                sorted_df = shuffled_df.sort_values(
                    by=[args.group_by, args.filter_by],
                    ascending=[False, False],
                    kind='mergesort',
                )
            else:
                raise ValueError(f"Invalid value for '--filter-mode': {args.filter_mode}")
    else:
        if args.filter_by is None:
            sorted_df = shuffled_df
        else:
            if args.filter_mode == 'ascending':
                sorted_df = shuffled_df.sort_values(
                    by=[args.filter_by],
                    ascending=True,
                    kind='mergesort',
                )
            elif args.filter_mode == 'descending':
                sorted_df = shuffled_df.sort_values(
                    by=[args.filter_by],
                    ascending=False,
                    kind='mergesort',
                )
            else:
                raise ValueError(f"Invalid value for '--filter-mode': {args.filter_mode}")
    sampled_df = sorted_df.groupby('group', group_keys=False).head(args.retain_per_clade)
    if 'group' in sampled_df.columns:
        sampled_df = sampled_df.drop(columns=['group'])
    sampled_df = sampled_df.merge(trait_df[['leaf_name', 'group']], on=['leaf_name'])
    sampled_df = sampled_df[trait_df.columns]
    return sampled_df

def iter_contrastive_nodes(root):
    if root.props.get('trait') != '_MIXED_':
        return
    stack = [root]
    while stack:
        node = stack.pop()
        flag_has_mixed_child = False
        for child in node.children:
            if child.props.get('trait') == '_MIXED_':
                stack.append(child)
                flag_has_mixed_child = True
        if not flag_has_mixed_child:
            yield node

def add_contrastive_clade_ids(trait_df, marked_tree):
    contrastive_clade_id = 1
    leafname2contrastive_clade = {}
    subtree_leaf_name_sets = get_subtree_leaf_name_sets(marked_tree)
    for node in iter_contrastive_nodes(root=marked_tree.root):
        for leaf_name in subtree_leaf_name_sets[node]:
            leafname2contrastive_clade[leaf_name] = contrastive_clade_id
        contrastive_clade_id += 1
    trait_df['contrastive_clade'] = trait_df['leaf_name'].map(leafname2contrastive_clade).astype('Int64')
    return trait_df

def skim_main(args):
    if args.retain_per_clade < 1:
        raise ValueError("'--retain-per-clade' must be a positive integer.")
    group_table_prefix = getattr(args, 'group_table_prefix', None)
    if getattr(args, 'output_groupfile', False) and group_table_prefix in ['', None]:
        if args.outfile == '-':
            raise ValueError("Legacy '--output-groupfile yes' requires '--outfile' to be a file path, not '-'.")
        group_table_prefix = args.outfile.removesuffix('.nwk')
    if group_table_prefix == '-':
        raise ValueError("'--group-table-prefix' must be a file path, not '-'.")
    tree = read_tree(args.infile, args.format, args.quoted_node_names)
    validate_unique_named_leaves(tree, option_name='--infile', context=" for 'skim'")
    trait_df = read_trait(args, tree)
    marked_tree = mark_traits_to_nodes(tree, trait_df, args)
    trait_df = add_group_ids(trait_df, marked_tree)
    if args.only_contrastive_clades:
        trait_df = add_contrastive_clade_ids(trait_df, marked_tree)
        sampled_trait_df = sample_from_groups(trait_df, args)
        sampled_trait_df = sampled_trait_df[~sampled_trait_df['contrastive_clade'].isna()]
    else:
        sampled_trait_df = sample_from_groups(trait_df, args)
    if sampled_trait_df.empty:
        raise ValueError('No leaves were selected for output. Adjust trait/grouping/sampling options.')
    if group_table_prefix not in ['', None]:
        trait_df.sort_values('group').to_csv(f'{group_table_prefix}.all.tsv', sep='\t', index=False)
        sampled_trait_df.sort_values('group').to_csv(f'{group_table_prefix}.sampled.tsv', sep='\t', index=False)
    tree.prune(sampled_trait_df['leaf_name'].tolist(), preserve_branch_length=True)
    write_tree(tree, args, format=args.outformat)

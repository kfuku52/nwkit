import pandas
import sys
from nwkit.util import *

def read_trait(args, tree):
    if args.trait is None:
        sys.stderr.write(f"'--trait' not specified. Sampling leaves at random.\n")
        trait_df = pandas.DataFrame({'leaf_name': list(tree.leaf_names())})
        return trait_df
    trait_df = pandas.read_csv(args.trait, sep='\t')
    leaf_names_list = list(tree.leaf_names())
    leaf_name_set = set(leaf_names_list)
    trait_df = trait_df[trait_df['leaf_name'].isin(leaf_name_set)]
    observed_leaf_names = set(trait_df['leaf_name'].tolist())
    missing_leaf_names = [leaf_name for leaf_name in leaf_names_list if leaf_name not in observed_leaf_names]
    if len(missing_leaf_names) > 0:
        log_txt = ''.join(
            f"'{leaf_name}' not found in '{args.trait}'. Treating its trait information as missing.\n"
            for leaf_name in missing_leaf_names
        )
        sys.stderr.write(log_txt)
        trait_df = pandas.concat(
            [trait_df, pandas.DataFrame({'leaf_name': missing_leaf_names})],
            ignore_index=True,
        )
    return trait_df

def mark_traits_to_nodes(tree, trait_df, args):
    if args.group_by is None:
        leafname2trait = {leaf_name: None for leaf_name in trait_df['leaf_name']}
    else:
        leafname2trait = {leaf_name: trait for leaf_name, trait in zip(trait_df['leaf_name'], trait_df[args.group_by])}
    for node in tree.traverse('postorder'):
        if node.is_leaf:
            node.add_props(trait=leafname2trait[node.name])
        else:
            trait = None
            for child in node.children:
                child_trait = child.props.get('trait')
                if pandas.isna(child_trait):
                    continue
                elif child_trait == '_MIXED_':
                    trait = '_MIXED_'
                    break
                elif pandas.isna(trait):
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
    for node in iter_frontier_nodes(root=marked_tree.root, stop_condition=lambda node: pandas.isna(node.props.get('trait')) or node.props.get('trait') != '_MIXED_'):
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
    shuffled_df = trait_df.sample(frac=1)
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
                raise ValueError(f"Invalid value for '--filter_mode': {args.filter_mode}")
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
                raise ValueError(f"Invalid value for '--filter_mode': {args.filter_mode}")
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
    tree = read_tree(args.infile, args.format, args.quoted_node_names)
    trait_df = read_trait(args, tree)
    marked_tree = mark_traits_to_nodes(tree, trait_df, args)
    trait_df = add_group_ids(trait_df, marked_tree)
    if args.only_contrastive_clades:
        trait_df = add_contrastive_clade_ids(trait_df, marked_tree)
        sampled_trait_df = sample_from_groups(trait_df, args)
        sampled_trait_df = sampled_trait_df[~sampled_trait_df['contrastive_clade'].isna()]
    else:
        sampled_trait_df = sample_from_groups(trait_df, args)
    if args.output_groupfile:
        filename = args.outfile.removesuffix('.nwk')
        trait_df.sort_values('group').to_csv(f'{filename}.all.tsv', sep='\t', index=False)
        sampled_trait_df.sort_values('group').to_csv(f'{filename}.sampled.tsv', sep='\t', index=False)
    tree.prune(sampled_trait_df['leaf_name'].tolist(), preserve_branch_length=True)
    write_tree(tree, args, format=args.outformat)

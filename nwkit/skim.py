import pandas
import sys
from nwkit.util import *

def read_trait(args, tree):
    if args.trait is None:
        sys.stderr.write(f"'--trait' not specified. Sampling leaves at random.\n")
        trait_df = pandas.DataFrame({'leaf_name': tree.get_leaf_names()})
        return trait_df
    trait_df = pandas.read_csv(args.trait, sep='\t')
    trait_df = trait_df[trait_df['leaf_name'].isin(tree.get_leaf_names())]
    for leaf_name in tree.get_leaf_names():
        if leaf_name not in trait_df['leaf_name'].values:
            sys.stderr.write(f"'{leaf_name}' not found in '{args.trait}'. Treating its trait information as missing.\n")
            trait_df = pandas.concat([trait_df, pandas.DataFrame({'leaf_name': [leaf_name]})], ignore_index=True)
    return trait_df

def mark_traits_to_nodes(tree, trait_df, args):
    if args.group_by is None:
        leafname2trait = {leaf_name: pandas.NA for leaf_name in trait_df['leaf_name']}
    else:
        leafname2trait = {leaf_name: trait for leaf_name, trait in zip(trait_df['leaf_name'], trait_df[args.group_by])}
    for node in tree.traverse('postorder'):
        if node.is_leaf():
            node.add_feature('trait', leafname2trait[node.name])
        else:
            trait = pandas.NA
            for child in node.children:
                if pandas.isna(child.trait):
                    continue
                elif child.trait == '_MIXED_':
                    trait = '_MIXED_'
                    break
                elif pandas.isna(trait):
                    trait = child.trait
                elif trait == child.trait:
                    continue
                else:
                    trait = '_MIXED_'
                    break
            node.add_feature('trait', trait)
    return tree

def iter_frontier_nodes(root, stop_condition):
    stack = [root]
    while stack:
        node = stack.pop()
        if stop_condition(node):
            yield node
            continue
        stack.extend(reversed(node.children))

def add_group_ids(trait_df, tree, args):
    group_id = 1
    leafname2group = {}
    marked_tree = mark_traits_to_nodes(tree, trait_df, args)
    for node in iter_frontier_nodes(root=marked_tree.get_tree_root(), stop_condition=lambda node: pandas.isna(node.trait) or node.trait != '_MIXED_'):
        if node.is_leaf():
            leafname2group[node.name] = group_id
            group_id += 1
        else:
            for leaf_name in node.get_leaf_names():
                leafname2group[leaf_name] = group_id
            group_id += 1
    trait_df['group'] = trait_df['leaf_name'].map(leafname2group)
    return trait_df

def sample_from_groups(trait_df, args):
    if args.prioritize_non_missing and args.group_by is not None:
        if args.filter_by is None:
            sampled_df = trait_df.groupby('group', group_keys=False).apply(lambda g: g.sample(frac=1).sort_values(by=[args.group_by]).head(args.retain_per_clade), include_groups=False)
        else:
            if args.filter_mode == 'ascending':
                sampled_df = trait_df.groupby('group', group_keys=False).apply(lambda g: g.sample(frac=1).sort_values(by=[args.group_by, args.filter_by], ascending=True).head(args.retain_per_clade), include_groups=False)
            elif args.filter_mode == 'descending':
                sampled_df = trait_df.groupby('group', group_keys=False).apply(lambda g: g.sample(frac=1).sort_values(by=[args.group_by, args.filter_by], ascending=False).head(args.retain_per_clade), include_groups=False)
            else:
                raise ValueError(f"Invalid value for '--filter_mode': {args.filter_mode}")
    else:
        if args.filter_by is None:
            sampled_df = trait_df.groupby('group', group_keys=False).apply(lambda g: g.sample(frac=1).head(args.retain_per_clade), include_groups=False)
        else:
            if args.filter_mode == 'ascending':
                sampled_df = trait_df.groupby('group', group_keys=False).apply(lambda g: g.sample(frac=1).sort_values(by=args.filter_by, ascending=True).head(args.retain_per_clade), include_groups=False)
            elif args.filter_mode == 'descending':
                sampled_df = trait_df.groupby('group', group_keys=False).apply(lambda g: g.sample(frac=1).sort_values(by=args.filter_by, ascending=False).head(args.retain_per_clade), include_groups=False)
            else:
                raise ValueError(f"Invalid value for '--filter_mode': {args.filter_mode}")
    return sampled_df.merge(trait_df[['leaf_name', 'group']], on=['leaf_name'])

def skim_main(args):
    tree = read_tree(args.infile, args.format, args.quoted_node_names)
    trait_df = read_trait(args, tree)
    trait_df = add_group_ids(trait_df, tree, args)
    sampled_trait_df = sample_from_groups(trait_df, args)
    tree.prune(sampled_trait_df['leaf_name'].tolist(), preserve_branch_length=True)
    if args.output_groupfile:
        filename = args.outfile.removesuffix('.nwk')
        trait_df.sort_values('group').to_csv(f'{filename}.all.tsv', sep='\t', index=False)
        sampled_trait_df.sort_values('group').to_csv(f'{filename}.sampled.tsv', sep='\t', index=False)
    write_tree(tree, args, format=args.outformat)
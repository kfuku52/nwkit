import os
from nwkit.util import *

def info_main(args):
    tree = read_tree(args.infile, args.format, args.quoted_node_names)
    tree_length = 0
    num_leaves = 0
    num_nodes = 0
    num_singleton_node = 0
    num_multifurcation_node = 0
    num_zero_branch_nodes = 0
    num_negative_branch_nodes = 0
    species_name_set = set()
    for node in tree.traverse():
        num_nodes += 1
        if node.is_leaf:
            num_leaves += 1
            species_name = extract_species_label(node.name, species_regex=args.species_regex)
            if species_name is None:
                species_name = str(node.name or '').strip()
            if species_name not in ['', None]:
                species_name_set.add(species_name)
        else:
            num_children = len(node.get_children())
            if num_children == 1:
                num_singleton_node += 1
            elif num_children > 2:
                num_multifurcation_node += 1
        if (not node.is_root) and (node.dist is not None):
            tree_length += node.dist
            if node.dist == 0:
                num_zero_branch_nodes += 1
        if (node.dist is not None) and (node.dist < 0):
            num_negative_branch_nodes += 1
    species_names = sorted(species_name_set)
    num_species = len(species_names)
    print(f'Tree file PATH: {os.path.realpath(args.infile)}')
    print(f'Tree length: {tree_length}')
    print(f'Number of leaves: {num_leaves}')
    print(f'Number of nodes: {num_nodes}')
    print(f'Number of singleton nodes: {num_singleton_node}')
    print(f'Number of multifurcation nodes: {num_multifurcation_node}')
    print(f'Number of nodes with zero branch length: {num_zero_branch_nodes}')
    print(f'Number of nodes with negative branch length: {num_negative_branch_nodes}')
    if args.species_regex == DEFAULT_SPECIES_REGEX:
        print(f'Number of species in the leaf name convention of GENUS_SPECIES_GENEID: {num_species}')
    else:
        print(f'Number of species parsed by --species_regex: {num_species}')
    print(f'Species names: {", ".join(species_names)}')

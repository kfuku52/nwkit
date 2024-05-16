import os
from nwkit.util import *

def info_main(args):
    tree = read_tree(args.infile, args.format, args.quoted_node_names)
    print(f'Tree file PATH: {os.path.realpath(args.infile)}')
    print(f'Tree length: {sum([ n.dist for n in tree.traverse() if not n.is_root() ])}')
    print(f'Number of leaves: {len(list(tree.iter_leaves()))}')
    print(f'Number of nodes: {len(list(tree.traverse()))}')
    num_singleton_node = 0
    for node in tree.traverse():
        if node.is_leaf():
            continue
        if len(node.get_children())==1:
            num_singleton_node += 1
    print(f'Number of singleton nodes: {num_singleton_node}')
    num_multifurcation_node = 0
    for node in tree.traverse():
        if node.is_leaf():
            continue
        if len(node.get_children())>2:
            num_multifurcation_node += 1
    print(f'Number of multifurcation nodes: {num_multifurcation_node}')
    print(f'Number of nodes with zero branch length: {len([ n.name for n in tree.traverse() if ((n.dist==0) & (not n.is_root())) ])}')
    print(f'Number of nodes with negative branch length: {len([ n for n in tree.traverse() if n.dist<0 ])}')
    species_names = sorted(list(set([ '_'.join(n.name.split('_')[:2]) for n in tree.iter_leaves() ])))
    num_species = len(species_names)
    print(f'Number of species in the leaf name convention of GENUS_SPECIES_GENEID: {num_species}')
    print(f'Species names: {", ".join(species_names)}')
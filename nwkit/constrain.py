import ete3

import re
import warnings

from nwkit.util import *

def initialize_tree(tree):
    for leaf in tree.get_leaves():
        leaf.has_taxon = False
        leaf.taxon_names = []
    return tree

def name_to_taxid(sp, ncbi):
    name2taxid = ncbi.get_name_translator([sp,])
    if (len(name2taxid)==0):
        genus = re.sub(' .*', '', sp)
        name2taxid = ncbi.get_name_translator([genus,])
    return name2taxid

def match_taxa(tree, splist):
    ncbi = ete3.NCBITaxa()
    leaf_names = [ ln.replace('_','') for ln in tree.get_leaf_names() ]
    leaf_name_set = set(leaf_names)
    for sp in splist:
        name2taxid = name_to_taxid(sp, ncbi)
        if (len(name2taxid)==0):
            warnings.warn('No genus-level match found in the NCBI database. Excluded from the output: {}'.format(sp))
            continue
        lineage = ncbi.get_lineage(name2taxid[sp][0])
        ancestor_names = ncbi.get_taxid_translator(lineage)
        ancestor = leaf_name_set.intersection(set(ancestor_names.values()))
        ancestor = list(ancestor)
        if (len(ancestor)>1):
            txt = 'Multiple hits. Taxon in the list = {}, Taxa in the tree = {} '
            warnings.warn(txt.format(sp, ','.join(list(ancestor))))
            continue
        elif (len(ancestor)==0):
            txt = 'No hit. Excluded from the output. Taxon = {}'
            warnings.warn(txt.format(sp, ','.join(list(ancestor))))
            continue
        for leaf in tree.get_leaves():
            if (leaf.name==ancestor[0]):
                leaf.has_taxon = True
                leaf.taxon_names.append(sp)
                break
    return tree

def delete_nomatch_leaves(tree):
    for leaf in tree.get_leaves():
        if not leaf.has_taxon:
            warnings.warn('Deleting taxon with no match: {}'.format(leaf.name))
            leaf.delete()
    return tree

def polytomize_one2many_matches(tree):
    for leaf in tree.get_leaves():
        if (len(leaf.taxon_names)==0):
            raise Exception('Leaf should have at least one match: {}'.format(leaf.name))
        if (len(leaf.taxon_names)==1):
            leaf.name = leaf.taxon_names[0]
        else:
            for i in range(len(leaf.taxon_names)):
                if (i==0):
                    leaf.name = leaf.taxon_names[i]
                else:
                    leaf.up.add_child(name=leaf.taxon_names[i])
    return tree

def constrain_main(args):
    tree = read_tree(args.infile, args.format)
    tree = initialize_tree(tree)
    with open(args.species_list, 'r') as f:
        splist = f.read().split('\n')
    splist = [ sp.replace('_', '') for sp in splist ]
    tree = match_taxa(tree, splist)
    tree = delete_nomatch_leaves(tree)
    tree = polytomize_one2many_matches(tree)
    write_tree(tree, args, format=9)
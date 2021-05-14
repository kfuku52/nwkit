import ete3
import numpy

import pkg_resources
import re
import sys
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

def get_lineage(sp, ncbi):
    name2taxid = name_to_taxid(sp, ncbi)
    if (len(name2taxid)==0):
        txt = 'Genus-level match was not found in the NCBI database. Excluded from the output: {}'
        warnings.warn(txt.format(sp))
        return []
    lineage = ncbi.get_lineage(name2taxid[sp][0])
    return lineage

def match_taxa(tree, splist):
    ncbi = ete3.NCBITaxa()
    leaf_names = [ ln.replace('_','') for ln in tree.get_leaf_names() ]
    leaf_name_set = set(leaf_names)
    for sp in splist:
        lineage = get_lineage(sp, ncbi)
        ancestor_names = ncbi.get_taxid_translator(lineage)
        ancestor = leaf_name_set.intersection(set(ancestor_names.values()))
        ancestor = list(ancestor)
        if (len(ancestor)>1):
            txt = 'Multiple hits. Excluded from the output. Taxon in the list = {}, Taxa in the tree = {} '
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
            sys.stderr.write('Deleting taxon in the tree with no match: {}\n'.format(leaf.name))
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
                leaf.add_child(name=leaf.taxon_names[i])
            leaf.name = ''
    return tree

def read_splist(species_list):
    with open(species_list, 'r') as f:
        splist = f.read().split('\n')
    splist = [ sp.replace('_', '') for sp in splist ]
    return splist

def get_lineages(splist):
    ncbi = ete3.NCBITaxa()
    lineages = dict()
    for sp in splist:
        lineages[sp] = get_lineage(sp, ncbi)
    return lineages

def get_taxid_counts(lineages):
    taxids = []
    for sp in lineages.keys():
        taxids += lineages[sp]
    uniq_taxids = list(set(taxids))
    counts = []
    for ut in uniq_taxids:
        counts.append(taxids.count(ut))
    taxid_counts = numpy.array([uniq_taxids,counts]).T
    count_order = numpy.argsort(taxid_counts[:,1])
    taxid_counts = taxid_counts[count_order,:]
    return taxid_counts

def get_mrca_taxid(multi_counts):
    ncbi = ete3.NCBITaxa()
    max_count = multi_counts[:,1].max()
    is_max_count = (multi_counts[:,1]==max_count)
    max_taxids = multi_counts[is_max_count,0]
    mrca_taxid = 1
    max_ancestor_num = 0
    for mt in max_taxids:
        ancestor_num = len(ncbi.get_lineage(mt))
        if (ancestor_num > max_ancestor_num):
            mrca_taxid = mt
            max_ancestor_num = ancestor_num
    return mrca_taxid

def get_max_ancestor_overlap_node(tree, ancestors):
    max_num_overlap = 0
    for node in tree.traverse():
        overlap_taxids = set(tree.ancestors).intersection(set(ancestors))
        num_overlap = len(overlap_taxids)
        if (num_overlap>max_num_overlap):
            max_num_overlap = num_overlap
            return_node = node
    return return_node

def populate_leaves(new_clade, taxid, lineages):
    for sp in lineages.keys():
        lineage = lineages[sp]
        if not taxid in lineage:
            continue
        new_leaf = ete3.PhyloNode(name=sp)
        new_leaf.ancestors = lineage
        new_clade.add_child(new_leaf)
    return new_clade

def add_new_clade(clades, new_clade):
    flag_append = True
    new_leaf_set = set(new_clade.get_leaf_names())
    delete_index = list()
    for j in range(len(clades)):
        clade = clades[j]
        clade_leaf_set = set(clade.get_leaf_names())
        intersection_set = clade_leaf_set.intersection(new_leaf_set)
        if (intersection_set==new_leaf_set):
            flag_append = False
            continue
        elif (intersection_set==clade_leaf_set):
            for nnode in new_clade.get_leaves():
                for int_name in intersection_set:
                    if nnode.name==int_name:
                        nnode.delete()
                        break
            new_clade.add_child(clade)
            delete_index.append(j)
    for di in sorted(delete_index)[::-1]:
        del clades[di]
    if flag_append:
        new_clade_leaves = new_clade.get_leaf_names()
        if (len(new_clade_leaves)!=len(set(new_clade_leaves))):
            txt = 'Redundant leaves in the new clade: {}'.format(', '.join(new_clade_leaves))
            raise Exception(txt)
        clades.append(new_clade)
    return clades

def get_polytomy_index(tree):
    num_polytomy = 0
    for node in tree.traverse():
        if node.is_leaf():
            continue
        num_polytomy += (len(node.get_children()) - 2)
    return num_polytomy

def taxid2tree(lineages, taxid_counts):
    ncbi = ete3.NCBITaxa()
    is_multiple = (taxid_counts[:,1]>1)
    multi_counts = taxid_counts[is_multiple,:]
    clades = list()
    for i in numpy.arange(multi_counts.shape[0]):
        taxid = multi_counts[i,0]
        count = multi_counts[i,1]
        ancestors = ncbi.get_lineage(taxid)
        new_clade = ete3.PhyloNode()
        new_clade.ancestors = ancestors
        new_clade = populate_leaves(new_clade, taxid, lineages)
        clades = add_new_clade(clades, new_clade)
    assert len(clades)==1, 'Failed to merge clades into a single tree.'
    tree = clades[0]
    return tree

def constrain_main(args):
    splist = read_splist(args.species_list)
    if (args.backbone=='ncbi'):
        lineages = get_lineages(splist)
        taxid_counts = get_taxid_counts(lineages)
        tree = taxid2tree(lineages, taxid_counts)
    else:
        if (args.backbone=='user'):
            tree = read_tree(args.infile, args.format)
        else:
            file_path = 'data_tree/'+args.backbone+'.nwk'
            nwk_string = pkg_resources.resource_string(__name__, file_path).decode("utf-8")
            tree = ete3.TreeNode(newick=nwk_string, format=0, quoted_node_names=True)
        tree = initialize_tree(tree)
        tree = match_taxa(tree, splist)
        tree = delete_nomatch_leaves(tree)
        tree = polytomize_one2many_matches(tree)
    tree = remove_singleton(tree, verbose=False, preserve_branch_length=False)
    for node in tree.traverse():
        node.name = node.name.replace(' ', '_')
    write_tree(tree, args, format=9)

import re
import requests
import sys
import time

from ete3 import NCBITaxa

from nwkit.util import *

def add_common_anc_constraint(tree, args):
    common_anc = tree.get_common_ancestor(args.left_species, args.right_species)
    if (args.lower_bound==args.upper_bound):
        constraint = '@' + args.lower_bound
    elif (args.lower_bound is not None) & (args.upper_bound is not None):
        constraint = 'B(' + ', '.join(
            [args.lower_bound, args.upper_bound, args.lower_tailProb, args.upper_tailProb]) + ')'
    elif (args.lower_bound is not None):
        constraint = 'L(' + ', '.join(
            [args.lower_bound, args.lower_offset, args.lower_scale, args.lower_tailProb]) + ')'
    elif (args.upper_bound is not None):
        constraint = 'U(' + ', '.join([args.upper_bound, args.upper_tailProb]) + ')'
    constraint = '\'' + constraint + '\''
    common_anc.name = constraint
    return tree

def check_leaf_taxid_availability(tree, ncbi):
    leaf_names = tree.get_leaf_names()
    leaf_names = [ ln.replace('_', ' ') for ln in leaf_names ]
    name2taxid = ncbi.get_name_translator(leaf_names)
    taxid_keys = name2taxid.keys()
    for ln in leaf_names:
        if not ln in taxid_keys:
            txt = 'NCBI Taxonomy ID was not found and thus not used as query to timetree.org: {}\n'
            sys.stderr.write(txt.format(ln))

def are_both_lineage_included(node, leaf_names):
    child_names_list = [ child.get_leaf_names() for child in node.get_children() ]
    flag_child_found = True
    for child_names in child_names_list:
        is_child_found = any([ child_name in leaf_names for child_name in child_names ])
        if not is_child_found:
            flag_child_found = False
            break
    return flag_child_found

def is_mrca_clade_root(node, timetree_result, ncbi):
    if 'missing_ids' not in timetree_result:
        return True
    missing_ids = timetree_result.replace('"', '').replace('\'', '').replace('\n', '')
    missing_ids = re.sub('.*missing_ids:\[', '', missing_ids)
    missing_ids = re.sub('\].*', '', missing_ids)
    if (len(missing_ids)==0):
        return True
    missing_ids = [ int(mid) for mid in missing_ids.split(',') ]
    taxid2name = ncbi.get_taxid_translator(missing_ids)
    missing_sci_names = list(taxid2name.values())
    missing_leaf_names = [ sn.replace(' ', '_') for sn in missing_sci_names ]
    return are_both_lineage_included(node=node, leaf_names=missing_leaf_names)

def are_two_lineage_rank_differentiated(node, taxids, ta_leaf_names):
    child_taxids = list()
    for child in node.get_children():
        child_leaf_names = child.get_leaf_names()
        ch_taxids = [ t for t,ln in zip(taxids,ta_leaf_names) if ln in child_leaf_names ]
        child_taxids.append(ch_taxids)
    assert len(child_taxids)==2, 'Non-bifurcation at the node containing: {}'.format(','.join(node.get_leaf_names()))
    if len(set(child_taxids[0]) - set(child_taxids[1]))==0:
        return False
    elif len(set(child_taxids[1]) - set(child_taxids[0]))==0:
        return False
    else:
        return True

def add_timetree_constraint(tree, args):
    endpoint_url = 'http://timetree.temple.edu/api'
    ncbi = NCBITaxa()
    check_leaf_taxid_availability(tree, ncbi)
    for node in tree.traverse():
        if node.is_leaf():
            continue
        node.name = ''
        leaf_names = node.get_leaf_names()
        sci_names = [ leaf_name.replace('_', ' ') for leaf_name in leaf_names ]
        name2taxid = ncbi.get_name_translator(sci_names)
        taxid_assigned_sci_names = list(name2taxid.keys())
        taxid_assigned_leaf_names = [ sci_name.replace(' ', '_') for sci_name in taxid_assigned_sci_names ]
        if not are_both_lineage_included(node=node, leaf_names=taxid_assigned_leaf_names):
            txt = "Skipping. Lack of NCBI Taxonomy information for the MRCA of {}\n"
            sys.stderr.write(txt.format(','.join(leaf_names)))
            continue
        species_taxids = [ t for ts in name2taxid.values() for t in ts ]
        lineage_taxids = list()
        for sp_taxid in species_taxids:
            lineages = ncbi.get_lineage(sp_taxid)
            ranks = ncbi.get_rank(lineages)
            lin_dict = dict()
            for taxid in ranks.keys():
                lin_dict[ranks[taxid]] = taxid
            lineage_taxids.append(lin_dict)
        search_ranks = ['species','genus','tribe','family','order','class','subphylum','phylum','kingdom','superkingdom']
        for search_rank in search_ranks:
            taxids = [ d[search_rank] for d in lineage_taxids if search_rank in list(d.keys()) ]
            ta_leaf_names = [ l for l,d in zip(taxid_assigned_leaf_names,lineage_taxids) if search_rank in list(d.keys()) ]
            if not are_both_lineage_included(node=node, leaf_names=ta_leaf_names):
                #txt = 'Rank annotation was not enough at {} for the node containing: {}\n'
                #sys.stderr.write(txt.format(search_rank, ','.join(leaf_names)))
                continue
            if not are_two_lineage_rank_differentiated(node=node, taxids=taxids, ta_leaf_names=ta_leaf_names):
                #txt = 'Taxonomic resolution was not enough at {} for the node containing: {}\n'
                #sys.stderr.write(txt.format(search_rank, ','.join(leaf_names)))
                continue
            if search_rank!='species':
                if not args.higher_rank_search:
                    continue
                sys.stderr.write('Searching higher taxonomic ranks to find MRCA at timetree.org: {}\n'.format(search_rank))
            request_url = '{}/mrca/id/{}'.format(endpoint_url, '+'.join([ str(t) for t in taxids ]))
            sys.stderr.write('Waiting for the REST API at timetree.org. ')
            start = time.time()
            response = requests.get(url=request_url)
            sys.stderr.write('Elapsed {:,} sec\n'.format(int(time.time() - start)))
            if response.status_code!=200:
                txt = 'Skipping. Failed to retrieve data from timetree.org for the MRCA of {}\n'
                sys.stderr.write(txt.format(','.join(leaf_names)))
            timetree_result = re.sub('.*;</script>', '', response.text)
            if "MRCA node not found" in timetree_result:
                txt = "Skipping. No MRCA found at timetree.org for the node containing: {}\n"
                sys.stderr.write(txt.format(','.join(leaf_names)))
                continue
            if "No study info found for node" in timetree_result:
                txt = "Skipping. No study info found at timetree.org for the node containing: {}\n"
                sys.stderr.write(txt.format(','.join(leaf_names)))
                continue
            if "No TimeTree study info available for this MRCA" in timetree_result:
                txt = "Skipping. No TimeTree study info available for this MRCA for the node containing: {}\n"
                sys.stderr.write(txt.format(','.join(leaf_names)))
                continue
            if not is_mrca_clade_root(node, timetree_result, ncbi):
                txt = "Skipping. Lack of timetree.org information for the MRCA of {}\n"
                sys.stderr.write(txt.format(','.join(leaf_names)))
                continue
            timetree_keys = re.sub('\r\n.*', '', re.sub('<br>.*', '', timetree_result)).split(',')
            timetree_values = re.sub('.*\r\n', '', re.sub('.*<br>', '', timetree_result)).split(',')
            timetree_dict = dict()
            for key,value in zip(timetree_keys, timetree_values):
                timetree_dict[key] = value
            if float(timetree_dict['precomputed_ci_high'])==0:
                txt = "Skipping. Upper age estimate at timetree.org is zero for the MRCA of {}\n"
                sys.stderr.write(txt.format(','.join(leaf_names)))
                continue
            if (args.timetree=='point'):
                constraint = '@' + timetree_dict['precomputed_age']
            elif (args.timetree=='ci'):
                constraint = 'B(' + ', '.join([timetree_dict['precomputed_ci_low'], timetree_dict['precomputed_ci_high'],
                                               args.lower_tailProb, args.upper_tailProb]) + ')'
            constraint = '\'' + constraint + '\''
            node.name = constraint
            break
    return tree

def remove_constraint_equal_upper(tree):
    removed_constraint_count = 0
    for node in tree.traverse(strategy='postorder'):
        if node.is_root():
            continue
        if node.is_leaf():
            continue
        if (node.name==''):
            continue
        if (node.up.name==''):
            continue
        if (node.name==node.up.name):
            node.name = ''
            removed_constraint_count += 1
    txt = 'Removed {:,} constraints that are equal to that of the parent node.\n'
    sys.stderr.write(txt.format(removed_constraint_count))
    return tree

def apply_min_clade_prop(tree, min_clade_prop):
    tree_size = len(tree.get_leaves())
    min_clade_size = min_clade_prop * tree_size
    removed_constraint_count = 0
    for node in tree.traverse():
        if node.is_root():
            continue
        if node.is_leaf():
            continue
        clade_size = len(node.get_leaves())
        if (clade_size < min_clade_size)&(node.name!=''):
            node.name = ''
            removed_constraint_count += 1
    txt = 'Removed {} constraints that are in clades smaller than {:,.1f}% ({:,} tips) of the tree size ({:,} tips).\n'
    sys.stderr.write(txt.format(removed_constraint_count, min_clade_prop*100, min_clade_size, tree_size))
    return tree

def mcmctree_main(args):
    tree = read_tree(args.infile, args.format, args.quoted_node_names)
    assert (len(list(tree.get_children()))==2), 'The input tree should be rooted.'
    for node in tree.traverse():
        if not node.is_leaf():
            if any([kw in node.name for kw in ['@', 'B(', 'L(', 'U(']]):
                node.name = '\'' + node.name + '\''
            else:
                node.name = ''
    if (args.timetree=='no'):
        tree = add_common_anc_constraint(tree, args)
    elif (args.timetree=='point'):
        tree = add_timetree_constraint(tree, args)
    elif (args.timetree=='ci'):
        tree = add_timetree_constraint(tree, args)
    tree = remove_constraint_equal_upper(tree)
    tree = apply_min_clade_prop(tree, min_clade_prop=args.min_clade_prop)
    nwk_text = tree.write(format=8, format_root_node=True, quoted_node_names=True)
    nwk_text = nwk_text.replace('NoName', '')
    nwk_text = nwk_text.replace('\"', '')
    if args.add_header:
        num_leaf = len(list(tree.get_leaf_names()))
        nwk_text = '{:} 1\n{}'.format(num_leaf, nwk_text)
    if args.outfile=='-':
        print(nwk_text)
    else:
        with open(args.outfile, mode='w') as f:
            f.write(nwk_text)

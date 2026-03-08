import re
import requests
import sys
import time

from nwkit.util import *

SEARCH_RANKS = [
    'species', 'genus', 'tribe', 'family', 'order',
    'class', 'subphylum', 'phylum', 'kingdom', 'superkingdom',
]

def add_common_anc_constraint(tree, args):
    if (args.left_species is None) or (args.right_species is None):
        raise ValueError("'--left_species' and '--right_species' are required when '--timetree no'.")
    if args.left_species == args.right_species:
        raise ValueError("'--left_species' and '--right_species' must be different species.")
    if (args.lower_bound is None) and (args.upper_bound is None):
        raise ValueError("Specify at least one of '--lower_bound' or '--upper_bound' when '--timetree no'.")
    leaf_name_set = set(tree.leaf_names())
    missing_species = [sp for sp in [args.left_species, args.right_species] if sp not in leaf_name_set]
    if missing_species:
        raise ValueError("Species not found in the input tree: {}".format(', '.join(missing_species)))
    common_anc = tree.common_ancestor([args.left_species, args.right_species])
    is_point_bound = False
    if (args.lower_bound is not None) and (args.upper_bound is not None):
        try:
            is_point_bound = (abs(float(args.lower_bound) - float(args.upper_bound)) < 10**-12)
        except ValueError:
            is_point_bound = (args.lower_bound == args.upper_bound)
    if (args.lower_bound is not None) and (args.upper_bound is not None) and is_point_bound:
        constraint = '@' + args.lower_bound
    elif (args.lower_bound is not None) and (args.upper_bound is not None):
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
    unnamed_leaves = [leaf for leaf in tree.leaves() if not leaf.name]
    if unnamed_leaves:
        raise ValueError('All leaves must have non-empty names when using "--timetree point/ci".')
    leaf_names = list(tree.leaf_names())
    leaf_names = [ ln.replace('_', ' ') for ln in leaf_names ]
    name2taxid = ncbi.get_name_translator(leaf_names)
    taxid_keys = set(name2taxid.keys())
    for ln in leaf_names:
        if not ln in taxid_keys:
            txt = 'NCBI Taxonomy ID was not found and thus not used as query to timetree.org: {}\n'
            sys.stderr.write(txt.format(ln))

def are_both_lineage_included(node, leaf_names, subtree_leaf_name_sets=None):
    if isinstance(leaf_names, set):
        leaf_name_set = leaf_names
    else:
        leaf_name_set = set(leaf_names)
    if len(leaf_name_set) == 0:
        return False
    for child in node.get_children():
        if subtree_leaf_name_sets is None:
            child_leaf_set = set(child.leaf_names())
        else:
            child_leaf_set = subtree_leaf_name_sets[child]
        if child_leaf_set.isdisjoint(leaf_name_set):
            return False
    return True

def is_mrca_clade_root(node, timetree_result, ncbi, subtree_leaf_name_sets=None):
    if 'missing_ids' not in timetree_result:
        return True
    missing_ids = timetree_result.replace('"', '').replace('\'', '').replace('\n', '')
    missing_ids = re.sub(r'.*missing_ids:\[', '', missing_ids)
    missing_ids = re.sub(r'\].*', '', missing_ids)
    if (len(missing_ids)==0):
        return True
    parsed_missing_ids = []
    for mid in missing_ids.split(','):
        mid = mid.strip()
        if mid == '':
            continue
        if mid.lower() == 'null':
            continue
        try:
            parsed_missing_ids.append(int(mid))
        except ValueError:
            continue
    if len(parsed_missing_ids) == 0:
        return True
    missing_ids = parsed_missing_ids
    taxid2name = ncbi.get_taxid_translator(missing_ids)
    missing_sci_names = list(taxid2name.values())
    missing_leaf_names = [ sn.replace(' ', '_') for sn in missing_sci_names ]
    return are_both_lineage_included(
        node=node,
        leaf_names=missing_leaf_names,
        subtree_leaf_name_sets=subtree_leaf_name_sets,
    )

def are_two_lineage_rank_differentiated(node, taxids, ta_leaf_names, subtree_leaf_name_sets=None):
    children = node.get_children()
    if len(children) != 2:
        return False
    if subtree_leaf_name_sets is None:
        child0_leaf_set = set(children[0].leaf_names())
        child1_leaf_set = set(children[1].leaf_names())
    else:
        child0_leaf_set = subtree_leaf_name_sets[children[0]]
        child1_leaf_set = subtree_leaf_name_sets[children[1]]
    child0_taxids = set()
    child1_taxids = set()
    for taxid, leaf_name in zip(taxids, ta_leaf_names):
        if leaf_name in child0_leaf_set:
            child0_taxids.add(taxid)
        elif leaf_name in child1_leaf_set:
            child1_taxids.add(taxid)
    if len(child0_taxids - child1_taxids)==0:
        return False
    elif len(child1_taxids - child0_taxids)==0:
        return False
    else:
        return True

def add_timetree_constraint(tree, args):
    endpoint_url = 'http://timetree.org/api'
    search_ranks = SEARCH_RANKS if args.higher_rank_search else SEARCH_RANKS[:1]
    unnamed_leaves = [leaf for leaf in tree.leaves() if not leaf.name]
    if unnamed_leaves:
        raise ValueError('All leaves must have non-empty names when using "--timetree point/ci".')
    ncbi = get_ete_ncbitaxa(args=args)
    try:
        check_leaf_taxid_availability(tree, ncbi)
        subtree_leaf_name_sets = get_subtree_leaf_name_sets(tree)
        leaf_name_to_sci_name = {leaf.name: leaf.name.replace('_', ' ') for leaf in tree.leaves()}
        taxid_lineage_rank_dict_cache = dict()
        for node in tree.traverse():
            if node.is_leaf:
                continue
            node.name = 'NoName'
            leaf_names = sorted(subtree_leaf_name_sets[node])
            sci_names = [leaf_name_to_sci_name[leaf_name] for leaf_name in leaf_names]
            name2taxid = ncbi.get_name_translator(sci_names)
            taxid_assigned_sci_names = list(name2taxid.keys())
            taxid_assigned_leaf_names = [ sci_name.replace(' ', '_') for sci_name in taxid_assigned_sci_names ]
            if not are_both_lineage_included(
                node=node,
                leaf_names=taxid_assigned_leaf_names,
                subtree_leaf_name_sets=subtree_leaf_name_sets,
            ):
                txt = "Skipping. Lack of NCBI Taxonomy information for the MRCA of {}\n"
                sys.stderr.write(txt.format(','.join(leaf_names)))
                continue
            species_taxids = [ t for ts in name2taxid.values() for t in ts ]
            lineage_taxids = list()
            for sp_taxid in species_taxids:
                if sp_taxid not in taxid_lineage_rank_dict_cache:
                    lineages = ncbi.get_lineage(sp_taxid)
                    ranks = ncbi.get_rank(lineages)
                    lin_dict = dict()
                    for taxid in ranks.keys():
                        lin_dict[ranks[taxid]] = taxid
                    taxid_lineage_rank_dict_cache[sp_taxid] = lin_dict
                lineage_taxids.append(taxid_lineage_rank_dict_cache[sp_taxid])
            leaf_rank_pairs = list(zip(taxid_assigned_leaf_names, lineage_taxids))
            for search_rank in search_ranks:
                taxids = [d[search_rank] for d in lineage_taxids if search_rank in d]
                ta_leaf_names = [l for l, d in leaf_rank_pairs if search_rank in d]
                ta_leaf_name_set = set(ta_leaf_names)
                if not are_both_lineage_included(
                    node=node,
                    leaf_names=ta_leaf_name_set,
                    subtree_leaf_name_sets=subtree_leaf_name_sets,
                ):
                    #txt = 'Rank annotation was not enough at {} for the node containing: {}\n'
                    #sys.stderr.write(txt.format(search_rank, ','.join(leaf_names)))
                    continue
                if not are_two_lineage_rank_differentiated(
                    node=node,
                    taxids=taxids,
                    ta_leaf_names=ta_leaf_names,
                    subtree_leaf_name_sets=subtree_leaf_name_sets,
                ):
                    #txt = 'Taxonomic resolution was not enough at {} for the node containing: {}\n'
                    #sys.stderr.write(txt.format(search_rank, ','.join(leaf_names)))
                    continue
                if search_rank!='species':
                    sys.stderr.write('Searching higher taxonomic ranks to find MRCA at timetree.org: {}\n'.format(search_rank))
                request_url = '{}/mrca/id/{}'.format(endpoint_url, '+'.join([ str(t) for t in taxids ]))
                sys.stderr.write('Waiting for the REST API at timetree.org. ')
                start = time.time()
                try:
                    response = requests.get(url=request_url, timeout=30)
                except requests.RequestException:
                    txt = 'Skipping. Failed to retrieve data from timetree.org for the MRCA of {}\n'
                    sys.stderr.write(txt.format(','.join(leaf_names)))
                    continue
                sys.stderr.write('Elapsed {:,} sec\n'.format(int(time.time() - start)))
                if response.status_code!=200:
                    txt = 'Skipping. Failed to retrieve data from timetree.org for the MRCA of {}\n'
                    sys.stderr.write(txt.format(','.join(leaf_names)))
                    continue
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
                if not is_mrca_clade_root(
                    node,
                    timetree_result,
                    ncbi,
                    subtree_leaf_name_sets=subtree_leaf_name_sets,
                ):
                    txt = "Skipping. Lack of timetree.org information for the MRCA of {}\n"
                    sys.stderr.write(txt.format(','.join(leaf_names)))
                    continue
                timetree_keys = re.sub('\r\n.*', '', re.sub('<br>.*', '', timetree_result)).split(',')
                timetree_values = re.sub('.*\r\n', '', re.sub('.*<br>', '', timetree_result)).split(',')
                timetree_dict = dict()
                for key,value in zip(timetree_keys, timetree_values):
                    timetree_dict[key] = value
                required_keys = ['precomputed_age', 'precomputed_ci_low', 'precomputed_ci_high']
                if any(key not in timetree_dict for key in required_keys):
                    txt = "Skipping. Unexpected response format from timetree.org for the node containing: {}\n"
                    sys.stderr.write(txt.format(','.join(leaf_names)))
                    continue
                try:
                    ci_high = float(timetree_dict['precomputed_ci_high'])
                except ValueError:
                    txt = "Skipping. Non-numeric age estimate from timetree.org for the node containing: {}\n"
                    sys.stderr.write(txt.format(','.join(leaf_names)))
                    continue
                if ci_high==0:
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
    finally:
        db = getattr(ncbi, 'db', None)
        if db is not None:
            db.close()

def remove_constraint_equal_upper(tree):
    removed_constraint_count = 0
    for node in tree.traverse(strategy='postorder'):
        if node.is_root:
            continue
        if node.is_leaf:
            continue
        if not node.name:
            continue
        if not node.up.name:
            continue
        if (node.name==node.up.name):
            node.name = 'NoName'
            removed_constraint_count += 1
    txt = 'Removed {:,} constraints that are equal to that of the parent node.\n'
    sys.stderr.write(txt.format(removed_constraint_count))
    return tree

def apply_min_clade_prop(tree, min_clade_prop):
    tree_size = len(list(tree.leaves()))
    min_clade_size = min_clade_prop * tree_size
    removed_constraint_count = 0
    subtree_leaf_name_sets = get_subtree_leaf_name_sets(tree)
    for node in tree.traverse():
        if node.is_root:
            continue
        if node.is_leaf:
            continue
        clade_size = len(subtree_leaf_name_sets[node])
        if (clade_size < min_clade_size) and node.name:
            node.name = 'NoName'
            removed_constraint_count += 1
    txt = 'Removed {} constraints that are in clades smaller than {:,.1f}% ({:,} tips) of the tree size ({:,} tips).\n'
    sys.stderr.write(txt.format(removed_constraint_count, min_clade_prop*100, min_clade_size, tree_size))
    return tree

def mcmctree_main(args):
    tree = read_tree(args.infile, args.format, args.quoted_node_names)
    if len(tree.get_children()) != 2:
        raise ValueError('The input tree should be rooted.')
    if (args.min_clade_prop < 0) or (args.min_clade_prop > 1):
        raise ValueError("'--min_clade_prop' must be between 0 and 1.")
    for node in tree.traverse():
        if not node.is_leaf:
            if any([kw in (node.name or '') for kw in ['@', 'B(', 'L(', 'U(']]):
                node.name = '\'' + node.name + '\''
            else:
                node.name = 'NoName'
    if (args.timetree=='no'):
        if (args.left_species is None) or (args.right_species is None):
            raise ValueError("'--left_species' and '--right_species' are required when '--timetree no'.")
        if (args.lower_bound is None) and (args.upper_bound is None):
            raise ValueError("Specify at least one of '--lower_bound' or '--upper_bound' when '--timetree no'.")
        tree = add_common_anc_constraint(tree, args)
    elif (args.timetree=='point'):
        tree = add_timetree_constraint(tree, args)
    elif (args.timetree=='ci'):
        tree = add_timetree_constraint(tree, args)
    else:
        raise ValueError("Unknown '--timetree' mode: {}. Choose from no/point/ci.".format(args.timetree))
    tree = remove_constraint_equal_upper(tree)
    tree = apply_min_clade_prop(tree, min_clade_prop=args.min_clade_prop)
    # Use parser=1 and post-process for MCMCtree format
    nwk_text = tree.write(parser=1, format_root_node=True)
    nwk_text = re.sub(r':[\d.eE+-]+', '', nwk_text)  # Remove branch lengths
    nwk_text = nwk_text.replace("'''", "'")  # Clean up triple quotes from ete4
    nwk_text = nwk_text.replace('NoName', '')
    nwk_text = nwk_text.replace('"', '')
    if args.add_header:
        num_leaf = len(list(tree.leaves()))
        nwk_text = '{:} 1\n{}'.format(num_leaf, nwk_text)
    if args.outfile=='-':
        print(nwk_text)
    else:
        with open(args.outfile, mode='w') as f:
            f.write(nwk_text)

import re
import requests
import sys
import time
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed

from nwkit.util import (
    get_ete_ncbitaxa,
    get_species_group_records,
    get_subtree_leaf_name_sets,
    get_subtree_sci_name_sets,
    read_tree,
    warn_cleanup_failure,
)

SEARCH_RANKS = [
    'species', 'genus', 'tribe', 'family', 'order',
    'class', 'subphylum', 'phylum', 'kingdom', 'superkingdom',
]


def _validate_threads(threads):
    try:
        threads = int(threads)
    except (TypeError, ValueError) as exc:
        raise ValueError("'--threads' must be an integer.") from exc
    if threads <= 0:
        raise ValueError("'--threads' must be positive.")
    return threads


def _fetch_timetree_url(request_url):
    start = time.time()
    try:
        response = requests.get(url=request_url, timeout=30)
    except requests.RequestException as exc:
        return {
            'url': request_url,
            'status_code': None,
            'text': '',
            'error': exc,
            'elapsed': int(time.time() - start),
        }
    return {
        'url': request_url,
        'status_code': response.status_code,
        'text': response.text,
        'error': None,
        'elapsed': int(time.time() - start),
    }


def _fetch_timetree_url_cached(request_url, response_cache):
    if request_url not in response_cache:
        sys.stderr.write('Waiting for the REST API at timetree.org. ')
        response_cache[request_url] = _fetch_timetree_url(request_url)
        sys.stderr.write('Elapsed {:,} sec\n'.format(response_cache[request_url]['elapsed']))
    return response_cache[request_url]


def _fetch_timetree_urls_parallel(request_urls, threads, response_cache):
    missing_urls = [
        request_url
        for request_url in dict.fromkeys(request_urls)
        if request_url not in response_cache
    ]
    if len(missing_urls) == 0:
        return {request_url: response_cache[request_url] for request_url in request_urls}
    max_workers = min(threads, len(missing_urls))
    if max_workers <= 1:
        for request_url in missing_urls:
            _fetch_timetree_url_cached(request_url, response_cache=response_cache)
        return {request_url: response_cache[request_url] for request_url in request_urls}
    sys.stderr.write(
        'Waiting for the REST API at timetree.org with {:,} parallel request(s).\n'.format(
            len(missing_urls)
        )
    )
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        future_to_url = {
            executor.submit(_fetch_timetree_url, request_url): request_url
            for request_url in missing_urls
        }
        for future in as_completed(future_to_url):
            request_url = future_to_url[future]
            response_cache[request_url] = future.result()
            sys.stderr.write(
                'Elapsed {:,} sec for TimeTree request: {}\n'.format(
                    response_cache[request_url]['elapsed'],
                    request_url,
                )
            )
    return {request_url: response_cache[request_url] for request_url in request_urls}


def _build_timetree_rank_attempt(context, search_rank, endpoint_url, subtree_species_label_sets):
    lineage_taxids = context['lineage_taxids']
    leaf_rank_pairs = context['leaf_rank_pairs']
    taxids = [d[search_rank] for d in lineage_taxids if search_rank in d]
    ta_leaf_names = [leaf_name for leaf_name, rank_by_name in leaf_rank_pairs if search_rank in rank_by_name]
    ta_leaf_name_set = set(ta_leaf_names)
    if not are_both_lineage_included(
        node=context['node'],
        leaf_names=ta_leaf_name_set,
        subtree_leaf_name_sets=subtree_species_label_sets,
    ):
        return None
    if not are_two_lineage_rank_differentiated(
        node=context['node'],
        taxids=taxids,
        ta_leaf_names=ta_leaf_names,
        subtree_leaf_name_sets=subtree_species_label_sets,
    ):
        return None
    if search_rank != 'species':
        sys.stderr.write('Searching higher taxonomic ranks to find MRCA at timetree.org: {}\n'.format(search_rank))
    request_url = '{}/mrca/id/{}'.format(endpoint_url, '+'.join([str(t) for t in taxids]))
    taxid_to_species_labels = defaultdict(set)
    for taxid, species_label in zip(taxids, ta_leaf_names):
        taxid_to_species_labels[int(taxid)].add(species_label)
    return {
        'context': context,
        'request_url': request_url,
        'taxid_to_species_labels': taxid_to_species_labels,
    }


def _constraint_from_timetree_response(attempt, response_record, ncbi, subtree_species_label_sets, args):
    context = attempt['context']
    node = context['node']
    leaf_names = context['leaf_names']
    if response_record['error'] is not None:
        txt = 'Skipping. Failed to retrieve data from timetree.org for the MRCA of {}\n'
        sys.stderr.write(txt.format(','.join(leaf_names)))
        return None
    if response_record['status_code'] != 200:
        txt = 'Skipping. Failed to retrieve data from timetree.org for the MRCA of {}\n'
        sys.stderr.write(txt.format(','.join(leaf_names)))
        return None
    timetree_result = re.sub('.*;</script>', '', response_record['text'])
    if "MRCA node not found" in timetree_result:
        txt = "Skipping. No MRCA found at timetree.org for the node containing: {}\n"
        sys.stderr.write(txt.format(','.join(leaf_names)))
        return None
    if "No study info found for node" in timetree_result:
        txt = "Skipping. No study info found at timetree.org for the node containing: {}\n"
        sys.stderr.write(txt.format(','.join(leaf_names)))
        return None
    if "No TimeTree study info available for this MRCA" in timetree_result:
        txt = "Skipping. No TimeTree study info available for this MRCA for the node containing: {}\n"
        sys.stderr.write(txt.format(','.join(leaf_names)))
        return None
    if not is_mrca_clade_root(
        node,
        timetree_result,
        ncbi,
        subtree_leaf_name_sets=subtree_species_label_sets,
        taxid_to_species_labels=attempt['taxid_to_species_labels'],
    ):
        txt = "Skipping. Lack of timetree.org information for the MRCA of {}\n"
        sys.stderr.write(txt.format(','.join(leaf_names)))
        return None
    timetree_keys = re.sub('\r\n.*', '', re.sub('<br>.*', '', timetree_result)).split(',')
    timetree_values = re.sub('.*\r\n', '', re.sub('.*<br>', '', timetree_result)).split(',')
    timetree_dict = dict()
    for key, value in zip(timetree_keys, timetree_values):
        timetree_dict[key] = value
    required_keys = ['precomputed_age', 'precomputed_ci_low', 'precomputed_ci_high']
    if any(key not in timetree_dict for key in required_keys):
        txt = "Skipping. Unexpected response format from timetree.org for the node containing: {}\n"
        sys.stderr.write(txt.format(','.join(leaf_names)))
        return None
    try:
        ci_high = float(timetree_dict['precomputed_ci_high'])
    except ValueError:
        txt = "Skipping. Non-numeric age estimate from timetree.org for the node containing: {}\n"
        sys.stderr.write(txt.format(','.join(leaf_names)))
        return None
    if ci_high == 0:
        txt = "Skipping. Upper age estimate at timetree.org is zero for the MRCA of {}\n"
        sys.stderr.write(txt.format(','.join(leaf_names)))
        return None
    if args.timetree == 'point':
        constraint = '@' + timetree_dict['precomputed_age']
    elif args.timetree == 'ci':
        constraint = 'B(' + ', '.join([
            timetree_dict['precomputed_ci_low'],
            timetree_dict['precomputed_ci_high'],
            _tail_probability(args, 'lower'),
            _tail_probability(args, 'upper'),
        ]) + ')'
    else:
        return None
    return '\'' + constraint + '\''


def _tail_probability(args, side):
    value = getattr(args, '{}_tail_prob'.format(side), None)
    if value is None:
        value = getattr(args, '{}_tailProb'.format(side), None)
    return '0.025' if value is None else value

def add_common_anc_constraint(tree, args):
    if (args.left_species is None) or (args.right_species is None):
        raise ValueError("'--left-species' and '--right-species' are required when '--timetree no'.")
    if args.left_species == args.right_species:
        raise ValueError("'--left-species' and '--right-species' must be different species.")
    if (args.lower_bound is None) and (args.upper_bound is None):
        raise ValueError("Specify at least one of '--lower-bound' or '--upper-bound' when '--timetree no'.")
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
        constraint = 'B(' + ', '.join([
            args.lower_bound,
            args.upper_bound,
            _tail_probability(args, 'lower'),
            _tail_probability(args, 'upper'),
        ]) + ')'
    elif (args.lower_bound is not None):
        constraint = 'L(' + ', '.join([
            args.lower_bound,
            args.lower_offset,
            args.lower_scale,
            _tail_probability(args, 'lower'),
        ]) + ')'
    elif (args.upper_bound is not None):
        constraint = 'U(' + ', '.join([
            args.upper_bound,
            _tail_probability(args, 'upper'),
        ]) + ')'
    constraint = '\'' + constraint + '\''
    common_anc.name = constraint
    return tree

def check_leaf_taxid_availability(species_names, ncbi):
    if isinstance(species_names, dict):
        species_label_to_taxonomy_query = species_names
    else:
        species_label_to_taxonomy_query = {
            species_name: str(species_name).replace('_', ' ')
            for species_name in species_names
        }
    query_names = sorted(set(species_label_to_taxonomy_query.values()))
    name2taxid = ncbi.get_name_translator(query_names)
    taxid_keys = set(name2taxid.keys())
    for ln in query_names:
        if ln not in taxid_keys:
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

def is_mrca_clade_root(node, timetree_result, ncbi, subtree_leaf_name_sets=None, taxid_to_species_labels=None):
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
    missing_leaf_names = list()
    if taxid_to_species_labels is not None:
        for missing_id in missing_ids:
            missing_leaf_names.extend(sorted(taxid_to_species_labels.get(int(missing_id), [])))
    else:
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
    endpoint_url = 'https://timetree.org/api'
    search_ranks = SEARCH_RANKS if args.higher_rank_search else SEARCH_RANKS[:1]
    threads = _validate_threads(getattr(args, 'threads', 1))
    unnamed_leaves = [leaf for leaf in tree.leaves() if not leaf.name]
    if unnamed_leaves:
        raise ValueError('All leaves must have non-empty names when using "--timetree point/ci".')
    leaf_name_to_species_label, species_to_leaf_names, species_label_to_taxonomy_query = get_species_group_records(
        tree,
        option_name='--infile',
        context=' for "--timetree point/ci"',
        args=args,
    )
    ncbi = get_ete_ncbitaxa(args=args)
    try:
        check_leaf_taxid_availability(species_label_to_taxonomy_query, ncbi)
        for leaf in tree.leaves():
            species_label = leaf_name_to_species_label[leaf.name]
            leaf.add_props(
                species_label=species_label,
                taxonomy_query=species_label_to_taxonomy_query[species_label],
                sci_name=species_label,
            )
        subtree_leaf_name_sets = get_subtree_leaf_name_sets(tree)
        subtree_species_label_sets = get_subtree_sci_name_sets(tree)
        name_to_taxid_cache = dict()
        timetree_response_cache = dict()
        taxid_lineage_rank_dict_cache = dict()
        node_contexts = list()
        for node in tree.traverse():
            if node.is_leaf:
                continue
            node.name = 'NoName'
            leaf_names = sorted(subtree_leaf_name_sets[node])
            species_labels = sorted(subtree_species_label_sets[node])
            query_name_to_species_labels = defaultdict(list)
            for species_label in species_labels:
                taxonomy_query = species_label_to_taxonomy_query.get(species_label)
                if taxonomy_query in ['', None]:
                    continue
                query_name_to_species_labels[taxonomy_query].append(species_label)
            query_names_key = tuple(sorted(query_name_to_species_labels.keys()))
            if query_names_key not in name_to_taxid_cache:
                name_to_taxid_cache[query_names_key] = ncbi.get_name_translator(list(query_names_key))
            name2taxid = name_to_taxid_cache[query_names_key]
            taxid_assigned_species_labels = list()
            for query_name, query_species_labels in query_name_to_species_labels.items():
                if query_name not in name2taxid:
                    continue
                taxid_assigned_species_labels.extend(query_species_labels)
            if not are_both_lineage_included(
                node=node,
                leaf_names=set(taxid_assigned_species_labels),
                subtree_leaf_name_sets=subtree_species_label_sets,
            ):
                txt = "Skipping. Lack of NCBI Taxonomy information for the MRCA of {}\n"
                sys.stderr.write(txt.format(','.join(leaf_names)))
                continue
            species_taxid_pairs = list()
            for query_name, taxids in name2taxid.items():
                for species_label in query_name_to_species_labels.get(query_name, []):
                    for taxid in taxids:
                        species_taxid_pairs.append((species_label, int(taxid)))
            lineage_taxids = list()
            for _, sp_taxid in species_taxid_pairs:
                if sp_taxid not in taxid_lineage_rank_dict_cache:
                    lineages = ncbi.get_lineage(sp_taxid)
                    ranks = ncbi.get_rank(lineages)
                    lin_dict = dict()
                    for taxid in ranks.keys():
                        lin_dict[ranks[taxid]] = taxid
                    taxid_lineage_rank_dict_cache[sp_taxid] = lin_dict
                lineage_taxids.append(taxid_lineage_rank_dict_cache[sp_taxid])
            leaf_rank_pairs = list(zip([species_label for species_label, _ in species_taxid_pairs], lineage_taxids))
            node_contexts.append(
                {
                    'node': node,
                    'leaf_names': leaf_names,
                    'lineage_taxids': lineage_taxids,
                    'leaf_rank_pairs': leaf_rank_pairs,
                }
            )
        pending_contexts = list(node_contexts)
        for search_rank in search_ranks:
            attempts = [
                attempt
                for attempt in (
                    _build_timetree_rank_attempt(
                        context=context,
                        search_rank=search_rank,
                        endpoint_url=endpoint_url,
                        subtree_species_label_sets=subtree_species_label_sets,
                    )
                    for context in pending_contexts
                )
                if attempt is not None
            ]
            if len(attempts) == 0:
                continue
            response_by_url = dict()
            if threads > 1:
                response_by_url = _fetch_timetree_urls_parallel(
                    request_urls=[attempt['request_url'] for attempt in attempts],
                    threads=threads,
                    response_cache=timetree_response_cache,
                )
            attempts_by_context_id = {
                id(attempt['context']): attempt
                for attempt in attempts
            }
            next_pending_contexts = list()
            for context in pending_contexts:
                attempt = attempts_by_context_id.get(id(context))
                if attempt is None:
                    next_pending_contexts.append(context)
                    continue
                request_url = attempt['request_url']
                response_record = response_by_url.get(request_url)
                if response_record is None:
                    response_record = _fetch_timetree_url_cached(
                        request_url=request_url,
                        response_cache=timetree_response_cache,
                    )
                constraint = _constraint_from_timetree_response(
                    attempt=attempt,
                    response_record=response_record,
                    ncbi=ncbi,
                    subtree_species_label_sets=subtree_species_label_sets,
                    args=args,
                )
                if constraint is None:
                    next_pending_contexts.append(context)
                    continue
                context['node'].name = constraint
            pending_contexts = next_pending_contexts
            if len(pending_contexts) == 0:
                break
        return tree
    finally:
        db = getattr(ncbi, 'db', None)
        if db is not None:
            try:
                db.close()
            except Exception as exc:
                warn_cleanup_failure('NCBI taxonomy database handle', exc)

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
        raise ValueError("'--min-clade-prop' must be between 0 and 1.")
    for node in tree.traverse():
        if not node.is_leaf:
            if any([kw in (node.name or '') for kw in ['@', 'B(', 'L(', 'U(']]):
                node.name = '\'' + node.name + '\''
            else:
                node.name = 'NoName'
    if (args.timetree=='no'):
        if (args.left_species is None) or (args.right_species is None):
            raise ValueError("'--left-species' and '--right-species' are required when '--timetree no'.")
        if (args.lower_bound is None) and (args.upper_bound is None):
            raise ValueError("Specify at least one of '--lower-bound' or '--upper-bound' when '--timetree no'.")
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

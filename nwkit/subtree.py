from nwkit.util import *

def subtree_main(args):
    tree = read_tree(args.infile, args.format, args.quoted_node_names)
    leaf_nodes = list(tree.leaves())
    leaf_name_set = set(leaf.name for leaf in leaf_nodes)
    leaf_by_name = {leaf.name: leaf for leaf in leaf_nodes}
    if args.leaves is not None:
        seed_leaf_names = args.leaves.split(',')
    else:
        seed_leaf_names = [args.left_leaf, args.right_leaf]
        seed_leaf_names = [ l for l in seed_leaf_names if l is not None ]
    for l in seed_leaf_names:
        if l not in leaf_name_set:
            raise Exception('Specified leaf not found in the tree: {}'.format(l))
    seed_leaf_nodes = [leaf_by_name[l] for l in seed_leaf_names]
    if len(seed_leaf_nodes)==1:
        mrca = seed_leaf_nodes[0]
    else:
        mrca = tree.common_ancestor(seed_leaf_nodes)
    if args.orthogroup:
        tree = annotate_scientific_names(tree)
        subtree_sci_name_sets = get_subtree_sci_name_sets(tree)
        tree = annotate_duplication_confidence_scores(tree, subtree_sci_name_sets=subtree_sci_name_sets)
        num_all_species = len(subtree_sci_name_sets[tree])
        sys.stderr.write('Number of species in input tree: {:,}\n'.format(num_all_species))
        num_species_last_updated = len(subtree_sci_name_sets[mrca])
        node_species_last_updated = mrca
        for current_node in mrca.ancestors():
            if current_node.is_root:
                break
            if num_species_last_updated == num_all_species:
                break
            current_num_species = len(subtree_sci_name_sets[current_node])
            is_species_num_increased = (current_num_species > num_species_last_updated)
            is_sufficiently_small_dup_conf_score = (current_node.props.get('dup_conf_score') <= args.dup_conf_score_threshold)
            if is_species_num_increased and is_sufficiently_small_dup_conf_score:
                node_species_last_updated = current_node
                num_species_last_updated = current_num_species
            #txt = 'i: {}, dup_conf_score: {}, num_species: {}'
            #print(txt.format(i, current_node.props.get('dup_conf_score'), current_num_species))
        subtree = node_species_last_updated
        num_species_subtree = len(subtree_sci_name_sets[subtree])
        sys.stderr.write('Number of species in output orthogroup tree: {:,}\n'.format(num_species_subtree))
    else:
        subtree = mrca
    write_tree(subtree, args, format=args.outformat)

from nwkit.util import *

def subtree_main(args):
    tree = read_tree(args.infile, args.format, args.quoted_node_names)
    leaf_names = tree.get_leaf_names()
    if args.leaves is not None:
        seed_leaf_names = args.leaves.split(',')
    else:
        seed_leaf_names = [args.left_leaf, args.right_leaf]
        seed_leaf_names = [ l for l in seed_leaf_names if l is not None ]
    for l in seed_leaf_names:
        if not l in leaf_names:
            raise Exception('Specified leaf not found in the tree: {}'.format(l))
    seed_leaf_nodes = [ l for l in tree.iter_leaves() if l.name in seed_leaf_names ]
    if len(seed_leaf_nodes)==1:
        mrca = seed_leaf_nodes[0]
    else:
        mrca = tree.get_common_ancestor(seed_leaf_nodes[0], *seed_leaf_nodes[1:len(seed_leaf_nodes)])
    if args.orthogroup:
        tree = annotate_scientific_names(tree)
        tree = annotate_duplication_confidence_scores(tree)
        all_sci_names = list(set([ l.sci_name for l in tree.iter_leaves() ]))
        num_all_species = len(all_sci_names)
        sys.stderr.write('Number of species in input tree: {:,}\n'.format(num_all_species))
        mrca_sci_names = list(set([ l.sci_name for l in mrca.iter_leaves() ]))
        num_species_last_updated = len(mrca_sci_names)
        node_species_last_updated = mrca
        for i,current_node in enumerate(mrca.iter_ancestors()):
            if current_node.is_root():
                break
            if num_species_last_updated == num_all_species:
                break
            current_sci_names = list(set([ l.sci_name for l in current_node.iter_leaves() ]))
            current_num_species = len(current_sci_names)
            is_species_num_increased = (current_num_species > num_species_last_updated)
            is_sufficiently_small_dup_conf_score = (current_node.dup_conf_score <= args.dup_conf_score_threshold)
            if is_species_num_increased & is_sufficiently_small_dup_conf_score:
                node_species_last_updated = current_node
                num_species_last_updated = current_num_species
            #txt = 'i: {}, dup_conf_score: {}, num_species: {}'
            #print(txt.format(i, current_node.dup_conf_score, current_num_species))
        subtree = node_species_last_updated
        num_species_subtree = len(list(set([ l.sci_name for l in subtree.iter_leaves() ])))
        sys.stderr.write('Number of species in output orthogroup tree: {:,}\n'.format(num_species_subtree))
    else:
        subtree = mrca
    write_tree(subtree, args, format=args.outformat)
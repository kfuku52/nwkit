from nwkit.util import *

def dist_main(args):
    if args.dist not in ['RF']:
        raise ValueError('`--dist {}` is not supported'.format(args.dist))
    if args.infile2 in ['', None]:
        raise ValueError("'--infile2' is required for 'dist'.")
    tree1 = read_tree(args.infile, args.format, args.quoted_node_names)
    tree2 = read_tree(args.infile2, args.format2, args.quoted_node_names)
    if (not is_rooted(tree1)) or (not is_rooted(tree2)):
        raise ValueError('RF distance currently requires rooted trees. Root the input trees first.')
    validate_unique_named_leaves(tree1, option_name='--infile', context=" for 'dist'")
    validate_unique_named_leaves(tree2, option_name='--infile2', context=" for 'dist'")
    leaf_names1 = list(tree1.leaf_names())
    leaf_names2 = list(tree2.leaf_names())
    if set(leaf_names1) != set(leaf_names2):
        raise ValueError('Leaf name(s) did not match.')
    if args.dist=='RF':
        out = tree1.robinson_foulds(t2=tree2, unrooted_trees=False)
        rf,rf_max,common_attrs,edges_t1,edges_t2,discarded_edges_t1,discarded_edges_t2 = out
        if args.outfile=='-':
            print('rf_dist\tmax_rf_dist')
            print("{}\t{}".format(rf, rf_max))
        else:
            with open(args.outfile, mode='w') as f:
                f.write('rf_dist\tmax_rf_dist\n')
                f.write("{}\t{}\n".format(rf, rf_max))

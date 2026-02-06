from nwkit.util import *

def dist_main(args):
    assert args.dist in ['RF',], '`--dist {}` is not supported'.format(args.dist)
    tree1 = read_tree(args.infile, args.format, args.quoted_node_names)
    tree2 = read_tree(args.infile2, args.format2, args.quoted_node_names)
    assert set(tree1.leaf_names())==set(tree2.leaf_names()), 'Leaf name(s) did not match.'
    if args.dist=='RF':
        out = tree1.robinson_foulds(t2=tree2, unrooted_trees=False)
        rf,rf_max,common_attrs,edges_t1,edges_t2,discarded_edges_t1,discarded_edges_t2 = out
        if args.outfile=='-':
            print('rf_dist\tmax_rf_dist')
            print("{}\t{}".format(rf, rf_max))
        else:
            with open(args.outfile, mode='w') as f:
                f.write('rf_dist\tmax_rf_dist')
                f.write("{}\t{}".format(rf, rf_max))

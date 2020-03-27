from nwkit.util import *


def dist_main(args):
    assert args.dist in ['RF',], '`--dist {}` is not supported'.format(args.dist)
    tree1 = read_tree(args.infile, args.format)
    tree2 = read_tree(args.infile2, args.format2)
    assert set(tree1.get_leaf_names())==set(tree2.get_leaf_names()), 'Leaf name(s) did not match.'
    if args.dist=='RF':
        out = tree1.robinson_foulds(t2=tree2, unrooted_trees=False)
        rf,rf_max,common_attrs,edges_t1,edges_t2,discarded_edges_t1,discarded_edges_t2 = out
        print('rf_dist\tmax_rf_dist')
        print("{}\t{}".format(rf, rf_max))




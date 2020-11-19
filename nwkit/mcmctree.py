from nwkit.util import *

def mcmctree_main(args):
    tree = read_tree(args.infile, args.format)
    assert (len(list(tree.get_children()))==2), 'The input tree should be rooted.'
    for node in tree.traverse():
        if not node.is_leaf():
            if any([kw in node.name for kw in ['@', 'B(', 'L(', 'U(']]):
                node.name = '\'' + node.name + '\''
            else:
                node.name = ''
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

    nwk_text = tree.write(format=8, format_root_node=True, quoted_node_names=True)
    nwk_text = nwk_text.replace('NoName', '')
    nwk_text = nwk_text.replace('\"', '')
    if args.add_header:
        num_leaf = len(list(tree.get_leaf_names()))
        print(num_leaf, '1')

    if args.outfile=='-':
        print(nwk_text)
    else:
        with open(args.outfile, mode='w') as f:
            f.write(nwk_text)

import sys
from ete3 import TreeNode

def read_tree(args):
    if args.infile=='-':
        nwk_string = sys.stdin.readlines()[0]
        tree = TreeNode(newick=nwk_string, format=args.format)
    else:
        tree = TreeNode(newick=args.infile, format=args.format)
    return tree

def write_tree(tree, args):
    tree_str = tree.write(format=args.format, format_root_node=True)
    tree_str = tree_str.replace(':123456789','').replace(':1.23457e+08','')
    tree_str = tree_str.replace('123456789','').replace('1.23457e+08','')
    if args.outfile=='-':
        print(tree_str)
    else:
        with open(args.outfile, mode='w') as f:
            f.write(tree_str)

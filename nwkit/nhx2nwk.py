from nwkit.util import *

def nhx2nwk_main(args):
    tree = read_tree(args.infile, args.format)
    write_tree(tree, args, format=args.format)

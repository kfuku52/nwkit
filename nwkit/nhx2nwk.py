from nwkit.util import *

def nhx2nwk_main(args):
    tree = read_tree(args)
    write_tree(tree, args)

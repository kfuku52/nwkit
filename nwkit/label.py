from nwkit.util import *

def label_main(args):
    tree = read_tree(args.infile, args.format, args.quoted_node_names)
    nodes = get_target_nodes(tree=tree, target=args.target)
    counter = 0
    for node in nodes:
        flag_label = False
        if not node.name:
            flag_label = True
        else:
            if args.force:
                flag_label = True
        if flag_label:
            node.name = args.prefix + str(counter)
            counter += 1
    sys.stderr.write(f'Number of labeled target nodes: {counter}/{len(nodes)}\n')
    write_tree(tree, args, format=1)

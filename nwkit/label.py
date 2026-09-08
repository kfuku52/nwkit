import sys

from nwkit.util import get_target_nodes, read_tree, write_tree


def _assign_unique_labels(tree, nodes, prefix, force, start=0):
    targets = {node for node in nodes if force or not node.name}
    reserved = {
        str(node.name) for node in tree.traverse() if node not in targets and node.name
    }
    counter = 0
    next_number = start
    for node in nodes:
        if node in targets:
            while prefix + str(next_number) in reserved:
                next_number += 1
            node.name = prefix + str(next_number)
            reserved.add(node.name)
            next_number += 1
            counter += 1
    return counter


def label_main(args):
    tree = read_tree(
        args.infile,
        args.format,
        args.quoted_node_names,
        rooted=getattr(args, "input_rooted", "auto"),
    )
    nodes = get_target_nodes(tree=tree, target=args.target)
    counter = _assign_unique_labels(
        tree, nodes, args.prefix, args.force, getattr(args, "start", 0)
    )
    sys.stderr.write(f"Number of labeled target nodes: {counter}/{len(nodes)}\n")
    outformat = args.outformat
    if outformat == "auto":
        outformat = 1
    write_tree(tree, args, format=outformat)

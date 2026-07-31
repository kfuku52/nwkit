import math

from nwkit.util import get_target_nodes, read_tree, write_tree


def rescale_main(args):
    if not math.isfinite(args.factor):
        raise ValueError("'--factor' must be a finite number.")
    tree = read_tree(args.infile, args.format, args.quoted_node_names)
    nodes = get_target_nodes(tree=tree, target=args.target)
    scaled_distances = list()
    for node in nodes:
        if node.dist is not None:
            scaled_distance = node.dist * args.factor
            if not math.isfinite(scaled_distance):
                raise ValueError('Rescaling produced a non-finite branch length.')
            scaled_distances.append((node, scaled_distance))
    for node, scaled_distance in scaled_distances:
        node.dist = scaled_distance
    write_tree(tree, args, format=args.outformat)

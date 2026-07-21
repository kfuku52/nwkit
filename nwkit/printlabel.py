import re
from contextlib import nullcontext

from nwkit.util import get_target_nodes, read_tree

def printlabel_main(args):
    tree = read_tree(args.infile, args.format, args.quoted_node_names)
    nodes = get_target_nodes(tree=tree, target=args.target)
    lines = list()
    for node in nodes:
        if re.fullmatch(args.pattern, node.name or ''):
            if args.sister:
                sister_names = [s.name for s in node.get_sisters() if s.name]
                lines.append(' '.join(sister_names))
            else:
                lines.append(node.name or '')
    outfile = getattr(args, 'outfile', '-')
    output_context = nullcontext(None) if outfile == '-' else open(outfile, mode='w')
    with output_context as handle:
        text = ''.join(f'{line}\n' for line in lines)
        if handle is None:
            print(text, end='')
        else:
            handle.write(text)

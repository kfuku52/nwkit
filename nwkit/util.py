import sys
import Bio.SeqIO
from ete3 import TreeNode

def read_tree(infile, format, quiet=False):
    if infile=='-':
        nwk_string = sys.stdin.readlines()[0]
        tree = TreeNode(newick=nwk_string, format=format)
    else:
        tree = TreeNode(newick=infile, format=format)
    if not quiet:
        num_leaves = len([ n for n in tree.traverse() if n.is_leaf() ])
        sys.stderr.write('number of leaves in input tree: {:,}\n'.format(num_leaves))
    return tree

def write_tree(tree, args, quiet=False):
    if not quiet:
        num_leaves = len([ n for n in tree.traverse() if n.is_leaf() ])
        sys.stderr.write('number of leaves in input tree: {:,}\n'.format(num_leaves))
    tree_str = tree.write(format=args.format, format_root_node=True)
    tree_str = tree_str.replace(':123456789','').replace(':1.23457e+08','')
    tree_str = tree_str.replace('123456789','').replace('1.23457e+08','')
    if args.outfile=='-':
        print(tree_str)
    else:
        with open(args.outfile, mode='w') as f:
            f.write(tree_str)

def read_seqs(seqfile, seqformat, quiet):
    if seqfile=='-':
        parsed = sys.stdin
    else:
        parsed = seqfile
    records = list(Bio.SeqIO.parse(parsed, seqformat))
    if not quiet:
        sys.stderr.write('number of input sequences: {}\n'.format(len(records)))
    return records

def write_seqs(records, outfile, seqformat='fasta', quiet=False):
    if not quiet:
        sys.stderr.write('number of output sequences: {}\n'.format(len(records)))
    if outfile=='-':
        Bio.SeqIO.write(records, sys.stdout, seqformat)
    else:
        Bio.SeqIO.write(records, outfile, seqformat)
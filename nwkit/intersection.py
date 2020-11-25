from nwkit.util import *
import itertools

def intersection_main(args):
    tree = read_tree(args.infile, args.format)
    seqs = read_seqs(seqfile=args.seqin, seqformat=args.seqformat, quiet=False)
    leaf_names = set(tree.get_leaf_names())
    seq_names = set([ s.name for s in seqs ])
    leaves_to_remove = leaf_names - seq_names
    seqs_to_remove = seq_names - leaf_names
    assert len(leaves_to_remove)!=len(leaf_names), 'No overlap was found. Check tree and sequence labels in the input files.'
    for leaf in tree.get_leaves():
        if leaf.name in leaves_to_remove:
            leaf.delete(preserve_branch_length=True)
    new_seqs = list()
    for seq in seqs:
        if not seq.name in seqs_to_remove:
            new_seqs.append(seq)
    write_tree(tree, args)
    write_seqs(records=new_seqs, outfile=args.seqout, seqformat=args.seqformat, quiet=False)


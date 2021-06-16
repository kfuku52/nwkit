from nwkit.util import *

def get_leaf_names(tree):
    leaf_names = tree.get_leaf_names()
    is_unique = (len(leaf_names)==len(set(leaf_names)))
    assert is_unique, 'Leaf names are not unique.'
    return leaf_names

def get_seq_names(seqs):
    seq_names = set([ s.name for s in seqs ])
    is_unique = (len(seq_names)==len(set(seq_names)))
    assert is_unique, 'Sequence names are not unique.'
    return seq_names

def get_new_seqs(seqs, seqs_to_remove):
    new_seqs = list()
    for seq in seqs:
        if not seq.name in seqs_to_remove:
            new_seqs.append(seq)
    return new_seqs

def match_complete(a1, a2):
    return a1==a2

def match_prefix(a1, a2):
    return a1.startswith(a2)|a2.startswith(a1)

def match_backward(a1, a2):
    return a1.endswith(a2)|a2.endswith(a1)

def get_remove_names(arr1, arr2, match):
    remove_names = list()
    if match=='complete':
        match_fun = match_complete
    elif match=='prefix':
        match_fun = match_prefix
    elif match=='backward':
        match_fun = match_backward
    for a1 in arr1:
        flag_remove = True
        for a2 in arr2:
            if match_fun(a1, a2):
                flag_remove = False
                break
        if flag_remove:
            remove_names.append(a1)
    return remove_names

def intersection_main(args):
    tree = read_tree(args.infile, args.format)
    leaf_names = get_leaf_names(tree)
    if (args.infile2 == ''):
        leaf_names2 = []
    else:
        tree2 = read_tree(args.infile2, args.format2)
        leaf_names2 = get_leaf_names(tree2)
    if (args.seqin == ''):
        seq_names = []
    else:
        seqs = read_seqs(seqfile=args.seqin, seqformat=args.seqformat, quiet=False)
        seq_names = get_seq_names(seqs)
        seqs_to_remove = get_remove_names(seq_names, leaf_names + leaf_names2, match=args.match)
        new_seqs = get_new_seqs(seqs, seqs_to_remove)
        write_seqs(records=new_seqs, outfile=args.seqout, seqformat=args.seqformat, quiet=False)
    leaves_to_remove = get_remove_names(leaf_names, seq_names + leaf_names2, match=args.match)
    assert len(leaves_to_remove)!=len(leaf_names), 'No overlap was found. Check tree and sequence labels in the input files.'
    for leaf in tree.get_leaves():
        if leaf.name in leaves_to_remove:
            leaf.delete(preserve_branch_length=True)
    write_tree(tree, args, format=args.format)


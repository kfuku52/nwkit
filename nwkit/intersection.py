from nwkit.util import *

TRIE_TERM = '__TERM__'

def get_leaf_names(tree):
    leaf_names = list(tree.leaf_names())
    is_unique = (len(leaf_names)==len(set(leaf_names)))
    if not is_unique:
        raise ValueError('Leaf names are not unique.')
    return leaf_names

def get_seq_names(seqs):
    seq_names = [s.name for s in seqs]
    is_unique = (len(seq_names)==len(set(seq_names)))
    if not is_unique:
        raise ValueError('Sequence names are not unique.')
    return seq_names

def get_new_seqs(seqs, seqs_to_remove):
    seqs_to_remove = set(seqs_to_remove)
    return [seq for seq in seqs if seq.name not in seqs_to_remove]

def match_complete(a1, a2):
    return a1==a2

def match_prefix(a1, a2):
    return a1.startswith(a2) or a2.startswith(a1)

def match_backward(a1, a2):
    return a1.endswith(a2) or a2.endswith(a1)

def build_prefix_trie(arr):
    trie = dict()
    for item in arr:
        node = trie
        for ch in item:
            node = node.setdefault(ch, dict())
        node[TRIE_TERM] = True
    return trie

def has_prefix_relation(item, trie):
    node = trie
    if TRIE_TERM in node:
        return True
    for ch in item:
        node = node.get(ch)
        if node is None:
            return False
        if TRIE_TERM in node:
            return True
    if TRIE_TERM in node:
        return len(node) > 1
    return len(node) > 0

def validate_labels_for_partial_match(arr, mode):
    has_missing = any(item is None for item in arr)
    if has_missing:
        raise ValueError("Missing labels are not supported with '--match {}'.".format(mode))
    return [str(item) for item in arr]

def get_remove_names(arr1, arr2, match):
    if match=='complete':
        arr2_set = set(arr2)
        return [a1 for a1 in arr1 if a1 not in arr2_set]
    if match=='prefix':
        arr1 = validate_labels_for_partial_match(arr1, mode=match)
        arr2 = validate_labels_for_partial_match(arr2, mode=match)
        trie = build_prefix_trie(arr2)
        return [a1 for a1 in arr1 if not has_prefix_relation(a1, trie)]
    elif match=='backward':
        arr1 = validate_labels_for_partial_match(arr1, mode=match)
        arr2 = validate_labels_for_partial_match(arr2, mode=match)
        reverse_trie = build_prefix_trie([a2[::-1] for a2 in arr2])
        return [a1 for a1 in arr1 if not has_prefix_relation(a1[::-1], reverse_trie)]
    else:
        raise ValueError('Unsupported match mode: {}'.format(match))

def intersection_main(args):
    tree = read_tree(args.infile, args.format, args.quoted_node_names)
    leaf_names = get_leaf_names(tree)
    if (args.infile2 == '') and (args.seqin == ''):
        raise ValueError("Specify at least one of '--infile2' or '--seqin'.")
    if args.seqin != '':
        seqout = args.seqout if args.seqout != '' else '-'
        if (args.outfile == '-') and (seqout == '-'):
            raise ValueError("Tree and sequence outputs cannot both be written to stdout. Set '--outfile' or '--seqout' to a file path.")
    if (args.infile2 == ''):
        leaf_names2 = []
    else:
        tree2 = read_tree(args.infile2, args.format2, args.quoted_node_names)
        leaf_names2 = get_leaf_names(tree2)
    if (args.seqin == ''):
        seq_names = []
    else:
        seqs = read_seqs(seqfile=args.seqin, seqformat=args.seqformat, quiet=False)
        seq_names = get_seq_names(seqs)
        if args.match == 'complete':
            reference_names = set(leaf_names)
            reference_names.update(leaf_names2)
            seqs_to_remove = [name for name in seq_names if name not in reference_names]
        else:
            seqs_to_remove = get_remove_names(seq_names, leaf_names + leaf_names2, match=args.match)
        new_seqs = get_new_seqs(seqs, seqs_to_remove)
        seqout = args.seqout if args.seqout != '' else '-'
        write_seqs(records=new_seqs, outfile=seqout, seqformat=args.seqformat, quiet=False)
    if args.match == 'complete':
        retained_names = set(seq_names)
        retained_names.update(leaf_names2)
        leaves_to_remove = {name for name in leaf_names if name not in retained_names}
    else:
        leaves_to_remove = set(get_remove_names(leaf_names, seq_names + leaf_names2, match=args.match))
    if len(leaves_to_remove) == len(leaf_names):
        raise ValueError('No overlap was found. Check tree and sequence labels in the input files.')
    for leaf in tree.leaves():
        if leaf.name in leaves_to_remove:
            leaf.delete(preserve_branch_length=True)
    write_tree(tree, args, format=args.outformat)

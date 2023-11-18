import sys
import Bio.SeqIO
from ete3 import TreeNode

def read_tree(infile, format, quoted_node_names, quiet=False):
    if infile=='-':
        nwk_string = sys.stdin.readlines()[0]
        tree = TreeNode(newick=nwk_string, format=format, quoted_node_names=quoted_node_names)
    else:
        tree = TreeNode(newick=infile, format=format, quoted_node_names=quoted_node_names)
    if not quiet:
        num_leaves = len([ n for n in tree.traverse() if n.is_leaf() ])
        sys.stderr.write('Number of leaves in input tree: {:,}\n'.format(num_leaves))
    return tree

def write_tree(tree, args, format, quiet=False):
    if not quiet:
        num_leaves = len([ n for n in tree.traverse() if n.is_leaf() ])
        sys.stderr.write('Number of leaves in output tree: {:,}\n'.format(num_leaves))
    tree_str = tree.write(format=format, format_root_node=True)
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
        sys.stderr.write('Number of input sequences: {:,}\n'.format(len(records)))
    return records

def write_seqs(records, outfile, seqformat='fasta', quiet=False):
    if not quiet:
        sys.stderr.write('Number of output sequences: {:,}\n'.format(len(records)))
    if outfile=='-':
        Bio.SeqIO.write(records, sys.stdout, seqformat)
    else:
        Bio.SeqIO.write(records, outfile, seqformat)

def remove_singleton(tree, verbose=False, preserve_branch_length=True):
    for node in tree.traverse():
        if node.is_leaf():
            continue
        num_children = len(node.get_children())
        if (num_children>1):
            continue
        if verbose:
            sys.stderr.write('Deleting a singleton node: {}\n'.format(node.name))
        node.delete(prevent_nondicotomic=False, preserve_branch_length=preserve_branch_length)
    return tree

def label2sciname(labels, in_delim='_', out_delim='_'):
    is_str_input = isinstance(labels, str)
    if is_str_input:
        labels = [labels,]
    scinames = list()
    for label in labels:
        splitted = label.split(in_delim)
        if len(splitted)>=2:
            sciname = splitted[0]+out_delim+splitted[1]
        else:
            sciname = None
        scinames.append(sciname)
    if is_str_input:
        scinames = scinames[0]
    return scinames

def read_item_per_line_file(file):
    with open(file, 'r') as f:
        out = f.read().split('\n')
    out = [ o for o in out if o!='' ]
    return out

def annotate_scientific_names(tree):
    for node in tree.iter_leaves():
        node.sci_name = label2sciname(node.name)
    return tree

def annotate_duplication_confidence_scores(tree):
    for node in tree.traverse():
        if node.is_leaf():
            continue
        sp_child1 = set([leaf.sci_name for leaf in node.children[0].iter_leaves()])
        sp_child2 = set([leaf.sci_name for leaf in node.children[1].iter_leaves()])
        num_union = len(sp_child1.union(sp_child2))
        num_intersection = len(sp_child1.intersection(sp_child2))
        node.dup_conf_score = num_intersection / num_union
    return tree

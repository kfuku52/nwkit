import sys
import Bio.SeqIO
from ete3 import TreeNode

def read_tree(infile, format, quoted_node_names, quiet=False):
    global INFILE_FORMAT
    if infile=='-':
        infile = sys.stdin.readlines()[0]
    if format=='auto':
        format_original = format
        for format in [0,1,2,3,4,5,6,7,8,9,100,'exception']:
            if format == 'exception':
                raise Exception('Failed to parse the input tree.')
            try:
                tree = TreeNode(newick=infile, format=format, quoted_node_names=True)
                INFILE_FORMAT = format
                break
            except:
                pass
            try:
                tree = TreeNode(newick=infile, format=format, quoted_node_names=False)
                INFILE_FORMAT = format
                break
            except:
                pass
    else:
        format_original = format
        format = int(format)
        tree = TreeNode(newick=infile, format=format, quoted_node_names=quoted_node_names)
    if format==0: # flexible with support values
        max_support = max([ node.support for node in tree.traverse() ])
        if max_support > 1.0:
            for node in tree.traverse():
                if (node.support - 1.0) < 10**-9: # 1.0 is the default value for missing support values
                    node.support = -999999
    if format==1: # flexible with internal node names
        for node in tree.traverse():
            node.support = -999999
    if not quiet:
        num_leaves = len(tree.get_leaf_names())
        txt = 'Number of leaves in input tree = {:,}, Input tree format = {}\n'
        sys.stderr.write(txt.format(num_leaves, format))
    return tree

def write_tree(tree, args, format, quiet=False):
    if format=='auto':
        format_original = format
        if sys.argv[1] == 'mark':
            format = 1
        else:
            format = INFILE_FORMAT
    else:
        format_original = format
        format = int(format)
    if not quiet:
        num_leaves = len(tree.get_leaf_names())
        txt = 'Number of leaves in output tree = {:,}, Output tree format = {}\n'
        sys.stderr.write(txt.format(num_leaves, format))
    tree_str = tree.write(format=format, format_root_node=True)
    tree_str = tree_str.replace(':-999999.0', '').replace(':-999999','')
    tree_str = tree_str.replace('-999999.0', '').replace('-999999','')
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

def is_all_leaf_names_identical(tree1, tree2, verbose=False):
    leaf_names1 = set(tree1.get_leaf_names())
    leaf_names2 = set(tree2.get_leaf_names())
    is_all_leaf_names_identical = (leaf_names1 == leaf_names2)
    if verbose:
        if not is_all_leaf_names_identical:
            unmatched_names = leaf_names1.symmetric_difference(leaf_names2)
            sys.stderr.write('Unmatched leaf labels: {}\n'.format(' '.join(unmatched_names)))
    return is_all_leaf_names_identical

def get_target_nodes(tree, target):
    if (target=='all'):
        nodes = list(tree.traverse())
    elif (target=='root'):
        nodes = [ node for node in tree.traverse() if node.is_root() ]
    elif (target=='leaf'):
        nodes = [ node for node in tree.traverse() if node.is_leaf() ]
    elif (target=='intnode'):
        nodes = [ node for node in tree.traverse() if not node.is_leaf() ]
    return nodes

def is_rooted(tree):
    num_subroot_nodes = len(tree.get_children())
    if num_subroot_nodes==2:
        return True
    elif num_subroot_nodes==3:
        return False
    else:
        raise Exception('Number of subroot nodes should be either 2 or 3.')

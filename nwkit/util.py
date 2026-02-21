import os
import re
import sys
from collections import Counter
import Bio.SeqIO
from ete4 import Tree

NODENAME_PLACEHOLDER_PATTERN = re.compile(r'NODENAME_PLACEHOLDER\d{10}')
TREE_FORMAT_PROP = '_nwkit_parser_format'

def read_tree(infile, format, quoted_node_names, quiet=False):
    global INFILE_FORMAT
    if infile=='-':
        infile = ''.join(sys.stdin.readlines()).strip()
        if infile == '':
            raise Exception('Failed to parse the input tree.')
    elif os.path.isfile(infile):
        with open(infile) as f:
            infile = f.read().strip()
    if format=='auto':
        format_original = format
        for format in [0,1,2,3,4,5,6,7,8,9,100,'exception']:
            if format == 'exception':
                raise Exception('Failed to parse the input tree.')
            try:
                tree = Tree(infile, parser=format)
                INFILE_FORMAT = format
                break
            except:
                pass
    else:
        format_original = format
        format = int(format)
        tree = Tree(infile, parser=format)
    INFILE_FORMAT = int(format)
    # Keep per-node format metadata so writing a subtree keeps the original parser.
    for node in tree.traverse():
        node.props[TREE_FORMAT_PROP] = INFILE_FORMAT
    if format==0: # flexible with support values
        # Single pass over nodes to detect max support and collect candidates.
        # If max support is > 1.0, treat non-root None/1.0 as missing support.
        max_support = None
        missing_support_candidates = list()
        for node in tree.traverse():
            support = node.support
            if support is not None:
                if (max_support is None) or (support > max_support):
                    max_support = support
            if not node.is_root:
                if support is None or (support - 1.0) < 10**-9: # 1.0 is default for missing support values
                    missing_support_candidates.append(node)
        if (max_support is not None) and (max_support > 1.0):
            for node in missing_support_candidates:
                node.support = -999999
    if format==1: # flexible with internal node names
        for node in tree.traverse():
            if not node.is_root:  # Don't set support on root (breaks ete4 set_outgroup)
                node.support = -999999
    if not quiet:
        num_leaves = len(list(tree.leaves()))
        txt = 'Number of leaves in input tree = {:,}, Input tree format = {}\n'
        sys.stderr.write(txt.format(num_leaves, format))
    return tree

def write_tree(tree, args, format, quiet=False):
    if format=='auto':
        format_original = format
        if len(sys.argv) > 1 and sys.argv[1] == 'mark':
            format = 1
        else:
            format = tree.props.get(TREE_FORMAT_PROP, globals().get('INFILE_FORMAT', 0))
        format = int(format)
    else:
        format_original = format
        format = int(format)
    if not quiet:
        num_leaves = len(list(tree.leaves()))
        txt = 'Number of leaves in output tree = {:,}, Output tree format = {}\n'
        sys.stderr.write(txt.format(num_leaves, format))
    node_name_dict = dict()
    i = 0
    for node in tree.traverse():
        if node.name and str(node.name) != '-999999':
            placeholder_name = 'NODENAME_PLACEHOLDER'+str(i).zfill(10)
            node_name_dict[placeholder_name] = node.name
            node.name = placeholder_name
            i += 1
    tree_str = tree.write(parser=format, format_root_node=True)
    tree_str = tree_str.replace('-999999.0', '').replace('-999999','')
    if tree_str.endswith(':;'):
        tree_str = tree_str[:-2]+';'
    if node_name_dict:
        tree_str = NODENAME_PLACEHOLDER_PATTERN.sub(lambda m: str(node_name_dict[m.group(0)]), tree_str)
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
        if node.is_leaf:
            continue
        num_children = len(node.get_children())
        if (num_children>1):
            continue
        if verbose:
            sys.stderr.write('Deleting a singleton node: {}\n'.format(node.name))
        node.delete(prevent_nondicotomic=False, preserve_branch_length=preserve_branch_length)
    return tree

def label2sciname(labels, in_delim='_', out_delim='_'):
    if labels is None:
        return None
    is_str_input = isinstance(labels, str)
    if is_str_input:
        labels = [labels,]
    scinames = list()
    for label in labels:
        if label is None:
            sciname = None
        else:
            label_str = str(label)
            splitted = label_str.split(in_delim)
            if len(splitted)>=2:
                sciname = splitted[0]+out_delim+splitted[1]
            else:
                sciname = None
        scinames.append(sciname)
    if is_str_input:
        scinames = scinames[0]
    return scinames

def get_subtree_leaf_name_sets(tree):
    subtree_leaf_name_sets = dict()
    for node in tree.traverse(strategy='postorder'):
        if node.is_leaf:
            subtree_leaf_name_sets[node] = {node.name}
            continue
        leaf_names = set()
        for child in node.get_children():
            leaf_names.update(subtree_leaf_name_sets[child])
        subtree_leaf_name_sets[node] = leaf_names
    return subtree_leaf_name_sets

def get_subtree_leaf_bitmasks(tree, leaf_name_to_bit):
    subtree_leaf_bitmasks = dict()
    for node in tree.traverse(strategy='postorder'):
        if node.is_leaf:
            if node.name not in leaf_name_to_bit:
                raise ValueError('Leaf label not found in reference mapping: {}'.format(node.name))
            subtree_leaf_bitmasks[node] = 1 << leaf_name_to_bit[node.name]
            continue
        bitmask = 0
        for child in node.get_children():
            bitmask |= subtree_leaf_bitmasks[child]
        subtree_leaf_bitmasks[node] = bitmask
    return subtree_leaf_bitmasks

def validate_unique_named_leaves(tree, option_name, context=''):
    leaf_names = list(tree.leaf_names())
    if len(leaf_names) != len(set(leaf_names)):
        raise ValueError("Duplicated leaf labels are not supported in '{}'{}.".format(option_name, context))
    has_missing = any((name is None) or (str(name) == '') for name in leaf_names)
    if has_missing:
        raise ValueError("Empty leaf labels are not supported in '{}'{}.".format(option_name, context))

def read_item_per_line_file(file):
    with open(file, 'r') as f:
        out = f.read().splitlines()
    out = [o.strip() for o in out if o.strip() != '']
    return out

def annotate_scientific_names(tree):
    for node in tree.leaves():
        node.add_props(sci_name=label2sciname(node.name))
    return tree

def get_subtree_sci_name_sets(tree):
    subtree_sci_name_sets = dict()
    for node in tree.traverse(strategy='postorder'):
        if node.is_leaf:
            subtree_sci_name_sets[node] = {node.props.get('sci_name')}
            continue
        sci_names = set()
        for child in node.get_children():
            sci_names.update(subtree_sci_name_sets[child])
        subtree_sci_name_sets[node] = sci_names
    return subtree_sci_name_sets

def annotate_duplication_confidence_scores(tree, subtree_sci_name_sets=None):
    if subtree_sci_name_sets is None:
        subtree_sci_name_sets = get_subtree_sci_name_sets(tree)
    for node in tree.traverse():
        if node.is_leaf:
            continue
        children = node.get_children()
        if len(children) != 2:
            node.add_props(dup_conf_score=0.0)
            continue
        sp_child1 = subtree_sci_name_sets[children[0]]
        sp_child2 = subtree_sci_name_sets[children[1]]
        num_union = len(sp_child1.union(sp_child2))
        num_intersection = len(sp_child1.intersection(sp_child2))
        if num_union == 0:
            node.add_props(dup_conf_score=0.0)
        else:
            node.add_props(dup_conf_score=num_intersection / num_union)
    return tree

def is_all_leaf_names_identical(tree1, tree2, verbose=False):
    leaf_names1 = list(tree1.leaf_names())
    leaf_names2 = list(tree2.leaf_names())
    counts1 = Counter(leaf_names1)
    counts2 = Counter(leaf_names2)
    is_all_leaf_names_identical = (counts1 == counts2)
    if verbose:
        if not is_all_leaf_names_identical:
            all_names = set(counts1.keys()).union(set(counts2.keys()))
            unmatched_names = [
                str(name) for name in sorted(all_names, key=lambda x: str(x))
                if counts1.get(name, 0) != counts2.get(name, 0)
            ]
            sys.stderr.write('Unmatched leaf labels: {}\n'.format(' '.join(unmatched_names)))
    return is_all_leaf_names_identical

def get_target_nodes(tree, target):
    if (target=='all'):
        nodes = list(tree.traverse())
    elif (target=='root'):
        nodes = [tree]
    elif (target=='leaf'):
        nodes = list(tree.leaves())
    elif (target=='intnode'):
        nodes = [ node for node in tree.traverse() if not node.is_leaf ]
    else:
        raise ValueError('Unknown target: {}'.format(target))
    return nodes

def is_rooted(tree):
    num_subroot_nodes = len(tree.get_children())
    if num_subroot_nodes <= 2:
        return True
    return False

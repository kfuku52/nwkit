"""Validated plain Newick exports for downstream phylogenetic programs."""

from pathlib import Path

import ete4


def prepare_species_tree(input_path, output_path):
    from ete4.parser import newick

    from nwkit.file_paths import validate_outputs_do_not_replace_inputs
    from nwkit.output_transaction import output_transaction
    from nwkit.rooting_state import get_rooting_info
    from nwkit.util import read_tree

    validate_outputs_do_not_replace_inputs(
        [("species tree", input_path)],
        [("GeneRax input", output_path)],
    )
    tree = read_tree(str(input_path), "auto", True, quiet=True)
    if get_rooting_info(tree).rooted is not True:
        raise ValueError("GeneRax species tree must be rooted.")
    tips = list(tree.leaf_names())
    if len(tips) < 2 or len(tips) != len(set(tips)):
        raise ValueError("GeneRax species tree requires at least two unique tips.")
    if any(len(node.children) != 2 for node in tree.traverse() if not node.is_leaf):
        raise ValueError("GeneRax requires a rooted, bifurcating species tree.")

    # ETE's default distance formatter rounds to six significant digits.
    parser = {
        kind: [dict(field) for field in fields]
        for kind, fields in newick.PARSERS[1].items()
    }
    for fields in parser.values():
        fields[1]["write"] = lambda value: format(value, ".17g")
    text = tree.write(parser=parser, props=[], format_root_node=True) + "\n"

    def signature(candidate):
        descendants = candidate.get_cached_content()
        return {
            frozenset(tip.name for tip in leaves): (
                node.name or "",
                node.props.get("dist"),
            )
            for node, leaves in descendants.items()
        }

    # Verify the root, every descendant clade, labels and branch lengths before
    # publishing the plain copy. Supports and other annotations are tool-irrelevant.
    from ete4 import Tree

    if signature(Tree(text, parser=1)) != signature(tree):
        raise ValueError("GeneRax serialization changed the rooted species tree.")
    with output_transaction([output_path], create_parents=True) as staged:
        Path(staged[output_path]).write_text(text, encoding="utf-8")


def prepare_target(input_path, output_path):
    from nwkit.file_paths import validate_outputs_do_not_replace_inputs
    from nwkit.output_transaction import output_transaction

    validate_outputs_do_not_replace_inputs(
        [("GeneRax tree", input_path)], [("UFBoot target", output_path)]
    )
    text = Path(input_path).read_text(encoding="utf-8").strip()
    if not text:
        raise ValueError(f"GeneRax tree is empty: {input_path}")
    tree = ete4.PhyloTree(text, parser=1)
    leaf_names = [str(leaf.name) for leaf in tree.leaves()]
    if len(leaf_names) < 4:
        raise ValueError("An unrooted UFBoot target requires at least four tips.")
    if len(leaf_names) != len(set(leaf_names)):
        raise ValueError("GeneRax UFBoot target contains duplicate leaf labels.")

    for node in tree.traverse():
        if node.is_leaf:
            continue
        node.name = ""
        node.props.pop("support", None)

    if len(tree.get_children()) == 2:
        tree.unroot()
    if len(tree.get_children()) < 3:
        raise ValueError(
            "Failed to suppress the GeneRax root for unrooted support mapping."
        )

    text = tree.write(parser=1)
    roundtrip = ete4.PhyloTree(text, parser=1)
    roundtrip_names = [str(leaf.name) for leaf in roundtrip.leaves()]
    if set(roundtrip_names) != set(leaf_names) or len(roundtrip_names) != len(
        leaf_names
    ):
        raise ValueError("Unrooting changed the GeneRax target tip set.")
    if len(roundtrip.get_children()) < 3:
        raise ValueError("Prepared GeneRax support target is still rooted.")
    with output_transaction([output_path], create_parents=True) as staged:
        Path(staged[output_path]).write_text(text, encoding="utf-8")


def topology_newick(text, parser=0):
    if not text.strip():
        raise ValueError("Input Newick is empty")
    return ete4.Tree(text, parser=parser).write(parser=9) + "\n"

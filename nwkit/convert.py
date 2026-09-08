"""Convert one tree between Newick, NHX, FigTree and MCMCtree containers."""

import sys
from decimal import Decimal
from pathlib import Path

from nwkit.output_transaction import output_transaction
from nwkit.rooting_state import extract_rooting_token, get_rooting_info
from nwkit.tree_formats import (
    _dated_tree,
    annotation_attributes,
    finite_decimal,
    read_container,
    select_mcmctree_tree,
    tokens,
    transform_annotations,
)
from nwkit.util import read_input_text, read_tree, validate_unique_named_leaves


def _select_statement(source, statements, tree_index):
    if tree_index is not None:
        if not isinstance(tree_index, int) or not 1 <= tree_index <= len(statements):
            raise ValueError(f"--tree-index must be between 1 and {len(statements)}.")
        statement = statements[tree_index - 1]
        if source == "mcmctree-output" and _dated_tree(statement) is None:
            raise ValueError(
                "Selected MCMCtree statement is a topology/index tree, not a dated tree."
            )
        return statement
    if source == "mcmctree-output":
        return select_mcmctree_tree(statements)
    if len(statements) != 1:
        raise ValueError(
            "Multiple input trees are ambiguous; select one with --tree-index."
        )
    return statements[0]


def _rewrite_rooting(statement, state):
    """Canonicalize explicit/overridden root declarations without NHX loss."""
    statement, _ = extract_rooting_token(statement)
    result = []
    for token in tokens(statement):
        if token.kind == "comment" and token.text.startswith("[&&NHX:"):
            attributes, _ = annotation_attributes(token.text)
            attributes.pop("nwkit_rooted", None)
            if attributes:
                result.append(
                    "[&&NHX:"
                    + ":".join(f"{key}={value}" for key, value in attributes.items())
                    + "]"
                )
        else:
            result.append(token.text)
    marker = "[&R]" if state == "rooted" else "[&U]"
    return marker + "".join(result)


def convert_tree_text(
    text,
    *,
    source="auto",
    target="nhx",
    time_factor=Decimal(1),
    age_ci="keep",
    tree_index=None,
    tree_format="auto",
    quoted_node_names=True,
    rooted="auto",
):
    """Serialize a validated result fully before the caller publishes any bytes."""
    if target not in {"newick", "nhx", "figtree"}:
        raise ValueError(f"Unsupported output format: {target}")
    factor = finite_decimal(time_factor, positive=True)
    detected, statements = read_container(text, source)
    statement = _select_statement(detected, statements, tree_index)
    tree = read_tree(
        statement, tree_format, quoted_node_names, quiet=True, rooted=rooted
    )
    validate_unique_named_leaves(tree, "--infile", context=" for convert")
    for node in tree.traverse():
        if node.dist is not None and node.dist < 0:
            raise ValueError("Tree branch lengths must be non-negative.")
    info = get_rooting_info(tree)
    if info.state != "unknown" and (
        info.source != "topology" or rooted != "auto" or target == "figtree"
    ):
        statement = _rewrite_rooting(statement, info.state)
    converted = transform_annotations(
        statement, output=target, factor=factor, age_ci=age_ci
    )
    # Verify the emitted numeric values and the exact downstream reader boundary.
    # Names/support are carried lexically, not renamed or rounded by an ETE writer.
    read_tree(converted, tree_format, quoted_node_names, quiet=True)
    if target == "figtree":
        return "#NEXUS\nBEGIN TREES;\n  UTREE 1 = " + converted + "\nEND;\n"
    return converted + "\n"


def convert_main(args):
    converted = convert_tree_text(
        read_input_text(args.infile),
        source=getattr(args, "from"),
        target=args.to,
        time_factor=args.time_factor,
        age_ci=args.age_ci,
        tree_index=args.tree_index,
        tree_format=args.format,
        quoted_node_names=args.quoted_node_names,
        rooted=getattr(args, "input_rooted", "auto"),
    )
    if args.outfile == "-":
        sys.stdout.write(converted)
    else:
        with output_transaction([args.outfile]) as staged:
            Path(staged[args.outfile]).write_text(converted, encoding="utf-8")

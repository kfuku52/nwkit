"""Convert one tree between Newick, NHX, FigTree and MCMCtree containers."""

import sys
from decimal import Decimal
from pathlib import Path

from nwkit.output_transaction import output_transaction
from nwkit.rooting_state import (
    ROOTED_PROP,
    extract_rooting_token,
    get_rooting_info,
    rooting_output_options,
    rooting_output_policy,
    topology_rooting,
)
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


def _strip_rooting(statement):
    """Remove only root declarations before applying the shared output policy."""
    statement, _ = extract_rooting_token(statement)
    result = []
    for token in tokens(statement):
        if token.kind == "comment" and token.text.startswith("[&&NHX:"):
            attributes, _ = annotation_attributes(token.text)
            attributes.pop(ROOTED_PROP, None)
            if attributes:
                result.append(
                    "[&&NHX:"
                    + ":".join(f"{key}={value}" for key, value in attributes.items())
                    + "]"
                )
        else:
            result.append(token.text)
    return "".join(result)


def _validate_reserved_properties(statement):
    """Reject attributes that ETE interprets as overrides of Newick fields."""
    for token in tokens(statement):
        if token.kind != "comment":
            continue
        attributes, _ = annotation_attributes(token.text)
        reserved = {"name", "dist", "support"} & attributes.keys()
        if reserved:
            raise ValueError(
                "Reserved NHX properties are ambiguous in convert: "
                + ", ".join(sorted(reserved))
                + ". Put node names, branch lengths and support in Newick fields, "
                "or rename these properties before conversion."
            )


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
    node_label="",
    properties="keep",
    rooting_token=False,
    rooting_nhx=False,
):
    """Serialize a validated result fully before the caller publishes any bytes."""
    if target not in {"newick", "nhx", "figtree"}:
        raise ValueError(f"Unsupported output format: {target}")
    factor = finite_decimal(time_factor, positive=True)
    detected, statements = read_container(text, source)
    statement = _select_statement(detected, statements, tree_index)
    _validate_reserved_properties(statement)
    tree = read_tree(
        statement, tree_format, quoted_node_names, quiet=True, rooted=rooted
    )
    validate_unique_named_leaves(tree, "--infile", context=" for convert")
    for node in tree.traverse():
        if node.dist is not None and node.dist < 0:
            raise ValueError("Tree branch lengths must be non-negative.")
    info = get_rooting_info(tree)
    prefix, nhx = rooting_output_policy(
        tree, rooting_token=rooting_token, rooting_nhx=rooting_nhx
    )
    if target == "newick":
        if rooting_nhx or (not prefix and info.rooted is not topology_rooting(tree)):
            raise ValueError(
                "Plain Newick without a rooting token cannot preserve this rooting "
                "state. Use --rooting-token yes (without --rooting-nhx yes), "
                "or --to nhx/figtree."
            )
        # Root NHX is a semantic declaration, not an ordinary property to drop.
        # Plain Newick uses the requested token or the unchanged binary root.
        nhx = False
    elif target == "figtree" and info.rooted is not None:
        # A NEXUS tree statement retains its format-specific declaration even
        # with standalone-Newick tokens disabled or root NHX requested.
        prefix = "[&R]" if info.rooted else "[&U]"
    if node_label:
        # Copy from the original attributes, including nwkit_rooted, before
        # canonicalizing rooting declarations or dropping/scaling properties.
        statement = transform_annotations(statement, node_label=node_label)
    statement = _strip_rooting(statement)
    converted = transform_annotations(
        statement,
        output=target,
        factor=factor,
        age_ci=age_ci,
        properties=properties,
    )
    if nhx:
        state = "unknown" if info.rooted is None else ("yes" if info.rooted else "no")
        converted = converted[:-1] + f"[&&NHX:{ROOTED_PROP}={state}];"
    converted = prefix + converted
    # Verify the emitted numeric values and the exact downstream reader boundary.
    # The input quoting policy was checked above. Generated labels are always
    # safely quoted, independent of whether input quotes were permitted.
    restored = read_tree(converted, 1 if node_label else tree_format, True, quiet=True)
    if get_rooting_info(restored).rooted is not info.rooted:
        raise ValueError("Converted tree lost its rooting interpretation.")
    if target == "figtree":
        kind = "UTREE" if info.rooted is False else "TREE"
        return f"#NEXUS\nBEGIN TREES;\n  {kind} 1 = " + converted + "\nEND;\n"
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
        node_label=getattr(args, "node_label", ""),
        properties=getattr(args, "properties", "keep"),
        **rooting_output_options(args),
    )
    if args.outfile == "-":
        sys.stdout.write(converted)
    else:
        with output_transaction([args.outfile]) as staged:
            Path(staged[args.outfile]).write_text(converted, encoding="utf-8")

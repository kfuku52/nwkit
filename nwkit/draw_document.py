"""N1/N2 document drawing entry point; no reconciliation or fitting is invoked."""

from __future__ import annotations

import sys

from nwkit.file_paths import validate_outputs_do_not_replace_inputs
from nwkit.output_transaction import output_transaction
from nwkit.reconciliation_document import (
    read_document,
    select_document_clade,
    validate_document,
)
from nwkit.reconciliation_export import document_text
from nwkit.reconciliation_import import import_saved_tsv
from nwkit.reconciliation_render import write_scene
from nwkit.reconciliation_scene import build_scene
from nwkit.result_plot import figure_format


def _validate_options(args):
    from nwkit.cli import parser

    defaults = vars(parser.parse_args(["draw"]))
    allowed = {
        "document",
        "document_out",
        "document_clade",
        "view",
        "pdf_pages",
        "pdf_page_size",
        "outfile",
        "image_format",
        "font_size",
        "font_family",
        "figure_width",
        "legend",
        "transparent",
        "debug",
        "audit",
    }
    if not args.document:
        allowed.update(
            {
                "reconciliation",
                "species_tree",
                "species_tree_format",
                "species_tree_rooted",
                "infile",
                "format",
                "quoted_node_names",
                "input_rooted",
            }
        )
    for key, value in vars(args).items():
        if (
            key not in allowed
            and not key.startswith("_nwkit_")
            and key in defaults
            and value != defaults[key]
        ):
            raise ValueError(
                f"--{key.replace('_', '-')} is not supported for document drawing."
            )
    if args.document == "-" or args.document_out == "-":
        raise ValueError("--document and --document-out require filesystem paths.")
    if args.document and args.infile != "-":
        raise ValueError("--document cannot be combined with --infile.")
    if not args.document and (not args.reconciliation or not args.species_tree):
        raise ValueError(
            "Embedded saved TSV drawing requires --reconciliation and --species-tree."
        )


def draw_document_main(args):
    _validate_options(args)
    inputs = (
        [("document", args.document)]
        if args.document
        else [
            ("infile", args.infile),
            ("species-tree", args.species_tree),
            ("reconciliation", args.reconciliation),
        ]
    )
    outputs = [("outfile", args.outfile), ("document-out", args.document_out)]
    validate_outputs_do_not_replace_inputs(inputs, outputs)
    document = (
        read_document(args.document)
        if args.document
        else validate_document(import_saved_tsv(args))
    )
    if args.document_clade:
        document = select_document_clade(document, args.document_clade)
    fmt = figure_format(args.outfile, args.image_format)
    width = 72 * (args.figure_width if args.figure_width is not None else 10)
    scene = build_scene(document, width=width, font_size=args.font_size)
    normalized = (
        document_text(document, args.document_out) if args.document_out else None
    )
    paths = [args.outfile] + ([args.document_out] if args.document_out else [])
    with output_transaction(paths, create_parents=True) as staged:
        write_scene(
            scene,
            staged[args.outfile],
            title=document.record.get("title", ""),
            image_format=fmt,
            font_size=args.font_size,
            font_family="DejaVu Sans"
            if args.font_family == "Helvetica"
            and not getattr(args, "_nwkit_draw_font_explicit", False)
            else args.font_family,
            transparent=args.transparent,
            pdf_pages=args.pdf_pages,
            pdf_page_size=args.pdf_page_size,
            legend=args.legend,
        )
        if args.document_out:
            staged.write_text(
                args.document_out, lambda stream: stream.write(normalized)
            )
    sys.stderr.write(f"Wrote document figure: {args.outfile}\n")

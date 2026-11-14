"""N1/N2 document drawing entry point; no reconciliation or fitting is invoked."""

from __future__ import annotations

import base64
import json
import sys
from pathlib import Path
from typing import Any

from nwkit.file_paths import validate_outputs_do_not_replace_inputs
from nwkit.output_transaction import output_transaction
from nwkit.reconciliation_composition import apply_view, comparison_scene
from nwkit.reconciliation_document import (
    read_document,
    select_document_clade,
    validate_document,
)
from nwkit.reconciliation_export import document_text
from nwkit.reconciliation_html import write_html
from nwkit.reconciliation_import import import_saved_tsv
from nwkit.reconciliation_render import write_scene
from nwkit.reconciliation_scene import build_scene
from nwkit.reconciliation_view import (
    default_view,
    read_comparison,
    read_view,
    validate_comparison_view,
    validate_view,
)
from nwkit.result_plot import figure_format


def _validate_options(args):
    from nwkit.cli import parser

    defaults = vars(parser.parse_args(["draw"]))
    allowed = {
        "document",
        "comparison",
        "view_config",
        "view_out",
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
    if not args.document and not args.comparison:
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
    if (args.document or args.comparison) and args.infile != "-":
        raise ValueError("--document cannot be combined with --infile.")
    if (
        not args.document
        and not args.comparison
        and (not args.reconciliation or not args.species_tree)
    ):
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
    if args.comparison:
        record, documents, panel_paths = read_comparison(args.comparison)
        inputs = [("comparison", args.comparison)] + [
            ("panel", path) for path in panel_paths
        ]
    else:
        record = documents = None
    if args.view_config:
        inputs.append(("view-config", args.view_config))
    outputs = [
        ("outfile", args.outfile),
        ("document-out", args.document_out),
        ("view-out", args.view_out),
    ]
    validate_outputs_do_not_replace_inputs(inputs, outputs)
    document = (
        None
        if args.comparison
        else (
            read_document(args.document)
            if args.document
            else validate_document(import_saved_tsv(args))
        )
    )
    if args.document_clade:
        if document is None:
            raise ValueError("Comparison clade selection belongs in --view-config.")
        document = select_document_clade(document, args.document_clade)
    view: dict[str, Any] = (
        read_view(args.view_config, document)
        if args.view_config
        else (record["view"] if record else default_view(document))
    )
    if args.view:
        view = dict(
            view,
            layout={
                "embedded": "topology_embedded",
                "observed-embedded": "observed_embedded",
            }.get(args.view, args.view),
        )
        if args.view == "observed-embedded":
            view["schema_version"] = "0.2.0"
    if document:
        validate_view(view, document)
    if documents:
        validate_comparison_view(view, documents)
    if record:
        record = dict(record, view=view)
    fmt = (
        "html"
        if Path(args.outfile).suffix.lower() == ".html" and args.image_format == "auto"
        else figure_format(args.outfile, args.image_format)
    )
    width = 72 * (args.figure_width if args.figure_width is not None else 10)
    font_size = view.get("font", {}).get("size", args.font_size)
    font_family = view.get("font", {}).get(
        "family",
        "DejaVu Sans"
        if args.font_family == "Helvetica"
        and not getattr(args, "_nwkit_draw_font_explicit", False)
        else args.font_family,
    )
    explicit: set[str] = getattr(args, "_nwkit_draw_explicit_options", set())
    if "--font-size" in explicit:
        font_size = args.font_size
    if "--font-family" in explicit:
        font_family = args.font_family
    pdf_pages = (
        args.pdf_pages
        if "--pdf-pages" in explicit
        else view.get("output", {}).get("pdf_pages", args.pdf_pages)
    )
    pdf_page_size = (
        args.pdf_page_size
        if "--pdf-page-size" in explicit
        else view.get("output", {}).get("pdf_page_size", args.pdf_page_size)
    )
    if record:
        scene = comparison_scene(record, documents, width=width, font_size=font_size)
        pdf_scene = scene
        if fmt == "html":
            scene = comparison_scene(
                record,
                documents,
                width=width,
                font_size=font_size,
                apply_display_view=False,
            )
        title = "Reconciliation comparison"
        captured = []
        for doc in documents.values():
            if doc.source_bytes is None:
                raise ValueError("Comparison document snapshot is missing.")
            captured.append(base64.b64encode(doc.source_bytes).decode())
        source = dict(
            record,
            documents=[d.record for d in documents.values()],
            document_bytes=captured,
        )
    else:
        assert document is not None
        import copy

        full_scene = build_scene(document, width=width, font_size=font_size, view=view)
        pdf_scene = apply_view(copy.deepcopy(full_scene), document, view)
        scene = (
            apply_view(full_scene, document, dict(view, collapsed_clades=[]))
            if fmt == "html"
            else pdf_scene
        )
        title = document.record.get("title", "")
        source = document.record
    if args.document_out and not document:
        raise ValueError(
            "Comparison sources stay separate; --document-out requires a single document."
        )
    normalized = (
        document_text(document, args.document_out) if args.document_out else None
    )
    paths = (
        [args.outfile]
        + ([args.document_out] if args.document_out else [])
        + ([args.view_out] if args.view_out else [])
    )
    view = dict(
        view,
        font={"family": font_family, "size": font_size},
        output={
            "format": fmt,
            "pdf_pages": pdf_pages,
            "pdf_page_size": pdf_page_size,
        },
    )
    with output_transaction(paths, create_parents=True) as staged:
        options = dict(
            title=title,
            font_size=font_size,
            font_family=font_family,
            transparent=args.transparent,
            pdf_pages=pdf_pages,
            pdf_page_size=pdf_page_size,
            legend=args.legend,
        )
        if fmt == "html":
            write_html(
                scene,
                staged[args.outfile],
                source=source,
                view=view,
                pdf_scene=pdf_scene,
                **options,
            )
        else:
            write_scene(scene, staged[args.outfile], image_format=fmt, **options)
        if args.view_out:
            staged.write_text(
                args.view_out,
                lambda stream: stream.write(
                    json.dumps(view, indent=2, allow_nan=False) + "\n"
                ),
            )
        if args.document_out:
            staged.write_text(
                args.document_out, lambda stream: stream.write(normalized)
            )
    sys.stderr.write(f"Wrote document figure: {args.outfile}\n")

"""Lightweight generic document options for nwkit draw."""


def register_document_options(draw):
    # Added to the existing mutually exclusive saved-result input group.
    group = next(
        g
        for g in draw._mutually_exclusive_groups
        if any(a.dest == "reconciliation" for a in g._group_actions)
    )
    group.add_argument(
        "--document",
        metavar="PATH",
        help="Versioned reconciliation JSON, recPhyloXML, phyloXML or Newick document; no inference is run.",
    )
    draw.add_argument(
        "--document-out",
        "--document_out",
        default=None,
        metavar="PATH",
        help="Write canonical JSON, or round-trip supported imported recPhyloXML (.xml/.recphyloxml), alongside the figure.",
    )
    draw.add_argument(
        "--document-clade",
        "--document_clade",
        default=None,
        metavar="TREE/NODE",
        help="Display only this lower or standalone clade; retain the full upper tree and record display selection.",
    )
    draw.add_argument(
        "--view",
        default=None,
        choices=["embedded"],
        help="Embed recorded lower-tree histories in the upper tree. Requires --document or --reconciliation.",
    )
    draw.add_argument(
        "--pdf-pages",
        "--pdf_pages",
        default="auto",
        choices=["auto", "single", "tiled"],
        help="default=auto: Paginate a tall document PDF or explicitly select single/tiled pages.",
    )
    draw.add_argument(
        "--pdf-page-size",
        "--pdf_page_size",
        default="A4",
        choices=["A4", "A3"],
        help="default=A4: Landscape paper size for tiled document PDF.",
    )

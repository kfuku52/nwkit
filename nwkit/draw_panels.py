"""Arrange ordinary NWKIT tree drawings in a labelled, row-major vector grid."""

import csv
import math
import re
from argparse import Namespace
from pathlib import Path
from tempfile import TemporaryDirectory
from xml.etree import ElementTree as ET

from defusedxml.ElementTree import parse
from matplotlib.font_manager import FontProperties
from matplotlib.textpath import TextToPath

from nwkit.draw import _validate_draw_paths, draw_main
from nwkit.draw_helpers import _resolve_image_format
from nwkit.image import load_cairosvg_module
from nwkit.output_transaction import output_transaction

SVG_NAMESPACE = "http://www.w3.org/2000/svg"
HREF = "{http://www.w3.org/1999/xlink}href"


def read_panel_table(path):
    source = Path(path)
    with source.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t", strict=True)
        if sorted(reader.fieldnames or []) != ["infile", "title"]:
            raise ValueError("--panel-table requires the TSV columns title and infile.")
        panels = []
        for row in reader:
            title = (row["title"] or "").strip()
            infile = (row["infile"] or "").strip()
            if None in row or not title or not infile or infile == "-":
                raise ValueError(
                    f"Invalid panel title or tree path on line {reader.line_num}."
                )
            if re.search(r"[\x00-\x08\x0b\x0c\x0e-\x1f\ufffe\uffff]", title):
                raise ValueError(
                    f"Invalid XML character in panel title on line {reader.line_num}."
                )
            tree_path = source.parent / infile
            if not tree_path.is_file():
                raise FileNotFoundError(
                    f"Panel tree on line {reader.line_num} is not a file: {tree_path}"
                )
            panels.append((title, tree_path))
    if not panels:
        raise ValueError("--panel-table contains no panels.")
    return panels


def _svg_panel(path, prefix):
    # Matplotlib uses fixed IDs for axes and clip paths across independent drawings.
    root = parse(path).getroot()
    for el in root.iter():
        if "id" in el.attrib:
            el.set("id", prefix + el.attrib["id"])
        for attr, value in list(el.attrib.items()):
            if attr in ("href", HREF) and value.startswith("#"):
                el.set(attr, "#" + prefix + value[1:])
            elif "url(#" in value:
                el.set(attr, re.sub(r"url\(#([^)]*)\)", rf"url(#{prefix}\1)", value))
    _, _, width, height = map(float, root.attrib["viewBox"].split())
    return root, width, height


def _title_lines(title, width, font, measure):
    lines = [""]
    for word in title.split():
        line = (lines[-1] + " " + word).strip()
        if (
            lines[-1]
            and measure.get_text_width_height_descent(line, font, False)[0] > width
        ):
            lines.append(word)
        else:
            lines[-1] = line
    return lines


def _sheet(panels, columns, font_size, family, transparent):
    columns = min(columns, len(panels))
    gap = font_size
    cell_width = max(width for _, _, width, _ in panels)
    width = columns * cell_width + (columns - 1) * gap
    root = ET.Element("svg", xmlns=SVG_NAMESPACE, version="1.1")
    if not transparent:
        ET.SubElement(root, "rect", width="100%", height="100%", fill="white")
    font = FontProperties(family=family, size=font_size, weight="bold")
    measure = TextToPath()
    y = 0.0
    for start in range(0, len(panels), columns):
        row = panels[start : start + columns]
        headings = [_title_lines(title, cell_width, font, measure) for title, *_ in row]
        header = (max(map(len, headings)) + 0.4) * font_size * 1.2
        for col, ((_, panel, panel_width, panel_height), lines) in enumerate(
            zip(row, headings, strict=True)
        ):
            x = col * (cell_width + gap)
            for n, line in enumerate(lines):
                el = ET.SubElement(
                    root, "text", x=str(x), y=str(y + font_size + n * 1.2 * font_size)
                )
                el.attrib.update(
                    {
                        "font-family": family,
                        "font-size": str(font_size),
                        "font-weight": "bold",
                    }
                )
                # Keep even an unbreakable identifier within its own column.
                line_width = measure.get_text_width_height_descent(line, font, False)[0]
                if line_width > cell_width:
                    el.set("x", "0")
                    el.set(
                        "transform",
                        f"translate({x:g}) scale({cell_width / line_width:g} 1)",
                    )
                el.text = line
            g = ET.SubElement(
                root,
                "g",
                transform=f"translate({x + (cell_width - panel_width) / 2:g} {y + header:g})",
            )
            # Preserve the single-tree viewport, including its clipping boundary.
            panel.attrib.update(
                width=str(panel_width), height=str(panel_height), overflow="hidden"
            )
            g.append(panel)
        y += header + max(h for _, _, _, h in row) + gap
    height = y - gap
    root.attrib.update(
        {
            "width": f"{width:g}pt",
            "height": f"{height:g}pt",
            "viewBox": f"0 0 {width:g} {height:g}",
        }
    )
    return ET.tostring(root, encoding="utf-8", xml_declaration=True), width


def draw_panels(args):
    document_options = (
        "document",
        "comparison",
        "view",
        "species_tree",
        "reconciliation",
        "radte_prefix",
        "view_config",
        "view_out",
        "document_out",
        "document_clade",
    )
    if (
        any(getattr(args, name, None) for name in document_options)
        or args.pdf_pages != "auto"
        or args.pdf_page_size != "A4"
    ):
        raise ValueError("--panel-table accepts ordinary tree drawings only.")
    columns = args.panel_columns if args.panel_columns is not None else 2
    if columns < 1:
        raise ValueError("--panel-columns must be positive.")
    if args.infile != "-" or args.outfile == "-" or args.layout_report:
        raise ValueError(
            "--panel-table requires a file output and cannot use --infile or --layout-report."
        )
    sources = read_panel_table(args.panel_table)
    inputs = {f"panel {i + 1}": str(path) for i, (_, path) in enumerate(sources)}
    inputs["--panel-table"] = args.panel_table
    for name in (
        "trait",
        "tip_image_manifest",
        "mcmctree_posterior",
        "densitree_trees",
        "species_map_tsv",
    ):
        value = getattr(args, name, None)
        if value == "-":
            raise ValueError(
                f"--{name.replace('_', '-')} cannot read stdin across panels."
            )
        if value:
            inputs[name] = value
    _validate_draw_paths(args.outfile, None, inputs)
    fmt = _resolve_image_format(args.outfile, args.image_format)
    cairo = load_cairosvg_module() if fmt != "svg" else None
    image_paths = set()
    with TemporaryDirectory() as temp:
        panels = []
        for i, (title, infile) in enumerate(sources):
            svg = Path(temp) / f"tree{i}.svg"
            single = Namespace(**vars(args))
            single.infile, single.outfile, single.image_format = (
                str(infile),
                str(svg),
                "svg",
            )
            draw_main(single)
            image_paths.update(single._nwkit_tip_image_paths)
            panels.append((title, *_svg_panel(svg, f"p{i + 1}-")))
    args._nwkit_tip_image_paths = sorted(image_paths)
    if image_paths:
        _validate_draw_paths(
            args.outfile,
            None,
            {f"tip image {i}": path for i, path in enumerate(image_paths)},
        )
    data, width = _sheet(
        panels, columns, args.font_size, args.font_family, args.transparent
    )
    with output_transaction([args.outfile], create_parents=True) as staged:
        output = staged[args.outfile]
        if cairo is None:
            Path(output).write_bytes(data)
        elif fmt == "pdf":
            cairo.svg2pdf(bytestring=data, write_to=str(output))
        else:
            cairo.svg2png(
                bytestring=data,
                write_to=str(output),
                output_width=math.ceil(width * 300 / 72),
            )

from xml.etree import ElementTree as ET

import pytest
from PIL import Image

from nwkit.cli import parser


def draw(manifest, output, *extra):
    args = parser.parse_args(
        [
            "draw",
            "--panel-table",
            str(manifest),
            "-o",
            str(output),
            "--tip-labels",
            "no",
            "--support-labels",
            "no",
            "--species-overlap-node-plot",
            "no",
            *extra,
        ]
    )
    return args.handler(args)


def test_grid_and_namespaced_references(tmp_path, monkeypatch):
    folder = tmp_path / "trees"
    folder.mkdir()
    (folder / "a.nwk").write_text("((A:1,B:2):1,C:2);")
    manifest = folder / "panels.tsv"
    manifest.write_text(
        "title\tinfile\nA & B\ta.nwk\nOther tree\ta.nwk\nLast panel\ta.nwk\n"
    )
    monkeypatch.chdir(tmp_path)
    result = tmp_path / "figure.svg"
    draw(manifest, result, "--panel-columns", "2")
    root = ET.parse(result).getroot()
    headings = [e.text for e in root.findall("{*}text")]
    assert headings == ["A & B", "Other tree", "Last panel"]
    gs = root.findall("{*}g")
    assert len(gs) == 3 and gs[0].get("transform") != gs[1].get("transform")
    ids = [n.get("id") for n in root.iter() if n.get("id")]
    assert len(ids) == len(set(ids))
    for n in root.iter():
        for value in n.attrib.values():
            if value.startswith("url(#"):
                assert value[5:-1] in ids
    assert not any(n.tag.endswith("}image") for n in root.iter())


@pytest.mark.parametrize("fmt", ["pdf", "png"])
def test_grid_other_formats(tmp_path, fmt):
    pytest.importorskip("cairosvg")
    (tmp_path / "t.nwk").write_text("(A:1,B:1);")
    manifest = tmp_path / "t.tsv"
    manifest.write_text("infile\ttitle\nt.nwk\tFirst\n")
    output = tmp_path / f"grid.{fmt}"
    draw(manifest, output, "--figure-width", "2", "--figure-height", "1")
    if fmt == "pdf":
        assert output.read_bytes().startswith(b"%PDF")
    else:
        with Image.open(output) as image:
            assert image.width == 600


def test_invalid_later_panel_preserves_output_and_inputs(tmp_path):
    t = tmp_path / "ok.nwk"
    t.write_text("(A,B);")
    manifest = tmp_path / "panels.tsv"
    manifest.write_text("title\tinfile\nGood\tok.nwk\nMissing\tfail.nwk\n")
    output = tmp_path / "out.svg"
    output.write_text("previous")
    with pytest.raises((OSError, ValueError)):
        draw(manifest, output)
    (tmp_path / "fail.nwk").write_text("(unbalanced")
    with pytest.raises(Exception, match="parse the input tree"):
        draw(manifest, output)
    assert output.read_text() == "previous"
    with pytest.raises(ValueError):
        draw(manifest, t)
    with pytest.raises(ValueError):
        draw(manifest, manifest)
    assert t.read_text() == "(A,B);"
    for extra in [
        ("--panel-columns", "0"),
        ("--infile", str(t)),
        ("--layout-report", str(t)),
        ("--view", "embedded"),
    ]:
        with pytest.raises(ValueError):
            draw(manifest, output, *extra)
    manifest.write_text("title\tinfile\n")
    with pytest.raises(ValueError, match="no panels"):
        draw(manifest, output)


def test_audit_protects_and_records_panel_inputs(tmp_path):
    import json

    from nwkit.cli import main

    tree = tmp_path / "input.nwk"
    tree.write_text("(A:1,B:1);")
    manifest = tmp_path / "grid.tsv"
    manifest.write_text("title\tinfile\nTest\tinput.nwk\n")
    output = tmp_path / "figure.svg"
    command = ["draw", "--panel-table", str(manifest), "-o", str(output), "--audit"]
    for path in (tree, manifest):
        with pytest.raises(ValueError):
            main([*command, str(path)])
    assert not output.exists()
    audit = tmp_path / "audit.json"
    main([*command, str(audit)])
    data = json.loads(audit.read_text())
    assert {str(tree.resolve()), str(manifest.resolve())} <= {
        r["path"] for r in data["inputs"]
    }


@pytest.mark.parametrize("scope,counts", [("branch", [1, 1, 3]), ("clade", [4, 1, 1])])
def test_nested_clade_colours_preserve_stroke_width(tmp_path, scope, counts):
    (tmp_path / "t.nwk").write_text(
        "(((A:1,B:1[&&NHX:fg=inner]):1,C:1):1[&&NHX:fg=outer],D:1);"
    )
    manifest = tmp_path / "m.tsv"
    manifest.write_text("title\tinfile\nNested\tt.nwk\n")
    output = tmp_path / "panel.svg"
    draw(
        manifest,
        output,
        "--legend",
        "no",
        "--branch-color-property",
        "fg",
        "--branch-color-scope",
        scope,
        "--property-color",
        "outer=#a8381b",
        "--property-color",
        "inner=#aabbcc",
        "--terminal-branch-color",
        "#010203",
        "--branch-width",
        "0.6",
    )
    paths = [
        n.get("style", "") for n in ET.parse(output).iter() if n.tag.endswith("}path")
    ]
    colours = ["#a8381b", "#aabbcc", "#010203"]
    # Four outer edges (stem, internal edge, A, C); B overrides and D is unmarked.
    for colour, count in zip(colours, counts, strict=True):
        strokes = [p for p in paths if f"stroke: {colour}" in p]
        assert len(strokes) == count
        assert all("stroke-width: 0.6;" in p for p in strokes)


def test_invalid_xml_panel_title_preserves_output(tmp_path):
    tree = tmp_path / "t.nwk"
    tree.write_text("(A,B);")
    manifest = tmp_path / "m.tsv"
    manifest.write_text("title\tinfile\nA\x01B\tt.nwk\n")
    out = tmp_path / "out.svg"
    out.write_text("previous")
    with pytest.raises(ValueError, match="Invalid XML character"):
        draw(manifest, out)
    assert out.read_text() == "previous"

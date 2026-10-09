"""Offline fossil import, placement, consumer round trips and failure recovery."""

import csv
import hashlib
import io
import json
import os
import shutil
import struct
import subprocess
from argparse import Namespace
from contextlib import contextmanager
from pathlib import Path

import pytest
import requests
from ete4 import Tree

import nwkit.angiocal as source
import nwkit.angiocal_constraints as placement
from nwkit.cli import main, parser
from nwkit.mcmctree import mcmctree_main
from nwkit.time_tree import annotate_mcmctree_calibrations, parse_mcmctree_calibration
from nwkit.util import read_tree
from nwkit.xls_compound import workbook_stream
from nwkit.xls_reader import read_xls

HEADERS = "fossil_id\tfossil_taxon\tminimum_age_ma\tplacement\tclade\n"


def dataset_file(tmp_path, text="1\tSynthetic fossil\t12.5\tcrown\tTestaceae\n"):
    path = tmp_path / "fossils.tsv"
    path.write_text(HEADERS + text, encoding="utf-8")
    return path


def args_for(tmp_path, dataset, **kwargs):
    args = parser.parse_args(
        [
            "mcmctree",
            "-i",
            "((A_a:1,B_b:1)Testaceae:1,C_c:2);",
            "--angiocal",
            "v1.0",
            "--angiocal-file",
            str(dataset),
            "--angiocal-taxonomy",
            "no",
            "--report",
            str(tmp_path / "report.tsv"),
            "-o",
            str(tmp_path / "result.tre"),
        ]
    )
    for key, value in kwargs.items():
        setattr(args, key, value)
    return args


def report_rows(path):
    with open(path, encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_real_xls_fixture_preserves_original_metadata():
    path = Path(__file__).parent / "data" / "angiocal-v1.0-synthetic.xls"
    dataset = source.read_angiocal(path)
    (record,) = dataset.records
    assert record.fossil_id == "1"
    assert record.minimum_age_ma == 12.5
    assert record.placement == "crown"
    assert record.clade == "Testaceae"
    assert record.source_row == 3
    assert record.relationship_reference == "synthetic metadata"
    assert dataset.sha256 == hashlib.sha256(path.read_bytes()).hexdigest()


def test_xls_format_is_detected_from_content(tmp_path):
    fixture = Path(__file__).parent / "data" / "angiocal-v1.0-synthetic.xls"
    path = tmp_path / "download-without-extension"
    path.write_bytes(fixture.read_bytes())
    assert source.read_angiocal(path).records == source.read_angiocal(fixture).records


def test_xls_symlink_does_not_change_source_format(tmp_path):
    fixture = Path(__file__).parent / "data" / "angiocal-v1.0-synthetic.xls"
    path = tmp_path / "source-data"
    path.write_bytes(fixture.read_bytes())
    link = tmp_path / "source.xls"
    try:
        link.symlink_to(path)
    except OSError:
        pytest.skip("Symlink creation is unavailable on this platform")
    dataset = source.read_angiocal(link)
    assert dataset.path == str(path.resolve())
    assert dataset.records == source.read_angiocal(fixture).records


def xls_with_invalid_cell(tmp_path, column, kind):
    fixture = Path(__file__).parent / "data" / "angiocal-v1.0-synthetic.xls"
    data = fixture.read_bytes()
    stream = workbook_stream(data)
    base, size = data.index(stream), len(stream)
    assert data.count(stream) == 1  # The fixture has one contiguous workbook stream.
    changed = bytearray(data)
    position = base
    while position < base + size:
        code, length = struct.unpack_from("<HH", data, position)
        payload = position + 4
        if kind == "date" and code == 0x00E0:  # First XF: built-in date format 14.
            struct.pack_into("<H", changed, payload + 2, 14)
            kind = "date_style_written"
        if code == 0x027E and struct.unpack_from("<HH", data, payload) == (2, column):
            if kind == "date_style_written":
                struct.pack_into("<H", changed, payload + 4, 0)
            else:  # BOOLERR has the same cell prefix; preserve the spare bytes.
                struct.pack_into("<H", changed, position, 0x0205)
                changed[payload + 6 : payload + 8] = (
                    b"\x07\x01" if kind == "error" else b"\x01\x00"
                )
        position += 4 + length
    path = tmp_path / "invalid-cell.xls"
    path.write_bytes(changed)
    return path


@pytest.mark.parametrize("column", [0, 2])
@pytest.mark.parametrize("kind", ["bool", "error", "date"])
def test_xls_numeric_cells_do_not_accept_boolean_date_or_error(tmp_path, column, kind):
    path = xls_with_invalid_cell(tmp_path, column, kind)
    expected_type = {
        "bool": "boolean",
        "error": "error",
        "date": "date",
    }
    sheet = read_xls(path.read_bytes())[0]
    assert sheet.cell(2, column).kind == expected_type[kind]
    with pytest.raises(
        ValueError, match="Excel error|numeric or decimal-text Excel cells"
    ):
        source.read_angiocal(path)


def test_numeric_zero_quality_metadata_is_preserved():
    # XLS numeric cells must not be mistaken for absent metadata.
    record = source._record(
        {
            "fossil_id": 1.0,
            "fossil_taxon": "Synthetic fossil",
            "minimum_age_ma": 12.5,
            "placement": "crown",
            "clade": "Testaceae",
            "age_quality_score": 0.0,
            "node_assignment_score": 0,
        },
        3,
    )
    assert record.age_quality_score == "0.0"
    assert record.node_assignment_score == "0"
    assert record.reconciliation_score == ""


@pytest.mark.parametrize(
    "change,error",
    [
        ("1\tF\tnan\tcrown\tTestaceae\n", "finite and positive"),
        ("1\tF\t0\tcrown\tTestaceae\n", "finite and positive"),
        ("1\tF\t-1\tcrown\tTestaceae\n", "finite and positive"),
        ("1\tF\tx\tcrown\tTestaceae\n", "numeric"),
        ("1.1\tF\t1\tcrown\tTestaceae\n", "positive integers"),
        ("1.00000000000000001\tF\t1\tcrown\tTestaceae\n", "positive integers"),
        ("9007199254740992.1\tF\t1\tcrown\tTestaceae\n", "positive integers"),
        ("1\tF\t1\ttip\tTestaceae\n", "placement"),
        ("1\tF\t1\tcrown\t\n", "Empty"),
        ("1\tF\t1\tcrown\tTestaceae\textra\n", "Malformed"),
        ("1\tF\t1\tcrown\n", "Malformed"),
        ("1\tF\t1\tcrown\tTestaceae\n1\tG\t2\tcrown\tTestaceae\n", "Duplicate"),
    ],
)
def test_invalid_source_rows(tmp_path, change, error):
    with pytest.raises(ValueError, match=error):
        source.read_angiocal(dataset_file(tmp_path, change))


@pytest.mark.parametrize(
    "header", ["fossil_id\tfossil_id\n", "fossil_id\t\n", "fossil_id\n"]
)
def test_invalid_headers(tmp_path, header):
    path = tmp_path / "bad.tsv"
    path.write_text(header, encoding="utf-8")
    with pytest.raises(ValueError, match="headers|columns"):
        source.read_angiocal(path)


def test_empty_and_oversized_input(tmp_path, monkeypatch):
    path = dataset_file(tmp_path, "")
    with pytest.raises(ValueError, match="no calibration"):
        source.read_angiocal(path)
    monkeypatch.setattr(source, "ANGIOCAL_MAX_BYTES", 10)
    with pytest.raises(ValueError, match="size limit"):
        source.read_angiocal(path)


def test_bom_and_fractional_age_round_trip(tmp_path):
    path = dataset_file(tmp_path, "1\tF\t12.3456789012345\tcrown\tTestaceae\n")
    path.write_bytes(b"\xef\xbb\xbf" + path.read_bytes())
    args = args_for(tmp_path, path, time_unit_ma=100, add_header=True)
    mcmctree_main(args)
    tree = read_tree(args.outfile, "auto", True)
    assert tree.common_ancestor(["A_a", "B_b"]).name.startswith("L(")
    calibration = parse_mcmctree_calibration(tree.common_ancestor(["A_a", "B_b"]).name)
    assert calibration["lower"] == 12.3456789012345 / 100
    assert ":" not in Path(args.outfile).read_text()
    assert annotate_mcmctree_calibrations(tree) == 1
    (row,) = report_rows(args.report)
    assert float(row["minimum_age_ma"]) == 12.3456789012345
    assert float(row["time_unit_ma"]) == 100
    assert row["status"] == "applied"


def test_large_fossil_ids_remain_distinct_through_explicit_mapping(tmp_path):
    path = dataset_file(
        tmp_path,
        "9007199254740992\tF\t10\tcrown\tUnknown\n"
        "9007199254740993\tG\t20\tstem\tUnknown\n",
    )
    mapping = tmp_path / "map.tsv"
    mapping.write_text(
        "fossil_id\tleft_species\tright_species\n"
        "9007199254740992\tA_a\tB_b\n9007199254740993\tA_a\tC_c\n"
    )
    args = args_for(tmp_path, path, calibration_map_tsv=str(mapping))
    mcmctree_main(args)
    rows = report_rows(args.report)
    assert [row["fossil_id"] for row in rows] == [
        "9007199254740992",
        "9007199254740993",
    ]
    tree = read_tree(args.outfile, "auto", True)
    assert parse_mcmctree_calibration(tree.name)["lower"] == 20
    assert parse_mcmctree_calibration(tree.children[0].name)["lower"] == 10


@pytest.mark.parametrize("tail", [0.9999999, 0.0123456789012345])
def test_prior_parameters_survive_export_without_rounding(tmp_path, tail):
    args = args_for(
        tmp_path,
        dataset_file(tmp_path),
        lower_tail_prob=tail,
        lower_offset=0.123456789012345,
        lower_scale=1.23456789012345,
    )
    mcmctree_main(args)
    tree = read_tree(args.outfile, "auto", True)
    calibration = parse_mcmctree_calibration(tree.children[0].name)
    assert calibration["lower_tail"] == tail
    assert calibration["offset"] == args.lower_offset
    assert calibration["scale"] == args.lower_scale
    assert report_rows(args.report)[0]["constraint"] == calibration["raw"]


def test_conflicting_source_placement_fails_before_output(tmp_path):
    path = tmp_path / "conflicting.tsv"
    path.write_text(
        HEADERS.strip() + "\tnode_calibrated\n"
        "1\tF\t12.5\tstem\tTestaceae\tcrown Testaceae\n"
    )
    args = args_for(tmp_path, path)
    with pytest.raises(ValueError, match="Conflicting.*crown/stem"):
        mcmctree_main(args)
    assert not Path(args.outfile).exists()


@pytest.mark.parametrize(
    "tree,error",
    [
        ("((A_a,B_b,C_c)Testaceae,D_d);", "fully bifurcating"),
        ("(((A_a),B_b)Testaceae,C_c);", "fully bifurcating"),
        ("((A_a,A_a)Testaceae,C_c);", "unique leaf names"),
    ],
)
def test_invalid_paml_topology_fails_before_loading_source(
    tmp_path, monkeypatch, tree, error
):
    args = args_for(tmp_path, dataset_file(tmp_path), infile=tree)
    monkeypatch.setattr(
        source, "load_angiocal", lambda args: pytest.fail("source must not load")
    )
    with pytest.raises(ValueError, match=error):
        mcmctree_main(args)
    assert not Path(args.outfile).exists()
    assert not Path(args.report).exists()


def test_map_header_cannot_have_an_empty_column(tmp_path):
    mapping = tmp_path / "map.tsv"
    mapping.write_text("fossil_id\tleft_species\tright_species\t\n1\tA_a\tB_b\t\n")
    args = args_for(tmp_path, dataset_file(tmp_path), calibration_map_tsv=str(mapping))
    with pytest.raises(ValueError, match="columns"):
        mcmctree_main(args)


@pytest.mark.parametrize("tip", ["NoName_species", "A_a-12"])
def test_fossil_export_preserves_tip_tokens(tmp_path, tip):
    args = args_for(
        tmp_path,
        dataset_file(tmp_path),
        infile=f"(('{tip}':1,B_b:1)Testaceae:1,C_c:2);",
    )
    mcmctree_main(args)
    result = read_tree(args.outfile, "auto", True)
    assert set(result.leaf_names()) == {tip, "B_b", "C_c"}


@pytest.mark.parametrize(
    "tip",
    [
        "A_a:12",
        "A a",
        "A_a#1",
        "Taxon_あ",
        "A_" + "a" * 99,
        "1",
        "01",
        "123456",
        "A_a\x00",
        "A_a\x07",
        "A_a\x7f",
    ],
)
def test_fossil_export_rejects_names_that_paml_cannot_preserve(tmp_path, tip):
    args = args_for(
        tmp_path,
        dataset_file(tmp_path),
        infile=f"(('{tip}',B_b)Testaceae,C_c);",
    )
    with pytest.raises(ValueError, match="tip identifiers"):
        mcmctree_main(args)
    assert not Path(args.outfile).exists()
    assert not Path(args.report).exists()


@pytest.mark.parametrize("tip", ["NoName_species", "A_a:12"])
def test_manual_export_preserves_tip_tokens(tmp_path, tip):
    args = parser.parse_args(
        [
            "mcmctree",
            "-i",
            f"(('{tip}':1,B_b:1):1,C_c:2);",
            "--left-species",
            tip,
            "--right-species",
            "B_b",
            "--lower-bound",
            "10",
            "-o",
            str(tmp_path / "manual.tre"),
        ]
    )
    mcmctree_main(args)
    result = read_tree(args.outfile, "auto", True)
    assert set(result.leaf_names()) == {tip, "B_b", "C_c"}


def test_same_node_fossils_select_oldest_and_keep_support(tmp_path):
    path = dataset_file(
        tmp_path, "1\tF\t10\tcrown\tTestaceae\n2\tG\t15\tcrown\tTestaceae\n"
    )
    args = args_for(tmp_path, path)
    mcmctree_main(args)
    rows = report_rows(args.report)
    assert [row["status"] for row in rows] == ["supporting", "applied"]
    assert {row["selected_minimum_age_ma"] for row in rows} == {"15.0"}
    assert len({row["node_clade_id"] for row in rows}) == 1
    assert all(
        parse_mcmctree_calibration(row["constraint"])["lower"] == 15 for row in rows
    )


def test_crown_is_not_guessed_from_sampled_mrca(tmp_path):
    path = dataset_file(tmp_path)
    args = args_for(tmp_path, path, infile="((A_a,B_b),C_c);")
    Path(args.outfile).write_text("old tree", encoding="utf-8")
    with pytest.raises(ValueError, match="No AngioCal"):
        mcmctree_main(args)
    assert Path(args.outfile).read_text() == "old tree"
    (row,) = report_rows(args.report)
    assert row["reason"] == "crown_requires_anchors"


def test_explicit_map_uses_final_mrca_for_stem_and_crown(tmp_path):
    path = dataset_file(
        tmp_path, "1\tF\t10\tstem\tUnknown clade\n2\tG\t5\tcrown\tUnknown crown\n"
    )
    mapping = tmp_path / "map.tsv"
    mapping.write_text(
        "fossil_id\tleft_species\tright_species\n1\tA_a\tC_c\n2\tA_a\tB_b\n",
        encoding="utf-8",
    )
    args = args_for(
        tmp_path, path, calibration_map_tsv=str(mapping), infile="((A_a,B_b),C_c);"
    )
    mcmctree_main(args)
    tree = read_tree(args.outfile, "auto", True)
    assert parse_mcmctree_calibration(tree.name)["lower"] == 10
    assert (
        parse_mcmctree_calibration(tree.common_ancestor(["A_a", "B_b"]).name)["lower"]
        == 5
    )
    assert {row["mapping_method"] for row in report_rows(args.report)} == {
        "explicit_map"
    }


@pytest.mark.parametrize(
    "mapping,error",
    [
        ("1\tA_a\tA_a\n", "different lineages"),
        ("1\tA_a\tMissing_sp\n", "not found"),
        ("9\tA_a\tB_b\n", "Unknown or duplicate"),
        ("1\tA_a\tB_b\n1\tA_a\tC_c\n", "Unknown or duplicate"),
        ("1\tA_a\n", "Malformed"),
    ],
)
def test_invalid_explicit_maps(tmp_path, mapping, error):
    path = dataset_file(tmp_path)
    m = tmp_path / "map.tsv"
    m.write_text("fossil_id\tleft_species\tright_species\n" + mapping, encoding="utf-8")
    args = args_for(tmp_path, path, calibration_map_tsv=str(m))
    with pytest.raises(ValueError, match=error):
        mcmctree_main(args)
    assert not Path(args.outfile).exists()
    assert not Path(args.report).exists()


def test_explicit_map_stdin_and_cli_ownership(tmp_path, monkeypatch):
    path = dataset_file(tmp_path)
    monkeypatch.setattr(
        "sys.stdin",
        io.StringIO("fossil_id\tleft_species\tright_species\n1\tA_a\tB_b\n"),
    )
    argv = [
        "mcmctree",
        "-i",
        "((A_a,B_b),C_c);",
        "--angiocal",
        "v1.0",
        "--angiocal-file",
        str(path),
        "--angiocal-taxonomy",
        "no",
        "--calibration-map-tsv",
        "-",
    ]
    main(argv)
    with pytest.raises(ValueError, match="STDIN"):
        main(argv[:1] + ["-i", "-"] + argv[3:])


@pytest.mark.parametrize(
    "changes",
    [
        {"timetree": "ci"},
        {"timetree": "point"},
        {"posterior": "posterior.txt"},
        {"left_species": "A_a"},
        {"lower_bound": "10"},
        {"report": "-"},
    ],
)
def test_incompatible_options_fail_before_input_or_network(
    tmp_path, changes, monkeypatch
):
    args = args_for(tmp_path, "missing.tsv", infile="missing.nwk", **changes)
    monkeypatch.setattr(
        source, "load_angiocal", lambda args: pytest.fail("source should not load")
    )
    with pytest.raises(ValueError, match="cannot be combined|file path"):
        mcmctree_main(args)


@pytest.mark.parametrize(
    "option", ["--angiocal-file", "--calibration-map-tsv", "--report"]
)
def test_auxiliary_options_require_angiocal(option):
    args = parser.parse_args(["mcmctree", option, "unused.tsv"])
    with pytest.raises(ValueError, match="require.*angiocal"):
        mcmctree_main(args)


@pytest.mark.parametrize("unit", [0, -1, float("inf"), float("nan")])
def test_invalid_units(tmp_path, unit):
    args = args_for(tmp_path, dataset_file(tmp_path), time_unit_ma=unit)
    with pytest.raises(ValueError, match="finite and positive"):
        mcmctree_main(args)


@pytest.mark.parametrize(
    "node_name,reason",
    [
        ("Testaceae", "applied"),
        ("'crown Testaceae'", "applied"),
    ],
)
def test_explicit_biological_node_labels(tmp_path, node_name, reason):
    args = args_for(
        tmp_path, dataset_file(tmp_path), infile=f"((A_a,B_b){node_name},C_c);"
    )
    mcmctree_main(args)
    assert report_rows(args.report)[0]["status"] == reason


def test_stem_named_clade_maps_to_parent_and_root_stem_is_skipped(tmp_path):
    path = dataset_file(tmp_path, "1\tF\t10\tstem\tTestaceae\n2\tG\t20\tstem\tAll\n")
    args = args_for(tmp_path, path, infile="((A_a,B_b)Testaceae,C_c)All;")
    mcmctree_main(args)
    rows = report_rows(args.report)
    assert rows[0]["node_tips"] == '["A_a","B_b","C_c"]'
    assert rows[1]["reason"] == "stem_outgroup_missing"


def test_tip_label_is_not_an_assertion_of_an_internal_clade(tmp_path):
    path = dataset_file(tmp_path, "1\tF\t20\tstem\tA_a\n2\tG\t12.5\tcrown\tTestaceae\n")
    args = args_for(tmp_path, path)
    mcmctree_main(args)
    rows = report_rows(args.report)
    assert rows[0]["status"] == "skipped"
    assert rows[0]["reason"] == "taxonomy_disabled"
    assert (
        parse_mcmctree_calibration(read_tree(args.outfile, "auto", True).name) is None
    )


def test_internal_clade_assertion_is_independent_of_tip_name(tmp_path):
    path = dataset_file(tmp_path, "1\tF\t20\tstem\tTestaceae\n")
    args = args_for(
        tmp_path,
        path,
        infile="((Testaceae,B_b)Testaceae,C_c);",
        species_regex="(.*)",
    )
    mcmctree_main(args)
    assert report_rows(args.report)[0]["node_tips"] == '["B_b","C_c","Testaceae"]'


def test_lowercase_existing_priors_are_exported_in_paml_case(tmp_path):
    path = dataset_file(tmp_path)
    mapping = tmp_path / "map.tsv"
    mapping.write_text("fossil_id\tleft_species\tright_species\n1\tA_a\tB_b\n")
    args = args_for(
        tmp_path,
        path,
        infile="((A_a,B_b)'l(20,0.1,1,1e-300)',C_c)'u(30,0.025)';",
        calibration_map_tsv=str(mapping),
    )
    mcmctree_main(args)
    tree = read_tree(args.outfile, "auto", True)
    assert tree.name == "U(30,0.025)"
    assert tree.children[0].name == "L(20,0.1,1,1e-300)"
    assert report_rows(args.report)[0]["constraint"] == tree.children[0].name


@pytest.mark.parametrize(
    "prior",
    [
        "G(20,2)",
        "SN(20,2,1)",
        "ST(20,2,1,3)",
        "S2N(0.5,20,2,1,25,2,1)",
        "U(30)",
        "L(20,0.1,1)",
        "B(10,30)",
        "L(nan,0.1,1,0.025)",
        "#1",
        "L(20,0.1,1,0.025) #1",
        "@unknown",
        "root U(30,0.025)",
        "calibrated G(20,2)",
    ],
)
def test_unsupported_existing_priors_are_not_silently_discarded(
    tmp_path, monkeypatch, prior
):
    args = args_for(
        tmp_path, dataset_file(tmp_path), infile=f"((A_a,B_b)Testaceae,C_c)'{prior}';"
    )
    Path(args.outfile).write_text("previous tree")
    Path(args.report).write_text("previous report")
    monkeypatch.setattr(
        source, "load_angiocal", lambda args: pytest.fail("source must not load")
    )
    with pytest.raises(ValueError, match="Unsupported.*calibration"):
        mcmctree_main(args)
    assert Path(args.outfile).read_text() == "previous tree"
    assert Path(args.report).read_text() == "previous report"


def test_explicit_map_rerun_preserves_tree_units_and_report(tmp_path):
    path = dataset_file(tmp_path)
    mapping = tmp_path / "map.tsv"
    mapping.write_text("fossil_id\tleft_species\tright_species\n1\tA_a\tB_b\n")
    args = args_for(
        tmp_path,
        path,
        infile="((A_a,B_b)Testaceae,C_c)'U(0.3,0.025)';",
        calibration_map_tsv=str(mapping),
        time_unit_ma=100,
        add_header=True,
    )
    mcmctree_main(args)
    first_tree = Path(args.outfile).read_text()
    first_report = report_rows(args.report)
    args.infile = args.outfile
    args.outfile = str(tmp_path / "rerun.tre")
    mcmctree_main(args)
    assert Path(args.outfile).read_text() == first_tree
    assert report_rows(args.report) == first_report


def test_no_crown_with_only_duplicate_species(tmp_path):
    path = dataset_file(tmp_path)
    args = args_for(tmp_path, path, infile="((A_a_1,A_a_2)Testaceae,C_c);")
    with pytest.raises(ValueError, match="No AngioCal"):
        mcmctree_main(args)
    assert report_rows(args.report)[0]["reason"] == "crown_not_represented"


@pytest.mark.parametrize(
    "label",
    ["'U(20, 0.01)'", "'U(20, 0.0123456789012345)'", "'L(20, 0.2, 2, 0.01)'"],
)
def test_compatible_existing_manual_calibrations(tmp_path, label):
    path = dataset_file(tmp_path)
    mapping = tmp_path / "map.tsv"
    mapping.write_text(
        "fossil_id\tleft_species\tright_species\n1\tA_a\tB_b\n", encoding="utf-8"
    )
    args = args_for(
        tmp_path,
        path,
        infile=f"((A_a,B_b){label},C_c);",
        calibration_map_tsv=str(mapping),
    )
    mcmctree_main(args)
    calibration = parse_mcmctree_calibration(
        read_tree(args.outfile, "auto", True).common_ancestor(["A_a", "B_b"]).name
    )
    if label.startswith("'U"):
        assert calibration["lower"] == 12.5 and calibration["upper"] == 20
        assert (
            calibration["upper_tail"] == parse_mcmctree_calibration(label)["upper_tail"]
        )
    else:
        assert calibration == parse_mcmctree_calibration(label)


@pytest.mark.parametrize(
    "label,error",
    [
        ("'U(10, 0.025)'", "younger"),
        ("'@10'", "does not treat them as fossil priors"),
        ("'L(10, 0.1, 1, 0.025)'", "separate prior"),
    ],
)
def test_conflicting_existing_calibrations_preserve_outputs(tmp_path, label, error):
    path = dataset_file(tmp_path, "1\tF\t12.5\tstem\tTestaceae\n")
    args = args_for(tmp_path, path, infile=f"((A_a,B_b)Testaceae,C_c){label};")
    Path(args.outfile).write_text("old tree", encoding="utf-8")
    Path(args.report).write_text("old report", encoding="utf-8")
    with pytest.raises(ValueError, match=error):
        mcmctree_main(args)
    assert Path(args.outfile).read_text() == "old tree"
    assert Path(args.report).read_text() == "old report"


def test_ancestor_upper_bound_cannot_contradict_descendant_fossil(tmp_path):
    path = dataset_file(tmp_path)
    args = args_for(tmp_path, path, infile="((A_a,B_b)Testaceae,C_c)'U(10, 0.025)';")
    with pytest.raises(ValueError, match="ancestor and its descendants"):
        mcmctree_main(args)
    assert not Path(args.outfile).exists()


def test_constraint_order_does_not_interpret_tip_names_as_calibrations():
    tree = Tree("(('@100',A_a),B_b)'U(20, 0.025)';", parser=1)
    placement.validate_calibration_order(tree)


@pytest.mark.parametrize(
    "label,lower_tail,error",
    [
        ("U(20, 0.6)", 0.6, "sum to less than 1"),
        ("U(20, 0.5)", 0.5, "sum to less than 1"),
        ("B(12.5, 12.5, 0.025, 0.025)", 0.025, "positive age range"),
        ("B(15, 20, 0.6, 0.5)", 0.025, "sum to less than 1"),
        ("@20", 0.025, "does not treat them as fossil priors"),
    ],
)
def test_invalid_existing_or_combined_fossil_priors_fail_without_output(
    tmp_path, label, lower_tail, error
):
    mapping = tmp_path / "map.tsv"
    mapping.write_text("fossil_id\tleft_species\tright_species\n1\tA_a\tB_b\n")
    args = args_for(
        tmp_path,
        dataset_file(tmp_path),
        infile=f"((A_a,B_b)'{label}',C_c)'U(30, 0.025)';",
        calibration_map_tsv=str(mapping),
        lower_tail_prob=lower_tail,
    )
    Path(args.outfile).write_text("old tree")
    Path(args.report).write_text("old report")
    with pytest.raises(ValueError, match=error):
        mcmctree_main(args)
    assert Path(args.outfile).read_text() == "old tree"
    assert Path(args.report).read_text() == "old report"


def test_point_annotation_on_unmapped_ancestor_is_also_rejected(tmp_path):
    args = args_for(
        tmp_path,
        dataset_file(tmp_path),
        infile="((A_a,B_b)Testaceae,C_c)'@30';",
    )
    with pytest.raises(ValueError, match="does not treat them as fossil priors"):
        mcmctree_main(args)


def test_zero_lower_prior_on_unmapped_ancestor_is_rejected(tmp_path):
    args = args_for(
        tmp_path,
        dataset_file(tmp_path),
        infile="((A_a,B_b)Testaceae,C_c)'L(0, 0.1, 1, 0.025)';",
    )
    with pytest.raises(ValueError, match="finite and positive"):
        mcmctree_main(args)


@pytest.mark.parametrize("parameter", ["lower_scale", "lower_offset"])
def test_finite_parameters_cannot_overflow_derived_lower_prior(tmp_path, parameter):
    args = args_for(tmp_path, dataset_file(tmp_path), **{parameter: 1e308})
    with pytest.raises(ValueError, match="finite and positive"):
        mcmctree_main(args)
    assert not Path(args.outfile).exists()
    assert not Path(args.report).exists()


def test_equal_ancestor_upper_and_descendant_minimum_has_no_positive_time_edge(
    tmp_path,
):
    args = args_for(
        tmp_path,
        dataset_file(tmp_path),
        infile="((A_a,B_b)Testaceae,C_c)'U(12.5, 0.025)';",
    )
    with pytest.raises(ValueError, match="ancestor and its descendants"):
        mcmctree_main(args)


def test_equal_upper_and_fossil_minimum_does_not_create_degenerate_bounded_prior(
    tmp_path,
):
    path = dataset_file(tmp_path, "1\tF\t12.5\tstem\tTestaceae\n")
    args = args_for(tmp_path, path, infile="((A_a,B_b)Testaceae,C_c)'U(12.5, 0.025)';")
    with pytest.raises(ValueError, match="must exceed"):
        mcmctree_main(args)


@pytest.mark.parametrize("unit,age", [(1e308, 1e-100), (1e-300, 1e100)])
def test_unit_conversion_never_silently_overflows_or_underflows(tmp_path, unit, age):
    args = args_for(
        tmp_path,
        dataset_file(tmp_path, f"1\tF\t{age}\tcrown\tTestaceae\n"),
        time_unit_ma=unit,
    )
    with pytest.raises(ValueError, match="finite and positive"):
        mcmctree_main(args)


def test_min_clade_filter_is_reported(tmp_path):
    path = dataset_file(
        tmp_path, "1\tF\t10\tcrown\tTestaceae\n2\tG\t20\tstem\tTestaceae\n"
    )
    args = args_for(tmp_path, path, min_clade_prop=0.9)
    mcmctree_main(args)
    assert report_rows(args.report)[0]["reason"] == "min_clade_prop"
    assert (
        parse_mcmctree_calibration(read_tree(args.outfile, "auto", True).name)["lower"]
        == 20
    )


@pytest.mark.parametrize("collision", ["source", "tree", "report", "map"])
def test_output_collisions(tmp_path, collision):
    path = dataset_file(tmp_path)
    tree_path = tmp_path / "tree.nwk"
    tree_path.write_text("((A_a,B_b)Testaceae,C_c);", encoding="utf-8")
    mapping = tmp_path / "map.tsv"
    mapping.write_text(
        "fossil_id\tleft_species\tright_species\n1\tA_a\tB_b\n", encoding="utf-8"
    )
    args = args_for(
        tmp_path, path, infile=str(tree_path), calibration_map_tsv=str(mapping)
    )
    args.outfile = {
        "source": str(path),
        "tree": str(tree_path),
        "report": args.report,
        "map": str(mapping),
    }[collision]
    originals = {p: p.read_bytes() for p in (path, tree_path, mapping)}
    with pytest.raises(ValueError, match="overwrite|distinct"):
        mcmctree_main(args)
    assert all(p.read_bytes() == data for p, data in originals.items())


def test_tree_and_report_rollback_on_commit_failure(tmp_path, monkeypatch):
    import nwkit.output_transaction as transaction

    args = args_for(tmp_path, dataset_file(tmp_path))
    Path(args.outfile).write_text("old tree", encoding="utf-8")
    Path(args.report).write_text("old report", encoding="utf-8")
    replace = transaction.replace_output

    def fail_second(source, target):
        if str(target) == args.report and ".stage." in str(source):
            raise OSError("injected installation failure")
        return replace(source, target)

    monkeypatch.setattr(transaction, "replace_output", fail_second)
    with pytest.raises(OSError, match="injected"):
        mcmctree_main(args)
    assert Path(args.outfile).read_text() == "old tree"
    assert Path(args.report).read_text() == "old report"


class DummyTaxonomy:
    db = None

    def get_name_translator(self, names):
        known = {
            "A a": [11],
            "B b": [12],
            "C c": [13],
            "Testaceae": [100],
            "Otheraceae": [200],
        }
        return {name: known[name] for name in names if name in known}

    def get_lineage(self, taxid):
        return {11: [1, 100, 11], 12: [1, 100, 12], 13: [1, 200, 13]}[taxid]


@pytest.mark.parametrize(
    "tree_text,clade,reason,target_size",
    [
        ("((A_a,B_b),C_c);", "Testaceae", "sampled_stem_ancestor", 3),
        ("((A_a,C_c),B_b);", "Testaceae", "nonmonophyletic_clade", 0),
        ("((A_a,B_b),C_c);", "Testaceae+Otheraceae", "stem_outgroup_missing", 0),
        ("((A_a,B_b),C_c);", "Unknown clade", "unresolved_clade_taxonomy", 0),
    ],
)
def test_taxonomy_stem_projection(
    tmp_path, monkeypatch, tree_text, clade, reason, target_size
):
    path = dataset_file(tmp_path, f"1\tF\t10\tstem\t{clade}\n")
    dataset = source.read_angiocal(path)
    args = args_for(tmp_path, path, angiocal_taxonomy=True)
    monkeypatch.setattr(placement, "get_ete_ncbitaxa", lambda args: DummyTaxonomy())
    placements, leaf_sets = placement.place_fossils(
        Tree(tree_text, parser=1), dataset, args
    )
    target, actual_reason, method = placements["1"]
    assert actual_reason == reason
    assert method == "taxonomy_stem"
    assert (len(leaf_sets[target]) if target is not None else 0) == target_size


def test_taxonomy_ambiguity_never_picks_arbitrary_taxid(tmp_path, monkeypatch):
    class Ambiguous(DummyTaxonomy):
        def get_name_translator(self, names):
            result = super().get_name_translator(names)
            if "A a" in names:
                result["A a"] = [11, 15]
            return result

    path = dataset_file(tmp_path, "1\tF\t10\tstem\tTestaceae\n")
    args = args_for(tmp_path, path, angiocal_taxonomy=True)
    monkeypatch.setattr(placement, "get_ete_ncbitaxa", lambda args: Ambiguous())
    placements, _ = placement.place_fossils(
        Tree("((A_a,B_b),C_c);", parser=1), source.read_angiocal(path), args
    )
    assert placements["1"][1] == "unresolved_tip_taxonomy"


def test_pinned_download_and_cache_reuse(tmp_path, monkeypatch):
    data = (Path(__file__).parent / "data" / "angiocal-v1.0-synthetic.xls").read_bytes()
    monkeypatch.setattr(source, "ANGIOCAL_SHA256", hashlib.sha256(data).hexdigest())
    calls = []

    class Response:
        def __enter__(self):
            return self

        def __exit__(self, *args):
            pass

        def raise_for_status(self):
            pass

        def iter_content(self, chunk_size):
            assert chunk_size == 64 * 1024
            yield data

    def get(url, **kwargs):
        calls.append((url, kwargs))
        return Response()

    monkeypatch.setattr(source.requests, "get", get)
    args = Namespace(angiocal="v1.0", angiocal_file=None, download_dir=str(tmp_path))
    first = source.load_angiocal(args)
    second = source.load_angiocal(args)
    assert first == second
    assert first.source == source.ANGIOCAL_URL
    assert len(calls) == 1
    assert calls[0][1] == {"stream": True, "timeout": 30}
    Path(first.path).write_bytes(b"bad cache")
    with pytest.raises(ValueError, match="checksum"):
        source.load_angiocal(args)
    assert len(calls) == 1


@pytest.mark.parametrize("role", ["--outfile", "--report", "--audit"])
@pytest.mark.parametrize("exists", [False, True])
def test_outputs_and_audit_cannot_overwrite_automatic_source_cache(
    tmp_path, monkeypatch, role, exists
):
    args = Namespace(angiocal_file=None, download_dir=str(tmp_path))
    cache = source.angiocal_input_path(args)
    original = b"protected source"
    if exists:
        cache.parent.mkdir(parents=True)
        cache.write_bytes(original)
    monkeypatch.setattr(
        source.requests, "get", lambda *args, **kwargs: pytest.fail("must not request")
    )
    with pytest.raises(ValueError, match="overwrite|distinct"):
        main(
            [
                "mcmctree",
                "-i",
                "((A_a,B_b)Testaceae,C_c);",
                "--angiocal",
                "v1.0",
                "--angiocal-taxonomy",
                "no",
                "--download-dir",
                str(tmp_path),
                role,
                str(cache),
            ]
        )
    assert cache.read_bytes() == original if exists else not cache.exists()


def test_first_download_is_included_in_audit_source_hashes(tmp_path, monkeypatch):
    data = (Path(__file__).parent / "data" / "angiocal-v1.0-synthetic.xls").read_bytes()
    digest = hashlib.sha256(data).hexdigest()
    monkeypatch.setattr(source, "ANGIOCAL_SHA256", digest)
    monkeypatch.setattr(source, "_download_angiocal", lambda: data)
    audit = tmp_path / "audit.jsonl"
    main(
        [
            "mcmctree",
            "-i",
            "((A_a,B_b)Testaceae,C_c);",
            "--angiocal",
            "v1.0",
            "--angiocal-taxonomy",
            "no",
            "--download-dir",
            str(tmp_path),
            "--audit",
            str(audit),
        ]
    )
    record = json.loads(audit.read_text())
    inputs = [item for item in record["inputs"] if item["argument"] == "angiocal_cache"]
    assert len(inputs) == 1
    assert inputs[0]["sha256"] == digest
    assert inputs[0]["bytes"] == len(data)


@pytest.mark.parametrize("role", ["--outfile", "--report", "--audit"])
def test_outputs_cannot_replace_ncbi_cache_assets(tmp_path, monkeypatch, role):
    cache = tmp_path / "ete4/taxa.sqlite"
    cache.parent.mkdir()
    cache.write_bytes(b"protected taxonomy")
    monkeypatch.setattr(
        source, "load_angiocal", lambda args: pytest.fail("source must not load")
    )
    with pytest.raises(ValueError, match="NCBI taxonomy cache"):
        main(
            [
                "mcmctree",
                "-i",
                "((A_a,B_b)Testaceae,C_c);",
                "--angiocal",
                "v1.0",
                "--download-dir",
                str(tmp_path),
                role,
                str(cache),
            ]
        )
    assert cache.read_bytes() == b"protected taxonomy"


@pytest.mark.parametrize(
    "replacement",
    [
        pytest.param(
            "symlink",
            marks=pytest.mark.skipif(
                os.name == "nt", reason="Windows symlinks may need privileges"
            ),
        ),
        "hardlink",
    ],
)
def test_cache_download_rejects_replaced_stage_before_writing(
    tmp_path, monkeypatch, replacement
):
    data = (Path(__file__).parent / "data" / "angiocal-v1.0-synthetic.xls").read_bytes()
    monkeypatch.setattr(source, "_download_angiocal", lambda: data)
    victim = tmp_path / "protected.bin"
    victim.write_bytes(b"untouched")
    original_transaction = source.output_transaction

    @contextmanager
    def replace_stage(paths, **kwargs):
        with original_transaction(paths, **kwargs) as staged:
            path = Path(staged[paths[0]])
            path.unlink()
            if replacement == "symlink":
                path.symlink_to(victim)
            else:
                path.hardlink_to(victim)
            yield staged

    monkeypatch.setattr(source, "output_transaction", replace_stage)
    args = Namespace(angiocal="v1.0", angiocal_file=None, download_dir=str(tmp_path))
    with pytest.raises((OSError, RuntimeError)):
        source.load_angiocal(args)
    assert victim.read_bytes() == b"untouched"
    assert not source.angiocal_input_path(args).exists()


@pytest.mark.parametrize("failure", ["checksum", "http", "size"])
def test_failed_download_never_installs_cache(tmp_path, monkeypatch, failure):
    class Response:
        def __enter__(self):
            return self

        def __exit__(self, *args):
            pass

        def raise_for_status(self):
            if failure == "http":
                raise requests.HTTPError("server error")

        def iter_content(self, chunk_size):
            yield b"invalid source"

    monkeypatch.setattr(source.requests, "get", lambda *args, **kwargs: Response())
    if failure == "size":
        monkeypatch.setattr(source, "ANGIOCAL_MAX_BYTES", 1)
    args = Namespace(angiocal="v1.0", angiocal_file=None, download_dir=str(tmp_path))
    with pytest.raises((ValueError, RuntimeError), match="checksum|download|size"):
        source.load_angiocal(args)
    assert not list(tmp_path.rglob("*.xls"))


@pytest.mark.integration
@pytest.mark.skipif(
    shutil.which("mcmctree") is None, reason="Optional PAML executable is not installed"
)
@pytest.mark.parametrize(
    "label,minimum",
    [
        ("Testaceae", 12.5),
        ("'l(20,0.1,1,1e-300)'", 20),
        ("'b(12.5,25,1e-300,0.025)'", 12.5),
    ],
)
def test_real_paml_consumes_angiocal_minimum_prior(tmp_path, label, minimum):
    mapping = tmp_path / "map.tsv"
    mapping.write_text("fossil_id\tleft_species\tright_species\n1\tA_a\tB_b\n")
    args = args_for(
        tmp_path,
        dataset_file(tmp_path),
        infile=f"((A_a,B_b){label},C_c)'u(30, 0.025)';",
        calibration_map_tsv=str(mapping),
        lower_tail_prob=1e-300,
        add_header=True,
    )
    mcmctree_main(args)
    (tmp_path / "input.phy").write_text(
        "3 20\nA_a  ACGTACGTACGTACGTACGT\n"
        "B_b  ACGTACGTACGTACGTACGT\nC_c  ACGTACGTACGTACGTACGT\n"
    )
    (tmp_path / "input.ctl").write_text(
        "seed = 1\nseqfile = input.phy\ntreefile = result.tre\noutfile = out.txt\n"
        "ndata = 1\nseqtype = 0\nusedata = 0\nclock = 1\nRootAge = <30\n"
        "BDparas = 1 1 0.1 M\nrgene_gamma = 2 20 1\nsigma2_gamma = 1 10 1\n"
        "burnin = 10\nsampfreq = 1\nnsample = 20\nprint = 1\n"
    )
    result = subprocess.run(
        [shutil.which("mcmctree"), "input.ctl"],
        cwd=tmp_path,
        stdin=subprocess.DEVNULL,
        capture_output=True,
        text=True,
        timeout=30,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    rows = report_rows(tmp_path / "mcmc.txt")
    assert len(rows) == 20
    # This short chain checks parsing and prior use, not MCMC convergence.
    assert all(float(row["t_n5"]) >= minimum for row in rows)
    assert all(float(row["t_n4"]) > float(row["t_n5"]) for row in rows)

import csv
import json
from pathlib import Path

import pytest
from ete4 import Tree

import nwkit.root_compare as root_compare_mod
from nwkit.cli import parser
from nwkit.root_compare import (
    evaluate_root_compare_methods,
    resolve_root_compare_methods,
    root_compare_main,
    root_compare_rows,
)
from nwkit.root_evaluation import RootingCandidate, RootingEvaluation
from tests.helpers import make_args


def _compare_args(**kwargs):
    defaults = {
        "methods": "all",
        "exclude_methods": None,
        "outgroup": None,
        "infile2": None,
        "format2": "auto",
        "taxon_mode": "exact",
        "species_tree": None,
        "species_tree_format": "auto",
        "duplication_cost": 1.0,
        "loss_cost": 1.0,
        "taxonomy_source": "ncbi,opentree,timetree",
        "taxid_tsv": None,
        "rank": "no",
        "figure_out": "comparison.pdf",
        "figure_width": None,
        "figure_height": None,
        "font_size": 8.0,
        "tip_labels": "auto",
        "unrooted_method": "auto",
        "layout_report": None,
    }
    defaults.update(kwargs)
    return make_args(**defaults)


def test_cli_defaults_to_all_and_does_not_expose_newick_output_format():
    args = parser.parse_args(["rootcompare", "--figure-out", "comparison.pdf"])

    assert args.methods == "all"
    assert args.figure_width is None
    assert args.unrooted_method == "auto"
    assert args.tip_labels == "auto"
    assert not hasattr(args, "outformat")


def test_all_adds_context_methods_and_each_taxonomy_source():
    tree = Tree(
        "(Homo_sapiens_g1:1,Pan_troglodytes_g1:1,Mus_musculus_g1:1);",
        parser=1,
    )
    args = _compare_args(
        outgroup="Mus_musculus_g1",
        infile2="reference.nwk",
        species_tree="species.nwk",
        taxonomy_source="opentree,ncbi,timetree",
    )

    methods, automatic = resolve_root_compare_methods(tree, args)

    assert automatic is True
    assert methods == [
        "midpoint",
        "mad",
        "mv",
        "outgroup",
        "transfer",
        "reconciliation",
        "taxonomy:opentree",
        "taxonomy:ncbi",
        "taxonomy:timetree",
    ]


def test_all_skips_taxonomy_when_species_names_cannot_be_parsed():
    tree = Tree("(A:1,B:1,C:1);", parser=1)

    methods, automatic = resolve_root_compare_methods(tree, _compare_args())

    assert automatic is True
    assert methods == ["midpoint", "mad", "mv"]


def test_taxid_table_enables_only_ncbi_when_names_cannot_be_parsed():
    tree = Tree("(A:1,B:1,C:1);", parser=1)

    methods, _ = resolve_root_compare_methods(
        tree,
        _compare_args(taxid_tsv="taxids.tsv"),
    )

    assert methods == ["midpoint", "mad", "mv", "taxonomy:ncbi"]


def test_taxonomy_exclusion_removes_every_source():
    tree = Tree(
        "(Homo_sapiens_g1:1,Pan_troglodytes_g1:1,Mus_musculus_g1:1);",
        parser=1,
    )

    methods, _ = resolve_root_compare_methods(
        tree,
        _compare_args(exclude_methods="taxonomy"),
    )

    assert methods == ["midpoint", "mad", "mv"]


def test_explicit_nontaxonomy_methods_do_not_validate_unused_taxonomy_sources():
    tree = Tree("(A:1,B:1,C:1);", parser=1)

    methods, automatic = resolve_root_compare_methods(
        tree,
        _compare_args(methods="midpoint", taxonomy_source="invalid"),
    )

    assert methods == ["midpoint"]
    assert automatic is False


def test_explicit_empty_method_list_is_rejected():
    tree = Tree("(A:1,B:1,C:1);", parser=1)

    with pytest.raises(ValueError, match="must name at least one method"):
        resolve_root_compare_methods(tree, _compare_args(methods=""))


def test_automatic_failures_are_reported_but_explicit_failures_raise():
    tree = Tree("(A:1,B:1);", parser=1)
    args = _compare_args()

    evaluations = evaluate_root_compare_methods(tree, ["mad"], args, automatic=True)

    assert evaluations[0].status == "failed"
    assert "at least 3 leaves" in evaluations[0].message
    with pytest.raises(ValueError, match="method 'mad' failed"):
        evaluate_root_compare_methods(tree, ["mad"], args, automatic=False)


def test_two_tip_mv_evaluation_reports_the_complete_root_edge():
    tree = Tree("(A:1,B:3):7;", parser=1)

    evaluations = evaluate_root_compare_methods(
        tree,
        ["mv"],
        _compare_args(),
        automatic=True,
    )

    evaluation = evaluations[0]
    assert evaluation.status == "ok"
    assert evaluation.evaluated_edges == 1
    assert len(evaluation.candidates) == 1
    candidate = evaluation.candidates[0]
    assert candidate.edge_length == pytest.approx(4.0)
    assert candidate.position_fraction_from_side_a == pytest.approx(0.5)
    assert candidate.score == pytest.approx(0.0)


def test_taxonomy_sources_are_evaluated_independently(monkeypatch):
    tree = Tree("(A:1,(B:1,C:1):1);", parser=1)
    calls = []

    def fake_taxonomy_rooting(tree, taxonomy_source, **kwargs):
        calls.append(taxonomy_source)
        if taxonomy_source == "ncbi":
            raise ValueError("NCBI unavailable")
        return root_compare_mod.outgroup_rooting(tree, "A")

    monkeypatch.setattr(
        root_compare_mod,
        "taxonomy_rooting",
        fake_taxonomy_rooting,
    )

    evaluations = evaluate_root_compare_methods(
        tree,
        ["taxonomy:ncbi", "taxonomy:opentree", "taxonomy:timetree"],
        _compare_args(),
        automatic=True,
    )

    assert calls == ["ncbi", "opentree", "timetree"]
    assert [evaluation.status for evaluation in evaluations] == [
        "failed",
        "ok",
        "ok",
    ]
    assert [evaluation.source for evaluation in evaluations] == [
        "ncbi",
        "opentree",
        "timetree",
    ]


def test_tsv_has_one_row_for_every_tied_best_position():
    split_a = (frozenset({"A"}), frozenset({"B", "C"}))
    split_b = (frozenset({"B"}), frozenset({"A", "C"}))
    evaluation = RootingEvaluation(
        method="mv",
        selection_basis="optimized",
        score_name="root_to_tip_variance",
        score_unit="branch_length_squared",
        candidates=(
            RootingCandidate(split_a, "exact", 0.25, 2.0, score=0.1),
            RootingCandidate(split_b, "exact", 0.75, 4.0, score=0.1),
        ),
        evaluated_edges=3,
        tie_rule="exact",
    )

    rows = root_compare_rows([evaluation])

    assert [row["tie_index"] for row in rows] == ["1", "2"]
    assert {row["num_best"] for row in rows} == {"2"}
    assert [row["is_canonical"] for row in rows] == ["yes", "no"]
    assert len({row["root_split_id"] for row in rows}) == 2


def test_zero_length_components_are_one_physical_root_position():
    tree = Tree("((A:1,B:1):0,C:1,D:1);", parser=1)
    all_taxa = frozenset(tree.leaf_names())
    splits = (
        root_compare_mod.canonical_split(frozenset({"A"}), all_taxa - frozenset({"A"})),
        root_compare_mod.canonical_split(frozenset({"B"}), all_taxa - frozenset({"B"})),
        root_compare_mod.canonical_split(
            frozenset({"A", "B"}), all_taxa - frozenset({"A", "B"})
        ),
    )
    evaluation = RootingEvaluation(
        method="mv",
        selection_basis="optimized",
        score_name="root_to_tip_variance",
        score_unit="branch_length_squared",
        candidates=(
            RootingCandidate(splits[0], "exact", 1.0, 1.0, score=0.0),
            RootingCandidate(splits[1], "exact", 1.0, 1.0, score=0.0),
            RootingCandidate(splits[2], "exact", 0.5, 0.0, score=0.0),
        ),
    )

    normalized = root_compare_mod._normalize_physical_root_positions(
        tree,
        evaluation,
    )

    assert len(normalized.candidates) == 1
    assert normalized.candidates[0].position_kind == "node"
    assert set(normalized.candidates[0].equivalent_splits) == set(splits)


def test_root_compare_main_writes_tsv_and_pdf(tmp_path):
    infile = tmp_path / "input.nwk"
    outfile = tmp_path / "summary.tsv"
    figure_out = tmp_path / "roots.pdf"
    layout_report = tmp_path / "layout.json"
    infile.write_text("(A:1,B:2,(C:3,D:4):1);\n")
    args = _compare_args(
        infile=str(infile),
        outfile=str(outfile),
        figure_out=str(figure_out),
        layout_report=str(layout_report),
    )

    evaluations = root_compare_main(args)

    with outfile.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert {evaluation.method for evaluation in evaluations} == {
        "midpoint",
        "mad",
        "mv",
    }
    assert {row["method"] for row in rows} == {"midpoint", "mad", "mv"}
    assert {row["status"] for row in rows} == {"ok"}
    assert figure_out.read_bytes().startswith(b"%PDF-")
    assert json.loads(layout_report.read_text())["final_collision_count"] == 0


def test_all_runs_every_supplied_context_method(tmp_path):
    data_dir = Path(__file__).with_name("data")
    outfile = tmp_path / "summary.tsv"
    figure_out = tmp_path / "roots.pdf"

    evaluations = root_compare_main(
        _compare_args(
            infile=str(data_dir / "rootcompare_gene.nwk"),
            infile2=str(data_dir / "rootcompare_reference.nwk"),
            species_tree=str(data_dir / "rootcompare_species.nwk"),
            outgroup="Gallus_gallus_g1",
            exclude_methods="taxonomy",
            outfile=str(outfile),
            figure_out=str(figure_out),
        )
    )

    assert [evaluation.method for evaluation in evaluations] == [
        "midpoint",
        "mad",
        "mv",
        "outgroup",
        "transfer",
        "reconciliation",
    ]
    assert {evaluation.status for evaluation in evaluations} == {"ok"}
    sources = {evaluation.method: evaluation.source for evaluation in evaluations}
    assert sources["transfer"] == str(data_dir / "rootcompare_reference.nwk")
    assert sources["reconciliation"] == str(data_dir / "rootcompare_species.nwk")


def test_transfer_preserves_exact_source_ratio_and_accepts_singleton_root(tmp_path):
    source = tmp_path / "source.nwk"
    source.write_text("(((A:1,B:1):0.4,(C:1,D:1):0.3):0.2);\n")
    tree = Tree("((A:2,B:2):2,(C:2,D:2):3);", parser=1)

    evaluation = evaluate_root_compare_methods(
        tree,
        ["transfer"],
        _compare_args(infile2=str(source)),
        automatic=False,
    )[0]

    assert evaluation.source == str(source)
    assert evaluation.candidates[0].position_kind == "exact"
    assert evaluation.candidates[0].position_fraction_from_side_a == pytest.approx(
        4.0 / 7.0
    )


@pytest.mark.parametrize(
    "source_newick",
    [
        "((A:1,B:1):0,(C:1,D:1):0);",
        "((A:1,B:1),(C:1,D:1):0.3);",
    ],
)
def test_transfer_reports_edge_only_when_source_ratio_is_unusable(
    tmp_path,
    source_newick,
):
    source = tmp_path / "source.nwk"
    source.write_text(source_newick + "\n")
    tree = Tree("((A:2,B:2):2,(C:2,D:2):3);", parser=1)

    evaluation = evaluate_root_compare_methods(
        tree,
        ["transfer"],
        _compare_args(infile2=str(source)),
        automatic=False,
    )[0]

    assert evaluation.candidates[0].position_kind == "edge_unspecified"
    assert evaluation.candidates[0].position_fraction_from_side_a is None


def test_automatic_context_failures_keep_their_source_paths(tmp_path):
    missing_transfer = tmp_path / "missing-reference.nwk"
    missing_species = tmp_path / "missing-species.nwk"
    tree = Tree("(A:1,(B:1,C:1):1);", parser=1)

    evaluations = evaluate_root_compare_methods(
        tree,
        ["transfer", "reconciliation"],
        _compare_args(
            infile2=str(missing_transfer),
            species_tree=str(missing_species),
        ),
        automatic=True,
    )

    assert [evaluation.source for evaluation in evaluations] == [
        str(missing_transfer),
        str(missing_species),
    ]


def test_root_compare_main_merges_endpoint_ties_at_one_physical_node(
    monkeypatch,
    tmp_path,
):
    infile = tmp_path / "input.nwk"
    outfile = tmp_path / "summary.tsv"
    figure_out = tmp_path / "roots.pdf"
    infile.write_text("(A:1,B:1,C:1,D:1);\n")
    captured = {}

    def fake_draw(tree, evaluations, args):
        captured["markers"] = root_compare_mod._rooting_markers(evaluations)
        figure_out.write_bytes(b"%PDF-1.4\n")

    monkeypatch.setattr(root_compare_mod, "_draw_root_comparison", fake_draw)

    root_compare_main(
        _compare_args(
            infile=str(infile),
            outfile=str(outfile),
            figure_out=str(figure_out),
            methods="mad",
        )
    )

    with outfile.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(captured["markers"]) == len(rows)
    assert len(rows) == 1
    assert rows[0]["position_kind"] == "node"
    assert rows[0]["num_equivalent_splits"] == "4"
    assert len(rows[0]["equivalent_root_split_ids"].split(",")) == 4


def test_auto_layout_switches_to_equal_angle_above_daylight_limit(monkeypatch):
    labels = ["T{}".format(index) for index in range(1002)]
    clades = ["{}:1".format(label) for label in labels]
    while len(clades) > 1:
        clades = [
            "({},{})".format(clades[index], clades[index + 1]) + ":1"
            if index + 1 < len(clades)
            else clades[index]
            for index in range(0, len(clades), 2)
        ]
    tree = Tree(clades[0].removesuffix(":1") + ";", parser=1)
    split = (frozenset({"T0"}), frozenset(labels[1:]))
    evaluation = RootingEvaluation(
        method="midpoint",
        selection_basis="diameter_midpoint",
        score_name="",
        score_unit="",
        candidates=(RootingCandidate(split, "exact", 0.5, 1.0),),
    )
    captured = {}

    def fake_draw_tree(**kwargs):
        captured.update(kwargs)

    monkeypatch.setattr(root_compare_mod, "_draw_tree", fake_draw_tree)

    root_compare_mod._draw_root_comparison(
        tree,
        [evaluation],
        _compare_args(unrooted_method="auto"),
    )

    assert captured["unrooted_method"] == "equal-angle"
    assert captured["tip_labels"] is False
    assert captured["figure_width"] == 7.2


def test_auto_figure_width_scales_with_tip_label_density():
    small_tree = Tree("(A:1,B:1,C:1);", parser=1)
    star_tree = Tree(
        "(" + ",".join("T{}:1".format(index) for index in range(100)) + ");",
        parser=1,
    )

    assert root_compare_mod._resolved_tip_labels(small_tree, "auto") is True
    assert root_compare_mod._resolved_tip_labels(star_tree, "auto") is True
    assert root_compare_mod._resolved_figure_width(small_tree, None, 8.0, True) == 7.2
    assert root_compare_mod._resolved_figure_width(star_tree, None, 8.0, True) == 12.0
    assert root_compare_mod._resolved_figure_width(star_tree, 9.0, 8.0, True) == 9.0


def test_all_failed_writes_failure_tsv_but_preserves_pdf(monkeypatch, tmp_path):
    infile = tmp_path / "input.nwk"
    outfile = tmp_path / "summary.tsv"
    figure_out = tmp_path / "roots.pdf"
    infile.write_text("(A:1,B:1,C:1);\n")
    figure_out.write_bytes(b"existing-pdf")

    monkeypatch.setattr(
        root_compare_mod,
        "evaluate_root_compare_methods",
        lambda *args, **kwargs: [
            RootingEvaluation.failed("midpoint", "synthetic failure")
        ],
    )

    with pytest.raises(ValueError, match="Every selected rooting method failed"):
        root_compare_main(
            _compare_args(
                infile=str(infile),
                outfile=str(outfile),
                figure_out=str(figure_out),
            )
        )

    with outfile.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert rows[0]["status"] == "failed"
    assert rows[0]["message"] == "synthetic failure"
    assert figure_out.read_bytes() == b"existing-pdf"


def test_draw_failure_preserves_existing_output_pair(monkeypatch, tmp_path):
    infile = tmp_path / "input.nwk"
    outfile = tmp_path / "summary.tsv"
    figure_out = tmp_path / "roots.pdf"
    infile.write_text("(A:1,B:2,C:3);\n")
    outfile.write_bytes(b"existing-tsv")
    figure_out.write_bytes(b"existing-pdf")

    def fail_draw(*args, **kwargs):
        raise RuntimeError("synthetic render failure")

    monkeypatch.setattr(root_compare_mod, "_draw_root_comparison", fail_draw)

    with pytest.raises(RuntimeError, match="synthetic render failure"):
        root_compare_main(
            _compare_args(
                infile=str(infile),
                outfile=str(outfile),
                figure_out=str(figure_out),
                methods="midpoint",
            )
        )

    assert outfile.read_bytes() == b"existing-tsv"
    assert figure_out.read_bytes() == b"existing-pdf"


def test_tsv_write_failure_preserves_existing_output_pair(monkeypatch, tmp_path):
    infile = tmp_path / "input.nwk"
    outfile = tmp_path / "summary.tsv"
    figure_out = tmp_path / "roots.pdf"
    infile.write_text("(A:1,B:2,C:3);\n")
    outfile.write_bytes(b"existing-tsv")
    figure_out.write_bytes(b"existing-pdf")

    def fail_write(path, rows):
        Path(path).write_text("partial")
        raise OSError("synthetic TSV failure")

    monkeypatch.setattr(root_compare_mod, "_write_root_compare_rows", fail_write)

    with pytest.raises(OSError, match="synthetic TSV failure"):
        root_compare_main(
            _compare_args(
                infile=str(infile),
                outfile=str(outfile),
                figure_out=str(figure_out),
                methods="midpoint",
            )
        )

    assert outfile.read_bytes() == b"existing-tsv"
    assert figure_out.read_bytes() == b"existing-pdf"


def test_commit_failure_rolls_back_both_existing_outputs(monkeypatch, tmp_path):
    infile = tmp_path / "input.nwk"
    outfile = tmp_path / "summary.tsv"
    figure_out = tmp_path / "roots.pdf"
    infile.write_text("(A:1,B:2,C:3);\n")
    outfile.write_bytes(b"existing-tsv")
    figure_out.write_bytes(b"existing-pdf")
    original_replace = root_compare_mod.os.replace
    failed = False

    def fail_figure_commit(source, destination):
        nonlocal failed
        source_name = Path(source).name
        if (
            not failed
            and Path(destination) == figure_out
            and source_name.startswith(".roots.pdf.")
        ):
            failed = True
            raise OSError("synthetic commit failure")
        return original_replace(source, destination)

    monkeypatch.setattr(root_compare_mod.os, "replace", fail_figure_commit)

    with pytest.raises(OSError, match="synthetic commit failure"):
        root_compare_main(
            _compare_args(
                infile=str(infile),
                outfile=str(outfile),
                figure_out=str(figure_out),
                methods="midpoint",
            )
        )

    assert failed is True
    assert outfile.read_bytes() == b"existing-tsv"
    assert figure_out.read_bytes() == b"existing-pdf"
    assert list(tmp_path.glob(".*.tmp")) == []


def test_invalid_tsv_parent_is_detected_before_drawing(monkeypatch, tmp_path):
    infile = tmp_path / "input.nwk"
    figure_out = tmp_path / "roots.pdf"
    infile.write_text("(A:1,B:2,C:3);\n")
    figure_out.write_bytes(b"existing-pdf")
    called = False

    def record_draw(*args, **kwargs):
        nonlocal called
        called = True

    monkeypatch.setattr(root_compare_mod, "_draw_root_comparison", record_draw)

    with pytest.raises(FileNotFoundError):
        root_compare_main(
            _compare_args(
                infile=str(infile),
                outfile=str(tmp_path / "missing" / "summary.tsv"),
                figure_out=str(figure_out),
                methods="midpoint",
            )
        )

    assert called is False
    assert figure_out.read_bytes() == b"existing-pdf"

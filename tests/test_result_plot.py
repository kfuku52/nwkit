"""Check figures against saved scientific results and exercise publication failures."""

import hashlib
import json
import shutil
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import pytest

from nwkit.cli import main
from nwkit.radte import radte_paths
from nwkit.result_plot import build_result_figure
from nwkit.result_plot_data import (
    RADTE_PLOT_SUFFIXES,
    dating_plot_data,
    load_dating_plot_data,
    read_result_table,
    reconciliation_plot_data,
)
from nwkit.util import read_tree


@pytest.fixture(scope="module")
def reference(tmp_path_factory):
    root = tmp_path_factory.mktemp("result-plot-reference")
    (root / "gene.nwk").write_text(
        "((A_1:0.1,B_1:0.1)S1:0.1,(A_2:0.2,B_2:0.2)S2:0.2)D;\n"
    )
    (root / "species.nwk").write_text("(A:10,B:10)AB;\n")
    (root / "map.tsv").write_text(
        "leaf_name\tspecies_label\nA_1\tA\nB_1\tB\nA_2\tA\nB_2\tB\n"
    )
    main(reconcile_args(root))
    main(radte_args(root))
    return root


@pytest.fixture
def result_dir(reference, tmp_path):
    shutil.copytree(reference, tmp_path, dirs_exist_ok=True)
    return tmp_path


def reconcile_args(root):
    return [
        "reconcile",
        "--infile",
        str(root / "gene.nwk"),
        "--species-tree",
        str(root / "species.nwk"),
        "--species-map-tsv",
        str(root / "map.tsv"),
        "--outfile",
        str(root / "events.tsv"),
    ]


def radte_args(root):
    return [
        "radte",
        "--gene-tree",
        str(root / "gene.nwk"),
        "--species-tree",
        str(root / "species.nwk"),
        "--species-map-tsv",
        str(root / "map.tsv"),
        "--reconcile",
        "lca",
        "--max-age",
        "30",
        "--rate-sd",
        "0.3",
        "--uncertainty",
        "profile",
        "--out-prefix",
        str(root / "family"),
    ]


def saved_args(root, mode, output):
    args = [
        "draw",
        "--species-tree",
        str(root / "species.nwk"),
        "--outfile",
        str(output),
    ]
    if mode == "radte":
        return args + ["--radte-prefix", str(root / "family")]
    return args + [
        "--infile",
        str(root / "gene.nwk"),
        "--reconciliation",
        str(root / "events.tsv"),
    ]


def load_data(root):
    return load_dating_plot_data(root / "family", root / "species.nwk")


@pytest.mark.parametrize("mode", ["reconcile", "radte"])
@pytest.mark.parametrize("extension", ["pdf", "svg"])
def test_direct_commands_publish_vector_report_with_results(
    result_dir, mode, extension
):
    output = result_dir / f"direct.{extension}"
    main(
        (reconcile_args(result_dir) if mode == "reconcile" else radte_args(result_dir))
        + ["--figure-out", str(output)]
    )
    content = output.read_bytes()
    assert content.startswith(b"%PDF") if extension == "pdf" else b"<svg" in content
    assert len(content) > 3000
    if mode == "radte":
        manifest = json.loads((result_dir / "family.manifest.json").read_text())
        assert (
            manifest["output_sha256"]["figure"] == hashlib.sha256(content).hexdigest()
        )


@pytest.mark.parametrize("mode", ["reconcile", "radte"])
def test_saved_results_render_without_running_reconciliation_or_dating(
    result_dir, mode, monkeypatch
):
    import nwkit.radte
    import nwkit.reconcile

    def forbidden(*args, **kwargs):
        raise AssertionError("Plotting must not perform inference")

    monkeypatch.setattr(nwkit.radte, "run_dating", forbidden)
    monkeypatch.setattr(nwkit.reconcile, "build_reconciliation_table", forbidden)
    output = result_dir / "saved.svg"
    before = {p: p.read_bytes() for p in result_dir.iterdir() if p.is_file()}
    main(saved_args(result_dir, mode, output))
    assert "A_1" in output.read_text()
    assert {p: p.read_bytes() for p in before} == before


def test_dated_coordinates_intervals_and_shared_age_guides_match_results(reference):
    data = load_data(reference)
    fig = build_result_figure(data)
    try:
        gene_ax, species_ax = fig.axes[:2]
        assert gene_ax.get_xlim() == species_ax.get_xlim()
        assert gene_ax.get_position().x0 == species_ax.get_position().x0
        assert gene_ax.get_position().width == species_ax.get_position().width
        guide_positions = [
            [
                line.get_xdata()[0]
                for line in ax.lines
                if line.get_gid() == "shared-speciation-age"
            ]
            for ax in (gene_ax, species_ax)
        ]
        assert guide_positions == [[10.0], [10.0]]
        nodes = {n.name: n for n in data.gene.traverse()}
        for name, age in [("S1", 10), ("S2", 10), ("D", 20)]:
            clade = data.gene_index.clade_id_for_node(nodes[name])
            marker = next(
                c
                for c in gene_ax.collections
                if (c.get_gid() or "").startswith("event:")
                and c.get_gid().endswith(clade)
            )
            assert float(marker.get_offsets()[0, 0]) == pytest.approx(age, abs=1e-4)
        root_id = data.gene_index.clade_id_for_node(data.gene)
        interval = next(
            c for c in gene_ax.collections if c.get_gid() == "age-interval:" + root_id
        )
        row = data.events[root_id]
        assert interval.get_segments()[0][:, 0].tolist() == pytest.approx(
            [row["interval_lower"], row["interval_upper"]]
        )
        assert gene_ax.get_xlim()[0] > row["interval_upper"]
        assert "Fixed calibration" in [t.get_text() for t in fig.legends[0].get_texts()]
    finally:
        plt.close(fig)


def test_unrequested_intervals_are_neither_invented_nor_labeled_posterior(reference):
    data = load_data(reference)
    data.manifest["uncertainty"] = "none"
    for row in data.events.values():
        row["interval_lower"] = row["interval_upper"] = None
    fig = build_result_figure(data)
    try:
        assert not any(
            (c.get_gid() or "").startswith("age-interval:")
            for c in fig.axes[0].collections
        )
        assert "Age interval" not in [t.get_text() for t in fig.legends[0].get_texts()]
        assert "Not requested; points only" in "\n".join(
            t.get_text() for ax in fig.axes for t in ax.texts
        )
    finally:
        plt.close(fig)


def test_soft_priors_are_not_labeled_fixed_calibrations(reference):
    data = load_data(reference)
    data.manifest["calibration_policy"] = "PAML-soft-root-and-speciation-priors"
    fig = build_result_figure(data)
    try:
        assert "Fixed calibration" not in [
            t.get_text() for t in fig.legends[0].get_texts()
        ]
        assert "original ranges" in fig.axes[2].get_title(loc="left")
    finally:
        plt.close(fig)


def test_reordered_species_and_gene_children_keep_clade_mapping(result_dir):
    (result_dir / "gene.nwk").write_text(
        "((B_2:0.2,A_2:0.2)second:0.2,(B_1:0.1,A_1:0.1)first:0.1)root;"
    )
    (result_dir / "species.nwk").write_text("(B:10,A:10)AB;")
    output = result_dir / "reordered.svg"
    main(saved_args(result_dir, "reconcile", output))
    assert "first / AB" in output.read_text()
    main(saved_args(result_dir, "radte", result_dir / "dated-reordered.svg"))


@pytest.mark.parametrize("table", ["nodes", "species", "tree"])
def test_tampered_bundle_fails_before_replacing_a_previous_plot(result_dir, table):
    source = Path(radte_paths(str(result_dir / "family"))[table])
    source.write_bytes(source.read_bytes() + b"\n")
    output = result_dir / "prior.pdf"
    output.write_bytes(b"previous figure")
    with pytest.raises(ValueError, match="does not match its manifest"):
        main(saved_args(result_dir, "radte", output))
    assert output.read_bytes() == b"previous figure"


@pytest.mark.parametrize(
    "failure",
    ["topology", "shared_age", "interval_endpoint", "interval_order", "missing_column"],
)
def test_inconsistent_scientific_tables_are_rejected(reference, failure):
    tree = read_tree(str(reference / "family.dated.nwk"), "auto", True)
    species = read_tree(str(reference / "species.nwk"), "auto", True)
    nodes = read_result_table(reference / "family.nodes.tsv")
    sp = read_result_table(reference / "family.species.tsv")
    if failure == "topology":
        nodes.loc[1, "parent_gene_clade_id"] = "wrong"
    elif failure == "shared_age":
        sp.loc[sp["node"] == "AB", "estimated_age"] = "11"
    elif failure == "interval_endpoint":
        nodes.loc[0, "interval_upper"] = "NA"
    elif failure == "interval_order":
        nodes.loc[0, ["interval_lower", "interval_upper"]] = ["40", "30"]
    else:
        nodes = nodes.drop(columns=["interval_upper"])
    with pytest.raises(ValueError):
        dating_plot_data(tree, species, nodes, sp, {})


@pytest.mark.parametrize("mode", ["reconcile", "radte"])
@pytest.mark.parametrize("link_type", ["hard", "symlink"])
def test_figure_cannot_replace_an_input_through_an_alias(result_dir, mode, link_type):
    output = result_dir / "alias.pdf"
    source = result_dir / "species.nwk"
    before = source.read_bytes()
    if link_type == "hard":
        output.hardlink_to(source)
    else:
        output.symlink_to(source)
    with pytest.raises(ValueError):
        main(saved_args(result_dir, mode, output))
    assert source.read_bytes() == before


def test_replot_cannot_overwrite_other_members_of_the_radte_bundle(result_dir):
    target = result_dir / "family.events.tsv"
    before = target.read_bytes()
    with pytest.raises(ValueError):
        main(saved_args(result_dir, "radte", target) + ["--image-format", "pdf"])
    assert target.read_bytes() == before


@pytest.mark.parametrize("mode", ["reconcile", "radte"])
def test_late_render_error_preserves_previous_numerical_results_and_figure(
    result_dir, mode, monkeypatch
):
    import nwkit.result_plot

    output = result_dir / "previous.pdf"
    output.write_bytes(b"previous figure")
    before = {p: p.read_bytes() for p in result_dir.iterdir() if p.is_file()}

    def fail_after_write(data, path, **kwargs):
        Path(path).write_bytes(b"partial figure")
        raise RuntimeError("injected render error")

    monkeypatch.setattr(nwkit.result_plot, "save_result_figure", fail_after_write)
    with pytest.raises(RuntimeError, match="injected"):
        main(
            (
                reconcile_args(result_dir)
                if mode == "reconcile"
                else radte_args(result_dir)
            )
            + ["--figure-out", str(output)]
        )
    assert {p: p.read_bytes() for p in before} == before


def test_reconcile_publication_error_rolls_back_the_table_already_installed(
    result_dir, monkeypatch
):
    import nwkit.output_transaction as transaction

    output = result_dir / "publication.pdf"
    output.write_bytes(b"previous figure")
    before = {p: p.read_bytes() for p in (result_dir / "events.tsv", output)}
    replace = transaction.os.replace
    triggered = []

    def fail_final_figure(source, destination):
        if (
            str(destination) == str(output)
            and ".stage." in str(source)
            and not triggered
        ):
            triggered.append(True)
            raise OSError("injected publication error")
        return replace(source, destination)

    monkeypatch.setattr(transaction.os, "replace", fail_final_figure)
    with pytest.raises(OSError, match="publication"):
        main(reconcile_args(result_dir) + ["--figure-out", str(output)])
    assert triggered
    assert {p: p.read_bytes() for p in before} == before


def test_saved_radte_audit_hashes_bundle_without_reading_stdin(result_dir, monkeypatch):
    class NoStdin:
        def read(self, *args):
            raise AssertionError("Saved result plotting must not consume STDIN")

    monkeypatch.setattr(sys, "stdin", NoStdin())
    audit = result_dir / "audit.jsonl"
    main(
        saved_args(result_dir, "radte", result_dir / "audit.svg")
        + ["--audit", str(audit)]
    )
    record = json.loads(audit.read_text().splitlines()[-1])
    assert "family.nodes.tsv" in json.dumps(record)
    assert "family.species.tsv" in json.dumps(record)


def test_report_rejects_single_tree_layout_options(result_dir):
    with pytest.raises(ValueError, match="single-tree option"):
        main(
            saved_args(result_dir, "radte", result_dir / "invalid.pdf")
            + ["--layout", "radial"]
        )


def test_reconciliation_requires_matching_species_clades(reference):
    gene = read_tree(str(reference / "gene.nwk"), "auto", True)
    species = read_tree("(A:1,C:1);", "auto", True)
    with pytest.raises(ValueError, match="species tree"):
        reconciliation_plot_data(
            gene, species, read_result_table(reference / "events.tsv")
        )


def test_plot_bundle_paths_follow_the_radte_output_contract():
    paths = radte_paths("example")
    assert all(
        paths[key] == "example" + suffix for key, suffix in RADTE_PLOT_SUFFIXES.items()
    )


def test_not_requested_intervals_are_points_only(result_dir):
    from nwkit.result_plot import _diagnostic_lines

    data = load_data(result_dir)
    data.manifest["uncertainty"] = "not-requested"
    lines = _diagnostic_lines(data)
    assert "Intervals: Not requested; points only" in lines
    assert not any("95%" in line for line in lines)


def test_paml_unestimated_species_nodes_use_explicit_input_positions():
    from ete4 import Tree

    from nwkit.radte import dated_newick, result_tables
    from nwkit.radte_inputs import build_chronology
    from nwkit.radte_model import fit_dates
    from nwkit.reconcile import build_reconciliation_table

    gene = Tree("(A_1:0.2,B_1:0.2)R;", parser=1)
    species = Tree("((A:10,C:10)AC:10,B:20)R;", parser=1)
    events = build_reconciliation_table(gene, species, {"A_1": "A", "B_1": "B"})
    chronology = build_chronology(gene, species, events, None, 30)
    fit, _ = fit_dates(chronology, rate_sd=0.3)
    tables = result_tables(chronology, fit, bound_policy="PAML-soft-prior")
    assert (
        tables["species"]
        .loc[tables["species"].node == "AC", "estimated_age"]
        .isna()
        .all()
    )
    data = dating_plot_data(
        read_tree(dated_newick(chronology, fit), "auto", True),
        species,
        tables["nodes"],
        tables["species"],
        {
            "calibration_policy": "PAML-soft-root-and-speciation-priors",
            "method": "mcmctree-posterior-mean",
            "interval_level": 0.95,
        },
    )
    fig = build_result_figure(data)
    try:
        ax = fig.axes[1]
        assert "Input species tree" in ax.get_title(loc="left")
        labels = [t.get_text() for t in ax.texts]
        assert any("AC" in text and "input only" in text for text in labels)
        ac = next(n for n in species.traverse() if n.name == "AC")
        sid = data.species_index.clade_id_for_node(ac)
        point = next(
            c for c in ax.collections if c.get_gid() == "event:speciation:" + sid
        )
        assert point.get_offsets()[0, 0] == 10
        assert data.species_rows[sid]["estimated_age"] is None
    finally:
        plt.close(fig)

"""External species evidence, conditional components, and three-mode reports."""

import hashlib
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest

from nwkit.cli import main, parser
from nwkit.radte_inputs import read_inputs
from nwkit.radte_species import attach_species_intervals
from nwkit.result_plot import build_result_figure
from nwkit.result_plot_data import load_dating_plot_data
from tests.test_radte import small_chronology


def setup_inputs(root):
    root.mkdir(exist_ok=True)
    (root / "gene.nwk").write_text(
        "(((A_1:0.1,B_1:0.1)S1:0.05,(A_2:0.1,B_2:0.1)S2:0.05)D:0.05,C_1:0.2)R;\n"
    )
    (root / "species.nwk").write_text("((A:10,B:10)AB:10,C:20)R;\n")
    (root / "map.tsv").write_text(
        "leaf_name\tspecies_label\nA_1\tA\nB_1\tB\nA_2\tA\nB_2\tB\nC_1\tC\n"
    )
    (root / "bounds.tsv").write_text("node\tage_min\tage_max\nAB\t8\t12\n")
    (root / "intervals.tsv").write_text(
        "node\tlower\tupper\tlevel\tkind\tsource\nAB\t8\t12\t0.95\tcredible\tSynthetic example\nR\t16\t24\t0.95\tcredible\tSynthetic example\n"
    )
    (root / "species-samples.nwk").write_text(
        "\n".join(
            f"((A:{10 * f},B:{10 * f})AB:{10 * f},C:{20 * f})R;"
            for f in np.linspace(0.9, 1.1, 20)
        )
    )
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
        "--species-node-intervals-tsv",
        str(root / "intervals.tsv"),
    ]


@pytest.fixture(scope="module")
def trio(tmp_path_factory):
    root = tmp_path_factory.mktemp("species-trio")
    base = setup_inputs(root)
    for mode in ("fixed", "bounded", "ensemble"):
        extra = ["--uncertainty", "profile"]
        if mode == "bounded":
            extra += ["--species-node-bounds-tsv", str(root / "bounds.tsv")]
        if mode == "ensemble":
            extra = [
                "--uncertainty",
                "input-ensemble",
                "--species-tree-ensemble",
                str(root / "species-samples.nwk"),
                "--ensemble-within-uncertainty",
                "profile",
            ]
        main(
            base
            + extra
            + [
                "--out-prefix",
                str(root / mode),
                "--figure-out",
                str(root / (mode + ".pdf")),
            ]
        )
    return root


def test_external_intervals_never_change_chronology_bounds(tmp_path):
    base = setup_inputs(tmp_path)
    args = parser.parse_args(base + ["--out-prefix", str(tmp_path / "x")])
    annotated = read_inputs(args)
    args.species_node_intervals_tsv = None
    plain = read_inputs(args)
    np.testing.assert_array_equal(annotated.lower, plain.lower)
    np.testing.assert_array_equal(annotated.upper, plain.upper)
    np.testing.assert_array_equal(annotated.initial, plain.initial)


@pytest.mark.parametrize(
    "row",
    [
        "AB\t12\t8\t0.95\tcredible\tx",
        "AB\t8\t12\t95\tcredible\tx",
        "AB\t8\t12\t0.95\thard\tx",
        "missing\t8\t12\t0.95\tcredible\tx",
        "AB\t8\t12\t0.95\tcredible\t",
    ],
)
def test_invalid_external_intervals_are_rejected(tmp_path, row):
    path = tmp_path / "interval.tsv"
    path.write_text("node\tlower\tupper\tlevel\tkind\tsource\n" + row + "\n")
    with pytest.raises(ValueError):
        attach_species_intervals(small_chronology(), path)


def test_species_plot_has_input_bounds_and_estimated_intervals(trio):
    data = load_dating_plot_data(str(trio / "bounded"), str(trio / "species.nwk"))
    fig = build_result_figure(data)
    try:
        gids = [a.get_gid() for a in fig.findobj() if a.get_gid()]
        assert any(g.startswith("species-input-interval:") for g in gids)
        assert any(g.startswith("species-hard-range:") for g in gids)
        assert any(g.startswith("species-estimated-interval:") for g in gids)
        row = next(r for r in data.species_rows.values() if r["node"] == "AB")
        external = next(
            a
            for a in fig.findobj()
            if a.get_gid() == "species-input-interval:" + row["species_event_id"]
        )
        np.testing.assert_allclose(
            external.get_segments()[0][:, 0],
            [row["input_interval_lower"], row["input_interval_upper"]],
        )

        for event in data.events.values():
            if (
                event.get("species_name") == "AB"
                and event["event_type"] == "speciation"
            ):
                assert event["interval_lower"] == row["interval_lower"]
                assert event["interval_upper"] == row["interval_upper"]
    finally:
        plt.close(fig)


def test_components_match_saved_conditional_fits(trio):
    components = pd.read_csv(trio / "ensemble.uncertainty-components.tsv", sep="\t")
    conditional = pd.read_csv(trio / "ensemble.conditional-intervals.tsv", sep="\t")
    samples = pd.read_csv(trio / "ensemble.age-samples.tsv", sep="\t")
    assert components.status.eq("separate-components").all()
    for row in components.itertuples():
        points = samples.loc[
            samples.shared_age_id == row.shared_age_id, "estimated_age"
        ]
        within = conditional.loc[conditional.shared_age_id == row.shared_age_id]
        assert row.input_refit_variance == pytest.approx(points.var(ddof=1))
        assert row.mean_conditional_interval_width == pytest.approx(
            (within.interval_upper - within.interval_lower).mean()
        )
        assert row.conditional_samples == 20
        if row.shared_age_id.startswith("S:"):
            assert row.mean_conditional_interval_width == 0
    species = pd.read_csv(trio / "ensemble.species.tsv", sep="\t").set_index("node")
    assert species.loc["AB", "interval_upper"] > species.loc["AB", "interval_lower"]
    assert species.loc["R", "interval_upper"] > species.loc["R", "interval_lower"]


def compare_args(root, output):
    return [
        "radte-compare",
        "--fixed-prefix",
        str(root / "fixed"),
        "--bounded-prefix",
        str(root / "bounded"),
        "--ensemble-prefix",
        str(root / "ensemble"),
        "--species-tree",
        str(root / "species.nwk"),
        "--out-prefix",
        str(output),
    ]


def test_saved_comparison_does_not_refit_and_has_verified_outputs(
    trio, tmp_path, monkeypatch
):
    import nwkit.radte

    monkeypatch.setattr(
        nwkit.radte, "run_dating", lambda *a, **kw: pytest.fail("must not refit")
    )
    prefix = tmp_path / "compare"
    main(compare_args(trio, prefix))
    manifest = json.loads(prefix.with_suffix(".manifest.json").read_text())
    assert (
        hashlib.sha256(prefix.with_suffix(".pdf").read_bytes()).hexdigest()
        == manifest["output_sha256"]["figure"]
    )
    table = pd.read_csv(prefix.with_suffix(".comparison.tsv"), sep="\t")
    assert set(table["mode"]) == {"fixed", "bounded", "ensemble"}
    assert set(table.kind) == {"species", "duplication"}


def test_conditional_failure_does_not_discard_input_point_fits(tmp_path, monkeypatch):
    from nwkit import radte_uncertainty

    def fail(*args, **kwargs):
        raise ValueError("injected conditional failure")

    base = setup_inputs(tmp_path)
    monkeypatch.setattr(radte_uncertainty, "profile_intervals", fail)
    main(
        base
        + [
            "--uncertainty",
            "input-ensemble",
            "--species-tree-ensemble",
            str(tmp_path / "species-samples.nwk"),
            "--ensemble-within-uncertainty",
            "profile",
            "--out-prefix",
            str(tmp_path / "failed-within"),
        ]
    )
    components = pd.read_csv(
        tmp_path / "failed-within.uncertainty-components.tsv", sep="\t"
    )
    conditional = pd.read_csv(
        tmp_path / "failed-within.conditional-intervals.tsv", sep="\t"
    )
    assert components.input_samples.eq(20).all()
    assert components.conditional_samples.eq(0).all()
    assert components.status.eq("unavailable-conditional-coverage").all()
    assert components.mean_conditional_interval_width.isna().all()
    assert conditional.error.eq("injected conditional failure").all()


def test_conditional_bootstrap_reports_variance_and_fixed_species_zeros(tmp_path):
    from nwkit.radte import run_dating
    from nwkit.radte_components import conditional_intervals

    base = setup_inputs(tmp_path)
    args = parser.parse_args(
        base
        + ["--out-prefix", str(tmp_path / "unused"), "--bootstrap-replicates", "20"]
    )
    chronology = read_inputs(args)
    fit, problem, _, _ = run_dating(chronology, args)
    args.ensemble_within_uncertainty = "bootstrap"
    rows = pd.DataFrame(
        conditional_intervals(chronology, fit, problem, args, 1, "branch-marginal")
    )
    assert rows.interval_lower.notna().all()
    assert (
        rows.loc[rows.shared_age_id.str.startswith("S:"), "conditional_variance"]
        .eq(0)
        .all()
    )
    assert (
        rows.loc[~rows.shared_age_id.str.startswith("S:"), "conditional_variance"] > 0
    ).all()


def test_comparison_refuses_overwriting_source_bundle(trio):
    path = trio / "ensemble.uncertainty-components.tsv"
    before = path.read_bytes()
    with pytest.raises(ValueError, match="replace|overwrite"):
        main(compare_args(trio, trio / "ensemble"))
    assert path.read_bytes() == before


def test_comparison_rejects_wrong_mode_without_replacing_output(trio, tmp_path):
    output = tmp_path / "compare"
    output.with_suffix(".pdf").write_bytes(b"keep previous")
    arguments = compare_args(trio, output)
    arguments[arguments.index("--fixed-prefix") + 1] = str(trio / "bounded")
    with pytest.raises(ValueError, match="fixed species ages"):
        main(arguments)
    assert output.with_suffix(".pdf").read_bytes() == b"keep previous"


def test_comparison_audit_cannot_overwrite_component_input(trio, tmp_path):
    path = trio / "ensemble.conditional-intervals.tsv"
    before = path.read_bytes()
    with pytest.raises(ValueError, match="replace|overwrite|collid|alias|distinct"):
        main(compare_args(trio, tmp_path / "compare") + ["--audit", str(path)])
    assert path.read_bytes() == before


def test_comparison_render_failure_preserves_all_outputs_and_closes_figures(
    trio, tmp_path, monkeypatch
):
    import nwkit.radte_compare as comparison

    prefix = tmp_path / "compare"
    paths = comparison.comparison_paths(str(prefix))
    for path in paths.values():
        Path(path).write_bytes(b"previous")
    existing = set(plt.get_fignums())

    def failed(*args):
        plt.figure()
        raise ValueError("injected renderer failure")

    monkeypatch.setattr(comparison, "_components_figure", failed)
    with pytest.raises(ValueError, match="injected renderer"):
        main(compare_args(trio, prefix))
    assert set(plt.get_fignums()) == existing
    assert all(Path(path).read_bytes() == b"previous" for path in paths.values())


def test_comparison_trees_share_age_scale(trio):
    from nwkit.radte_compare_plot import comparison_figure

    runs = {
        mode: load_dating_plot_data(str(trio / mode), str(trio / "species.nwk"))
        for mode in ("fixed", "bounded", "ensemble")
    }
    fig = comparison_figure(runs, 0.95)
    try:
        assert len(fig.axes) == 6
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()
        legend_box = fig.legends[0].get_window_extent(renderer)
        for ax in fig.axes[:3]:
            title = next(
                a
                for a in ax.get_children()
                if hasattr(a, "get_text") and a.get_text() == ax.get_title(loc="left")
            )
            assert not title.get_window_extent(renderer).overlaps(legend_box)
        extent = fig.axes[0].get_xlim()
        assert extent[0] > extent[1]  # Older ages appear on the left.
        for ax in fig.axes:
            np.testing.assert_allclose(ax.get_xlim(), extent)
            assert ax.get_position().width == pytest.approx(
                fig.axes[0].get_position().width
            )
    finally:
        plt.close(fig)


def test_component_points_align_with_tree_nodes(trio):
    from nwkit.radte_compare_plot import components_figure
    from nwkit.result_plot import _layout, _species_tip_labels

    reference = load_dating_plot_data(str(trio / "fixed"), str(trio / "species.nwk"))
    components = pd.read_csv(trio / "ensemble.uncertainty-components.tsv", sep="\t")
    fig = components_figure(components, reference)
    try:
        fig.canvas.draw()
        for row, gene in enumerate((True, False)):
            tree = reference.gene if gene else reference.species
            index = reference.gene_index if gene else reference.species_index
            y = _layout(tree, None if gene else _species_tip_labels(tree))[1]
            tree_ax = fig.axes[row * 3]
            for col, field in enumerate(
                ("input_refit_sd", "mean_conditional_interval_width"), 1
            ):
                ax = fig.axes[row * 3 + col]
                np.testing.assert_allclose(ax.get_xlim(), fig.axes[col].get_xlim())
                artists = {a.get_gid(): a for a in ax.collections}
                for node in tree.traverse():
                    if node.is_leaf:
                        continue
                    sid = index.clade_id_for_node(node)
                    point = artists[f"component:{field}:{sid}"].get_offsets()[0]
                    # Screen coordinates verify both data layout and panel placement.
                    assert ax.transData.transform(point)[1] == pytest.approx(
                        tree_ax.transData.transform((0, y[node]))[1]
                    )
    finally:
        plt.close(fig)


def test_comparison_caption_reports_missing_intervals(trio):
    from nwkit.radte_compare_plot import _interval_caption

    data = load_dating_plot_data(str(trio / "fixed"), str(trio / "species.nwk"))
    for row in data.events.values():
        row["interval_lower"] = row["interval_upper"] = None
    assert _interval_caption(data) == "profile; gene intervals available: 0/4"


def test_comparison_honors_species_rooting_override(trio, tmp_path):
    with pytest.raises(ValueError, match="must be rooted"):
        main(
            compare_args(trio, tmp_path / "unrooted") + ["--species-tree-rooted", "no"]
        )
    assert not (tmp_path / "unrooted.pdf").exists()

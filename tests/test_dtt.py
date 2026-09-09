"""DTT reference, independent distance/GLS oracles and CLI regression checks."""

import json
from io import StringIO
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from scipy.spatial.distance import pdist

from nwkit.cli import main
from nwkit.disparity import (
    clade_disparities,
    curve_area,
    dtt_curve,
    fit_brownian,
    make_design,
    simulate_curves,
    transform_traits,
    validate_work,
)
from nwkit.dtt import _ultrametric_depths
from nwkit.evolution import build_evolutionary_process
from nwkit.util import read_tree

EXAMPLE = Path(__file__).resolve().parents[1] / "examples" / "dtt"
TREE = (EXAMPLE / "tree.nwk").read_text().strip()
DATA = pd.read_csv(EXAMPLE / "traits.tsv", sep="\t").iloc[:, 1:3].to_numpy()


def design_for(text=TREE):
    tree = read_tree(text, 0, False)
    _, height = _ultrametric_depths(tree)
    for node in tree.traverse():
        if not node.is_root:
            node.dist /= height
    depths, _ = _ultrametric_depths(tree)
    names = sorted(tree.leaf_names())
    return tree, make_design(tree, names, depths)


def command(tmp_path, extra=(), data=DATA, tree=TREE):
    table = tmp_path / "traits.tsv"
    pd.DataFrame(
        {"leaf_name": list("ABCDEFGH")[: len(data)], "x": data[:, 0], "y": data[:, 1]}
    ).to_csv(table, sep="\t", index=False)
    return [
        "dtt",
        "-i",
        tree,
        "--trait",
        str(table),
        "--columns",
        "x,y",
        "--n-sim",
        "9",
        *extra,
    ]


def test_geiger_2011_reference_and_independent_pairwise_oracle():
    tree, design = design_for()
    observed = dtt_curve(design, DATA)
    assert design.times == pytest.approx([0, 0, 0.25, 0.45, 0.55, 0.6, 0.65, 1])
    assert observed == pytest.approx(
        [
            1,
            0.54238574826661268,
            0.78328862796700227,
            0.83229583143546393,
            1.17475580748013542,
            1.31079508072270823,
            0,
            0,
        ],
        rel=1e-13,
    )
    disparities, _ = clade_disparities(design, DATA)
    lookup = dict(zip("ABCDEFGH", DATA, strict=True))
    dense = []
    for node in tree.traverse("preorder"):
        tips = np.array([lookup[name] for name in node.leaf_names()])
        dense.append(pdist(tips, metric="sqeuclidean").mean() if len(tips) > 1 else 0)
    assert disparities == pytest.approx(dense, abs=1e-13)
    assert observed == pytest.approx(
        design.weights.toarray() @ (np.array(dense) / dense[0])
    )


def test_bm_rate_geiger_and_gls_oracle():
    tree, _ = design_for()
    covariance = build_evolutionary_process(tree).tip_covariance(list("ABCDEFGH"))
    center, rate, factor = fit_brownian(covariance, DATA)
    inverse = np.linalg.inv(covariance)
    expected_center = np.ones(8) @ inverse @ DATA / inverse.sum()
    residual = DATA - expected_center
    assert center == pytest.approx(expected_center)
    assert rate == pytest.approx(residual.T @ inverse @ residual / 7)
    assert rate == pytest.approx(
        2
        * np.array(
            [
                [0.53835301552053938, 0.33303418946191948],
                [0.33303418946191948, 1.91578728670767195],
            ]
        ),
        rel=1e-13,
    )
    assert factor @ factor.T == pytest.approx(rate)


def test_parallel_seed_and_empirical_curve_reference():
    tree, design = design_for()
    covariance = build_evolutionary_process(tree).tip_covariance(list("ABCDEFGH"))
    _, _, factor = fit_brownian(covariance, DATA)
    serial = simulate_curves(design, factor, 12, seed=42)
    parallel = simulate_curves(design, factor, 12, seed=42, threads=2)
    assert np.array_equal(serial, parallel)
    assert np.array_equal(serial[:4], simulate_curves(design, factor, 4, seed=42))
    assert (serial[:, 0] == 1).all()
    assert (serial[:, -2:] == 0).all()


@pytest.mark.parametrize("scale", ["raw", "standardize"])
def test_units_and_translation_invariance(scale):
    _, design = design_for()
    transformed, _, _ = transform_traits(DATA, scale)
    factors = np.array([1e90, 2e90]) if scale == "standardize" else 1e90
    other, _, _ = transform_traits((DATA + 100) * factors, scale)
    assert dtt_curve(design, other) == pytest.approx(
        dtt_curve(design, transformed), abs=1e-13
    )


def test_area_duplicate_root_and_interpolated_range():
    times = np.array([0, 0, 0.5, 1])
    curves = np.array([[1, 2, 0, 0], [1, 4, 0, 0]])
    assert curve_area(times, curves) == pytest.approx([0.5, 1])
    assert curve_area(times, curves, (0.1, 0.4)) == pytest.approx([0.3, 0.6])


@pytest.mark.parametrize(
    "format_name,magic", [("png", b"\x89PNG"), ("pdf", b"%PDF"), ("svg", b"<?xml")]
)
def test_cli_artifacts_and_mdi_reconstruction(tmp_path, capsys, format_name, magic):
    paths = {
        "summary-out": tmp_path / "summary.tsv",
        "model-out": tmp_path / "model.json",
        "simulations-out": tmp_path / "simulations.tsv",
        "clades-out": tmp_path / "clades.tsv",
        "figure-out": tmp_path / ("figure." + format_name),
    }
    flags = [item for role, path in paths.items() for item in ("--" + role, str(path))]
    main(command(tmp_path, flags))
    curve = pd.read_csv(StringIO(capsys.readouterr().out), sep="\t")
    summary = pd.read_csv(paths["summary-out"], sep="\t").iloc[0]
    simulations = (
        pd.read_csv(paths["simulations-out"], sep="\t")
        .relative_disparity.to_numpy()
        .reshape(9, -1)
    )
    assert curve.bm_median.to_numpy() == pytest.approx(np.median(simulations, axis=0))
    assert summary.mdi == pytest.approx(
        curve_area(
            curve.relative_time.to_numpy(),
            curve.relative_disparity.to_numpy() - np.median(simulations, axis=0),
        ),
        abs=1e-11,
    )
    model = json.loads(paths["model-out"].read_text())
    assert model["used_taxa"] == list("ABCDEFGH")
    assert model["mdi"] == pytest.approx(summary.mdi)
    assert len(pd.read_csv(paths["clades-out"], sep="\t")) == 15
    assert paths["figure-out"].read_bytes().startswith(magic)


def test_observed_only_rank_deficiency_and_json(tmp_path, capsys):
    path = tmp_path / "model.json"
    data = np.column_stack([DATA[:, 0], DATA[:, 0] * 2])
    main(
        command(
            tmp_path,
            [
                "--n-sim",
                "0",
                "--model-out",
                str(path),
                "--figure-out",
                str(tmp_path / "observed.png"),
            ],
            data=data,
        )
    )
    curve = pd.read_csv(StringIO(capsys.readouterr().out), sep="\t")
    assert curve.bm_median.isna().all()
    assert json.loads(path.read_text())["mdi"] is None
    with pytest.raises(ValueError, match="independent"):
        main(command(tmp_path, data=data))


def test_missing_drop_rebases_crown_and_keeps_original_ids(tmp_path, capsys):
    data = DATA.copy()
    data[4:, :] = np.nan
    model, clades = tmp_path / "model.json", tmp_path / "clades.tsv"
    with pytest.raises(ValueError, match="complete"):
        main(command(tmp_path, data=data))
    main(
        command(
            tmp_path,
            [
                "--missing",
                "drop",
                "--model-out",
                str(model),
                "--clades-out",
                str(clades),
            ],
            data=data,
        )
    )
    capsys.readouterr()
    metadata = json.loads(model.read_text())
    assert metadata["used_taxa"] == list("ABCD")
    assert metadata["excluded_taxa"] == list("EFGH")
    assert metadata["crown_age"] == 1.5
    assert metadata["crown_root_branch_id"] != 0
    nodes = pd.read_csv(clades, sep="\t")
    assert len(nodes) == 7
    assert nodes.iloc[0].parent_branch_id == -1
    assert nodes.iloc[0].branch_id == metadata["crown_root_branch_id"]


@pytest.mark.parametrize(
    "extra,message",
    [
        (["--n-sim", "1"], "n-sim"),
        (["--threads", "0"], "threads"),
        (["--seed", "-1"], "seed"),
        (["--ci-level", "nan"], "ci-level"),
        (["--mdi-range", "0.9,0.2"], "mdi-range"),
        (["--n-sim", "0", "--simulations-out", "unused.tsv"], "requires BM"),
        (["--figure-out", "unused.jpg"], "PNG"),
        (["--summary-out", "-"], "file path"),
    ],
)
def test_invalid_options(tmp_path, extra, message):
    with pytest.raises(ValueError, match=message):
        main(command(tmp_path, extra))


@pytest.mark.parametrize(
    "tree,message",
    [
        (TREE.replace("A:0.7", "A:0.8"), "ultrametric"),
        (TREE.replace("A:0.7", "A:-0.7"), "non-negative"),
        ("[&R](A:0,B:0,C:0);", "positive crown"),
    ],
)
def test_invalid_trees(tmp_path, tree, message):
    with pytest.raises(ValueError, match=message):
        main(command(tmp_path, tree=tree))


def test_star_is_descriptive_but_uninformative(tmp_path, capsys):
    summary = tmp_path / "summary.tsv"
    main(
        command(
            tmp_path,
            ["--summary-out", str(summary)],
            tree="[&R](A:1,B:1,C:1,D:1,E:1,F:1,G:1,H:1);",
        )
    )
    curve = pd.read_csv(StringIO(capsys.readouterr().out), sep="\t")
    assert curve.relative_disparity.tolist() == [1, 0, 0]
    row = pd.read_csv(summary, sep="\t").iloc[0]
    assert row.status == "uninformative_topology"
    assert row.mdi == 0


def test_constants_and_singular_bm():
    with pytest.raises(ValueError, match="nonconstant"):
        transform_traits(np.ones((4, 1)), "raw")
    with pytest.raises(ValueError, match="nonsingular"):
        fit_brownian(np.ones((8, 8)), DATA)
    with pytest.raises(ValueError, match="more tips"):
        fit_brownian(np.eye(2), DATA[:2])


def test_transaction_rolls_back_without_stdout(tmp_path, capsys, monkeypatch):
    import nwkit.dtt_figure

    output, figure = tmp_path / "curve.tsv", tmp_path / "figure.png"
    output.write_text("old table")
    figure.write_bytes(b"old figure")

    def fail(*args, **kwargs):
        raise RuntimeError("render failure")

    monkeypatch.setattr(nwkit.dtt_figure, "draw_dtt", fail)
    for extra in (["-o", str(output)], []):
        with pytest.raises(RuntimeError, match="render failure"):
            main(command(tmp_path, [*extra, "--figure-out", str(figure)]))
        assert capsys.readouterr().out == ""
        assert output.read_text() == "old table"
        assert figure.read_bytes() == b"old figure"


def test_path_alias_and_work_limits(tmp_path):
    with pytest.raises(ValueError, match="input"):
        main(command(tmp_path, ["-o", str(tmp_path / "traits.tsv")]))
    _, design = design_for()
    with pytest.raises(ValueError, match="retained"):
        validate_work(design, 1, 300000)
    with pytest.raises(ValueError, match="500,000"):
        validate_work(design, 1, 70000, True)
    with pytest.raises(ValueError, match="operations"):
        validate_work(design, 64, 10000)


def test_scalar_zero_length_polytomy_and_time_units(tmp_path, capsys):
    tree = "[&R](((A:1,B:1):0,C:1):1,(D:1,E:1,F:1):1,G:2,H:2);"
    output = tmp_path / "scalar.tsv"
    main(command(tmp_path, ["--columns", "x", "-o", str(output)], tree=tree))
    first = pd.read_csv(output, sep="\t")
    assert first.relative_time.tolist() == [0, 0, 0.5, 1]
    assert first.num_clades.tolist() == [1, 2, 0, 0]
    main(
        command(
            tmp_path,
            ["--columns", "x", "-o", str(output)],
            tree=tree.replace(":1", ":1e-90").replace(":2", ":2e-90"),
        )
    )
    second = pd.read_csv(output, sep="\t")
    assert first.relative_disparity.to_numpy() == pytest.approx(
        second.relative_disparity.to_numpy()
    )
    assert first.bm_median.to_numpy() == pytest.approx(second.bm_median.to_numpy())
    assert capsys.readouterr().out == ""


def test_bm_simulations_match_independent_dense_tip_distribution():
    tree, design = design_for()
    covariance = build_evolutionary_process(tree).tip_covariance(list("ABCDEFGH"))
    factor = np.array([[1.0]])
    branch_curves = simulate_curves(design, factor, 1000, seed=15)
    rng = np.random.default_rng(26)
    tips = rng.multivariate_normal(np.zeros(8), covariance, size=1000)
    dense_curves = np.array([dtt_curve(design, row[:, None]) for row in tips])
    # Two independent Monte Carlo samples must agree within their sampling error.
    standard_error = np.sqrt(
        (branch_curves.var(axis=0) + dense_curves.var(axis=0)) / 1000
    )
    assert np.all(
        np.abs(branch_curves.mean(axis=0) - dense_curves.mean(axis=0))
        <= 5 * standard_error + 1e-12
    )


def test_explicit_unrooted_is_rejected(tmp_path):
    with pytest.raises(ValueError, match="rooted tree"):
        main(command(tmp_path, ["--input-rooted", "no"]))


def test_trait_tree_panels_use_retained_original_values(tmp_path, monkeypatch, capsys):
    from matplotlib.figure import Figure

    figures = []
    save = Figure.savefig

    def capture(figure, *args, **kwargs):
        figures.append(figure)
        return save(figure, *args, **kwargs)

    monkeypatch.setattr(Figure, "savefig", capture)
    data = DATA.copy()
    data[4:] = np.nan
    main(
        command(
            tmp_path,
            [
                "--figure-layout",
                "trees",
                "--missing",
                "drop",
                "--scale",
                "standardize",
                "--figure-out",
                str(tmp_path / "trees.png"),
            ],
            data=data,
        )
    )
    capsys.readouterr()
    panels = [ax for ax in figures[0].axes if ax.get_title(loc="left") in {"x", "y"}]
    assert len(panels) == 2
    dtt_axis = next(
        ax
        for ax in figures[0].axes
        if ax.get_title(loc="left") == "Disparity through time"
    )
    for panel in panels:
        endpoints = np.array([[0, 0], [1, 0]])
        assert panel.transData.transform(endpoints)[:, 0] == pytest.approx(
            dtt_axis.transData.transform(endpoints)[:, 0], abs=1e-6
        )
    for index, panel in enumerate(panels):
        assert np.asarray(panel.collections[0].get_array()) == pytest.approx(
            (DATA[:4, index] - DATA[:4, index].min()) / np.ptp(DATA[:4, index])
        )
        assert [text.get_text().split()[0] for text in panel.texts] == list("ABCD")
        assert np.asarray(panel.collections[0].get_offsets())[:, 0] == pytest.approx(
            np.ones(4)
        )


@pytest.mark.parametrize("scale", ["raw", "standardize"])
@pytest.mark.parametrize("suffix", ["png", "pdf", "svg"])
def test_heatmap_order_scaling_alignment_and_unchanged_analysis(
    tmp_path, capsys, monkeypatch, scale, suffix
):
    from matplotlib.figure import Figure

    figures = []
    save = Figure.savefig

    def capture(figure, *args, **kwargs):
        result = save(figure, *args, **kwargs)
        figures.append(figure)
        return result

    monkeypatch.setattr(Figure, "savefig", capture)
    base = [
        "dtt",
        "-i",
        TREE,
        "--trait",
        str(EXAMPLE / "traits.tsv"),
        "--columns",
        "size,shape,performance",
        "--n-sim",
        "0",
    ]
    main(base)
    baseline = capsys.readouterr().out
    main(
        [
            *base,
            "--figure-layout",
            "heatmap",
            "--figure-columns",
            "performance,size",
            "--figure-scale",
            scale,
            "--figure-out",
            str(tmp_path / ("heat." + suffix)),
        ]
    )
    assert capsys.readouterr().out == baseline
    axes = {ax.get_title(loc="left"): ax for ax in figures[0].axes}
    heat = axes["Observed traits"]
    actual = np.asarray(heat.images[0].get_array())
    source = pd.read_csv(EXAMPLE / "traits.tsv", sep="\t")[
        ["performance", "size"]
    ].to_numpy()
    expected = (
        (source - source.min()) / np.ptp(source)
        if scale == "raw"
        else (source - source.mean(axis=0)) / source.std(axis=0, ddof=1)
    )
    assert actual == pytest.approx(expected)
    assert [tick.get_text() for tick in heat.get_xticklabels()] == [
        "performance",
        "size",
    ]
    tree = axes["Retained tree"]
    endpoints = np.array([[0, 0], [1, 0]])
    assert tree.transData.transform(endpoints)[:, 0] == pytest.approx(
        axes["Disparity through time"].transData.transform(endpoints)[:, 0], abs=1e-6
    )
    # First input row is drawn at the top, on the same line as tip A.
    assert heat.images[0].origin == "upper"
    rows = np.column_stack([np.zeros(8), np.arange(7, -1, -1)])
    assert heat.transData.transform(rows)[:, 1] == pytest.approx(
        tree.transData.transform(rows)[:, 1], abs=1e-6
    )
    assert [label.get_text() for label in tree.texts] == list("ABCDEFGH")


@pytest.mark.parametrize(
    "displayed,expected_heatmap", [(None, True), ("size,shape", True), ("size", True)]
)
def test_default_heatmap_for_any_displayed_count(
    tmp_path, capsys, monkeypatch, displayed, expected_heatmap
):
    from matplotlib.figure import Figure

    figures = []
    original = Figure.savefig

    def capture(figure, *args, **kwargs):
        figures.append(figure)
        return original(figure, *args, **kwargs)

    monkeypatch.setattr(Figure, "savefig", capture)
    flags = ["--figure-columns", displayed] if displayed else []
    main(
        [
            "dtt",
            "-i",
            TREE,
            "--trait",
            str(EXAMPLE / "traits.tsv"),
            "--columns",
            "size,shape,performance",
            "--n-sim",
            "0",
            "--figure-out",
            str(tmp_path / "auto.png"),
            *flags,
        ]
    )
    capsys.readouterr()
    assert any(ax.images for ax in figures[0].axes) == expected_heatmap
    assert "3 traits analyzed" in figures[0]._suptitle.get_text()


@pytest.mark.parametrize("columns", ["absent", "x,x", "", "x,"])
def test_invalid_display_columns(tmp_path, columns):
    with pytest.raises(ValueError):
        main(command(tmp_path, ["--figure-columns", columns]))

"""Independent R, algebraic and CLI verification of phylogenetic PCA."""

import json
from io import StringIO

import numpy as np
import pandas as pd
import pytest
from scipy.optimize import minimize

from nwkit.cli import main
from nwkit.evolution import build_evolutionary_covariance
from nwkit.phylogenetic_pca import fit_pca
from nwkit.util import read_tree
from tests.test_signal import TREE, VALUES

MIXED = np.array([1, 4, 2, 3, 1.5, 4.2, 1.8, 3.8])
DATA = np.column_stack([VALUES, MIXED])


def covariance():
    return build_evolutionary_covariance(read_tree(TREE, 0, False), list("ABCDEFGH"))


@pytest.mark.parametrize(
    "mode,expected",
    [
        ("cov", [1.9920819577, 0.462058344532]),
        ("corr", [1.32793034219, 0.672069657807]),
    ],
)
def test_phytools_230_reference(mode, expected):
    # R: phyl.pca(read.tree(TREE), DATA, method="BM", mode=mode).
    fit = fit_pca(covariance(), DATA, mode=mode)
    assert fit.eigenvalues == pytest.approx(expected, rel=1e-10)
    if mode == "cov":
        assert fit.loadings[:, 0] == pytest.approx([0.429554290824, 0.993968472466])
        assert fit.scores[0] == pytest.approx([-2.042969, -1.4631407])
    else:
        assert np.abs(fit.rotation) == pytest.approx(np.full((2, 2), np.sqrt(0.5)))
        assert np.abs(fit.loadings[:, 0]) == pytest.approx([0.8148406, 0.8148406])


def test_phytools_fixed_lambda_likelihood():
    fit = fit_pca(covariance(), DATA, model="LAMBDA", lambda_value=0.4)
    # phytools 2.3.0 plugs n-1 covariance into its likelihood. NWKIT profiles
    # the ML covariance (divisor n); the difference is constant in lambda.
    r_likelihood = -25.36392
    correction = 8 * 2 / 2 * np.log(8 / 7) - 2 / 2
    assert fit.log_likelihood == pytest.approx(r_likelihood + correction, abs=1e-5)
    assert fit.eigenvalues == pytest.approx([1.1359193, 0.5874151], abs=1e-7)


@pytest.mark.parametrize("mode", ["cov", "corr"])
def test_gls_orthogonality_and_reconstruction(mode):
    fit = fit_pca(covariance(), DATA, mode=mode)
    inverse = np.linalg.inv(fit.covariance)
    assert np.ones(8) @ inverse @ fit.scores == pytest.approx(np.zeros(2), abs=1e-12)
    assert fit.scores.T @ inverse @ fit.scores / 7 == pytest.approx(
        np.diag(fit.eigenvalues), abs=1e-12
    )
    restored = (fit.scores @ fit.rotation.T) * fit.scale + fit.center
    assert restored == pytest.approx(DATA, abs=1e-12)
    assert fit.rotation.T @ fit.rotation == pytest.approx(np.eye(2), abs=1e-12)
    assert np.abs(fit.loadings).max() <= 1 + 1e-12


def test_lambda_joint_fit_against_independent_matrix_normal():
    c = covariance()
    fit = fit_pca(c, DATA, model="LAMBDA")

    def objective(params):
        lam = params[0]
        matrix = lam * c + (1 - lam) * np.diag(np.diag(c))
        inverse = np.linalg.inv(matrix)
        center = (np.ones(8) @ inverse @ DATA) / inverse.sum()
        residual = DATA - center
        trait = residual.T @ inverse @ residual / 8
        return 0.5 * (
            16 * (np.log(2 * np.pi) + 1)
            + 2 * np.linalg.slogdet(matrix)[1]
            + 8 * np.linalg.slogdet(trait)[1]
        )

    candidates = [
        minimize(objective, [start], bounds=[(0, 1)], method="L-BFGS-B")
        for start in [0, 0.25, 0.5, 0.75, 1]
    ]
    best = min(candidates, key=lambda item: item.fun)
    assert fit.log_likelihood == pytest.approx(-best.fun, abs=1e-7)
    assert fit.lambda_value == pytest.approx(best.x[0], abs=1e-5)


def test_correlation_pca_and_lambda_are_invariant_to_units():
    first = fit_pca(covariance(), DATA, model="LAMBDA", mode="corr")
    second = fit_pca(
        covariance(),
        DATA * np.array([1e10, 1e-8]) + [1e11, 1e-7],
        model="LAMBDA",
        mode="corr",
    )
    assert second.lambda_value == pytest.approx(first.lambda_value, abs=1e-5)
    assert second.eigenvalues == pytest.approx(first.eigenvalues, abs=1e-5)
    assert np.abs(second.scores) == pytest.approx(np.abs(first.scores), abs=1e-5)


def test_rank_deficiency_is_supported_only_for_fixed_tree_model():
    data = np.column_stack([VALUES, VALUES * 2])
    fit = fit_pca(covariance(), data)
    assert fit.status == "rank_deficient"
    assert fit.scores.shape == (8, 1)
    assert (fit.scores @ fit.rotation.T) + fit.center == pytest.approx(data)
    with pytest.raises(ValueError, match="full-rank"):
        fit_pca(covariance(), data, model="LAMBDA")


def command(tmp_path, extra=(), data=None):
    path = tmp_path / "traits.tsv"
    if data is None:
        data = pd.DataFrame({"leaf_name": list("ABCDEFGH"), "x": VALUES, "y": MIXED})
    data.to_csv(path, sep="\t", index=False)
    return ["pca", "-i", TREE, "--trait", str(path), "--columns", "x,y", *extra]


def test_cli_all_artifacts_and_asr(tmp_path, capsys):
    roles = {
        "loadings-out": "loadings.tsv",
        "eigenvalues-out": "eigenvalues.tsv",
        "model-out": "model.json",
        "ancestral-out": "ancestors.tsv",
        "figure-out": "figure.png",
    }
    flags = [
        token
        for role, name in roles.items()
        for token in ["--" + role, str(tmp_path / name)]
    ]
    main(command(tmp_path, flags))
    scores = pd.read_csv(StringIO(capsys.readouterr().out), sep="\t")
    assert list(scores.columns) == ["leaf_name", "PC1", "PC2"]
    model = json.loads((tmp_path / "model.json").read_text())
    assert model["used_taxa"] == list("ABCDEFGH")
    assert (tmp_path / "figure.png").read_bytes().startswith(b"\x89PNG")
    ancestors = pd.read_csv(tmp_path / "ancestors.tsv", sep="\t")
    assert len(ancestors) == 30
    tips = ancestors[ancestors.node_class == "leaf"].pivot(
        index="name", columns="component", values="mean"
    )
    assert tips.to_numpy() == pytest.approx(scores.set_index("leaf_name").to_numpy())
    assert (ancestors[ancestors.node_class == "leaf"].variance == 0).all()
    assert ancestors.loc[ancestors.node_class == "root", "mean"].abs().max() < 1e-10
    assert ancestors.variance.min() >= 0


def test_ancestral_covariance_matches_gaussian_conditioning(tmp_path, capsys):
    from nwkit.cli import parser
    from nwkit.evolution import build_evolutionary_process
    from nwkit.pca import _read_data, ancestral_scores

    args = parser.parse_args(command(tmp_path))
    tree, _, names, _, values = _read_data(args)
    process = build_evolutionary_process(tree)
    leaves = {n.name: n for n in tree.leaves()}
    nodes = list(tree.traverse())
    joint = process.covariance(nodes + [leaves[name] for name in names])
    c = joint[len(nodes) :, len(nodes) :]
    cross = joint[: len(nodes), len(nodes) :]
    inverse = np.linalg.inv(c)
    fit = fit_pca(c, values)
    table = ancestral_scores(tree, names, fit, 0.95)
    expected_mean = cross @ inverse @ fit.scores
    expected_var = (
        np.diag(joint[: len(nodes), : len(nodes)] - cross @ inverse @ cross.T)
        + (1 - cross @ inverse @ np.ones(8)) ** 2 / inverse.sum()
    )
    for j in range(2):
        selected = table[table.component == f"PC{j + 1}"].set_index("branch_id")
        for i, node in enumerate(nodes):
            row = selected.loc[node.props["_nwkit_pca_input_branch_id"]]
            assert row["mean"] == pytest.approx(expected_mean[i, j], abs=1e-10)
            assert row.variance == pytest.approx(
                expected_var[i] * fit.eigenvalues[j], abs=1e-10
            )


def test_drop_preserves_input_ids_and_omits_incomplete_tips(tmp_path, capsys):
    data = pd.DataFrame({"leaf_name": list("ABCDEFGH"), "x": VALUES, "y": MIXED})
    data.loc[4:, "y"] = np.nan
    with pytest.raises(ValueError, match="complete trait"):
        main(command(tmp_path, data=data))
    main(
        command(
            tmp_path,
            [
                "--missing",
                "drop",
                "--model-out",
                str(tmp_path / "model.json"),
                "--ancestral-out",
                str(tmp_path / "nodes.tsv"),
            ],
            data,
        )
    )
    result = pd.read_csv(StringIO(capsys.readouterr().out), sep="\t")
    assert result.leaf_name.tolist() == list("ABCD")
    model = json.loads((tmp_path / "model.json").read_text())
    assert model["root_branch_id"] == 1
    assert model["excluded_taxa"] == list("EFGH")
    nodes = pd.read_csv(tmp_path / "nodes.tsv", sep="\t")
    assert set(nodes[nodes.node_class == "leaf"]["name"]) == set("ABCD")
    assert set(nodes.branch_id) == {1, 3, 4, 7, 8, 9, 10}


def test_render_failure_rolls_back_all_outputs(tmp_path, monkeypatch):
    import nwkit.pca_figure

    score, model, figure = [
        tmp_path / name for name in ["scores.tsv", "model.json", "figure.png"]
    ]
    for path in (score, model, figure):
        path.write_text("original")

    def fail(*args):
        raise ValueError("render failed")

    monkeypatch.setattr(nwkit.pca_figure, "draw_pca", fail)
    with pytest.raises(ValueError, match="render failed"):
        main(
            command(
                tmp_path,
                [
                    "-o",
                    str(score),
                    "--model-out",
                    str(model),
                    "--figure-out",
                    str(figure),
                ],
            )
        )
    assert all(path.read_text() == "original" for path in (score, model, figure))


@pytest.mark.parametrize(
    "extra",
    [
        ["--lambda-value", ".5"],
        ["--model", "LAMBDA", "--lambda-value", "nan"],
        ["--ci-level", "0"],
        ["--figure-out", "bad.jpg"],
        ["--model-out", "-"],
        ["--columns", "x,x"],
        ["--input-rooted", "no"],
    ],
)
def test_invalid_options(tmp_path, extra):
    with pytest.raises(ValueError):
        main(command(tmp_path, extra))


def test_output_aliases_and_trait_stdin(tmp_path, monkeypatch, capsys):
    output = tmp_path / "both.tsv"
    with pytest.raises(ValueError, match="distinct"):
        main(command(tmp_path, ["-o", str(output), "--loadings-out", str(output)]))
    monkeypatch.setattr(
        "sys.stdin", StringIO("leaf_name\tx\ty\nA\t1\t2\nB\t2\t4\nC\t3\t1\n")
    )
    main(["pca", "-i", "((A:1,B:1):1,C:2);", "--trait", "-", "--columns", "x,y"])
    assert len(pd.read_csv(StringIO(capsys.readouterr().out), sep="\t")) == 3


@pytest.mark.parametrize("mode", ["cov", "corr"])
def test_branch_unit_rescaling(mode):
    first = fit_pca(covariance(), DATA, mode=mode)
    second = fit_pca(covariance() * 10, DATA, mode=mode)
    assert second.eigenvalues == pytest.approx(
        first.eigenvalues / (10 if mode == "cov" else 1)
    )
    assert second.scores == pytest.approx(
        first.scores * (1 if mode == "cov" else np.sqrt(10))
    )


@pytest.mark.parametrize("extension,magic", [("pdf", b"%PDF-"), ("svg", b"<?xml")])
def test_figure_formats_and_zero_lambda_ancestors(tmp_path, capsys, extension, magic):
    picture = tmp_path / f"figure.{extension}"
    nodes = tmp_path / "ancestors.tsv"
    main(
        command(
            tmp_path,
            [
                "--model",
                "LAMBDA",
                "--lambda-value",
                "0",
                "--figure-out",
                str(picture),
                "--ancestral-out",
                str(nodes),
                "--figure-tip-labels",
                "no",
            ],
        )
    )
    capsys.readouterr()
    assert picture.read_bytes().startswith(magic)
    table = pd.read_csv(nodes, sep="\t")
    ancestors = table[table.node_class != "leaf"]
    assert ancestors["mean"].abs().max() < 1e-10
    for _, group in ancestors.groupby("component"):
        assert group.variance.max() == pytest.approx(group.variance.min())


def test_more_traits_than_tips_and_star_model():
    rng = np.random.default_rng(55)
    x = rng.normal(size=(4, 7))
    fit = fit_pca(np.eye(4), x, mode="corr")
    assert fit.scores.shape == (4, 3)
    assert fit.status == "rank_deficient"
    assert (fit.scores @ fit.rotation.T) * fit.scale + fit.center == pytest.approx(x)
    with pytest.raises(ValueError, match="unidentifiable"):
        fit_pca(np.eye(4), x, model="LAMBDA")


def test_input_file_and_audit_aliases_are_protected(tmp_path):
    arguments = command(tmp_path)
    data = tmp_path / "traits.tsv"
    original = data.read_bytes()
    with pytest.raises(ValueError, match="overwrite input"):
        main([*arguments, "--loadings-out", str(data)])
    assert data.read_bytes() == original
    alias = tmp_path / "alias.tsv"
    with pytest.raises(ValueError):
        main([*arguments, "--loadings-out", str(alias), "--audit", str(alias)])

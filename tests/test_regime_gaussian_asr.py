import math

import numpy as np
import pandas as pd
import pytest

from nwkit.asr_regimes import RegimeAssignment
from nwkit.cli import main
from nwkit.gaussian_inference import _relative_parameter_rank
from nwkit.ou_asr import compute_ou_marginals
from nwkit.regime_gaussian_asr import (
    build_regime_ou_process,
    compute_regime_ou_marginals,
    regime_parameter_columns,
)
from nwkit.util import read_tree


def tree_from(source):
    return read_tree(source, "1", True, quiet=True, rooted="yes")


def assignment(tree, values):
    nodes = list(tree.traverse())
    return RegimeAssignment(
        tuple(dict.fromkeys(values)), dict(zip(nodes, values, strict=True)), "memory"
    )


@pytest.mark.parametrize(
    "model, columns",
    [
        ("OUM", ("theta",)),
        ("OUMA", ("theta", "alpha")),
        ("OUMV", ("theta", "sigma2")),
        ("OUMVA", ("theta", "alpha", "sigma2")),
    ],
)
def test_regime_parameter_columns(model, columns):
    assert regime_parameter_columns(model) == columns


def test_fixed_equal_regime_process_reduces_to_single_optimum_ou():
    tree = tree_from("((A:0.4,B:1.2)I:0.7,C:1.5,D:0.8)R;")
    values = {"A": -1.0, "B": 2.0, "C": 4.0, "D": None}
    errors = {"A": 0.2, "B": 0.3, "C": 0.4}
    regimes = assignment(tree, ["cold", "cold", "warm", "cold", "warm", "warm"])
    expected, expected_fit = compute_ou_marginals(
        tree,
        values,
        alpha=0.7,
        sigma2=1.4,
        theta=0.3,
        standard_errors=errors,
    )
    fixed = {
        regime: {"theta": 0.3, "alpha": 0.7, "sigma2": 1.4}
        for regime in regimes.regimes
    }
    observed, fit = compute_regime_ou_marginals(
        tree,
        values,
        regimes,
        model="OUMVA",
        regime_parameters=fixed,
        standard_errors=errors,
    )
    assert fit.log_likelihood == pytest.approx(expected_fit.log_likelihood, abs=2e-12)
    for node in tree.traverse():
        assert observed[node].mean == pytest.approx(expected[node].mean, abs=2e-12)
        assert observed[node].variance == pytest.approx(
            expected[node].variance, abs=2e-12
        )


@pytest.mark.parametrize("model", ["OUM", "OUMA", "OUMV", "OUMVA"])
def test_single_regime_ou_variants_use_the_canonical_ou_fit(model):
    tree = tree_from("((A:0.4,B:1.2)I:0.7,C:1.5,D:0.8,E:2)R;")
    values = {"A": -1.0, "B": 2.0, "C": 4.0, "D": 1.0, "E": 3.0}
    regimes = assignment(tree, ["shared"] * len(list(tree.traverse())))
    expected_posterior, expected_fit = compute_ou_marginals(tree, values, theta=0.3)
    observed_posterior, observed_fit = compute_regime_ou_marginals(
        tree,
        values,
        regimes,
        model=model,
        theta=0.3,
    )
    assert observed_fit.alpha_by_regime["shared"] == expected_fit.alpha
    assert observed_fit.sigma2_by_regime["shared"] == expected_fit.sigma2
    assert observed_fit.theta_by_regime["shared"] == expected_fit.theta
    assert observed_fit.log_likelihood == expected_fit.log_likelihood
    assert observed_fit.fit_status == expected_fit.fit_status
    assert observed_fit.optimizer_starts == expected_fit.optimizer_starts
    assert "single-regime reduction" in observed_fit.optimizer_message
    for node in tree.traverse():
        assert observed_posterior[node].mean == expected_posterior[node].mean
        assert observed_posterior[node].variance == expected_posterior[node].variance


def test_parameter_rank_tolerance_does_not_grow_with_summary_row_count():
    normalized = np.zeros((100_000, 2), dtype=float)
    normalized[0, 0] = 1.0
    normalized[0, 1] = 1.0
    normalized[1, 1] = 1e-4
    normalized[:, 1] /= np.linalg.norm(normalized[:, 1])
    assert _relative_parameter_rank(normalized) == 2


def test_regime_process_uses_incoming_edge_alpha_sigma_and_theta():
    tree = tree_from("(A:2,B:3)R;")
    regimes = assignment(tree, ["root", "slow", "fast"])
    process = build_regime_ou_process(
        tree,
        regimes,
        alpha_by_regime={"root": 0.5, "slow": 0.25, "fast": 1.5},
        sigma2_by_regime={"root": 2.0, "slow": 1.0, "fast": 4.0},
        theta_by_regime={"root": 0.0, "slow": -2.0, "fast": 3.0},
    )
    slow = process.transitions[tree["A"]]
    fast = process.transitions[tree["B"]]
    assert slow.slope == pytest.approx(math.exp(-0.25 * 2.0))
    assert slow.intercept == pytest.approx(-2.0 * (1.0 - slow.slope))
    assert slow.variance == pytest.approx(
        1.0 / (2.0 * 0.25) * (1.0 - math.exp(-2.0 * 0.25 * 2.0))
    )
    assert fast.slope == pytest.approx(math.exp(-1.5 * 3.0))
    assert fast.intercept == pytest.approx(3.0 * (1.0 - fast.slope))
    assert process.root.variance == pytest.approx(2.0)


def test_free_regime_parameter_rejects_regime_outside_observed_subtrees():
    tree = tree_from("(A:1,B:1,C:1,D:1,E:1)R;")
    regimes = RegimeAssignment(
        ("root", "background", "ghost"),
        {
            tree: "root",
            tree["A"]: "background",
            tree["B"]: "background",
            tree["C"]: "background",
            tree["D"]: "background",
            tree["E"]: "ghost",
        },
        "memory",
    )
    with pytest.raises(ValueError, match="unrepresented regime.*ghost"):
        compute_regime_ou_marginals(
            tree,
            {"A": 0.0, "B": 1.0, "C": 2.0, "D": 3.0, "E": None},
            regimes,
            model="OUMV",
            alpha=0.7,
            theta=0.0,
        )


def test_oumva_rejects_root_only_alpha_sigma_confounding():
    tree = tree_from("(A:1,B:1,C:1,D:1,E:1,F:1)R;")
    regimes = RegimeAssignment(
        ("root", "background"),
        {tree: "root", **{leaf: "background" for leaf in tree.leaves()}},
        "memory",
    )
    with pytest.raises(ValueError, match="Root-regime alpha and sigma2 are confounded"):
        compute_regime_ou_marginals(
            tree,
            {name: float(index) for index, name in enumerate("ABCDEF")},
            regimes,
            model="OUMVA",
            theta=0.0,
        )


def test_oum_jointly_estimates_shared_process_and_regime_optima():
    tree = tree_from(
        "(((A:.4,B:.5)I:.3,(C:.6,D:.7)J:.2)K:.5,"
        "((E:.4,F:.5)L:.3,(G:.6,H:.7)M:.2)N:.5)R;"
    )
    warm_tips = set("EFGH")
    by_node = {
        node: (
            "warm"
            if node is not tree
            and {str(leaf.name) for leaf in node.leaves()} <= warm_tips
            else "cold"
        )
        for node in tree.traverse()
    }
    regimes = RegimeAssignment(("cold", "warm"), by_node, "memory")
    _, fit = compute_regime_ou_marginals(
        tree,
        dict(
            zip(
                "ABCDEFGH",
                (-0.2, 0.1, 0.3, 0.5, 2.8, 3.1, 3.4, 3.6),
                strict=True,
            )
        ),
        regimes,
        model="OUM",
    )
    assert fit.alpha > 0.0
    assert fit.sigma2 > 0.0
    assert fit.theta_by_regime["warm"] > fit.theta_by_regime["cold"]
    assert fit.optimizer_converged_starts > 0


@pytest.mark.integration
@pytest.mark.parametrize("model", ["OUM", "OUMA", "OUMV", "OUMVA"])
def test_regime_ou_cli_fixed_parameters(model, tmp_path):
    traits = tmp_path / "traits.tsv"
    regime_map = tmp_path / "regimes.tsv"
    parameters = tmp_path / "parameters.tsv"
    output = tmp_path / "asr.tsv"
    metadata = tmp_path / "model.tsv"
    traits.write_text("leaf_name\tvalue\nA\t0\nB\t1\nC\t3\nD\t4\n")
    regime_map.write_text(
        "branch_id\tregime\n0\tlow\n1\tlow\n2\thigh\n3\tlow\n4\tlow\n5\thigh\n6\thigh\n"
    )
    columns = regime_parameter_columns(model)
    rows = []
    for regime, theta, alpha, sigma2 in (
        ("low", 0.0, 0.7, 1.2),
        ("high", 3.0, 1.1, 2.0),
    ):
        values = {"theta": theta, "alpha": alpha, "sigma2": sigma2}
        rows.append("\t".join([regime, *(str(values[column]) for column in columns)]))
    parameters.write_text(
        "regime\t" + "\t".join(columns) + "\n" + "\n".join(rows) + "\n"
    )
    options = []
    if "alpha" not in columns:
        options += ["--alpha", "0.7"]
    if "sigma2" not in columns:
        options += ["--sigma2", "1.2"]
    main(
        [
            "asr",
            "-i",
            "((A:1,B:1):1,(C:1,D:1):1)R;",
            "--input-rooted",
            "yes",
            "--trait",
            str(traits),
            "--state-column",
            "value",
            "--model",
            model,
            "--regime-map",
            str(regime_map),
            "--regime-parameters",
            str(parameters),
            "--model-out",
            str(metadata),
            "-o",
            str(output),
            *options,
        ]
    )
    assert len(pd.read_csv(output, sep="\t")) == 7
    row = pd.read_csv(metadata, sep="\t").iloc[0]
    assert row["model"] == model
    assert row["root_prior"] == "stationary"
    assert row["regimes"] == "low,high"

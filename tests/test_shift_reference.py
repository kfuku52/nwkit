"""Same-parameter checks against the separately installed kfl1ou engine."""

import os
import subprocess

import numpy as np
import pandas as pd
import pytest

from nwkit.shift_backend import R_SCRIPT
from nwkit.shift_reference import evaluate_shift_model
from nwkit.util import read_tree

TREE = "(((t0:1,t1:1):1,(t2:1,t3:1):1):1,((t4:1,t5:1):1,(t6:1,t7:1):1):1);"
VALUES = [1, 1.2, 0.9, 1.4, 6, 6.2, 5.8, 6.1]


def test_brownian_boundary_has_finite_mean_effects():
    tree = read_tree("((A:1,B:1):1,C:2);", "auto", True, quiet=True)
    node = next(n for n in tree.traverse() if set(n.leaf_names()) == {"A", "B"})
    mean, covariance, score = evaluate_shift_model(
        tree,
        observations=[2, 2.5, 1],
        alpha=0,
        sigma2=1,
        intercept=1,
        mean_effects={node: 1},
        root_model="OUfixedRoot",
    )
    np.testing.assert_allclose(mean, [2, 2, 1])
    np.testing.assert_allclose(covariance, [[2, 1, 0], [1, 2, 0], [0, 0, 2]])
    assert np.isfinite(score)


@pytest.mark.integration
@pytest.mark.skipif(
    not os.environ.get("NWKIT_TEST_RSCRIPT"),
    reason="Set NWKIT_TEST_RSCRIPT to Rscript with kfl1ou >= 3.0.9",
)
@pytest.mark.parametrize(
    "alpha,root_model",
    [
        (0.0, "OUfixedRoot"),
        (1e-7, "OUfixedRoot"),
        (1e-7, "OUrandomRoot"),
        (0.5, "OUfixedRoot"),
        (0.5, "OUrandomRoot"),
    ],
)
@pytest.mark.parametrize("nested", [False, True])
@pytest.mark.parametrize("with_error", [False, True])
def test_kfl1ou_fixed_parameter_equivalence(
    tmp_path, root_model, alpha, nested, with_error
):
    # R fits means/variance at a fixed alpha and configuration, then NWKIT
    # evaluates the very same parameter snapshot, with no Python optimization.
    (tmp_path / "tree.nwk").write_text(TREE)
    pd.DataFrame({"leaf_name": [f"t{i}" for i in range(8)], "value": VALUES}).to_csv(
        tmp_path / "trait.tsv", sep="\t", index=False
    )
    errors = (
        np.array([0.01, 0.05, 0, 0.07, 0.02, 0.1, 0.03, 0.08])
        if with_error
        else np.zeros(8)
    )
    table = pd.read_csv(tmp_path / "trait.tsv", sep="\t")
    if with_error:
        table["observation_variance"] = errors**2
    table.to_csv(tmp_path / "trait.tsv", sep="\t", index=False)
    start = R_SCRIPT.index("fit <- kfl1ou::estimate_shift_configuration(")
    end = R_SCRIPT.index("tr <- fit$tree", start)
    fit_code = """
tr <- dat$tree
wanted <- if (args[3] == "TRUE") c("t0/t1/t2/t3", "t0/t1") else character()
edge_keys <- vapply(tr$edge[,2], function(node) {
    tips <- if (node <= length(tr$tip.label)) tr$tip.label[node] else
        ape::extract.clade(tr, node)$tip.label
    paste(sort(tips), collapse="/")
}, "")
fit <- kfl1ou::fit_OU(tr, dat$Y,
    shift.configuration=which(edge_keys %in% wanted), criterion="BIC",
    root.model=args[2], alpha.lower=as.numeric(args[1]),
    alpha.upper=as.numeric(args[1]), alpha.starting.value=as.numeric(args[1]),
    compute.hessian=FALSE, input_error=input.error, measurement_error=FALSE)
fit$profile <- list(scores=fit$score, configurations=list(fit$shift.configuration))
fit$search.diagnostics <- list(strategy="fixed")
"""
    script = (
        R_SCRIPT[:start]
        + fit_code
        + R_SCRIPT[end:]
        + """
S <- kfl1ou::sqrt_OU_covariance(fit$tree, alpha=fit$alpha,
    root.model=args[2])$sqrtSigma
covariance <- fit$sigma2 * tcrossprod(S)
if (!is.null(input.error)) covariance <- covariance + diag(input.error[fit$tree$tip.label,1])
write_tsv(as.data.frame(covariance), "covariance.tsv")
"""
    )
    (tmp_path / "reference.R").write_text(script)
    result = subprocess.run(
        [
            os.environ["NWKIT_TEST_RSCRIPT"],
            "--vanilla",
            str(tmp_path / "reference.R"),
            str(alpha),
            root_model,
            "TRUE" if nested else "FALSE",
        ],
        cwd=tmp_path,
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0, result.stderr
    tree = read_tree(TREE, "auto", True, quiet=True)
    nodes = {"/".join(sorted(n.leaf_names())): n for n in tree.traverse()}
    effects = pd.read_csv(tmp_path / "shifts.tsv", sep="\t")
    parameters = pd.read_csv(tmp_path / "model.tsv", sep="\t").iloc[0]
    assert parameters.alpha == pytest.approx(alpha, abs=1e-14)
    tips = pd.read_csv(tmp_path / "tips.tsv", sep="\t")
    order = [list(tips.token).index(n.name) for n in tree.leaves()]
    mean, covariance, score = evaluate_shift_model(
        tree,
        alpha=float(parameters.alpha),
        sigma2=float(parameters.sigma2),
        intercept=float(parameters.intercept),
        mean_effects={nodes[r.clade]: r.mean_effect for r in effects.itertuples()},
        root_model=root_model,
        observations=VALUES,
        standard_errors=errors,
    )
    expected_cov = pd.read_csv(tmp_path / "covariance.tsv", sep="\t").to_numpy()[
        np.ix_(order, order)
    ]
    np.testing.assert_allclose(
        mean, tips.predicted.to_numpy()[order], rtol=1e-7, atol=1e-8
    )
    np.testing.assert_allclose(covariance, expected_cov, rtol=1e-7, atol=1e-8)
    # Stationary-root covariance is ill-conditioned at tiny alpha; retain a
    # separate absolute likelihood tolerance without weakening other cases.
    tolerance = 1e-5 if root_model == "OUrandomRoot" and 0 < alpha < 1e-6 else 1e-7
    assert score == pytest.approx(parameters.log_likelihood, abs=tolerance)
    # Changing the time unit must preserve all observables when alpha/sigma2
    # are transformed inversely. No tree-height normalization is involved.
    for node in tree.traverse():
        if not node.is_root:
            node.dist *= 1000
    scaled = evaluate_shift_model(
        tree,
        alpha=float(parameters.alpha) / 1000,
        sigma2=float(parameters.sigma2) / 1000,
        intercept=float(parameters.intercept),
        mean_effects={nodes[r.clade]: r.mean_effect for r in effects.itertuples()},
        root_model=root_model,
        observations=VALUES,
        standard_errors=errors,
    )
    np.testing.assert_allclose(scaled[0], mean, rtol=1e-10, atol=1e-10)
    np.testing.assert_allclose(scaled[1], covariance, rtol=1e-10, atol=1e-10)
    assert scaled[2] == pytest.approx(score, abs=tolerance)


def test_random_root_boundary_is_not_silently_replaced():
    tree = read_tree(TREE, "auto", True, quiet=True)
    with pytest.raises(ValueError, match="undefined at alpha=0"):
        evaluate_shift_model(
            tree,
            observations=VALUES,
            alpha=0,
            sigma2=1,
            intercept=1,
            mean_effects={},
            root_model="OUrandomRoot",
        )

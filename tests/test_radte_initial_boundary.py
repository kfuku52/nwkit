"""Regression for a feasible family whose branch initializer collapses an edge."""

from pathlib import Path

import numpy as np
import pytest

from nwkit.cli import parser
from nwkit.radte import run_dating
from nwkit.radte_inputs import read_inputs
from nwkit.radte_model import solve_problem


def test_sequence_fit_starts_inside_chronology_when_branch_initializer_collapses():
    root = Path(__file__).parent / "data/radte-initial-boundary"
    args = parser.parse_args(
        [
            "radte",
            "--gene-tree",
            str(root / "gene.nwk"),
            "--species-tree",
            str(root / "species.nwk"),
            "--species-map-tsv",
            str(root / "mapping.tsv"),
            "--reconcile",
            "lca",
            "--alignment",
            str(root / "alignment.fa"),
            "--substitution-model",
            "gy94",
            "--gamma-categories",
            "4",
            "--max-age",
            "1000",
            "--maxiter",
            "1000",
            "--seed",
            "20830099",
            "--uncertainty",
            "none",
            "--out-prefix",
            "unused",
        ]
    )
    chronology = read_inputs(args)
    fit, problem, _, _ = run_dating(chronology, args)
    assert "sequence_initial_ages_reset_from_duration_boundary" in fit.diagnostics
    assert problem.feasible(fit.parameters)
    assert np.all(chronology.durations(fit.ages) >= chronology.min_duration * 0.99)
    assert np.all(fit.ages >= chronology.lower)
    assert np.all(fit.ages <= chronology.upper)
    # Independently restart the unchanged exact objective from the generic
    # chronology interior, with different seeds. This is not a truth-age check.
    reference, _ = solve_problem(problem, starts=5, seed=817, maxiter=1000)
    assert problem.value_gradient(reference)[0] == pytest.approx(
        fit.objective, abs=1e-5
    )
    np.testing.assert_allclose(problem.unpack_ages(reference), fit.ages, atol=1e-5)

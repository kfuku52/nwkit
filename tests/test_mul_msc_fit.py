import json
import math
from dataclasses import replace

import numpy as np
import pandas as pd
import pytest

from nwkit.cli import main
from nwkit.mul_msc_fit import (
    FitSettings,
    MscObjective,
    fit_candidate,
    local_diagnostics,
    optimize_coordinates,
)
from nwkit.mul_msc_model import dated_candidates, score_gene
from tests.test_mul_coalescent import all_topologies, oracle_distribution
from tests.test_mul_msc import SPECIES, invoke
from tests.test_mul_reconcile import parser, tree


def bounded_candidate(h2="A", bounds=(0.1, 1.9)):
    return next(
        c
        for c in dated_candidates(tree(SPECIES), "X", h2, None, age_bounds=bounds)[0]
        if c.status == "evaluated"
    )


def test_analytic_three_tip_ridge_is_locally_unidentified():
    candidate = bounded_candidate()
    genes = [tree(t + ";") for t in all_topologies(("a_A", "x1_X", "x2_X"))]
    settings = FitSettings("joint", (0.1, 1.9), (0.1, 10.0))
    objective = MscObjective(candidate, genes, parser(), settings)
    p = np.exp(objective.evaluate([0.5, 0.5])[1])[objective.indices]
    assert sorted(p) == pytest.approx(
        sorted(
            [math.exp(-0.5) / 3, 0.5 - math.exp(-0.5) / 6, 0.5 - math.exp(-0.5) / 6]
        ),
        abs=2e-13,
    )
    alternative = [(0.5 - 0.1) / 1.8, math.log(1.5 / 0.1) / math.log(100)]
    assert objective.evaluate(alternative)[1][objective.indices] == pytest.approx(
        np.log(p), abs=2e-13
    )
    diagnostic = local_diagnostics(objective, np.array([0.5, 0.5]))
    assert diagnostic["rank"] == 1
    assert diagnostic["status"] == "locally-unidentified"
    assert diagnostic["step_ranks"] == [1, 1]


def test_global_time_and_ne_scaling_does_not_change_likelihood():
    gene = tree("((a_A,x1_X),(b_B,x2_X));")
    original = next(
        c
        for c in dated_candidates(tree(SPECIES), "X", "B", 1.0)[0]
        if c.status == "evaluated"
    )
    species = tree(SPECIES)
    for node in species.traverse():
        if node.dist is not None:
            node.dist *= 17
    scaled = next(
        c for c in dated_candidates(species, "X", "B", 17)[0] if c.status == "evaluated"
    )
    assert score_gene(gene, original, parser(), time_scale=2)[0] == pytest.approx(
        score_gene(gene, scaled, parser(), time_scale=34)[0], abs=2e-13
    )


def test_age_fit_matches_independent_three_tip_mle():
    genes = [tree("((a_A,x1_X),x2_X);"), tree("((x1_X,x2_X),a_A);")]
    settings = FitSettings("age", (0.1, 1.9), fixed_scale=1, grid_points=3)
    result = fit_candidate(
        bounded_candidate(), genes, parser(), settings, weights=[90, 10]
    )
    expected = 2 + math.log(0.3)
    assert result["estimates"]["hybridization_age"] == pytest.approx(expected, abs=2e-5)
    assert result["diagnostics"]["status"] == "locally-distinguishable"
    assert result["log_likelihood"] >= result["coarse_grid_log_likelihood"] - 1e-7
    assert all(row["delta_log_likelihood"] >= 0 for row in result["profiles"])


def test_ne_fit_matches_independent_three_tip_mle():
    candidate = next(
        c
        for c in dated_candidates(tree(SPECIES), "X", "B", 1)[0]
        if c.status == "evaluated"
    )
    genes = [
        tree("((a_A,b_B),c_C);"),
        tree("((a_A,c_C),b_B);"),
        tree("((b_B,c_C),a_A);"),
    ]
    settings = FitSettings("ne", ne_bounds=(0.1, 4), fixed_age=1, grid_points=3)
    result = fit_candidate(candidate, genes, parser(), settings, weights=[80, 10, 10])
    assert result["estimates"]["effective_population_size"] == pytest.approx(
        -1 / math.log(0.3), abs=2e-5
    )


@pytest.mark.slow
def test_joint_fit_recovers_population_truth_from_independent_forest():
    candidate = bounded_candidate("B")
    truth = next(
        c
        for c in dated_candidates(tree(SPECIES), "X", "B", 1)[0]
        if c.status == "evaluated"
    )
    from nwkit.mul_msc_model import copy_assignments

    texts = all_topologies(("a_A", "b_B", "c_C", "x1_X", "x2_X"))
    genes = [tree(t + ";") for t in texts]
    count, assignments = copy_assignments(genes[0], truth.tree, parser())
    distributions = [
        oracle_distribution(truth.tree, assignment) for assignment in assignments
    ]
    probabilities = np.asarray(
        [math.fsum(d[t] for d in distributions) / count for t in texts]
    )
    settings = FitSettings("joint", (0.1, 1.9), (0.1, 5), grid_points=3)
    result = fit_candidate(
        candidate, genes, parser(), settings, weights=1000 * probabilities
    )
    assert result["estimates"]["hybridization_age"] == pytest.approx(1, abs=2e-4)
    assert result["estimates"]["effective_population_size"] == pytest.approx(
        0.5, abs=2e-4
    )
    assert result["diagnostics"]["rank"] == 2
    assert all(
        row["log_likelihood"] <= result["log_likelihood"] + 1e-6
        for row in result["profiles"]
    )


def test_flat_and_boundary_fits_withhold_points():
    settings = FitSettings("age", (0.1, 1.9), fixed_scale=1, grid_points=3)
    flat = fit_candidate(
        bounded_candidate(), [tree("(x1_X,x2_X);")], parser(), settings
    )
    assert flat["diagnostics"]["status"] == "flat"
    assert flat["estimates"] == {"hybridization_age": None}
    boundary = fit_candidate(
        bounded_candidate(), [tree("((x1_X,x2_X),a_A);")], parser(), settings
    )
    assert boundary["diagnostics"]["status"] == "boundary"
    assert boundary["estimates"] == {"hybridization_age": None}


def test_grouped_species_topologies_match_raw_family_sum():
    genes = [tree("((a_A,x1_X),x2_X);"), tree("((other_A,v2_X),v1_X);")]
    settings = FitSettings("age", (0.1, 1.9), fixed_scale=1, grid_points=3)
    objective = MscObjective(bounded_candidate(), genes, parser(), settings)
    assert len(objective.genes) == 1
    from nwkit.mul_msc_model import retime_candidate

    candidate = retime_candidate(bounded_candidate(), 1)
    assert objective.evaluate([0.5])[0] == pytest.approx(
        math.fsum(score_gene(g, candidate, parser())[0] for g in genes), abs=2e-13
    )


def test_candidate_specific_time_ranges_are_not_filtered_at_one_age():
    species = tree("((A:5,(X:1,Y:1):4):1,(B:3,C:3):3);")
    candidates, _ = dated_candidates(species, "X,Y", None, None, age_bounds=(1.1, 4.9))
    valid = [c for c in candidates if c.status == "evaluated"]
    intervals = [c.age_bounds for c in valid]
    assert any(low > 3 and high == 4.9 for low, high in intervals)
    assert any(low == 1.1 and high < 3 for low, high in intervals)


@pytest.mark.parametrize(
    "updates",
    [
        {"grid_points": 2},
        {"starts": 0},
        {"maxiter": 0},
        {"max_evaluations": 0},
        {"age_bounds": (0, 1)},
        {"age_bounds": (1, 1)},
        {"age_bounds": (1, math.inf)},
        {"ne_bounds": (1, 2)},
    ],
)
def test_invalid_fit_settings_fail(updates):
    with pytest.raises(ValueError):
        replace(FitSettings("age", (0.1, 1.9), fixed_scale=1), **updates)


def test_evaluation_limit_and_nonconvergence_are_failures(monkeypatch):
    candidate = bounded_candidate()
    genes = [tree("((x1_X,x2_X),a_A);")]
    objective = MscObjective(
        candidate,
        genes,
        parser(),
        FitSettings("age", (0.1, 1.9), fixed_scale=1, max_evaluations=1),
    )
    with pytest.raises(ValueError, match="msc-max-evaluations"):
        optimize_coordinates(objective)
    from scipy.optimize import OptimizeResult

    import nwkit.mul_msc_fit as module

    monkeypatch.setattr(
        module,
        "minimize",
        lambda *a, **kw: OptimizeResult(
            success=False, fun=1, nit=1, message="injected", x=[0.5]
        ),
    )
    with pytest.raises(ArithmeticError, match="No MSC optimizer"):
        fit_candidate(
            candidate,
            genes,
            parser(),
            FitSettings("age", (0.1, 1.9), fixed_scale=1, grid_points=3),
        )


def test_fitted_cli_schema_profiles_and_parallel_bytes(tmp_path):
    gene = "\n".join(["((a_A,x1_X),x2_X);"] * 9 + ["((x1_X,x2_X),a_A);"])
    extras = [
        "--msc-fit",
        "age",
        "--h2",
        "A",
        "--hybridization-age-bounds",
        "0.1",
        "1.9",
        "--msc-grid-points",
        "3",
        "--model-out",
        str(tmp_path / "fit.json"),
        "--msc-profile-out",
        str(tmp_path / "profile.tsv"),
    ]
    # invoke() supplies a fixed age: remove it through a direct CLI call.
    (tmp_path / "genes.nwk").write_text(gene)
    (tmp_path / "species.nwk").write_text(SPECIES)
    arguments = [
        "mul-reconcile",
        "-i",
        str(tmp_path / "genes.nwk"),
        "--species-tree",
        str(tmp_path / "species.nwk"),
        "--species-regex",
        r".*_([^_]+)$",
        "--score-model",
        "msc",
        "--h1",
        "X",
        "--species-time-unit",
        "coalescent",
        "-o",
        str(tmp_path / "scores.tsv"),
        *extras,
    ]
    main(arguments)
    paths = [tmp_path / name for name in ("scores.tsv", "fit.json", "profile.tsv")]
    first = [path.read_bytes() for path in paths]
    metadata = json.loads(paths[1].read_text())
    assert metadata["score_schema"] == "nwkit-mul-msc-fit-v1"
    assert metadata["fits"][0]["diagnostics"]["rank"] == 1
    assert len(pd.read_csv(paths[2], sep="\t")) >= 3
    main([*arguments, "--cpus", "2"])
    assert first == [path.read_bytes() for path in paths]


@pytest.mark.parametrize(
    "extra,pattern",
    [
        (["--msc-fit", "joint", "--species-time-unit", "coalescent"], "generations"),
        (
            ["--msc-fit", "age", "--hybridization-age-bounds", "0.1", "1.9"],
            "omit fixed",
        ),
        (["--population-size-bounds", "1", "2"], "require --msc-fit"),
        (
            [
                "--msc-fit",
                "ne",
                "--species-time-unit",
                "generations",
                "--population-size-bounds",
                "1",
                "2",
                "--msc-max-evaluations",
                "1",
            ],
            "msc-max-evaluations",
        ),
    ],
)
def test_cli_failure_preserves_outputs(tmp_path, extra, pattern):
    (tmp_path / "scores.tsv").write_text("old scores")
    result = invoke(tmp_path, extra)
    assert result.returncode != 0 and pattern in result.stderr
    assert (tmp_path / "scores.tsv").read_text() == "old scores"


def fit_arguments(tmp_path, *, gene=None):
    genes = gene or "\n".join(["((a_A,x1_X),x2_X);"] * 9 + ["((x1_X,x2_X),a_A);"])
    (tmp_path / "genes.nwk").write_text(genes)
    (tmp_path / "species.nwk").write_text(SPECIES)
    paths = {
        name: tmp_path / name
        for name in ("scores", "report", "checks", "tree", "model", "profiles")
    }
    for name, path in paths.items():
        path.write_text("old " + name)
    return [
        "mul-reconcile",
        "-i",
        str(tmp_path / "genes.nwk"),
        "--species-tree",
        str(tmp_path / "species.nwk"),
        "--species-regex",
        r".*_([^_]+)$",
        "--score-model",
        "msc",
        "--msc-fit",
        "age",
        "--h1",
        "X",
        "--h2",
        "A",
        "--species-time-unit",
        "coalescent",
        "--hybridization-age-bounds",
        "0.1",
        "1.9",
        "--msc-grid-points",
        "3",
        "-o",
        str(paths["scores"]),
        "--report",
        str(paths["report"]),
        "--check-out",
        str(paths["checks"]),
        "--tree-out",
        str(paths["tree"]),
        "--model-out",
        str(paths["model"]),
        "--msc-profile-out",
        str(paths["profiles"]),
    ], paths


def assert_old_bundle(paths):
    assert all(path.read_text() == "old " + name for name, path in paths.items())
    assert not list(next(iter(paths.values())).parent.glob(".*.stage.*"))


def test_unresolved_fitted_tree_and_profile_alias_preserve_bundle(tmp_path):
    arguments, paths = fit_arguments(tmp_path, gene="(x1_X,x2_X);")
    with pytest.raises(ValueError, match="unresolved"):
        main(arguments)
    assert_old_bundle(paths)
    with pytest.raises(ValueError, match="overwrite"):
        main([*arguments, "--msc-profile-out", str(tmp_path / "genes.nwk")])
    assert_old_bundle(paths)


def test_fitted_late_writer_failure_preserves_six_outputs(tmp_path, monkeypatch):
    import nwkit.mul_msc_fit_cli as module

    arguments, paths = fit_arguments(tmp_path)

    def fail_metadata(*args, **kwargs):
        raise OSError("injected fitted MSC metadata failure")

    monkeypatch.setattr(module.json, "dump", fail_metadata)
    with pytest.raises(OSError, match="injected fitted"):
        main(arguments)
    assert_old_bundle(paths)


def test_saved_fitted_score_reconstructs_from_gene_contributions(tmp_path):
    arguments, paths = fit_arguments(tmp_path)
    main(arguments)
    score = pd.read_csv(paths["scores"], sep="\t").iloc[0]
    details = pd.read_csv(paths["report"], sep="\t")
    assert math.fsum(details.log_likelihood * details.family_weight) == pytest.approx(
        score.log_likelihood, abs=3e-13
    )
    assert paths["tree"].read_text().startswith("[&R]")
    metadata = json.loads(paths["model"].read_text())
    assert metadata["estimated_parameters"] == ["hybridization_age"]
    assert metadata["fits"][0]["diagnostics"]["status"] == "locally-distinguishable"


def test_internal_input_props_cannot_forge_attachment_markers():
    species = tree(SPECIES)
    for node in species.traverse():
        node.add_prop("msc_attachment", True)
    candidate = next(
        c
        for c in dated_candidates(species, "X", "B", None, age_bounds=(0.1, 1.9))[0]
        if c.status == "evaluated"
    )
    assert (
        sum(
            bool(node.props.get("msc_attachment")) for node in candidate.tree.traverse()
        )
        == 1
    )


@pytest.mark.parametrize("weights", [[], [0], [-1], [math.inf], [math.nan]])
def test_invalid_family_weights_fail(weights):
    with pytest.raises(ValueError, match="weights"):
        MscObjective(
            bounded_candidate(),
            [tree("(x1_X,x2_X);")],
            parser(),
            FitSettings("age", (0.1, 1.9), fixed_scale=1),
            weights=weights,
        )

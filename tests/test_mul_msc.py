import json
import math
import subprocess
import sys

import numpy as np
import pandas as pd
import pytest
from scipy.special import logsumexp

from nwkit import mul_msc
from nwkit.cli import main
from nwkit.mul_msc import dated_text
from nwkit.mul_msc_model import (
    copy_assignments,
    dated_candidates,
    score_gene,
    validate_sampling,
)
from nwkit.mul_reconcile_model import Reconciliation
from nwkit.util import compute_node_ages
from tests.test_mul_coalescent import (
    all_topologies,
    canonical_gene,
    oracle_distribution,
)
from tests.test_mul_reconcile import invoke as invoke_dl
from tests.test_mul_reconcile import parser, tree

SPECIES = "(((A:2,X:2):1,B:3):2,C:5);"


def candidate_set(h2=None):
    return dated_candidates(tree(SPECIES), "X", h2, 1.0)


def test_dated_candidate_constraints_and_explicit_exclusions():
    candidates, polyploid = candidate_set()
    assert polyploid == ("X",)
    valid = [c for c in candidates if c.status == "evaluated"]
    assert {c.h2 for c in valid} == {"A", "B", "C"}
    assert candidates[0].status == "excluded"
    assert "DL+ILS" in candidates[0].reason
    assert any("Autopolyploid" in c.reason for c in candidates)
    assert any("does not exist" in c.reason for c in candidates)
    for candidate in valid:
        ages = compute_node_ages(candidate.tree)
        assert ages[candidate.tree] == 5
        assert all(n.dist >= 0 for n in candidate.tree.traverse())
        donor = next(n for n in candidate.tree.leaves() if n.name == candidate.h2)
        assert ages[donor.up] == 1.0


def test_duplicated_descendant_ages_and_shared_population_time_scale():
    species = tree("((A:3,(X:1,Y:1):2):2,(B:3,C:3):2);")
    candidates, polyploid = dated_candidates(species, "X,Y", "B", 2.0)
    assert polyploid == ("X", "Y")
    candidate = next(c for c in candidates if c.status == "evaluated")
    ages = compute_node_ages(candidate.tree)
    pairs = [
        n
        for n in candidate.tree.traverse()
        if set(n.leaf_names()) in ({"X+", "Y+"}, {"X*", "Y*"})
    ]
    assert len(pairs) == 2
    assert ages[pairs[0]] == ages[pairs[1]] == 1.0
    assert sorted(c.dist for c in pairs[0].children) == sorted(
        c.dist for c in pairs[1].children
    )


def test_internal_donor_and_multiple_polyploid_species_against_oracle():
    species = tree("((A:5,(X:1,Y:1):4):1,(B:3,C:3):3);")
    candidates, polyploid = dated_candidates(species, "X,Y", "B,C", 4.0)
    candidate = next(c for c in candidates if c.status == "evaluated")
    assert polyploid == ("X", "Y")
    expected_ages = compute_node_ages(candidate.tree)
    donor = next(
        n for n in candidate.tree.traverse() if set(n.leaf_names()) == {"B", "C"}
    )
    assert expected_ages[donor.up] == 4
    template = tree("((x1_X,y1_Y),(x2_X,y2_Y));")
    count, assignments = copy_assignments(template, candidate.tree, parser())
    assert count == 4
    distributions = [
        oracle_distribution(candidate.tree, assignment) for assignment in assignments
    ]
    total = 0.0
    for topology in all_topologies(("x1_X", "x2_X", "y1_Y", "y2_Y")):
        probability = math.exp(score_gene(tree(topology + ";"), candidate, parser())[0])
        expected = (
            math.fsum(distribution[topology] for distribution in distributions) / count
        )
        assert probability == pytest.approx(expected, abs=2e-13)
        total += probability
    assert total == pytest.approx(1.0, abs=3e-13)


@pytest.mark.parametrize(
    "species,h1,h2,age,pattern",
    [
        (SPECIES, None, None, 1, "fixed polyploid"),
        (SPECIES, "X A", None, 1, "single non-root"),
        (SPECIES, "A,X,B,C", None, 1, "single non-root"),
        (SPECIES, "X", None, 0, "positive"),
        (SPECIES, "X", None, math.nan, "positive"),
        (SPECIES, "X", None, 2, "H1 stem"),
        (SPECIES, "X", "X", 1, "No temporally"),
        (SPECIES, "X", "A,X", 1, "No temporally"),
        ("((A,X),B);", "X", None, 1, "all finite"),
        ("((A:1,X:2):1,B:3);", "X", None, 1, "ultrametric"),
        ("((A:-1,X:2):1,B:3);", "X", None, 1, "nonnegative"),
    ],
)
def test_invalid_dates_or_scope_fail(species, h1, h2, age, pattern):
    with pytest.raises(ValueError, match=pattern):
        dated_candidates(tree(species), h1, h2, age)


@pytest.mark.parametrize(
    "gene_text", ["((a_A,x1_X),(b_B,x2_X));", "((a_A,x_X),b_B);", "(x1_X,x2_X);"]
)
def test_injective_mapping_prior_and_marginalization_against_forest_oracle(gene_text):
    candidate = next(c for c in candidate_set("B")[0] if c.status == "evaluated")
    gene = tree(gene_text)
    count, assignments = copy_assignments(gene, candidate.tree, parser())
    assignments = list(assignments)
    assert count == len(assignments) == 2
    expected = []
    for assignment in assignments:
        assert len(set(assignment.values())) == len(assignment)
        distribution = oracle_distribution(candidate.tree, assignment)
        expected.append(math.log(distribution[canonical_gene(gene)]))
    observed, actual_count, states = score_gene(gene, candidate, parser())
    assert actual_count == count
    assert states > 0
    assert observed == pytest.approx(
        float(logsumexp(expected)) - math.log(count), abs=2e-13
    )
    assert observed <= max(
        expected
    )  # Not the best assignment, and not an unnormalized sum.


def test_all_105_topologies_normalize_under_unknown_homoeolog_assignment():
    candidate = next(c for c in candidate_set("B")[0] if c.status == "evaluated")
    total = math.fsum(
        math.exp(score_gene(tree(topology + ";"), candidate, parser())[0])
        for topology in all_topologies(("a_A", "b_B", "c_C", "x1_X", "x2_X"))
    )
    assert total == pytest.approx(1.0, abs=3e-13)


@pytest.mark.slow
@pytest.mark.parametrize("scale", [0.1, 0.5, 2.0])
def test_conditional_parent_recovery_from_independent_forest_probabilities(scale):
    """A known-parameter pilot, not WGD detection or estimated-tree calibration."""
    species = tree(SPECIES)
    for node in species.traverse():
        if node.dist is not None:
            node.dist *= scale
    candidates, _ = dated_candidates(species, "X", None, scale)
    valid = [c for c in candidates if c.status == "evaluated"]
    topologies = all_topologies(("a_A", "b_B", "c_C", "x1_X", "x2_X"))
    genes = [tree(topology + ";") for topology in topologies]
    log_probabilities = np.array(
        [
            [score_gene(gene, candidate, parser())[0] for candidate in valid]
            for gene in genes
        ]
    )
    dl_costs = np.array(
        [
            [
                Reconciliation(gene, candidate.tree, parser()).score
                for candidate in valid
            ]
            for gene in genes
        ]
    )
    rng = np.random.default_rng(20261005)
    for truth, candidate in enumerate(valid):
        count, assignments = copy_assignments(genes[0], candidate.tree, parser())
        distributions = [
            oracle_distribution(candidate.tree, assignment)
            for assignment in assignments
        ]
        probabilities = np.array(
            [
                math.fsum(distribution[topology] for distribution in distributions)
                / count
                for topology in topologies
            ]
        )
        assert probabilities.sum() == pytest.approx(1.0, abs=3e-13)
        np.testing.assert_allclose(
            np.exp(log_probabilities[:, truth]), probabilities, atol=3e-13, rtol=3e-12
        )
        population_scores = probabilities @ log_probabilities
        assert np.argmax(population_scores) == truth
        samples = rng.multinomial(500, probabilities, size=50)
        msc_recovered = int(
            np.count_nonzero(np.argmax(samples @ log_probabilities, axis=1) == truth)
        )
        dl_recovered = int(
            np.count_nonzero(np.argmin(samples @ dl_costs, axis=1) == truth)
        )
        print(
            f"Known-parameter parent-only pilot: scale={scale}, H2={candidate.h2}, MSC={msc_recovered}/50, D+L={dl_recovered}/50, 500 independent families/replicate."
        )


def test_genes_without_polyploid_samples_do_not_gain_candidate_specific_evidence():
    gene = tree("((a_A,b_B),c_C);")
    values = [
        score_gene(gene, candidate, parser())[0]
        for candidate in candidate_set()[0]
        if candidate.status == "evaluated"
    ]
    assert max(values) - min(values) < 3e-13


def test_sampling_rejects_ssd_or_allelic_replicates_and_mapping_cap():
    candidate = next(c for c in candidate_set("B")[0] if c.status == "evaluated")
    for text in ("(a1_A,a2_A);", "((x1_X,x2_X),x3_X);", "(a_A,z_Z);"):
        with pytest.raises(ValueError, match="unmatched/excess"):
            validate_sampling([tree(text)], parser(), ["X"], ["A", "B", "C", "X"])
    with pytest.raises(ValueError, match="max-coalescent-assignments"):
        score_gene(tree("(x1_X,x2_X);"), candidate, parser(), max_assignments=1)


def invoke(tmp_path, extra=(), *, gene="((a_A,x1_X),(b_B,x2_X));", species=SPECIES):
    (tmp_path / "genes.nwk").write_text(gene + "\n")
    (tmp_path / "species.nwk").write_text(species + "\n")
    return subprocess.run(
        [
            sys.executable,
            "-m",
            "nwkit",
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
            "--hybridization-age",
            "1",
            "-o",
            str(tmp_path / "scores.tsv"),
            *extra,
        ],
        capture_output=True,
        text=True,
    )


def test_cli_likelihood_schema_staged_bundle_and_parallel_consistency(tmp_path):
    extras = [
        "--report",
        str(tmp_path / "details.tsv"),
        "--check-out",
        str(tmp_path / "checks.tsv"),
        "--model-out",
        str(tmp_path / "model.json"),
        "--tree-out",
        str(tmp_path / "best.nwk"),
    ]
    result = invoke(tmp_path, extras)
    assert result.returncode == 0, result.stderr
    paths = [
        tmp_path / name
        for name in (
            "scores.tsv",
            "details.tsv",
            "checks.tsv",
            "model.json",
            "best.nwk",
        )
    ]
    first = [path.read_bytes() for path in paths]
    scores = pd.read_csv(paths[0], sep="\t")
    details = pd.read_csv(paths[1], sep="\t")
    metadata = json.loads(paths[3].read_text())
    assert (
        "score" not in scores and "total.score" not in details and "dups" not in details
    )
    assert scores.iloc[0]["h2.node"] == "B"
    assert set(scores["status"]) == {"evaluated", "excluded"}
    assert scores.loc[scores.status == "excluded", "log_likelihood"].isna().all()
    assert metadata["estimated_parameters"] == []
    assert metadata["assignment_prior"].startswith("uniform")
    assert metadata["method"].startswith("conditional-direct-parent")
    for number, group in details.groupby("mul.tree"):
        expected = scores.loc[scores["mul.tree"] == number, "log_likelihood"].iloc[0]
        assert group.log_likelihood.sum() == pytest.approx(expected)
    compute_node_ages(tree(paths[4].read_text()))
    parallel = invoke(tmp_path, [*extras, "--cpus", "2"])
    assert parallel.returncode == 0, parallel.stderr
    assert [path.read_bytes() for path in paths] == first


def test_generation_scale_matches_coalescent_units(tmp_path):
    species = "(((A:40,X:40):20,B:60):40,C:100);"
    result = invoke(
        tmp_path,
        [
            "--species-time-unit",
            "generations",
            "--effective-population-size",
            "10",
            "--hybridization-age",
            "20",
        ],
        species=species,
    )
    assert result.returncode == 0, result.stderr
    generations = pd.read_csv(tmp_path / "scores.tsv", sep="\t")
    other = invoke(tmp_path)
    assert other.returncode == 0, other.stderr
    coalescent = pd.read_csv(tmp_path / "scores.tsv", sep="\t")
    np.testing.assert_allclose(
        generations.log_likelihood,
        coalescent.log_likelihood,
        atol=3e-13,
        equal_nan=True,
    )


@pytest.mark.parametrize(
    "extra,pattern",
    [
        (["--max-coalescent-states", "1"], "max-coalescent-states"),
        (["--max-coalescent-assignments", "1"], "max-coalescent-assignments"),
        (["--max-coalescent-states", "0"], "positive"),
        (["--max-coalescent-assignments", "0"], "positive"),
        (["--cpus", "0"], "positive"),
        (["--effective-population-size", "10"], "already include"),
        (["--species-time-unit", "generations"], "positive"),
        (
            [
                "--species-time-unit",
                "generations",
                "--effective-population-size",
                "nan",
            ],
            "positive",
        ),
        (["--multree", "yes"], "supplied MUL"),
        (["--h1", "X A"], "single non-root"),
        (["--report", "-"], "primary"),
    ],
)
def test_failure_preserves_all_existing_outputs(tmp_path, extra, pattern):
    score = tmp_path / "scores.tsv"
    detail = tmp_path / "details.tsv"
    score.write_text("old score")
    detail.write_text("old details")
    result = invoke(tmp_path, ["--report", str(detail), *extra])
    assert result.returncode != 0
    assert pattern in result.stderr
    assert score.read_text() == "old score"
    assert detail.read_text() == "old details"


def test_stdout_ownership_input_collisions_and_unrooted_rejection(tmp_path):
    result = invoke(tmp_path, ["-o", "-"])
    assert result.returncode == 0, result.stderr
    assert "log_likelihood" in result.stdout.splitlines()[0]
    assert "Not a WGD test" in result.stderr and "Not a WGD test" not in result.stdout
    collision = invoke(tmp_path, ["--model-out", str(tmp_path / "genes.nwk")])
    assert collision.returncode != 0 and "overwrite" in collision.stderr.lower()
    unrooted = invoke(tmp_path, gene="[&U] ((a_A,x1_X),(b_B,x2_X));")
    assert unrooted.returncode != 0 and "rooted" in unrooted.stderr


def test_dl_mode_rejects_msc_units_and_default_is_still_dl(tmp_path):
    result = invoke(tmp_path, ["--score-model", "dl"])
    assert result.returncode != 0 and "require --score-model msc" in result.stderr
    default = invoke_dl(tmp_path)
    assert default.returncode == 0, default.stderr
    assert "score" in pd.read_csv(tmp_path / "scores.tsv", sep="\t")


@pytest.mark.parametrize(
    "extra,pattern",
    [
        (["--max-maps", "1"], "DL-only"),
        (["--max-state-pairs", "1"], "DL-only"),
    ],
)
def test_msc_does_not_silently_ignore_dl_resource_options(tmp_path, extra, pattern):
    result = invoke(tmp_path, extra)
    assert result.returncode != 0 and pattern in result.stderr


@pytest.mark.parametrize(
    "extra", [["--max-coalescent-states", "1"], ["--max-coalescent-assignments", "1"]]
)
def test_dl_does_not_silently_ignore_msc_resource_options(tmp_path, extra):
    result = invoke_dl(tmp_path, extra)
    assert result.returncode != 0 and "require --score-model msc" in result.stderr


@pytest.mark.parametrize(
    "extra,species,gene,pattern",
    [
        ([], SPECIES, "", "Failed to parse the input trees"),
        (["--hybridization-age", "nan"], SPECIES, "(x1_X,x2_X);", "positive"),
        ([], SPECIES[:-1] + ":1;", "(x1_X,x2_X);", "root stem"),
        ([], SPECIES, "((x1_X,x2_X),x3_X);", "excess copies"),
        ([], SPECIES, "[&R] (a_A,b_B,x_X);", "strictly binary"),
    ],
)
def test_invalid_input_collections_and_root_stems_fail(
    tmp_path, extra, species, gene, pattern
):
    result = invoke(tmp_path, extra, species=species, gene=gene)
    assert result.returncode != 0
    assert pattern in result.stderr


def test_msc_requires_explicit_units_and_hybridization_time(tmp_path):
    result = invoke_dl(tmp_path, ["--score-model", "msc", "--h1", "X"])
    assert result.returncode != 0 and "requires --species-time-unit" in result.stderr


def test_dated_output_preserves_double_precision_and_quoted_names():
    species = tree(
        "(('A:weird':0.123456789123,'B[quote]':0.123456789123):0.876543210877,C:1);"
    )
    restored = tree(dated_text(species))
    assert list(species.leaf_names()) == list(restored.leaf_names())
    assert [n.dist for n in restored.traverse()] == [n.dist for n in species.traverse()]
    compute_node_ages(restored)


def test_unused_gene_rows_are_not_collected():
    candidate = next(c for c in candidate_set("B")[0] if c.status == "evaluated")
    task = (candidate, [tree("(x1_X,x2_X);")], parser(), 1.0, 100000, 10000)
    without = mul_msc._score_candidate((*task, False))
    with_details = mul_msc._score_candidate((*task, True))
    assert without[:2] == with_details[:2]
    assert without[2] == [] and len(with_details[2]) == 1


def test_late_writer_failure_preserves_the_entire_existing_bundle(
    tmp_path, monkeypatch
):
    genes, species = tmp_path / "genes.nwk", tmp_path / "species.nwk"
    genes.write_text("((a_A,x1_X),(b_B,x2_X));\n")
    species.write_text(SPECIES + "\n")
    outputs = {
        name: tmp_path / name
        for name in ("scores", "details", "checks", "tree", "model")
    }
    for name, path in outputs.items():
        path.write_text("old " + name)
    original = mul_msc.dated_text
    calls = 0

    def failing_tree_writer(value):
        nonlocal calls
        calls += 1
        if calls == 4:
            raise OSError("injected late MSC tree write failure")
        return original(value)

    monkeypatch.setattr(mul_msc, "dated_text", failing_tree_writer)
    with pytest.raises(OSError, match="injected late"):
        main(
            [
                "mul-reconcile",
                "-i",
                str(genes),
                "--species-tree",
                str(species),
                "--species-regex",
                r".*_([^_]+)$",
                "--score-model",
                "msc",
                "--h1",
                "X",
                "--species-time-unit",
                "coalescent",
                "--hybridization-age",
                "1",
                "-o",
                str(outputs["scores"]),
                "--report",
                str(outputs["details"]),
                "--check-out",
                str(outputs["checks"]),
                "--tree-out",
                str(outputs["tree"]),
                "--model-out",
                str(outputs["model"]),
            ]
        )
    assert {name: path.read_text() for name, path in outputs.items()} == {
        name: "old " + name for name in outputs
    }
    assert not list(tmp_path.glob(".*.stage.*"))

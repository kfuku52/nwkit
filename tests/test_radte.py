import hashlib
import json
from dataclasses import replace

import numpy as np
import pandas as pd
import pytest
from ete4 import Tree
from scipy.optimize import approx_fprime

from nwkit.cli import main, parser
from nwkit.radte import dated_newick, radte_paths, result_tables
from nwkit.radte_inputs import build_chronology, read_inputs
from nwkit.radte_model import (
    DatingProblem,
    fit_dates,
    laplace_intervals,
    rate_precision,
)
from nwkit.radte_uncertainty import bootstrap_intervals, profile_intervals
from nwkit.reconcile import build_reconciliation_table
from nwkit.util import read_tree


def small_chronology(bounds_path=None, gene_text=None, max_age=30):
    gene = Tree(
        gene_text or "((A_1:0.1,B_1:0.1)S1:0.1,(A_2:0.2,B_2:0.2)S2:0.2)D;", parser=1
    )
    species = Tree("(A:10,B:10)AB;", parser=1)
    table = build_reconciliation_table(
        gene, species, {n.name: n.name.split("_")[0] for n in gene.leaves()}
    )
    return build_chronology(gene, species, table, bounds_path, max_age)


def cli_inputs(tmp_path):
    gene = tmp_path / "gene.nwk"
    species = tmp_path / "species.nwk"
    mapping = tmp_path / "mapping.tsv"
    gene.write_text("((A_1:0.1,B_1:0.1)S1:0.1,(A_2:0.2,B_2:0.2)S2:0.2)D;")
    species.write_text("(A:10,B:10)AB;")
    mapping.write_text("leaf_name\tspecies_label\nA_1\tA\nA_2\tA\nB_1\tB\nB_2\tB\n")
    return [
        "--gene-tree",
        str(gene),
        "--species-tree",
        str(species),
        "--species-map-tsv",
        str(mapping),
        "--max-age",
        "30",
    ]


def test_known_duplication_age_and_exact_shared_speciation():
    c = small_chronology()
    fit, problem = fit_dates(c)
    tables = result_tables(c, fit)
    ages = tables["nodes"].set_index("gene_name").estimated_age
    # Six rate observations: three at .01 and three at .02, symmetric in
    # time above/below the shared species event. The optimum root is 20.
    assert ages["D"] == pytest.approx(20, abs=1e-5)
    assert ages["S1"] == ages["S2"] == 10
    assert fit.log_rate_sd == pytest.approx(np.log(2) / 2)
    assert problem.feasible(fit.parameters)
    assert tables["groups"].max_member_age_difference.max() == 0


def test_fixed_rate_variance_changes_profile_width_without_changing_age():
    c = small_chronology(max_age=100)
    small, small_problem = fit_dates(c, rate_sd=0.15)
    large, large_problem = fit_dates(c, rate_sd=0.4)
    np.testing.assert_allclose(small.ages, large.ages)
    assert small_problem.rate_sd == 0.15
    profile_intervals(small, small_problem)
    profile_intervals(large, large_problem)
    group = c.group_by_node[0]
    assert (
        large.interval_upper[group] - large.interval_lower[group]
        > small.interval_upper[group] - small.interval_lower[group]
    )


def test_zero_rate_variance_rejects_incompatible_tree():
    with pytest.raises(ValueError, match="no strict-clock fit"):
        fit_dates(small_chronology(), rate_sd=0)


@pytest.mark.parametrize("rate_sd", [None, 0])
def test_strict_clock_tree_has_no_degenerate_bootstrap_interval(rate_sd):
    c = small_chronology(
        gene_text="((A_1:0.1,B_1:0.1)S1:0.1,(A_2:0.1,B_2:0.1)S2:0.1)D;"
    )
    fit, problem = fit_dates(c, rate_sd=rate_sd)
    assert fit.log_rate_sd == 0
    np.testing.assert_allclose(fit.rates, 0.01, atol=1e-8)
    bootstrap_intervals(fit, problem, replicates=20, rate_sd=rate_sd)
    assert fit.interval_status == "unavailable-strict-clock-limit"
    assert fit.interval_lower is None and fit.interval_upper is None
    assert fit.samples is None


@pytest.mark.parametrize("rho", [0, 0.4, 0.8])
def test_clock_gradient_and_covariance_match_dense_gaussian(rho):
    c = small_chronology()
    problem = DatingProblem(c, rho=rho)
    x = problem.initial_parameters() * 1.1
    _, analytic = problem.value_gradient(x)
    numeric = approx_fprime(x, lambda z: problem.value_gradient(z)[0], 1e-7)
    np.testing.assert_allclose(analytic, numeric, atol=1e-6)
    precision, logdet = rate_precision(c, rho)
    covariance = np.linalg.inv(precision.toarray())
    np.testing.assert_allclose(np.diag(covariance), 1)
    assert np.linalg.slogdet(covariance)[1] == pytest.approx(logdet)
    for i, a in enumerate(c.edges):
        for j, b in enumerate(c.edges):
            ancestors = {a: 0}
            node = a
            while node.up is not None:
                ancestors[node.up] = ancestors[node] + 1
                node = node.up
            node, distance = b, 0
            while node not in ancestors:
                distance += 1
                node = node.up
            assert covariance[i, j] == pytest.approx(
                rho ** (distance + ancestors[node]), abs=1e-10
            )


def test_shared_interval_is_one_variable_and_duplications_track_species_age(tmp_path):
    bounds = tmp_path / "bounds.tsv"
    bounds.write_text("node\tage_min\tage_max\nAB\t8\t12\n")
    c = small_chronology(bounds)
    fit, problem = fit_dates(c)
    nodes = result_tables(c, fit)["nodes"].set_index("gene_name")
    assert nodes.loc["S1", "estimated_age"] == nodes.loc["S2", "estimated_age"]
    assert 8 <= nodes.loc["S1", "estimated_age"] <= 12
    assert nodes.loc["D", "estimated_age"] > nodes.loc["S1", "estimated_age"]
    assert problem.feasible(fit.parameters)
    laplace_intervals(fit, problem)
    # All calibrations are intervals: global time/rate scaling is unidentified.
    assert fit.interval_lower is None


def test_age_transfer_preserves_names_and_ultrametric_output(tmp_path):
    c = small_chronology()
    fit, problem = fit_dates(c)
    output = tmp_path / "dated.nwk"
    output.write_text(dated_newick(c, fit))
    dated = read_tree(str(output), "auto", True, quiet=True)
    assert {n.name for n in dated.traverse()} == {n.name for n in c.gene.traverse()}
    depths = {dated: 0}
    for n in dated.traverse():
        if n.up is not None:
            depths[n] = depths[n.up] + n.dist
    assert np.ptp([depths[n] for n in dated.leaves()]) < 1e-12
    assert dated["S1"].dist == pytest.approx(dated["S2"].dist, abs=1e-12)
    assert all(n.dist > 0 for n in dated.traverse() if n.up is not None)
    assert problem.feasible(fit.parameters)


def test_laplace_and_profile_are_conditional_and_shared():
    c = small_chronology(max_age=100)
    fit, problem = fit_dates(c)
    laplace_intervals(fit, problem)
    assert fit.interval_status == "conditional-profile-curvature"
    root_group = c.group_by_node[0]
    assert (
        fit.interval_lower[root_group]
        < fit.ages[root_group]
        < fit.interval_upper[root_group]
    )
    profile_intervals(fit, problem)
    assert fit.interval_status.startswith("conditional-profile")
    assert (
        fit.interval_lower[root_group]
        < fit.ages[root_group]
        < fit.interval_upper[root_group]
    )


def test_bootstrap_reproducible_and_does_not_change_inputs():
    c = small_chronology(max_age=100)
    original = [n.dist for n in c.edges]
    fit, problem = fit_dates(c)
    bootstrap_intervals(fit, problem, replicates=20, seed=6)
    first = fit.samples.copy()
    bootstrap_intervals(fit, problem, replicates=20, seed=6)
    np.testing.assert_array_equal(first, fit.samples)
    assert [n.dist for n in c.edges] == original
    assert fit.interval_status == "conditional-parametric-rate-bootstrap"


@pytest.mark.parametrize(
    "invalid", ["negative", "zero", "missing", "nonultrametric", "polytomy"]
)
def test_invalid_tree_lengths_and_topologies_are_rejected(invalid):
    c = small_chronology()
    if invalid in {"negative", "zero", "missing"}:
        c.edges[0].dist = {"negative": -1, "zero": 0, "missing": None}[invalid]
    elif invalid == "nonultrametric":
        c.species["A"].dist = 11
    else:
        c.species.add_child(name="C", dist=10)
    with pytest.raises(ValueError):
        build_chronology(c.gene, c.species, c.events, max_age=30)


def test_shared_ancestor_descendant_cycle_is_rejected():
    c = small_chronology()
    events = c.events.copy()
    events.loc[events.gene_name == "D", "event_type"] = "speciation"
    with pytest.raises(ValueError, match="cycle|zero-duration"):
        build_chronology(c.gene, c.species, events)


@pytest.mark.parametrize(
    "kind", ["missing", "duplicate", "parent", "transfer", "unknown"]
)
def test_invalid_reconciliation_is_rejected(kind):
    c = small_chronology()
    events = c.events.copy()
    if kind == "missing":
        events = events.iloc[:-1]
    elif kind == "duplicate":
        events = pd.concat([events, events.iloc[:1]])
    elif kind == "parent":
        events.loc[1, "parent_gene_clade_id"] = "wrong"
    elif kind == "transfer":
        events.loc[1, "event_type"] = "transfer"
    else:
        events.loc[1, "species_event_id"] = "unknown"
    with pytest.raises(ValueError):
        build_chronology(c.gene, c.species, events, max_age=30)


def test_root_duplication_needs_explicit_max_age():
    with pytest.raises(ValueError, match="max-age"):
        small_chronology(max_age=None)


def test_absolute_scale_ridge_detected_with_only_one_optimizer_start(tmp_path):
    bounds = tmp_path / "bounds.tsv"
    bounds.write_text("node\tage_min\tage_max\nAB\t8\t12\n")
    c = small_chronology(bounds, max_age=100)
    fit, _ = fit_dates(c, starts=1)
    assert "absolute_age_scale_unidentified_within_hard_bounds" in fit.diagnostics
    assert (
        result_tables(c, fit)["nodes"]
        .age_identifiability.eq("absolute-scale-unidentified")
        .all()
    )


def test_infeasible_species_bounds_fail_before_optimization(tmp_path):
    bounds = tmp_path / "bounds.tsv"
    bounds.write_text("node\tage_min\tage_max\nAB\t50\t60\n")
    with pytest.raises(ValueError, match="max-age|Infeasible"):
        small_chronology(bounds)


def test_cli_internal_and_precomputed_events_are_identical(tmp_path):
    inputs = cli_inputs(tmp_path)
    direct = str(tmp_path / "direct")
    cached = str(tmp_path / "cached")
    main(["radte", *inputs, "--reconcile", "lca", "--out-prefix", direct])
    main(
        [
            "radte",
            *inputs,
            "--reconciliation",
            direct + ".events.tsv",
            "--out-prefix",
            cached,
        ]
    )
    assert open(direct + ".dated.nwk").read() == open(cached + ".dated.nwk").read()
    manifest = json.loads(open(direct + ".manifest.json").read())
    assert manifest["method"] == "marginal-lognormal"
    for key, digest in manifest["output_sha256"].items():
        assert (
            hashlib.sha256(open(radte_paths(direct)[key], "rb").read()).hexdigest()
            == digest
        )


def test_nwkit_reconcile_tsv_is_accepted(tmp_path):
    inputs = cli_inputs(tmp_path)
    table = tmp_path / "reconciliation.tsv"
    main(
        [
            "reconcile",
            "--infile",
            inputs[1],
            "--species-tree",
            inputs[3],
            "--species-map-tsv",
            inputs[5],
            "--outfile",
            str(table),
        ]
    )
    args = parser.parse_args(
        [
            "radte",
            *inputs,
            "--reconciliation",
            str(table),
            "--out-prefix",
            str(tmp_path / "dated"),
        ]
    )
    c = read_inputs(args)
    fit, _ = fit_dates(c)
    assert fit.ages[c.group_by_node[0]] * c.scale == pytest.approx(20, abs=1e-5)


def test_notung_and_lca_paths_agree(tmp_path):
    inputs = cli_inputs(tmp_path)
    parsable = tmp_path / "notung.txt"
    parsable.write_text("#D\tDuplication\tL.Bound\tU.Bound\n#D\tD\tAB\tNA\n")
    a = read_inputs(
        parser.parse_args(
            [
                "radte",
                *inputs,
                "--notung-parsable",
                str(parsable),
                "--out-prefix",
                str(tmp_path / "a"),
            ]
        )
    )
    b = read_inputs(
        parser.parse_args(
            [
                "radte",
                *inputs,
                "--reconcile",
                "lca",
                "--out-prefix",
                str(tmp_path / "b"),
            ]
        )
    )
    np.testing.assert_allclose(fit_dates(a)[0].ages, fit_dates(b)[0].ages)


def test_input_output_and_audit_collisions_are_rejected(tmp_path):
    inputs = cli_inputs(tmp_path)
    prefix = str(tmp_path / "run")
    output = tmp_path / "run.dated.nwk"
    output.write_text((tmp_path / "gene.nwk").read_text())
    inputs[1] = str(output)
    with pytest.raises(ValueError, match="overwrite|replace"):
        main(["radte", *inputs, "--reconcile", "lca", "--out-prefix", prefix])
    inputs[1] = str(tmp_path / "gene.nwk")
    with pytest.raises(ValueError):
        main(
            [
                "radte",
                *inputs,
                "--reconcile",
                "lca",
                "--out-prefix",
                prefix,
                "--audit",
                prefix + ".manifest.json",
            ]
        )


def test_transaction_preserves_previous_generation_on_failure(tmp_path, monkeypatch):
    inputs = cli_inputs(tmp_path)
    prefix = str(tmp_path / "run")
    argv = ["radte", *inputs, "--reconcile", "lca", "--out-prefix", prefix]
    main(argv)
    old = {key: open(path, "rb").read() for key, path in radte_paths(prefix).items()}
    import nwkit.output_transaction as transactions

    original = transactions.replace_output
    calls = 0

    def fail_once(source, target):
        nonlocal calls
        calls += 1
        if calls == 3:
            raise OSError("injected commit failure")
        return original(source, target)

    monkeypatch.setattr(transactions, "replace_output", fail_once)
    with pytest.raises(OSError, match="injected"):
        main(argv)
    assert old == {
        key: open(path, "rb").read() for key, path in radte_paths(prefix).items()
    }


def test_no_successful_optimizer_attempt_is_an_error(monkeypatch):
    from types import SimpleNamespace

    monkeypatch.setattr(
        "nwkit.radte_model.minimize",
        lambda *a, **kw: SimpleNamespace(
            success=False, x=a[1], fun=1, nit=1, message="test failure"
        ),
    )
    with pytest.raises(ValueError, match="no constraints were dropped"):
        fit_dates(small_chronology(), starts=1)


def test_scale_invariance_of_time_and_substitution_units():
    c = small_chronology()
    before, _ = fit_dates(c)
    time_scaled = replace(c, scale=c.scale * 100)
    after, _ = fit_dates(time_scaled)
    np.testing.assert_allclose(before.ages, after.ages)
    np.testing.assert_allclose(before.rates / 100, after.rates)
    for node in c.edges:
        node.dist *= 100
    after, _ = fit_dates(c)
    np.testing.assert_allclose(before.ages, after.ages, atol=1e-6)
    np.testing.assert_allclose(before.rates * 100, after.rates, rtol=1e-6)

import hashlib
import itertools
import json
import math
from collections import Counter
from dataclasses import replace

import numpy as np
import pandas as pd
import pytest

from nwkit.cli import main
from nwkit.mul_locus import (
    Locus,
    LocusCoalescent,
    LocusParameters,
    detected_signature,
    sample_locus_tree,
    sample_selected,
    signature_size,
)
from nwkit.mul_locus_cli import finite_json, json_pairs
from nwkit.mul_locus_mc import (
    LocusBank,
    build_bank,
    calibrate,
    category_bound,
    make_tasks,
    pattern_probability,
    probability_interval,
    search_banks,
    validate_model,
)
from nwkit.mul_msc_fit import species_topology_signature
from nwkit.wgd_tree_model import lineage_log_parameters
from tests.test_mul_coalescent import (
    oracle_ancestral,
    oracle_branch,
    oracle_distribution,
)
from tests.test_mul_reconcile import parser, tree

SPECIES = "((A:1,X:1):1,B:2);"


def model(samples=1000):
    return {
        "schema": "nwkit-mul-locus-mc-model-v1",
        "copy_role": "distinct-loci",
        "root_locus_count": 1,
        "species_time_unit": "generations",
        "ancestral_stem": 0.5,
        "detection": {"A": 1, "B": 1, "X": 1},
        "max_observed_tips": 4,
        "parameter_grid": [
            {"duplication": 0.06, "loss": 0.03, "ne": 0.5, "hybridization_age": 0.5}
        ],
        "samples": samples,
        "seed": 20261020,
        "confidence": 0.99,
        "max_attempts": 100000,
        "max_locus_nodes": 1000,
        "max_coalescent_states": 10000,
    }


def color(text):
    return species_topology_signature(tree(text + ";"), parser())


def strip_ids(gene):
    return detected_signature(
        gene, {name: 1 for name in ("A", "B", "C", "D", "X")}, np.random.default_rng(1)
    )


def test_count_dp_matches_independent_full_forest_and_daughter_normalization():
    daughter = Locus(
        1,
        "speciation",
        children=[Locus(0, "tip", "A"), Locus(0, "tip", "B")],
        daughter=True,
    )
    root = Locus(1.5, "duplication", children=[Locus(0, "tip", "C"), daughter])
    solver = LocusCoalescent(root, 0.5)
    assert solver.base[daughter] == {2: 0.0}
    oracle = oracle_branch(("A", "B"), 0.5)
    assert sum(p for state, p in oracle.items() if len(state) == 1) == pytest.approx(
        1 - math.exp(-0.5), abs=2e-14
    )
    assert solver.output[daughter] == {1: 0.0}
    assert math.exp(solver.bound_log_normalizers[daughter]) == pytest.approx(
        1 - math.exp(-0.5), abs=2e-14
    )
    rng = np.random.default_rng(55)
    assert {strip_ids(solver.sample(rng)) for _ in range(100)} == {
        color("((a_A,b_B),c_C)")
    }


@pytest.mark.slow
def test_zero_dl_sample_distribution_matches_full_forest_oracle():
    population = tree("((A:1,B:1):0.5,C:1.5);")
    locus = sample_locus_tree(
        population, LocusParameters(0, 0, 0.5, 0.5), 0, np.random.default_rng(5)
    )
    sampler = LocusCoalescent(locus, 0.5)
    expected = oracle_distribution(population, {n.name: n for n in population.leaves()})
    rng = np.random.default_rng(20261021)
    counts = Counter(strip_ids(sampler.sample(rng)) for _ in range(8000))
    assert sum(counts.values()) == 8000
    for text, p in expected.items():
        lo, hi = probability_interval(
            counts[
                color(text.replace("A", "a_A").replace("B", "b_B").replace("C", "c_C"))
            ],
            8000,
            0.001,
        )
        assert lo <= p <= hi


@pytest.mark.slow
@pytest.mark.parametrize("candidate", [1, 2])
def test_zero_dl_allopolyploid_distribution_matches_full_colored_forest(candidate):
    point = LocusParameters(0, 0, 0.7, 0.5)
    task = next(
        t
        for t in make_tasks(tree(SPECIES), "X", "A B", [point], 100)[0]
        if t[0] == candidate
    )
    population = task[3]
    oracle_population = population.copy()
    assignment = {}
    for i, leaf in enumerate(oracle_population.leaves()):
        assignment[f"g{i}_{leaf.props.get('mul_species', leaf.name)}"] = leaf
    for node in oracle_population.traverse():
        node.dist /= 2 * point.ne
    expected = Counter()
    for topology, probability in oracle_distribution(
        oracle_population, assignment
    ).items():
        expected[color(topology)] += probability
    assert sum(expected.values()) == pytest.approx(1, abs=3e-13)
    locus = sample_locus_tree(population, point, 0.5, np.random.default_rng(1))
    sampler = LocusCoalescent(locus, point.ne)
    rng = np.random.default_rng(20261028 + candidate)
    counts = Counter(strip_ids(sampler.sample(rng)) for _ in range(8000))
    for signature, probability in expected.items():
        lo, hi = probability_interval(counts[signature], 8000, 0.0001)
        assert lo <= probability <= hi


@pytest.mark.slow
def test_three_tip_daughter_distribution_matches_independent_bounded_forest():
    ab = Locus(0.5, "speciation", children=[Locus(0, "tip", "A"), Locus(0, "tip", "B")])
    daughter = Locus(
        1, "speciation", children=[ab, Locus(0, "tip", "C")], daughter=True
    )
    root = Locus(2.5, "duplication", children=[Locus(0, "tip", "D"), daughter])
    sampler = LocusCoalescent(root, 0.5)
    expected = Counter()
    for forest, p in oracle_branch(("A", "B"), 0.5).items():
        for final, q in oracle_branch(tuple(sorted((*forest, "C"))), 1.5).items():
            if len(final) == 1:
                expected[final[0]] += p * q
    total = sum(expected.values())
    assert math.exp(sampler.bound_log_normalizers[daughter]) == pytest.approx(
        total, abs=2e-13
    )
    rng = np.random.default_rng(20261024)
    counts = Counter(strip_ids(sampler.sample(rng)) for _ in range(8000))
    for topology, p in expected.items():
        names = topology.replace("A", "a_A").replace("B", "b_B").replace("C", "c_C")
        signature = color("(" + names + ",d_D)")
        lo, hi = probability_interval(counts[signature], 8000, 0.001)
        assert lo <= p / total <= hi


def joint_locus_forest(node, parent_age, ne, labels):
    """Unnormalized all-pair histories: impose all daughter bounds jointly."""
    current = Counter({(labels[node],) if node.kind == "tip" else (): 1.0})
    for child in node.children:
        joined = Counter()
        for (first, p), (second, q) in itertools.product(
            current.items(), joint_locus_forest(child, node.age, ne, labels).items()
        ):
            joined[tuple(sorted(first + second))] += p * q
        current = joined
    if parent_age is None:
        return current
    result = Counter()
    for forest, p in current.items():
        for final, q in oracle_branch(
            forest, (parent_age - node.age) / (2 * ne)
        ).items():
            if not node.daughter or len(final) <= 1:
                result[final] += p * q
    return result


def nested_locus_with_loss():
    a, b, x1, x2, c = (Locus(0, "tip", s) for s in ("A", "B", "X", "X", "C"))
    surviving = Locus(0.2, "speciation", children=[c, Locus(0.1, "loss")])
    inner = Locus(0.4, "speciation", children=[x2, surviving], daughter=True)
    mother = Locus(0.7, "speciation", children=[b, x1])
    outer = Locus(1.4, "duplication", children=[mother, inner], daughter=True)
    original = Locus(1, "speciation", children=[a, Locus(0.3, "loss")])
    root = Locus(
        3, "origin", children=[Locus(2, "duplication", children=[original, outer])]
    )
    labels = {tip: f"g{i}_{tip.species}" for i, tip in enumerate((a, b, x1, x2, c))}
    return root, labels


@pytest.mark.slow
def test_nested_daughters_loss_and_detection_match_joint_forest_oracle():
    root, labels = nested_locus_with_loss()
    ne = 0.7
    forest = joint_locus_forest(root, None, ne, labels)
    joint_bound = math.fsum(forest.values())
    sampler = LocusCoalescent(root, ne)
    assert len(sampler.bound_log_normalizers) == 2
    assert math.exp(math.fsum(sampler.bound_log_normalizers.values())) == pytest.approx(
        joint_bound, abs=2e-13
    )
    topologies = Counter()
    for state, p in forest.items():
        for topology, q in oracle_ancestral(state).items():
            topologies[topology] += p * q / joint_bound
    assert math.fsum(topologies.values()) == pytest.approx(1, abs=2e-13)
    detection = {"A": 0.85, "B": 0.8, "X": 0.65, "C": 0.4}
    expected = Counter()
    for topology, p in topologies.items():
        gene = tree(topology + ";")
        tips = list(gene.leaf_names())
        for mask in itertools.product((False, True), repeat=len(tips)):
            kept = [name for name, present in zip(tips, mask, strict=True) if present]
            if not 2 <= len(kept) <= 4:
                continue
            probability = math.prod(
                detection[name.rsplit("_", 1)[1]]
                if present
                else 1 - detection[name.rsplit("_", 1)[1]]
                for name, present in zip(tips, mask, strict=True)
            )
            pruned = gene.copy()
            pruned.prune(kept)
            expected[species_topology_signature(pruned, parser())] += p * probability
    selection = math.fsum(expected.values())
    rng = np.random.default_rng(20261030)
    counts = Counter()
    for _ in range(12000):
        signature = detected_signature(sampler.sample(rng), detection, rng)
        if 2 <= signature_size(signature) <= 4:
            counts[signature] += 1
    total = sum(counts.values())
    lo, hi = probability_interval(total, 12000, 0.0001)
    assert lo <= selection <= hi
    assert counts.keys() <= expected.keys()
    for signature, p in expected.items():
        lo, hi = probability_interval(counts[signature], total, 0.0001)
        assert lo <= p / selection <= hi


@pytest.mark.slow
@pytest.mark.parametrize("duplication,loss", [(0.2, 0), (0.2, 0.2), (0.1, 0.3)])
def test_ancestral_linear_bd_counts_match_independent_closed_form(duplication, loss):
    population = tree("(A:0,B:0);")
    rng = np.random.default_rng(20261022)
    counts = Counter()
    for _ in range(5000):
        locus = sample_locus_tree(
            population, LocusParameters(duplication, loss, 1, 0.5), 1, rng
        )
        pending, a = [locus], 0
        while pending:
            node = pending.pop()
            a += node.kind == "tip" and node.species == "A"
            pending.extend(node.children)
        counts[a] += 1
    extinction, geometric, survival, success = lineage_log_parameters(
        duplication, loss, 1
    )
    for n in range(5):
        p = math.exp(extinction if n == 0 else survival + success + (n - 1) * geometric)
        lo, hi = probability_interval(counts[n], 5000, 0.001)
        assert lo <= p <= hi


def test_detection_is_applied_after_full_bounded_genealogy():
    gene = ("node", ("tip", "A", 0), ("node", ("tip", "X", 1), ("tip", "X", 2)))
    assert detected_signature(gene, {"A": 1, "X": 0}, np.random.default_rng(1)) == (
        "tip",
        "A",
    )
    assert signature_size(gene) == 3


def test_singleton_sampling_preserves_numpy_choice_rng_stream():
    actual = np.random.default_rng(53)
    expected = np.random.default_rng(53)
    assert LocusCoalescent.choose([("value", -1000)], actual) == "value"
    assert expected.choice(1, p=[1.0]) == 0
    assert actual.bit_generator.state == expected.bit_generator.state
    with pytest.raises(ArithmeticError, match="weights"):
        LocusCoalescent.choose([("value", -math.inf)], actual)


def test_resource_caps_fail_without_discarding_histories():
    population = tree(SPECIES)
    point = LocusParameters(0, 0, 0.5, 0.5)
    with pytest.raises(ValueError, match="node cap"):
        sample_locus_tree(population, point, 0, np.random.default_rng(1), max_nodes=1)
    locus = sample_locus_tree(population, point, 0, np.random.default_rng(1))
    with pytest.raises(ValueError, match="state cap"):
        LocusCoalescent(locus, 0.5, max_states=1)
    config = model()
    config["detection"] = {name: 0 for name in config["detection"]}
    config["max_attempts"] = 2
    with pytest.raises(ValueError, match="no partial bank"):
        sample_selected(population, point, config, np.random.default_rng(1), count=1)


def test_state_cap_is_checked_before_dense_transition_allocation(monkeypatch):
    import nwkit.mul_locus as module

    daughters = [Locus(0, "tip", "A"), Locus(0, "tip", "B")]
    locus = Locus(1, "origin", children=[Locus(0, "speciation", children=daughters)])
    original = module.transition

    def guarded(k, j, time):
        assert k - j + 1 < 2, "Dense transition called before work cap"
        return original(k, j, time)

    monkeypatch.setattr(module, "transition", guarded)
    with pytest.raises(ValueError, match="state cap"):
        LocusCoalescent(locus, 0.5, max_states=7)


@pytest.mark.parametrize(
    "update",
    [
        {"copy_role": "alleles"},
        {"species_time_unit": "years"},
        {"root_locus_count": 2},
        {"root_locus_count": True},
        {"samples": 0},
        {"samples": True},
        {"seed": -1},
        {"confidence": 1},
        {"detection": {"A": 0}},
        {"parameter_grid": []},
        {"ancestral_stem": math.inf},
        {"unknown": 1},
    ],
)
def test_explicit_model_validation(update):
    with pytest.raises(ValueError):
        validate_model(model() | update, tree(SPECIES))


@pytest.mark.parametrize("field", ["duplication", "loss", "ne", "hybridization_age"])
@pytest.mark.parametrize("value", [True, "1", None, 10**500])
def test_scientific_grid_parameters_require_finite_numbers(field, value):
    config = model()
    config["parameter_grid"][0][field] = value
    with pytest.raises(ValueError):
        validate_model(config, tree(SPECIES))


@pytest.mark.parametrize(
    "update",
    [
        {"ancestral_stem": True},
        {"detection": {"A": True, "B": 1, "X": 1}},
        {"detection": ["A", "B", "X"]},
        {"detection": None},
        {"confidence": "0.99"},
        {"parameter_grid": [None]},
        {"parameter_grid": "grid"},
        {"ancestral_stem": 10**500},
        {"confidence": True},
        {"integration": {}},
    ],
)
def test_scientific_json_types_fail_with_validation_error(update):
    with pytest.raises(ValueError):
        validate_model(model() | update, tree(SPECIES))


@pytest.mark.parametrize("config", [None, [], "model", 1, True])
def test_scientific_model_requires_json_object(config):
    with pytest.raises(ValueError):
        validate_model(config, tree(SPECIES))


@pytest.mark.parametrize("rate", [1e308, 10**308])
def test_combined_scientific_rate_overflow_fails_validation(rate):
    with pytest.raises(ValueError, match="Combined locus rate"):
        LocusParameters(rate, rate, 1, 0.5)


@pytest.mark.parametrize("ne", [1e308, 10**308])
def test_diploid_scale_overflow_fails_validation(ne):
    with pytest.raises(ValueError, match="2\\*Ne"):
        LocusParameters(0, 0, ne, 0.5)
    with pytest.raises(ValueError, match="Ne"):
        LocusCoalescent(Locus(0, "origin"), ne)


def test_mc_zero_support_has_no_pseudocount_fallback():
    signature = color("(a_A,x_X)")
    other = color("(a_A,b_B)")
    point = LocusParameters(0, 0, 0.5, 0.5)
    banks = [
        LocusBank(i, "A", 0, None, point, Counter({other: 100}), 100, 100)
        for i in (0, 1)
    ]
    with pytest.raises(ValueError, match="no pseudocount"):
        search_banks(banks, [signature], 0.001)
    assert probability_interval(0, 100, 0.001)[0] == 0
    assert finite_json({"log_likelihood": -math.inf}) == {"log_likelihood": None}
    with pytest.raises(ValueError, match="Duplicate"):
        json.loads('{"a":1,"a":2}', object_pairs_hook=json_pairs)


def test_stratified_prior_underflow_fails_instead_of_dropping_a_stratum():
    config = model(10) | {
        "integration": "ancestral-stratified",
        "ancestral_stem": 1e-300,
    }
    point = LocusParameters(1e-300, 1e300, 0.5, 0.5)
    with pytest.raises(ArithmeticError, match="birth stratum prior underflowed"):
        build_bank((0, "NA", 0, tree(SPECIES), point), config)


def test_stratum_selection_mass_and_prior_weights_are_normalized_together():
    a, b = color("(a_A,x_X)"), color("(a_A,b_B)")
    strata = [
        {
            "weight": 0.9,
            "samples": 100,
            "selected": 50,
            "counts": Counter({a: 30, b: 20}),
        },
        {
            "weight": 0.05,
            "samples": 100,
            "selected": 80,
            "counts": Counter({a: 20, b: 60}),
        },
    ]
    bank = LocusBank(
        0,
        "NA",
        0,
        None,
        LocusParameters(0.1, 0.1, 1, 0.5),
        Counter({a: 50, b: 80}),
        200,
        200,
        strata,
    )
    p, lo, hi = pattern_probability(bank, a, 0.001)
    assert p == pytest.approx((0.9 * 0.3 + 0.05 * 0.2) / (0.9 * 0.5 + 0.05 * 0.8))
    assert lo <= p <= hi
    assert p + pattern_probability(bank, b, 0.001)[0] == pytest.approx(1)


@pytest.mark.slow
def test_stratified_ancestral_dl_selection_matches_closed_form():
    population = tree("(A:0,B:0);")
    config = model(1600)
    config.update(
        ancestral_stem=1, detection={"A": 1, "B": 1}, integration="ancestral-stratified"
    )
    point = LocusParameters(0.2, 0.1, 0.5, 0.5)
    bank = build_bank((0, "NA", 0, population, point), config)
    _, geometric, survival, success = lineage_log_parameters(0.2, 0.1, 1)
    p1 = math.exp(survival + success)
    p2 = math.exp(survival + success + geometric)
    estimated, lo, hi = pattern_probability(bank, color("(a_A,b_B)"), 0.001)
    assert lo <= p1 / (p1 + p2) <= hi
    assert 0 < estimated < 1
    assert sum(
        pattern_probability(bank, signature, 0.001)[0] for signature in bank.counts
    ) == pytest.approx(1)
    assert bank.attempts == 1600


def test_bootstrap_refits_full_candidate_grid_and_finite_sample_correction(monkeypatch):
    import nwkit.mul_locus_mc as module

    a, b = color("(a_A,x_X)"), color("(a_A,b_B)")
    point = LocusParameters(0.1, 0.1, 1, 0.5)
    banks = [
        LocusBank(0, "NA", 0, None, point, Counter({a: 90, b: 10}), 100, 100),
        LocusBank(
            0, "NA", 1, None, replace(point, ne=2), Counter({a: 10, b: 90}), 100, 100
        ),
        LocusBank(1, "A", 0, None, point, Counter({a: 95, b: 5}), 100, 100),
        LocusBank(2, "B", 1, None, point, Counter({a: 5, b: 95}), 100, 100),
    ]
    draws = iter([[a], [b]])
    monkeypatch.setattr(
        module, "sample_selected", lambda *args, **kwargs: (next(draws), 1)
    )
    fit, calibration = calibrate(banks, [a], model(), 0.001, 2)
    assert fit["null"]["grid"] == 0
    assert [r["null_grid"] for r in calibration["rows"]] == [0, 1]
    assert [r["alternative_candidate"] for r in calibration["rows"]] == [1, 2]
    assert calibration["p_value"] == 1
    assert (
        calibration["mc_p_lower"] <= calibration["p_value"] <= calibration["mc_p_upper"]
    )
    assert category_bound(2, 3) == 28


def arguments(tmp_path, samples=400):
    (tmp_path / "genes.nwk").write_text("((a_A,x_X),b_B);\n")
    (tmp_path / "species.nwk").write_text(SPECIES)
    (tmp_path / "config.json").write_text(json.dumps(model(samples)))
    paths = {
        name: tmp_path / name
        for name in ("scores", "report", "checks", "model", "calibration")
    }
    args = [
        "mul-reconcile",
        "--score-model",
        "locus-mc",
        "-i",
        str(tmp_path / "genes.nwk"),
        "--species-tree",
        str(tmp_path / "species.nwk"),
        "--species-regex",
        r".*_([^_]+)$",
        "--h1",
        "X",
        "--locus-model",
        str(tmp_path / "config.json"),
        "--locus-bootstrap",
        "2",
        "-o",
        str(paths["scores"]),
        "--report",
        str(paths["report"]),
        "--check-out",
        str(paths["checks"]),
        "--model-out",
        str(paths["model"]),
        "--locus-calibration-out",
        str(paths["calibration"]),
    ]
    return args, paths


@pytest.mark.parametrize(
    "integration",
    ["selected-histogram", "ancestral-stratified", "detection-rb", "hybrid-rb"],
)
def test_cli_serial_parallel_schema_and_full_bundle(tmp_path, integration, capsys):
    args, paths = arguments(tmp_path)
    config = json.loads((tmp_path / "config.json").read_text())
    config["integration"] = integration
    if integration in ("detection-rb", "hybrid-rb"):
        # Raw draws are not the selected-family budget of the legacy fixture.
        config["samples"] = 2000
    (tmp_path / "config.json").write_text(json.dumps(config))
    main(args)
    first = [p.read_bytes() for p in paths.values()]
    main([*args, "--cpus", "2"])
    assert [p.read_bytes() for p in paths.values()] == first
    saved = json.loads(paths["model"].read_text())
    captured = capsys.readouterr()
    assert captured.out == ""
    assert f"P={saved['calibration']['p_value']:.8g}" in captured.err
    assert "MC-overlap interval" in captured.err
    assert saved["schema"] == "nwkit-mul-locus-mc-v1"
    assert saved["scores"]["null"]["mul.tree"] == 0
    assert len(saved["calibration"]["rows"]) == 2
    assert 0 in pd.read_csv(paths["scores"], sep="\t")["mul.tree"].values
    if integration in ("detection-rb", "hybrid-rb"):
        from nwkit.mul_locus_cli import finite_json
        from nwkit.mul_locus_integral import (
            IntegratedBank,
            Moments,
            score_integrated_bank,
        )
        from nwkit.mul_locus_mc import search_banks

        def freeze(value):
            return tuple(freeze(v) for v in value) if isinstance(value, list) else value

        restored = []
        for bank in saved["banks"]:
            population = tree(bank["population_tree"].removeprefix("[&R]"))
            for leaf in population.leaves():
                leaf.props["mul_species"] = bank["population_leaf_species"][leaf.name]
            strata = [
                {
                    **s,
                    "selected": Moments(**s["selected"]),
                    "patterns": {
                        freeze(p["signature"]): Moments(p["n"], p["mean"], p["m2"])
                        for p in s["patterns"]
                    },
                }
                for s in bank["strata"]
            ]
            restored.append(
                IntegratedBank(
                    bank["candidate"],
                    bank["h2"],
                    bank["grid"],
                    population,
                    LocusParameters(**saved["model"]["parameter_grid"][bank["grid"]]),
                    strata,
                    bank["method"],
                    bank["interval_method"],
                )
            )
        replay = search_banks(
            restored,
            [freeze(s) for s in saved["observations"]],
            saved["per_probability_alpha"],
            scorer=score_integrated_bank,
        )
        assert finite_json(replay) == saved["scores"]


@pytest.mark.parametrize("integration", ["selected-histogram", "ancestral-stratified"])
def test_saved_population_and_observations_are_self_contained(tmp_path, integration):
    args, paths = arguments(tmp_path)
    config_path = tmp_path / "config.json"
    config = json.loads(config_path.read_text()) | {"integration": integration}
    config_path.write_text(json.dumps(config))
    main(args)
    saved = json.loads(paths["model"].read_text())
    assert saved["hypothesis_scope"]["h1"] == "X"
    assert saved["hypothesis_scope"]["h2"] is None
    original = tree(saved["hypothesis_scope"]["species_tree"].removeprefix("[&R]"))
    assert list(original.leaf_names()) == ["A", "X", "B"]
    observed = saved["observations"]
    assert len(observed) == saved["num_gene_trees"] == 1
    assert observed[0] == json.loads(json.dumps(color("((a_A,x_X),b_B)")))

    def restore(value):
        return tuple(restore(v) if isinstance(v, list) else v for v in value)

    def counts(values):
        return Counter({restore(v["signature"]): v["hits"] for v in values})

    rebuilt = []
    for bank in saved["banks"]:
        population = tree(bank["population_tree"].removeprefix("[&R]"))
        assert set(population.leaf_names()) == set(bank["population_leaf_species"])
        assert set(bank["population_leaf_species"].values()) <= {"A", "B", "X"}
        for leaf in population.leaves():
            leaf.props["mul_species"] = bank["population_leaf_species"][leaf.name]
        if bank["candidate"] == 0:
            assert bank["parameters"]["hybridization_age"] is None
        else:
            assert bank["parameters"]["hybridization_age"] == 0.5
            assert sum(v == "X" for v in bank["population_leaf_species"].values()) == 2
        strata = bank["strata"]
        if strata is not None:
            strata = [{**s, "counts": counts(s["counts"])} for s in strata]
        restored_bank = LocusBank(
            bank["candidate"],
            bank["h2"],
            bank["grid"],
            population,
            LocusParameters(**saved["model"]["parameter_grid"][bank["grid"]]),
            counts(bank["counts"]),
            bank["samples"],
            bank["attempts"],
            strata,
        )
        regenerated = build_bank(
            (
                restored_bank.candidate,
                restored_bank.h2,
                restored_bank.grid,
                restored_bank.population,
                restored_bank.parameters,
            ),
            saved["model"],
        )
        assert regenerated.counts == restored_bank.counts
        assert regenerated.samples == restored_bank.samples
        assert regenerated.attempts == restored_bank.attempts
        assert regenerated.strata == restored_bank.strata
        rebuilt.append(restored_bank)
    rescored = search_banks(
        rebuilt,
        [restore(v) for v in saved["observations"]],
        saved["per_probability_alpha"],
    )
    assert json.loads(json.dumps(finite_json(rescored))) == saved["scores"]
    assert set(pd.read_csv(paths["report"], sep="\t")["integration"]) == {integration}


def test_cli_failures_and_aliases_preserve_all_files(tmp_path, monkeypatch):
    import nwkit.mul_locus_cli as module

    args, paths = arguments(tmp_path)
    for name, path in paths.items():
        path.write_text("old " + name)
    with pytest.raises(ValueError, match="overwrite"):
        main([*args, "--locus-model", str(paths["model"])])
    assert all(p.read_text() == "old " + name for name, p in paths.items())

    def fail(*a, **kw):
        raise OSError("injected locus writer")

    monkeypatch.setattr(module.json, "dump", fail)
    with pytest.raises(OSError, match="injected"):
        main(args)
    assert all(p.read_text() == "old " + name for name, p in paths.items())
    assert not list(tmp_path.glob(".*.stage.*"))


def test_audit_hashes_locus_model_and_all_five_outputs(tmp_path):
    args, paths = arguments(tmp_path)
    audit = tmp_path / "audit.jsonl"
    main([*args, "--audit", str(audit)])
    saved = json.loads(audit.read_text())
    inputs = {r["argument"]: r for r in saved["inputs"]}
    assert (
        inputs["locus_model"]["sha256"]
        == hashlib.sha256((tmp_path / "config.json").read_bytes()).hexdigest()
    )
    outputs = {r["argument"]: r for r in saved["outputs"]}
    for argument, path in zip(
        ("outfile", "report", "check_out", "model_out", "locus_calibration_out"),
        paths.values(),
        strict=True,
    ):
        assert (
            outputs[argument]["sha256"] == hashlib.sha256(path.read_bytes()).hexdigest()
        )


@pytest.mark.parametrize("max_tips", [100, 1000, 10**9])
def test_unrepresentable_confidence_precision_fails_before_sampling(
    tmp_path, monkeypatch, max_tips
):
    import nwkit.mul_locus_cli as module

    args, paths = arguments(tmp_path)
    config = model(10) | {"max_observed_tips": max_tips}
    (tmp_path / "config.json").write_text(json.dumps(config))

    def fail(*a, **kw):
        pytest.fail("Precision failure should precede bank simulation")

    monkeypatch.setattr(module, "build_bank", fail)
    with pytest.raises(ValueError, match="confidence precision"):
        main(args)
    assert all(not p.exists() for p in paths.values())


@pytest.mark.parametrize("target", ["input", "checks", "calibration"])
def test_audit_alias_cannot_modify_model_or_companion_outputs(tmp_path, target):
    args, paths = arguments(tmp_path)
    original = (tmp_path / "config.json").read_bytes()
    for name, path in paths.items():
        path.write_text("old " + name)
    audit = tmp_path / "config.json" if target == "input" else paths[target]
    with pytest.raises(ValueError):
        main([*args, "--audit", str(audit)])
    assert (tmp_path / "config.json").read_bytes() == original
    assert all(p.read_text() == "old " + name for name, p in paths.items())


@pytest.mark.parametrize("score", ["dl", "msc"])
def test_locus_options_in_other_modes_fail(tmp_path, score):
    args, paths = arguments(tmp_path)
    with pytest.raises(ValueError, match="require --score-model locus-mc"):
        main([*args, "--score-model", score])
    assert not paths["scores"].exists()


def test_candidate_grid_keeps_null_and_temporal_exclusions():
    config = model(10)
    points = validate_model(config, tree(SPECIES))
    tasks, excluded = make_tasks(tree(SPECIES), "X", None, points, 100)
    assert tasks[0][0] == 0
    assert any(r["reason"].startswith("Autopolyploid") for r in excluded)
    bank = build_bank(tasks[0], config)
    assert sum(bank.counts.values()) == 10 and bank.attempts >= 10

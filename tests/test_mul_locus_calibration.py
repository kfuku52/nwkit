import itertools
import json
import math
from collections import Counter

import numpy as np
import pandas as pd
import pytest

from nwkit.cli import main
from nwkit.mul_locus import LocusParameters
from nwkit.mul_locus_mc import LocusBank, calibrate, search_banks
from tests.test_mul_locus import arguments, model


def toy_banks():
    return [
        LocusBank(0, "NA", grid, None, LocusParameters(0, 0, ne, 0.5), Counter(), 1, 1)
        for grid, ne in ((3, 1), (7, 2))
    ]


def toy_search(banks, data, alpha):
    value = data[0]
    return {
        "null": {"grid": 3},
        "alternative": {"grid": 0, "mul.tree": 1},
        "contrast": value,
        "contrast_lower": value,
        "contrast_upper": value,
    }


def scored_banks():
    signature = ("node", ("tip", "A"), ("tip", "B"))
    other = ("node", ("tip", "A"), ("tip", "A"))
    banks = toy_banks()
    banks[1].candidate = 1
    for bank, hits in zip(banks, (3, 1), strict=True):
        bank.counts = Counter({signature: hits, other: 4 - hits})
        bank.samples = bank.attempts = 4
    return banks, signature


def test_one_shot_observations_are_identical_for_every_bank():
    banks, signature = scored_banks()
    observations = [signature, signature]
    assert search_banks(banks, iter(observations), 0.001) == search_banks(
        banks, observations, 0.001
    )


@pytest.mark.parametrize("mode", ["plug-in", "grid-supremum"])
def test_calibration_accepts_one_shot_complete_family_sets(mode):
    banks, signature = scored_banks()

    def sampler(pop, point, config, rng, *, count):
        return iter([signature] * count), count

    def list_sampler(pop, point, config, rng, *, count):
        return [signature] * count, count

    assert calibrate(
        banks,
        iter([signature, signature]),
        model(),
        0.001,
        2,
        sampler=sampler,
        null_calibration=mode,
    ) == calibrate(
        banks,
        [signature, signature],
        model(),
        0.001,
        2,
        sampler=list_sampler,
        null_calibration=mode,
    )


@pytest.mark.parametrize("mode", ["plug-in", "grid-supremum"])
def test_calibration_reuses_one_shot_banks(mode):
    banks, signature = scored_banks()

    def sampler(pop, point, config, rng, *, count):
        return [signature] * count, count

    assert calibrate(
        iter(banks),
        [signature],
        model(),
        0.001,
        2,
        sampler=sampler,
        null_calibration=mode,
    ) == calibrate(
        banks, [signature], model(), 0.001, 2, sampler=sampler, null_calibration=mode
    )


@pytest.mark.parametrize("mode", ["plug-in", "grid-supremum"])
@pytest.mark.parametrize("size", [0, 2])
def test_incomplete_or_extra_null_families_fail_before_refit(monkeypatch, mode, size):
    import nwkit.mul_locus_mc as module

    searches = []

    def search(banks, data, alpha):
        searches.append(tuple(data))
        return toy_search(banks, data, alpha)

    def sampler(pop, point, config, rng, *, count):
        return [0] * size, max(1, size)

    monkeypatch.setattr(module, "search_banks", search)
    with pytest.raises(ValueError, match=f"expected 1 families, got {size}"):
        calibrate(
            toy_banks(),
            [1],
            model(),
            0.001,
            2,
            sampler=sampler,
            null_calibration=mode,
        )
    assert searches == [(1,)]


@pytest.mark.parametrize("mode", ["plug-in", "grid-supremum"])
def test_numerical_failure_preserves_type_and_replicate_context(monkeypatch, mode):
    import nwkit.mul_locus_mc as module

    monkeypatch.setattr(module, "search_banks", toy_search)

    def sampler(*a, **kw):
        raise ArithmeticError("unrepresentable bound")

    context = " generating grid 3" if mode == "grid-supremum" else ""
    with pytest.raises(
        ArithmeticError,
        match=f"Null bootstrap{context} replicate 1: unrepresentable bound",
    ):
        calibrate(
            toy_banks(),
            [1],
            model(),
            0.001,
            2,
            sampler=sampler,
            null_calibration=mode,
        )


def test_grid_supremum_evaluates_every_null_and_refits_search(monkeypatch):
    import nwkit.mul_locus_mc as module

    monkeypatch.setattr(module, "search_banks", toy_search)
    draws = {1: iter((0, 0, 0)), 2: iter((1, 1, 0))}
    called = []

    def sampler(population, point, config, rng, *, count):
        called.append(point.ne)
        return [next(draws[point.ne])], 1

    observed, result = calibrate(
        toy_banks()[::-1],
        [1],
        model(),
        0.001,
        3,
        sampler=sampler,
        null_calibration="grid-supremum",
    )
    assert observed["null"]["grid"] == 3
    assert called == [1, 1, 1, 2, 2, 2]
    assert result["p_value"] == result["mc_p_lower"] == result["mc_p_upper"] == 0.75
    assert [s["p_value"] for s in result["null_grid_calibrations"]] == [0.25, 0.75]
    assert result["least_favorable_grid"] == 7
    assert result["num_null_grid_points"] == 2
    assert result["total_replicates"] == 6
    assert [r["generating_null_grid"] for r in result["rows"]] == [3] * 3 + [7] * 3
    assert all(r["null_grid"] == 3 for r in result["rows"])
    assert all(
        r["null_generating_parameters"]["hybridization_age"] is None
        for r in result["null_grid_calibrations"]
    )


def test_exhaustive_discrete_null_rank_size_with_supremum_and_ties(monkeypatch):
    import nwkit.mul_locus_mc as module

    monkeypatch.setattr(module, "search_banks", toy_search)
    probabilities = (0.2, 0.8)
    rejection = [0.0, 0.0]
    mass = 0.0
    for first, second in itertools.product(
        itertools.product((0, 1), repeat=3), repeat=2
    ):
        weight = math.prod(
            p if draw else 1 - p
            for p, draws in zip(probabilities, (first, second), strict=True)
            for draw in draws
        )
        mass += weight
        for observed in (0, 1):
            pools = {1: iter(first), 2: iter(second)}

            def sampler(pop, point, config, rng, *, count, pools=pools):
                return [next(pools[point.ne])], 1

            _, result = calibrate(
                toy_banks(),
                [observed],
                model(),
                0.001,
                3,
                sampler=sampler,
                null_calibration="grid-supremum",
            )
            if result["p_value"] <= 0.25:
                for i, p in enumerate(probabilities):
                    rejection[i] += weight * (p if observed else 1 - p)
    assert mass == pytest.approx(1, abs=2e-14)
    assert rejection == pytest.approx([p * 0.8**3 * 0.2**3 for p in probabilities])
    assert all(p <= 0.25 for p in rejection)


@pytest.mark.parametrize("mode", ["plug-in", "grid-supremum"])
def test_calibration_seed_namespaces_and_reordered_banks(monkeypatch, mode):
    import nwkit.mul_locus_mc as module

    monkeypatch.setattr(module, "search_banks", toy_search)
    seen = []

    def sampler(pop, point, config, rng, *, count):
        seen.append(int(rng.integers(2**31)))
        return [0], 1

    calibrate(
        toy_banks(), [1], model(), 0.001, 3, sampler=sampler, null_calibration=mode
    )
    expected = [
        int(np.random.default_rng(np.random.SeedSequence(seed)).integers(2**31))
        for grid in ((3, 7) if mode == "grid-supremum" else (3,))
        for replicate in range(3)
        for seed in (
            [model()["seed"], 2, grid, replicate]
            if mode == "grid-supremum"
            else [model()["seed"], 1, replicate],
        )
    ]
    assert seen == expected
    seen.clear()
    calibrate(
        toy_banks()[::-1],
        [1],
        model(),
        0.001,
        3,
        sampler=sampler,
        null_calibration=mode,
    )
    assert seen == expected


def test_point_and_mc_upper_can_have_different_least_favorable_grids(monkeypatch):
    import nwkit.mul_locus_mc as module

    def uncertain_search(banks, data, alpha):
        fit = toy_search(banks, data, alpha)
        if data[0] == 0:
            fit["contrast_upper"] = 2
        return fit

    monkeypatch.setattr(module, "search_banks", uncertain_search)

    def sampler(pop, point, config, rng, *, count):
        return [int(point.ne == 2)], 1

    _, result = calibrate(
        toy_banks(),
        [1],
        model(),
        0.001,
        3,
        sampler=sampler,
        null_calibration="grid-supremum",
    )
    assert result["least_favorable_grid"] == 7
    assert result["least_favorable_mc_upper_grid"] == 3
    assert result["p_value"] == result["mc_p_upper"] == 1


@pytest.mark.parametrize(
    "mode,replicates",
    [
        ("other", 3),
        ("grid-supremum", 0),
        ("grid-supremum", -1),
        ("grid-supremum", True),
    ],
)
def test_invalid_calibration_settings_fail_before_search(monkeypatch, mode, replicates):
    import nwkit.mul_locus_mc as module

    monkeypatch.setattr(
        module, "search_banks", lambda *a: pytest.fail("search started")
    )
    with pytest.raises(ValueError):
        calibrate(toy_banks(), [1], model(), 0.001, replicates, null_calibration=mode)


def test_grid_failure_identifies_generating_point_and_replicate(monkeypatch):
    import nwkit.mul_locus_mc as module

    monkeypatch.setattr(module, "search_banks", toy_search)

    def sampler(pop, point, config, rng, *, count):
        if point.ne == 2:
            raise ValueError("selection cap reached")
        return [0], 1

    with pytest.raises(
        ValueError, match="generating grid 7 replicate 1.*selection cap"
    ):
        calibrate(
            toy_banks(),
            [1],
            model(),
            0.001,
            2,
            sampler=sampler,
            null_calibration="grid-supremum",
        )


def test_cli_grid_calibration_full_bundle_and_process_equality(tmp_path, capsys):
    args, paths = arguments(tmp_path)
    config = model(400)
    config["parameter_grid"].append({**config["parameter_grid"][0], "ne": 1})
    (tmp_path / "config.json").write_text(json.dumps(config))
    args += ["--locus-null-calibration", "grid-supremum"]
    main(args)
    first = [p.read_bytes() for p in paths.values()]
    main([*args, "--cpus", "2"])
    assert [p.read_bytes() for p in paths.values()] == first
    saved = json.loads(paths["model"].read_text())["calibration"]
    assert saved["method"] == "finite-grid-supremum-Monte-Carlo-test"
    assert saved["total_replicates"] == 4
    rows = pd.read_csv(paths["calibration"], sep="\t")
    assert len(rows) == 4
    assert set(rows["generating_null_grid"]) == {0, 1}
    assert "Finite-grid supremum MC" in capsys.readouterr().err


@pytest.mark.parametrize("score", ["dl", "msc"])
def test_other_models_reject_null_calibration(tmp_path, score):
    args, _ = arguments(tmp_path)
    args[args.index("locus-mc")] = score
    for flag in ("--locus-model", "--locus-bootstrap", "--locus-calibration-out"):
        position = args.index(flag)
        del args[position : position + 2]
    with pytest.raises(ValueError, match="Locus options"):
        main([*args, "--locus-null-calibration", "grid-supremum"])


def test_cli_grid_failure_preserves_prior_bundle(tmp_path, monkeypatch):
    import nwkit.mul_locus_cli as module

    args, paths = arguments(tmp_path)
    for path in paths.values():
        path.write_text("prior result")

    def fail(*a, **kw):
        raise ValueError("Null bootstrap generating grid 0 replicate 1: injected")

    monkeypatch.setattr(module, "calibrate", fail)
    with pytest.raises(ValueError, match="generating grid 0"):
        main([*args, "--locus-null-calibration", "grid-supremum"])
    assert all(p.read_text() == "prior result" for p in paths.values())
    assert not list(tmp_path.glob(".*.stage.*"))

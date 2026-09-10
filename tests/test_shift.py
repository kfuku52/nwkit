"""Shift adapter tests; the optional real R boundary is tested when configured."""

import json
import os

import pandas as pd
import pytest

from nwkit.cli import main


def fake_backend(directory, args):
    """A fixed nested-shift result in backend clade coordinates."""
    from nwkit.util import read_tree

    tree = read_tree(str(directory / "tree.nwk"), "auto", True, quiet=True)
    # t0/t1 is not always a root child; use its incoming-branch start time.
    clade = next(n for n in tree.traverse() if set(n.leaf_names()) == {"t0", "t1"})
    import math

    delta = 0.4 / -math.expm1(
        -0.5
        * tree.get_distance(clade.up, next(n for n in tree.leaves() if n.name == "t0"))
    )
    leaf_delta = -0.1 / -math.expm1(
        -0.5 * next(n for n in tree.leaves() if n.name == "t1").dist
    )
    (directory / "shifts.tsv").write_text(
        f"clade\tmean_effect\toptimum_effect\nt0/t1\t0.4\t{delta}\nt1\t-0.1\t{leaf_delta}\n"
    )
    data = pd.read_csv(directory / "trait.tsv", sep="\t")
    rows = []
    for row in data.itertuples():
        mean = (
            2
            + (0.4 if row.leaf_name in {"t0", "t1"} else 0)
            - (0.1 if row.leaf_name == "t1" else 0)
        )
        optimum = (
            2
            + (delta if row.leaf_name in {"t0", "t1"} else 0)
            + (leaf_delta if row.leaf_name == "t1" else 0)
        )
        rows.append([row.leaf_name, row.value, mean, row.value - mean, optimum])
    pd.DataFrame(
        rows, columns=["token", "observed", "predicted", "residual", "optimum"]
    ).to_csv(directory / "tips.tsv", sep="\t", index=False)
    (directory / "search.tsv").write_text(
        "strategy\tconfiguration_space_size\tevaluated_configurations\tcoverage\tglobally_optimal\tensemble_attempted\tensemble_successful\tensemble_failed\talpha_lower\talpha_upper\nexhaustive\t2\t2\t1\tTRUE\t0\t0\t0\t0\t10\n"
    )
    (directory / "model.tsv").write_text(
        "backend_version\talpha\tsigma2\tintercept\tlog_likelihood\tscore\n"
        "3.0.9\t0.5\t1.0\t2.0\t-5\t20\n"
    )
    (directory / "candidates.tsv").write_text("score\tclades\n20\tt0/t1;t1\n25\t\n")
    (directory / "diagnostics.txt").write_text('list(strategy="exhaustive")')
    (directory / "fit.rds").write_bytes(b"test fit")
    return {"stdout": "", "stderr": "", "rscript": "fake"}


@pytest.fixture
def shift_inputs(tmp_path):
    tree = tmp_path / "tree.nwk"
    tree.write_text("((A:1,B:1):1,(C:1,D:1):1);")
    trait = tmp_path / "trait.tsv"
    trait.write_text("leaf_name\tx\nD\t4\nB\t2\nA\t1\nC\t3\n")
    return [
        "shift",
        "--selection",
        "ic",
        "-i",
        str(tree),
        "--trait",
        str(trait),
        "--state-column",
        "x",
        "--model-out",
        str(tmp_path / "model.json"),
    ]


def test_nested_shifts_and_tip_alignment(shift_inputs, tmp_path, monkeypatch):
    def backend(directory, args):
        data = pd.read_csv(directory / "trait.tsv", sep="\t")
        assert data.value.tolist() == [1, 2, 3, 4]
        assert ":1" in (directory / "tree.nwk").read_text()
        return fake_backend(directory, args)

    monkeypatch.setattr("nwkit.shift.run_backend", backend)
    output = tmp_path / "regimes.tsv"
    main([*shift_inputs, "-o", str(output), "--fit-out", str(tmp_path / "fit.rds")])
    rows = pd.read_csv(output, sep="\t").set_index("branch_id").regime.to_dict()
    assert rows == {
        0: "baseline",
        1: "shift_1",
        2: "baseline",
        3: "shift_1",
        4: "shift_4",
        5: "baseline",
        6: "baseline",
    }
    model = json.loads((tmp_path / "model.json").read_text())
    assert model["shift_branch_ids"] == [1, 4]
    assert model["tip_tokens"] == {"t0": "A", "t1": "B", "t2": "C", "t3": "D"}
    assert model["root_model"] == "OUfixedRoot"
    assert not model["tree_normalized"]
    assert model["candidates"][1]["shift_branch_ids"] == []


@pytest.mark.parametrize(
    "tree,match",
    [
        ("((A:1,B:2):1,(C:1,D:1):1);", "ultrametric"),
        ("((A:0,B:1):1,(C:1,D:1):1);", "positive"),
        ("(A:2,B:2,C:2,D:2);", "bifurcating"),
        ("((A:1,A:1):1,(C:1,D:1):1);", "Duplicated"),
    ],
)
def test_invalid_trees_fail_before_backend(shift_inputs, tmp_path, tree, match):
    (tmp_path / "tree.nwk").write_text(tree)
    with pytest.raises(ValueError, match=match):
        main([*shift_inputs, "--input-rooted", "yes"])
    assert not (tmp_path / "model.json").exists()


@pytest.mark.parametrize(
    "values,match",
    [
        ("A\t1\nB\tNA\nC\t3\nD\t4\n", "finite"),
        ("A\t1\nB\t1\nC\t1\nD\t1\n", "invariant"),
        ("A\t1\nB\t2\nC\t3\n", "tips differ"),
    ],
)
def test_invalid_traits(shift_inputs, tmp_path, values, match):
    (tmp_path / "trait.tsv").write_text("leaf_name\tx\n" + values)
    with pytest.raises(ValueError, match=match):
        main(shift_inputs)


def test_failure_preserves_outputs(shift_inputs, tmp_path, monkeypatch):
    def invalid(directory, args):
        result = fake_backend(directory, args)
        (directory / "shifts.tsv").write_text(
            "clade\tmean_effect\toptimum_effect\nt0/t2\t0.4\t1\n"
        )
        return result

    monkeypatch.setattr("nwkit.shift.run_backend", invalid)
    model = tmp_path / "model.json"
    model.write_text("previous")
    with pytest.raises(ValueError, match="invalid shift clade"):
        main(shift_inputs)
    assert model.read_text() == "previous"


def test_protect_inputs(shift_inputs, tmp_path):
    with pytest.raises(ValueError, match="input"):
        main([*shift_inputs, "-o", str(tmp_path / "tree.nwk")])


def test_missing_rscript(shift_inputs):
    with pytest.raises(ValueError, match="Rscript was not found"):
        main([*shift_inputs, "--rscript", "/does-not-exist/Rscript"])


def test_stdout_is_table(shift_inputs, monkeypatch, capsys):
    monkeypatch.setattr("nwkit.shift.run_backend", fake_backend)
    main(shift_inputs)
    assert capsys.readouterr().out.startswith("branch_id\tregime\n")


@pytest.mark.integration
@pytest.mark.skipif(
    not os.environ.get("NWKIT_TEST_RSCRIPT"),
    reason="Set NWKIT_TEST_RSCRIPT to Rscript with kfl1ou >= 3.0.9",
)
@pytest.mark.parametrize("strategy", ["exhaustive", "lasso", "ensemble", "auto"])
def test_real_kfl1ou(tmp_path, strategy):
    tree = tmp_path / "tree.nwk"
    tree.write_text("(((A:1,B:1):1,(C:1,D:1):1):1,((E:1,F:1):1,(G:1,H:1):1):1);")
    trait = tmp_path / "trait.tsv"
    trait.write_text(
        "leaf_name\tx\n"
        + "".join(
            f"{n}\t{x}\n"
            for n, x in zip(
                "ABCDEFGH", [1, 1.2, 0.9, 1.4, 6, 6.2, 5.8, 6.1], strict=True
            )
        )
    )
    output = tmp_path / "regimes.tsv"
    model = tmp_path / "model.json"
    main(
        [
            "shift",
            "--selection",
            "ic",
            "-i",
            str(tree),
            "--trait",
            str(trait),
            "--state-column",
            "x",
            "--model-out",
            str(model),
            "-o",
            str(output),
            "--rscript",
            os.environ["NWKIT_TEST_RSCRIPT"],
            "--max-shifts",
            "1",
            "--criterion",
            "BIC",
            "--search-strategy",
            strategy,
        ]
    )
    saved = json.loads(model.read_text())
    assert saved["backend_version"]
    assert len(pd.read_csv(output, sep="\t")) == 15
    assert len(saved["shift_branch_ids"]) == 1
    # ASR must accept the complete per-branch map (a new stationary-root refit).
    main(
        [
            "asr",
            "-i",
            str(tree),
            "--trait",
            str(trait),
            "--state-column",
            "x",
            "--trait-type",
            "continuous",
            "--model",
            "OUM",
            "--regime-map",
            str(output),
            "-o",
            str(tmp_path / "asr.tsv"),
        ]
    )
    assert (tmp_path / "asr.tsv").is_file()


def test_zero_shifts_and_audit(shift_inputs, tmp_path, monkeypatch):
    def zero(directory, args):
        result = fake_backend(directory, args)
        (directory / "shifts.tsv").write_text("clade\tmean_effect\toptimum_effect\n")
        tips = pd.read_csv(directory / "tips.tsv", sep="\t")
        tips["predicted"] = 2.0
        tips["residual"] = tips["observed"] - 2.0
        tips["optimum"] = 2.0
        tips.to_csv(directory / "tips.tsv", sep="\t", index=False)
        (directory / "candidates.tsv").write_text("score\tclades\n20\t\n")
        return result

    monkeypatch.setattr("nwkit.shift.run_backend", zero)
    audit = tmp_path / "audit.jsonl"
    fit = tmp_path / "fit.rds"
    output = tmp_path / "regimes.tsv"
    main(
        [
            *shift_inputs,
            "--max-shifts",
            "0",
            "--audit",
            str(audit),
            "--fit-out",
            str(fit),
            "-o",
            str(output),
        ]
    )
    assert set(pd.read_csv(output, sep="\t").regime) == {"baseline"}
    record = json.loads(audit.read_text())
    assert record["status"] == "ok"
    assert any(row["path"] == str(fit) for row in record["outputs"])


def test_backend_failure(shift_inputs, tmp_path, monkeypatch):
    from types import SimpleNamespace

    monkeypatch.setattr("nwkit.shift_backend.shutil.which", lambda _: "/mock/Rscript")
    monkeypatch.setattr(
        "nwkit.shift_backend.subprocess.run",
        lambda *a, **kw: SimpleNamespace(
            returncode=1, stderr="missing package", stdout=""
        ),
    )
    with pytest.raises(RuntimeError, match="missing package"):
        main(shift_inputs)
    assert not (tmp_path / "model.json").exists()


def test_stdout_input_ownership(shift_inputs):
    with pytest.raises(ValueError, match="only one input"):
        main([*shift_inputs, "-i", "-", "--trait", "-"])


def test_effects_predictions_and_search_outputs(shift_inputs, tmp_path, monkeypatch):
    monkeypatch.setattr("nwkit.shift.run_backend", fake_backend)
    paths = [tmp_path / name for name in ("effects.tsv", "parameters.tsv", "tips.tsv")]
    main(
        [
            *shift_inputs,
            "--effects-out",
            str(paths[0]),
            "--regime-parameters-out",
            str(paths[1]),
            "--tip-summary-out",
            str(paths[2]),
        ]
    )
    effects, regimes, tips = [pd.read_csv(p, sep="\t") for p in paths]
    assert effects.mean_effect.tolist() == [0.4, -0.1]
    assert tips.leaf_name.tolist() == ["A", "B", "C", "D"]
    assert tips.predicted.tolist() == pytest.approx([2.4, 2.3, 2, 2])
    assert tips.residual.tolist() == pytest.approx([-1.4, -0.3, 1, 2])
    assert regimes.optimum.tolist()[0] == 2
    model = json.loads((tmp_path / "model.json").read_text())
    assert model["schema_version"] == 5
    assert model["search"]["coverage"] == 1
    assert model["search"]["globally_optimal"] is True
    assert model["search"]["alpha_at_upper_bound"] is False
    assert model["tip_predictions"][1]["optimum"] == pytest.approx(
        regimes.optimum.iloc[2]
    )


@pytest.mark.parametrize("corruption", ["tip", "effect", "residual", "optimum"])
def test_inconsistent_backend_results_preserve_outputs(
    shift_inputs, tmp_path, monkeypatch, corruption
):
    def corrupt(directory, args):
        result = fake_backend(directory, args)
        path = directory / ("shifts.tsv" if corruption == "effect" else "tips.tsv")
        frame = pd.read_csv(path, sep="\t")
        column = {
            "tip": "token",
            "effect": "optimum_effect",
            "residual": "residual",
            "optimum": "optimum",
        }[corruption]
        frame.loc[0, column] = "unknown" if corruption == "tip" else 100.0
        frame.to_csv(path, sep="\t", index=False)
        return result

    monkeypatch.setattr("nwkit.shift.run_backend", corrupt)
    output = tmp_path / "effects.tsv"
    output.write_text("previous")
    with pytest.raises(ValueError):
        main([*shift_inputs, "--effects-out", str(output)])
    assert output.read_text() == "previous"
    assert not (tmp_path / "model.json").exists()


def test_bm_optima_are_missing(shift_inputs, tmp_path, monkeypatch):
    def boundary(directory, args):
        result = fake_backend(directory, args)
        model = pd.read_csv(directory / "model.tsv", sep="\t")
        model["alpha"] = 0.0
        model.to_csv(directory / "model.tsv", sep="\t", index=False)
        for name, column in [("shifts.tsv", "optimum_effect"), ("tips.tsv", "optimum")]:
            table = pd.read_csv(directory / name, sep="\t")
            table[column] = "NA"
            table.to_csv(directory / name, sep="\t", index=False)
        return result

    monkeypatch.setattr("nwkit.shift.run_backend", boundary)
    main(shift_inputs)
    model = json.loads((tmp_path / "model.json").read_text())
    assert all(
        row["optimum"] is None and not row["optimum_identifiable"]
        for row in model["regime_parameters"]
    )
    assert model["shift_effects"][0]["mean_effect"] == 0.4
    assert model["search"]["alpha_at_lower_bound"]

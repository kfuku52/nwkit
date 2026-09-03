import numpy as np
import pandas as pd
import pytest

from nwkit.cli import main
from nwkit.discrete_asr_models import (
    build_rate_matrix,
    model_equivalence_family,
    read_rate_design,
)


def test_rate_design_builds_shared_direct_rates(tmp_path):
    path = tmp_path / "design.tsv"
    path.write_text(
        "from_state\tto_state\trate_class\n"
        "low\tmid\tup\n"
        "mid\thigh\tup\n"
        "high\tmid\tdown\n"
        "mid\tlow\tdown\n"
    )
    design = read_rate_design(path, ("low", "mid", "high"))
    assert design.class_names == ("up", "down")
    matrix = build_rate_matrix(
        "MK-DESIGN",
        ("low", "mid", "high"),
        (0.2, 0.7),
        rate_design=design,
    )
    np.testing.assert_allclose(
        matrix,
        [[-0.2, 0.2, 0.0], [0.7, -0.9, 0.2], [0.0, 0.7, -0.7]],
    )


def test_rate_design_rejects_duplicate_edges_and_unknown_states(tmp_path):
    duplicate = tmp_path / "duplicate.tsv"
    duplicate.write_text("from_state\tto_state\trate_class\na\tb\tx\na\tb\ty\n")
    with pytest.raises(ValueError, match="duplicated edge"):
        read_rate_design(duplicate, ("a", "b"))

    unknown = tmp_path / "unknown.tsv"
    unknown.write_text("from_state\tto_state\trate_class\na\tc\tx\n")
    with pytest.raises(ValueError, match="absent from the model state space"):
        read_rate_design(unknown, ("a", "b"))


@pytest.mark.integration
def test_mk_design_cli_reports_fitted_rate_classes(tmp_path):
    trait = tmp_path / "traits.tsv"
    design = tmp_path / "design.tsv"
    output = tmp_path / "asr.tsv"
    model = tmp_path / "model.tsv"
    trait.write_text(
        "leaf_name\tstate\nA\tlow\nB\tmid\nC\thigh\nD\tlow\nE\tmid\nF\thigh\n"
    )
    design.write_text(
        "from_state\tto_state\trate_class\n"
        "low\tmid\tup\nmid\thigh\tup\n"
        "high\tmid\tdown\nmid\tlow\tdown\n"
    )
    main(
        [
            "asr",
            "-i",
            "[&R]((A:.2,B:.4):.3,(C:.7,D:1.1):.5,E:.6,F:1.3)R;",
            "--trait",
            str(trait),
            "--state-column",
            "state",
            "--model",
            "MK-DESIGN",
            "--rate-design",
            str(design),
            "--model-out",
            str(model),
            "-o",
            str(output),
        ]
    )
    result = pd.read_csv(output, sep="\t")
    fit = pd.read_csv(model, sep="\t").iloc[0]
    assert np.allclose(result[["p_low", "p_mid", "p_high"]].sum(axis=1), 1.0)
    assert fit["model"] == "MK-DESIGN"
    assert fit["rate_classes"] == "up,down"
    assert fit["rate_class_up"] > 0.0
    assert fit["rate_class_down"] > 0.0


@pytest.mark.integration
def test_asrcompare_excludes_rate_design_equivalent_to_er(tmp_path):
    trait = tmp_path / "traits.tsv"
    design_path = tmp_path / "design.tsv"
    comparison = tmp_path / "comparison.tsv"
    trait.write_text("leaf_name\tstate\nA\ta\nB\ta\nC\tb\nD\tb\nE\ta\n")
    design_path.write_text("from_state\tto_state\trate_class\na\tb\tq\nb\ta\tq\n")
    design = read_rate_design(design_path, ("a", "b"))
    assert model_equivalence_family(
        "MK-DESIGN", ("a", "b"), rate_design=design
    ) == model_equivalence_family("ER", ("a", "b"))
    main(
        [
            "asrcompare",
            "-i",
            "[&R](((A:1,B:1):1,C:2):1,(D:1,E:1):2)R;",
            "--trait",
            str(trait),
            "--state-column",
            "state",
            "--models",
            "ER,MK-DESIGN",
            "--rate-design",
            str(design_path),
            "-o",
            str(comparison),
        ]
    )
    table = pd.read_csv(comparison, sep="\t").set_index("model")
    assert table.loc["MK-DESIGN", "status"] == "equivalent"
    assert table.loc["MK-DESIGN", "equivalent_to"] == "ER"

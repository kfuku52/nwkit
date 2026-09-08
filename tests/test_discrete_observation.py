from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from nwkit.cli import main
from nwkit.discrete_observation import apply_discrete_observation_model


def test_misclassification_orientation_and_missing(tmp_path):
    path = tmp_path / "matrix.tsv"
    path.write_text("state\t0\t1\n0\t0.9\t0.1\n1\t0.2\t0.8\n")
    actual = apply_discrete_observation_model(
        ["0", "1"],
        {"a": [1, 0], "b": [0, 1], "c": [1, 1]},
        SimpleNamespace(misclassification_matrix=path),
    )
    assert actual["a"] == pytest.approx([0.9, 0.2])
    assert actual["b"] == pytest.approx([0.1, 0.8])
    assert actual["c"] == pytest.approx([1, 1])


@pytest.mark.parametrize(
    "row", ["a\t0\t0", "a\tNaN\t1", "a\t-0.1\t1", "a\t0.2\t1.2", "unknown\t0.2\t0.8"]
)
def test_invalid_tip_likelihoods_rejected(tmp_path, row):
    path = tmp_path / "likelihoods.tsv"
    path.write_text("leaf_name\t0\t1\n" + row + "\n")
    with pytest.raises(ValueError):
        apply_discrete_observation_model(
            ["0", "1"], {"a": [1, 0]}, SimpleNamespace(tip_likelihoods=path)
        )


@pytest.mark.integration
def test_discrete_observation_cli_and_comparison_consistent(tmp_path):
    traits = tmp_path / "traits.tsv"
    traits.write_text("leaf_name\tstate\nA\t0\nB\t1\nC\t0\n")
    likelihoods = tmp_path / "likelihoods.tsv"
    likelihoods.write_text("leaf_name\t0\t1\nA\t0.6\t0.4\nB\t0.3\t0.7\nC\t0.8\t0.2\n")
    options = [
        "-i",
        "(A:1,B:1,C:1)R;",
        "--input-rooted",
        "yes",
        "--trait",
        str(traits),
        "--trait-type",
        "discrete",
        "--state-column",
        "state",
        "--states",
        "0,1",
        "--tip-likelihoods",
        str(likelihoods),
        "--rate",
        "0.4",
    ]
    main(
        [
            "asr",
            *options,
            "--model",
            "ER",
            "--model-out",
            str(tmp_path / "model.tsv"),
            "-o",
            str(tmp_path / "asr.tsv"),
        ]
    )
    main(
        ["asrcompare", *options, "--models", "ER", "-o", str(tmp_path / "compare.tsv")]
    )
    model = pd.read_csv(tmp_path / "model.tsv", sep="\t")
    comparison = pd.read_csv(tmp_path / "compare.tsv", sep="\t")
    assert comparison.log_likelihood.iloc[0] == pytest.approx(
        model.log_likelihood.iloc[0]
    )
    output = pd.read_csv(tmp_path / "asr.tsv", sep="\t")
    probabilities = output.filter(regex="^p_").to_numpy()
    assert probabilities.size
    assert np.all((probabilities > 0) & (probabilities < 1))

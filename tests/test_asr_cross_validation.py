import numpy as np
import pandas as pd
import pytest
from scipy.stats import norm

from nwkit.asr_cross_validation import gaussian_cross_validation, holdout_groups
from nwkit.cli import main
from nwkit.evolution import build_evolutionary_process
from nwkit.util import read_tree


def fixture():
    tree = read_tree("((A:1,B:1):1,(C:1,D:1):1)R;", "1", True, quiet=True, rooted="yes")
    values = {"A": 0.0, "B": 1.0, "C": 3.0, "D": 4.0}
    return tree, values


def test_holdout_clades_partition_observed_tips():
    tree, values = fixture()
    assert holdout_groups(tree, values, "clade") == [("A", "B"), ("C", "D")]
    values["A"] = None
    assert holdout_groups(tree, values, "clade") == [("B",), ("C", "D")]
    with pytest.raises(ValueError, match="mode"):
        holdout_groups(tree, values, "wrong")
    with pytest.raises(ValueError, match="two nonempty"):
        holdout_groups(tree, {"A": 1.0})


@pytest.mark.parametrize("mode", ["tip", "clade"])
def test_cross_validation_matches_dense_conditioning_without_leakage(mode):
    tree, values = fixture()
    process = build_evolutionary_process(tree, model="brownian", root_mode="flat")
    training_sets = []

    def refit(training):
        training_sets.append(dict(training))
        return process

    errors = dict.fromkeys(values, 0.2)
    result = gaussian_cross_validation(tree, values, refit, errors=errors, mode=mode)
    covariance = (
        np.array([[2, 1, 0, 0], [1, 2, 0, 0], [0, 0, 2, 1], [0, 0, 1, 2]])
        + np.eye(4) * 0.04
    )
    names = list(values)
    for row in result.itertuples():
        training = training_sets[row.fold]
        assert training[row.name] is None
        indices = [i for i, name in enumerate(names) if training[name] is not None]
        inverse = np.linalg.inv(covariance[np.ix_(indices, indices)])
        ones = np.ones(len(indices))
        data = np.array([values[names[i]] for i in indices])
        root_var = 1 / (ones @ inverse @ ones)
        root_mean = root_var * (ones @ inverse @ data)
        cross = covariance[names.index(row.name), indices]
        mean = root_mean + cross @ inverse @ (data - root_mean)
        variance = (
            covariance[names.index(row.name), names.index(row.name)]
            - cross @ inverse @ cross
            + root_var * (1 - cross @ inverse @ ones) ** 2
        )
        assert row.predicted_mean == pytest.approx(mean)
        assert row.predicted_sd**2 == pytest.approx(variance)
        assert row.log_score == pytest.approx(
            norm.logpdf(values[row.name], mean, np.sqrt(variance))
        )
    assert values == {"A": 0.0, "B": 1.0, "C": 3.0, "D": 4.0}


@pytest.mark.integration
def test_cross_validation_cli_refits(tmp_path):
    traits = tmp_path / "traits.tsv"
    traits.write_text("leaf_name\tvalue\nA\t0\nB\t1\nC\t3\nD\t4\n")
    cv = tmp_path / "cv.tsv"
    main(
        [
            "asr",
            "-i",
            "((A:1,B:1):1,(C:1,D:1):1)R;",
            "--input-rooted",
            "yes",
            "--trait",
            str(traits),
            "--state-column",
            "value",
            "--cross-validation-out",
            str(cv),
            "-o",
            str(tmp_path / "asr.tsv"),
        ]
    )
    result = pd.read_csv(cv, sep="\t")
    assert len(result) == 4
    assert set(result.num_training) == {3}
    assert np.isfinite(result.log_score).all()

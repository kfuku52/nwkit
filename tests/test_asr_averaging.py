from types import SimpleNamespace

import numpy as np
import pytest
from scipy.stats import norm

from nwkit.asr_averaging import (
    average_state_probabilities,
    compatible_model_weights,
    gaussian_mixture_summary,
)
from nwkit.asr_tree_ensemble import (
    align_ensemble_nodes,
    summarize_tree_ensemble,
    summarize_vector_tree_ensemble,
)
from nwkit.util import read_tree


@pytest.mark.integration
@pytest.mark.parametrize(
    "trait_type,model,values,extra",
    [
        ("continuous", "BM", [1, 2, 4, 5], ["--sigma2", "1"]),
        ("discrete", "ER", ["a", "a", "b", "b"], ["--rate", "0.4"]),
    ],
)
def test_tree_ensemble_cli(tmp_path, trait_type, model, values, extra):
    import pandas as pd

    from nwkit.cli import main

    reference = "((A:1,B:1)AB:1,(C:1,D:1)CD:1)R;"
    sample = tmp_path / "trees.nwk"
    sample.write_text(reference + "\n((A:1,C:1)AC:1,(B:1,D:1)BD:1)R;\n")
    traits = tmp_path / "traits.tsv"
    traits.write_text(
        "leaf_name\tstate\n"
        + "".join(
            f"{name}\t{value}\n" for name, value in zip("ABCD", values, strict=True)
        )
    )
    output = tmp_path / "ensemble.tsv"
    main(
        [
            "asr",
            "-i",
            reference,
            "--input-rooted",
            "yes",
            "--trait",
            str(traits),
            "--trait-type",
            trait_type,
            "--state-column",
            "state",
            "--model",
            model,
            *extra,
            "--tree-ensemble",
            str(sample),
            "--tree-ensemble-out",
            str(output),
            "-o",
            str(tmp_path / "reference.tsv"),
        ]
    )
    table = pd.read_csv(output, sep="\t")
    assert table.loc[table.name == "AB", "matched_tree_weight"].iloc[0] == 0.5
    assert table.loc[table.name == "R", "matched_tree_weight"].iloc[0] == 1
    if trait_type == "discrete":
        assert table.filter(regex="^p_").sum(axis=1).to_numpy() == pytest.approx(
            np.ones(len(table))
        )


def tree_from(text):
    return read_tree(text, "1", True, quiet=True, rooted="yes")


def test_vector_tree_ensemble_marginals():
    tree = tree_from("(A:1,B:1)R;")
    posterior = {
        node: SimpleNamespace(
            mean=np.array([1, 2]), covariance=np.array([[3, 1], [1, 4]])
        )
        for node in tree.traverse()
    }
    table = summarize_vector_tree_ensemble(
        tree, [(tree, posterior)], trait_names=("x", "y")
    )
    assert set(table.trait) == {"x", "y"}
    assert table.loc[table.trait == "x", "variance"].to_numpy() == pytest.approx(
        [3, 3, 3]
    )
    assert table.loc[table.trait == "y", "variance"].to_numpy() == pytest.approx(
        [4, 4, 4]
    )


@pytest.mark.integration
def test_model_average_cli(tmp_path):
    import pandas as pd

    from nwkit.cli import main

    traits = tmp_path / "traits.tsv"
    traits.write_text("leaf_name\tstate\nA\t0\nB\t1\nC\t2\nD\t0\nE\t1\nF\t2\n")
    output = tmp_path / "average.tsv"
    main(
        [
            "asrcompare",
            "-i",
            "((A:1,B:2):1,(C:1,D:2):1,(E:2,F:1):1);",
            "--input-rooted",
            "yes",
            "--trait",
            str(traits),
            "--state-column",
            "state",
            "--trait-type",
            "discrete",
            "--models",
            "ER,SYM",
            "--model-average-out",
            str(output),
            "-o",
            str(tmp_path / "comparison.tsv"),
        ]
    )
    table = pd.read_csv(output, sep="\t")
    assert table.filter(regex="^p_").sum(axis=1).to_numpy() == pytest.approx(
        np.ones(len(table))
    )
    assert table.num_models.min() >= 1


@pytest.mark.integration
def test_continuous_model_average_cli(tmp_path):
    import pandas as pd

    from nwkit.cli import main

    traits = tmp_path / "traits.tsv"
    traits.write_text("leaf_name\tstate\nA\t0\nB\t1\nC\t3\nD\t2\nE\t5\nF\t4\n")
    output = tmp_path / "average.tsv"
    main(
        [
            "asrcompare",
            "-i",
            "((A:1,B:2):1,(C:1,D:2):1,(E:2,F:1):1);",
            "--input-rooted",
            "yes",
            "--trait",
            str(traits),
            "--state-column",
            "state",
            "--trait-type",
            "continuous",
            "--models",
            "BM,LAMBDA",
            "--evolution-parameter",
            "0.5",
            "--model-average-out",
            str(output),
            "-o",
            str(tmp_path / "comparison.tsv"),
        ]
    )
    table = pd.read_csv(output, sep="\t")
    assert table.num_models.min() == 2
    assert table.variance.to_numpy() == pytest.approx(
        (table.within_model_variance + table.between_model_variance).to_numpy()
    )


def test_gaussian_mixture_uses_total_variance_and_mixture_quantiles():
    result = gaussian_mixture_summary([-4, 4], [1, 1], [1, 1])
    assert result.mean == 0
    assert result.variance == 17
    assert result.within_variance == 1
    assert result.between_variance == 16
    assert (
        norm.cdf(result.lower + 4) + norm.cdf(result.lower - 4)
    ) / 2 == pytest.approx(0.025)
    assert result.upper == pytest.approx(-result.lower)


def test_gaussian_mixture_handles_point_masses():
    result = gaussian_mixture_summary([0, 10], [0, 0], [0.99, 0.01])
    assert result.lower == pytest.approx(0, abs=1e-50)
    assert result.upper == pytest.approx(0, abs=1e-50)
    assert result.mean == pytest.approx(0.1)
    assert result.variance == pytest.approx(0.99)


def test_gaussian_mixture_extreme_confidence_uses_stable_survival_tail():
    level = np.nextafter(1.0, 0.0)
    result = gaussian_mixture_summary([0], [1], [1], level=level)
    expected = norm.isf((1 - level) / 2)
    assert result.upper == pytest.approx(expected)
    assert result.lower == pytest.approx(-expected)


def test_compatible_weights_reject_different_root_groups():
    assert compatible_model_weights([10, 12], ["BM:flat", "BM:flat"]) == pytest.approx(
        [1 / (1 + np.exp(-1)), 1 / (1 + np.exp(1))]
    )
    with pytest.raises(ValueError, match="incompatible"):
        compatible_model_weights([10, 12], ["BM:flat", "OU:stationary"])
    assert average_state_probabilities(
        [[0.2, 0.8], [0.6, 0.4]], [1, 3]
    ) == pytest.approx([0.5, 0.5])


def test_tree_mapping_uses_clades_not_branch_order():
    reference = tree_from("((A:1,B:1)AB:1,(C:1,D:1)CD:1)R;")
    reordered = tree_from("((D:1,C:1)otherCD:1,(B:1,A:1)otherAB:1)R;")
    different = tree_from("((A:1,C:1)AC:1,(B:1,D:1)BD:1)R;")
    node = next(node for node in reference.traverse() if node.name == "AB")
    assert align_ensemble_nodes(reference, reordered)[node].name == "otherAB"
    assert node not in align_ensemble_nodes(reference, different)
    assert align_ensemble_nodes(reference, different, mapping="mrca")[node] is different
    with pytest.raises(ValueError, match="exactly the reference tip set"):
        align_ensemble_nodes(reference, tree_from("(A:1,B:1)R;"))


def test_tree_ensemble_reports_clade_support_and_between_tree_uncertainty():
    reference = tree_from("((A:1,B:1)AB:1,(C:1,D:1)CD:1)R;")
    different = tree_from("((A:1,C:1)AC:1,(B:1,D:1)BD:1)R;")
    posteriors = [
        (
            tree,
            {node: SimpleNamespace(mean=mean, variance=1) for node in tree.traverse()},
        )
        for tree, mean in [(reference, 0), (different, 4)]
    ]
    result = summarize_tree_ensemble(reference, posteriors).set_index("name")
    assert result.loc["AB", "matched_tree_weight"] == 0.5
    assert result.loc["AB", "mean"] == 0
    assert result.loc["R", "mean"] == 2
    assert result.loc["R", "variance"] == 5
    mrca = summarize_tree_ensemble(reference, posteriors, mapping="mrca").set_index(
        "name"
    )
    assert mrca.loc["AB", "matched_tree_weight"] == 1
    assert mrca.loc["AB", "mean"] == 2


def test_multivariate_model_average_cli(tmp_path):
    import pandas as pd

    from nwkit.cli import main

    traits = tmp_path / "traits.tsv"
    traits.write_text("leaf_name\tx\ty\nA\t0\t1\nB\t2\t0\nC\t1\t3\nD\t4\t2\n")
    output = tmp_path / "average.tsv"
    main(
        [
            "asrcompare",
            "-i",
            "(A:1,B:2,C:3,D:4)R;",
            "--input-rooted",
            "yes",
            "--trait",
            str(traits),
            "--state-column",
            "x,y",
            "--models",
            "MV-BM",
            "--model-average-out",
            str(output),
            "-o",
            str(tmp_path / "models.tsv"),
        ]
    )
    result = pd.read_csv(output, sep="\t")
    assert set(result.trait) == {"x", "y"}
    assert len(result) == 10
    assert result.num_models.min() == 1
    assert result.between_model_variance.max() == 0

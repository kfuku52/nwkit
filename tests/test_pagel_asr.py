import numpy as np
import pandas as pd
import pytest

from nwkit.asr import _is_missing_observation
from nwkit.cli import main
from nwkit.pagel_asr import (
    compute_pagel_marginals,
    pagel_rate_design,
    prepare_pagel_data,
)
from nwkit.util import read_tree

TREE = "[&R](((A:.2,B:.3):.4,(C:.2,D:.5):.3):.6,((E:.2,F:.4):.5,(G:.3,H:.7):.2):.9)R;"


def _trait_file(tmp_path):
    path = tmp_path / "traits.tsv"
    path.write_text(
        "leaf_name\tfirst\tsecond\n"
        "A\t0\tx\nB\t0\ty\nC\t1\tx\nD\t1\ty\n"
        "E\t0\tx\nF\t1\ty\nG\t1\tx\nH\t0\ty\n"
    )
    return path


def test_pagel_designs_have_four_and_eight_rates_without_double_changes(tmp_path):
    trait = _trait_file(tmp_path)
    tree = read_tree(TREE, "1", True, quiet=True, rooted="yes")
    data = prepare_pagel_data(trait, ("first", "second"), list(tree.leaf_names()))
    independent = pagel_rate_design(data, "PAGEL-INDEPENDENT")
    dependent = pagel_rate_design(data, "PAGEL-DEPENDENT")
    assert len(independent.class_names) == 4
    assert len(dependent.class_names) == 8
    assert independent.graph.sum() == dependent.graph.sum() == 8
    for source in range(4):
        for target in range(4):
            differing_bits = (source // 2 != target // 2) + (source % 2 != target % 2)
            if differing_bits == 2:
                assert not independent.graph[source, target]
                assert not dependent.graph[source, target]


def test_pagel_dependent_nests_independent_likelihood(tmp_path):
    trait = _trait_file(tmp_path)
    tree = read_tree(TREE, "1", True, quiet=True, rooted="yes")
    data = prepare_pagel_data(trait, ("first", "second"), list(tree.leaf_names()))
    _, independent = compute_pagel_marginals(tree, data, model="PAGEL-INDEPENDENT")
    _, dependent = compute_pagel_marginals(tree, data, model="PAGEL-DEPENDENT")
    assert dependent["log_likelihood"] >= independent["log_likelihood"] - 1e-7
    assert len(independent["rates"]) == 4
    assert len(dependent["rates"]) == 8
    assert independent["likelihood_kind"] == dependent["likelihood_kind"]
    assert independent["optimizer_starts"] == dependent["optimizer_starts"] == 6


def test_pagel_partial_missing_trait_retains_partial_likelihood(tmp_path):
    trait = _trait_file(tmp_path)
    frame = pd.read_csv(trait, sep="\t", dtype=str)
    frame.loc[frame["leaf_name"] == "A", "second"] = "NA"
    frame.to_csv(trait, sep="\t", index=False)
    tree = read_tree(TREE, "1", True, quiet=True, rooted="yes")
    data = prepare_pagel_data(
        trait,
        ("first", "second"),
        list(tree.leaf_names()),
        missing_values_arg="NA",
    )
    likelihood = data.likelihood_by_leaf["A"]
    assert _is_missing_observation(data.observed_state_by_leaf["A"])
    assert "null" in data.observed_state_by_leaf["A"]
    assert likelihood.sum() == pytest.approx(2.0)
    assert np.count_nonzero(likelihood) == 2
    _, fit = compute_pagel_marginals(
        tree, data, model="PAGEL-INDEPENDENT", root_prior_mode="empirical"
    )
    expected_counts = np.zeros(4)
    for _leaf_name, leaf_likelihood in data.likelihood_by_leaf.items():
        expected_counts += leaf_likelihood / leaf_likelihood.sum()
    assert fit["root_prior"] == pytest.approx(expected_counts / expected_counts.sum())


def test_pagel_explicit_states_allow_an_unobserved_binary_state(tmp_path):
    trait = _trait_file(tmp_path)
    frame = pd.read_csv(trait, sep="\t", dtype=str)
    frame["second"] = "x"
    frame.to_csv(trait, sep="\t", index=False)
    tree = read_tree(TREE, "1", True, quiet=True, rooted="yes")
    with pytest.raises(ValueError, match="at least two states"):
        prepare_pagel_data(trait, ("first", "second"), list(tree.leaf_names()))
    data = prepare_pagel_data(
        trait,
        ("first", "second"),
        list(tree.leaf_names()),
        states_arg="0,1;x,y",
    )
    assert data.trait_states == (("0", "1"), ("x", "y"))
    assert len(data.joint_states) == 4


@pytest.mark.integration
def test_pagel_cli_and_batch_comparison(tmp_path):
    trait = _trait_file(tmp_path)
    output = tmp_path / "asr.tsv"
    model = tmp_path / "model.tsv"
    comparison = tmp_path / "comparison.tsv"
    common = [
        "-i",
        TREE,
        "--trait",
        str(trait),
        "--state-column",
        "first,second",
        "--trait-type",
        "discrete",
    ]
    main(
        [
            "asr",
            *common,
            "--model",
            "PAGEL-INDEPENDENT",
            "--model-out",
            str(model),
            "-o",
            str(output),
        ]
    )
    result = pd.read_csv(output, sep="\t")
    fit = pd.read_csv(model, sep="\t").iloc[0]
    assert len([column for column in result if column.startswith("p_")]) == 4
    assert fit["model"] == "PAGEL-INDEPENDENT"
    assert fit["pagel_trait_columns"] == "first,second"
    assert len(fit["rate_classes"].split(",")) == 4

    main(
        [
            "asrcompare",
            *common,
            "--models",
            "PAGEL-INDEPENDENT,PAGEL-DEPENDENT",
            "-o",
            str(comparison),
        ]
    )
    table = pd.read_csv(comparison, sep="\t")
    assert set(table["model"]) == {"PAGEL-INDEPENDENT", "PAGEL-DEPENDENT"}
    assert table["comparison_group"].nunique() == 1
    assert set(table["likelihood_kind"]) == {"pagel_joint_ml"}
    assert table["aic_weight"].sum() == pytest.approx(1.0)


@pytest.mark.integration
def test_asrcompare_all_routes_two_discrete_columns_to_mixture_and_pagel(tmp_path):
    trait = _trait_file(tmp_path)
    comparison = tmp_path / "all.tsv"
    main(
        [
            "asrcompare",
            "-i",
            TREE,
            "--trait",
            str(trait),
            "--state-column",
            "first,second",
            "--trait-type",
            "discrete",
            "--models",
            "all",
            "--rate-categories",
            "2",
            "-o",
            str(comparison),
        ]
    )
    table = pd.read_csv(comparison, sep="\t").set_index("model")
    assert table.loc["PAGEL-INDEPENDENT", "likelihood_kind"] == "pagel_joint_ml"
    assert table.loc["PAGEL-DEPENDENT", "likelihood_kind"] == "pagel_joint_ml"
    assert table.loc["MK-MIXTURE", "likelihood_kind"] == "discrete_ml"
    assert (
        table.loc["PAGEL-INDEPENDENT", "comparison_group"]
        == table.loc["PAGEL-DEPENDENT", "comparison_group"]
    )
    assert (
        table.loc["MK-MIXTURE", "comparison_group"]
        != table.loc["PAGEL-INDEPENDENT", "comparison_group"]
    )
    assert table.loc["ER", "status"] == "not_applicable"

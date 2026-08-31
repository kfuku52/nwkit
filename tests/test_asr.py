import math

import numpy as np
import pandas as pd
import pytest
from ete4 import Tree

import nwkit.asr as asr
from nwkit import rsc_diagnostics
from nwkit.asr import _er_transition_matrix, asr_main
from nwkit.clade_index import CladeIndex
from nwkit.cli import main
from nwkit.util import MISSING_SUPPORT_VALUE, read_tree
from tests.helpers import make_args


def _write_trait(tmp_path, rows, name="traits.tsv"):
    path = tmp_path / name
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False)
    return str(path)


def test_initial_rate_stays_scaled_for_huge_finite_branch_lengths():
    tree = Tree(
        "((A:1e308,B:1e308):1e308,C:1e308);",
        parser=1,
    )

    rate = asr._initial_rate_value(
        tree,
        rate=None,
        rate_bounds=(1e-320, 1e3),
    )

    assert rate == pytest.approx(1e-308, rel=1e-12, abs=0.0)


def test_cli_defaults_and_stochastic_maps_do_not_depend_on_time_units(tmp_path):
    trait = tmp_path / "trait.tsv"
    trait.write_text("leaf_name\tstate\nA\tx\nB\ty\n")
    maps, posteriors = [], []
    for index, scale in enumerate((1.0, 1e9)):
        tree = tmp_path / f"tree{index}.nwk"
        tree.write_text(f"(A:{scale},B:{scale});")
        table = tmp_path / f"posterior{index}.tsv"
        mapping = tmp_path / f"map{index}.tsv"
        main(
            [
                "asr",
                "-i",
                str(tree),
                "--trait",
                str(trait),
                "--state-column",
                "state",
                "--rate",
                str(0.2 / scale),
                "--n-sim",
                "20",
                "--seed",
                "9",
                "--stochastic-map-out",
                str(mapping),
                "-o",
                str(table),
            ]
        )
        maps.append(pd.read_csv(mapping, sep="\t"))
        posteriors.append(pd.read_csv(table, sep="\t"))
    pd.testing.assert_frame_equal(maps[0], maps[1])
    pd.testing.assert_frame_equal(posteriors[0], posteriors[1])
    assert maps[0]["total_count"].sum() >= 20


def test_tiny_asymmetric_rates_are_not_misclassified_as_er():
    matrix = np.asarray([[-1.0, 1.0], [2.0, -2.0]])
    assert asr._get_er_rate_from_matrix(matrix * 1e-18) is None
    assert asr._transition_matrix(matrix, 1.0) == pytest.approx(
        asr._transition_matrix(matrix * 1e-18, 1e18)
    )


def test_rsc_origin_mapping_preserves_a_fixed_process_across_time_units(monkeypatch):
    original_fit = asr.compute_mk_marginals
    results, omissions = [], []
    for scale in (1.0, 1e9):
        tree = Tree("((A:1,B:1):1,(C:1,D:1):1);", parser=1)
        for node in tree.traverse():
            if not node.is_root:
                node.dist *= scale

        def fixed_process(*args, _scale=scale, **kwargs):
            return original_fit(*args, rate=0.2 / _scale, **kwargs)

        # Hold Q*t fixed at the fitting boundary; exercise the real mapping
        # consumer, origin labels, and descendant-event sensitivity sets.
        monkeypatch.setattr(rsc_diagnostics, "compute_mk_marginals", fixed_process)
        clades = CladeIndex(tree)
        frame, omitted = rsc_diagnostics.build_categorical_origin_diagnostics(
            tree,
            {"state": {"A": "x", "B": "y", "C": "x", "D": "y"}},
            ["state"],
            [
                clades.clade_id_for_node(node)
                for node in tree.traverse()
                if not node.is_leaf
            ],
            num_simulations=20,
            seed=9,
        )
        results.append(frame.drop(columns=["mk_rates"]))
        omissions.append(omitted)
    pd.testing.assert_frame_equal(results[0], results[1])
    assert omissions[0] == omissions[1]
    assert results[0]["total_count"].sum() >= 40


@pytest.mark.parametrize(
    "length, matrix", [(0.0, [[-1, 1], [1, -1]]), (1.0, [[0, 0], [0, 0]])]
)
def test_impossible_zero_event_bridge_is_rejected(length, matrix):
    with pytest.raises(ValueError, match="cannot change state"):
        asr._sample_bridge_event_count(
            0, 1, length, np.asarray(matrix), np.random.default_rng(0)
        )


def test_unreachable_stochastic_bridge_is_rejected():
    matrix = np.asarray([[-1.0, 1.0, 0.0], [1.0, -1.0, 0.0], [0.0, 0.0, 0.0]])
    with pytest.raises(ValueError, match="No feasible stochastic bridge"):
        asr._sample_bridge_event_count(0, 2, 1.0, matrix, np.random.default_rng(0))


@pytest.mark.parametrize(
    "tree_text, expected_names",
    [
        ("((A:1,B:1)'95':1,C:2)'20';", {"95", "20"}),
        ("((A:1,B:1):1,C:2);", set()),
    ],
)
def test_annotated_tree_preserves_numeric_names_and_missing_support(
    tmp_path, tree_text, expected_names
):
    source = tmp_path / "source.nwk"
    source.write_text(tree_text)
    tree = read_tree(str(source), "auto", True)
    posterior = {node: np.asarray([0.25, 0.75]) for node in tree.traverse()}
    destination = tmp_path / "annotated.nwk"
    asr._write_annotated_tree(
        tree,
        ["x", "y"],
        posterior,
        {},
        make_args(
            tree_out=str(destination), tree_outformat="auto", tree_annotation="all"
        ),
    )
    assert str(int(MISSING_SUPPORT_VALUE)) not in destination.read_text()
    restored = read_tree(str(destination), "auto", True)
    assert {
        node.name for node in restored.traverse() if not node.is_leaf and node.name
    } == expected_names
    assert all(node.props["asr_state"] == "y" for node in restored.traverse())


class TestAsrMain:
    def test_rejects_colliding_primary_and_model_outputs(self, tmp_path):
        output = tmp_path / "same.tsv"
        args = make_args(
            outfile=str(output),
            model_out=str(output),
            tree_out=None,
            stochastic_map_out=None,
        )
        with pytest.raises(ValueError, match="Output paths must be distinct"):
            asr_main(args)

    def test_probabilities_report_internal_nodes_and_missing_tips(
        self, tmp_nwk, tmp_path
    ):
        infile = tmp_nwk("((A:1,B:1):1,(C:1,D:1):1);", "tree.nwk")
        trait = _write_trait(
            tmp_path,
            [
                {"leaf_name": "A", "habitat": "x"},
                {"leaf_name": "B", "habitat": "x"},
                {"leaf_name": "C", "habitat": "y"},
                {"leaf_name": "D", "habitat": ""},
            ],
        )
        outfile = tmp_path / "asr.tsv"
        args = make_args(
            infile=infile,
            outfile=str(outfile),
            trait=trait,
            state_column="habitat",
            states="x,y",
            missing_values=None,
            model="ER",
            rate=0.1,
            rate_bounds=None,
            root_prior="equal",
            target="intnode,missing_tip",
            output="probabilities",
        )
        asr_main(args)
        table = pd.read_csv(outfile, sep="\t")
        assert {
            "branch_id",
            "parent",
            "node_class",
            "name",
            "map_state",
            "map_probability",
            "p_x",
            "p_y",
        }.issubset(table.columns)
        assert len(table.index) == 4
        assert set(table["node_class"]) == {"root", "intnode", "leaf"}
        for _, row in table.iterrows():
            assert abs((row["p_x"] + row["p_y"]) - 1.0) < 1e-9
            assert row["map_state"] in ["x", "y"]
            assert abs(row["map_probability"] - max(row["p_x"], row["p_y"])) < 1e-9
        d_row = table.loc[table["name"] == "D"].iloc[0]
        assert bool(d_row["is_imputed"]) is True
        assert d_row["p_y"] > d_row["p_x"]

    def test_map_output_can_report_all_nodes(self, tmp_nwk, tmp_path):
        infile = tmp_nwk("((A:1,B:1):1,C:2);", "tree.nwk")
        trait = _write_trait(
            tmp_path,
            [
                {"leaf_name": "A", "state": "red"},
                {"leaf_name": "B", "state": "red"},
                {"leaf_name": "C", "state": "blue"},
            ],
        )
        outfile = tmp_path / "asr_map.tsv"
        args = make_args(
            infile=infile,
            outfile=str(outfile),
            trait=trait,
            state_column="state",
            states="red,blue",
            missing_values=None,
            model="ER",
            rate=0.2,
            rate_bounds=None,
            root_prior="equal",
            target="all",
            output="map",
        )
        asr_main(args)
        table = pd.read_csv(outfile, sep="\t")
        assert list(table.columns) == [
            "branch_id",
            "parent",
            "node_class",
            "name",
            "observed_state",
            "is_imputed",
            "state",
            "probability",
        ]
        assert len(table.index) == 5
        a_row = table.loc[table["name"] == "A"].iloc[0]
        assert a_row["state"] == "red"
        assert abs(a_row["probability"] - 1.0) < 1e-9

    def test_default_output_includes_observed_tips(self, tmp_nwk, tmp_path):
        infile = tmp_nwk("((A:1,B:1):1,C:2);", "tree.nwk")
        trait = _write_trait(
            tmp_path,
            [
                {"leaf_name": "A", "state": "red"},
                {"leaf_name": "B", "state": "red"},
                {"leaf_name": "C", "state": "blue"},
            ],
        )
        outfile = tmp_path / "asr_default.tsv"
        args = make_args(
            infile=infile,
            outfile=str(outfile),
            trait=trait,
            state_column="state",
            states="red,blue",
            missing_values=None,
            model="ER",
            rate=0.2,
            rate_bounds=None,
            root_prior="equal",
            output="probabilities",
        )
        asr_main(args)
        table = pd.read_csv(outfile, sep="\t")
        assert len(table.index) == 5
        tip_rows = table.loc[table["node_class"] == "leaf"].set_index("name")
        assert set(tip_rows.index) == {"A", "B", "C"}
        assert tip_rows.loc["A", "observed_state"] == "red"
        assert bool(tip_rows.loc["A", "is_imputed"]) is False
        assert tip_rows.loc["A", "p_red"] == pytest.approx(1.0)
        assert tip_rows.loc["C", "observed_state"] == "blue"
        assert tip_rows.loc["C", "p_blue"] == pytest.approx(1.0)

    def test_omitted_trait_leaf_is_imputed(self, tmp_nwk, tmp_path):
        infile = tmp_nwk("((A:1,B:1):1,C:2);", "tree.nwk")
        trait = _write_trait(
            tmp_path,
            [
                {"leaf_name": "A", "state": "red"},
                {"leaf_name": "B", "state": "red"},
            ],
            name="partial_traits.tsv",
        )
        outfile = tmp_path / "imputed.tsv"
        args = make_args(
            infile=infile,
            outfile=str(outfile),
            trait=trait,
            state_column="state",
            states="red,blue",
            missing_values=None,
            model="ER",
            rate=0.2,
            rate_bounds=None,
            root_prior="equal",
            target="missing_tip",
            output="probabilities",
        )
        asr_main(args)
        table = pd.read_csv(outfile, sep="\t")
        assert table["name"].tolist() == ["C"]
        assert bool(table.loc[0, "is_imputed"]) is True
        assert abs((table.loc[0, "p_red"] + table.loc[0, "p_blue"]) - 1.0) < 1e-9

    def test_rejects_observed_state_not_listed_in_states(self, tmp_nwk, tmp_path):
        infile = tmp_nwk("(A:1,B:1);", "tree.nwk")
        trait = _write_trait(
            tmp_path,
            [
                {"leaf_name": "A", "state": "red"},
                {"leaf_name": "B", "state": "green"},
            ],
        )
        args = make_args(
            infile=infile,
            outfile=str(tmp_path / "out.tsv"),
            trait=trait,
            state_column="state",
            states="red,blue",
            missing_values=None,
            model="ER",
            rate=0.2,
            rate_bounds=None,
            root_prior="equal",
            target="all",
            output="probabilities",
        )
        with pytest.raises(ValueError, match="--states"):
            asr_main(args)

    def test_zero_rate_rejects_conflicting_tip_states(self, tmp_nwk, tmp_path):
        infile = tmp_nwk("(A:1,B:1);", "tree.nwk")
        trait = _write_trait(
            tmp_path,
            [
                {"leaf_name": "A", "state": "red"},
                {"leaf_name": "B", "state": "blue"},
            ],
        )
        args = make_args(
            infile=infile,
            outfile=str(tmp_path / "out.tsv"),
            trait=trait,
            state_column="state",
            states="red,blue",
            missing_values=None,
            model="ER",
            rate=0.0,
            rate_bounds=None,
            root_prior="equal",
            target="all",
            output="probabilities",
        )
        with pytest.raises(ValueError, match="zero likelihood"):
            asr_main(args)

    def test_er_transition_matrix_matches_two_state_formula(self):
        matrix = _er_transition_matrix(branch_length=1.5, rate=0.2, num_states=2)
        decay = math.exp(-2.0 * 0.2 * 1.5)
        assert matrix[0, 0] == pytest.approx(0.5 + 0.5 * decay)
        assert matrix[0, 1] == pytest.approx(0.5 - 0.5 * decay)
        assert matrix[1, 0] == pytest.approx(0.5 - 0.5 * decay)
        assert matrix[1, 1] == pytest.approx(0.5 + 0.5 * decay)

    def test_inside_likelihood_reuses_transition_matrix_for_equal_branch_lengths(
        self, monkeypatch
    ):
        tree = Tree("((A:1,B:1):1,C:1);", parser=1)
        states = ["red", "blue"]
        likelihood_by_leaf = {
            "A": pd.Series([1.0, 0.0]).to_numpy(),
            "B": pd.Series([1.0, 0.0]).to_numpy(),
            "C": pd.Series([0.0, 1.0]).to_numpy(),
        }
        rate_matrix = asr._build_rate_matrix("ER", states, [0.2])
        call_count = {"value": 0}
        original_transition_matrix = asr._transition_matrix

        def counted_transition_matrix(rate_matrix_arg, branch_length_arg):
            call_count["value"] += 1
            return original_transition_matrix(rate_matrix_arg, branch_length_arg)

        monkeypatch.setattr(asr, "_transition_matrix", counted_transition_matrix)
        asr._compute_inside_likelihoods(tree, likelihood_by_leaf, rate_matrix)
        assert call_count["value"] == 1

    def test_model_out_reports_fixed_er_metadata(self, tmp_nwk, tmp_path):
        infile = tmp_nwk("((A:1,B:1):1,C:2);", "tree.nwk")
        trait = _write_trait(
            tmp_path,
            [
                {"leaf_name": "A", "state": "red"},
                {"leaf_name": "B", "state": "red"},
                {"leaf_name": "C", "state": "blue"},
            ],
        )
        model_out = tmp_path / "model.tsv"
        args = make_args(
            infile=infile,
            outfile=str(tmp_path / "asr.tsv"),
            trait=trait,
            state_column="state",
            states="red,blue",
            missing_values=None,
            model="ER",
            rate=0.2,
            rate_bounds=None,
            root_prior="equal",
            target="all",
            output="map",
            model_out=str(model_out),
        )
        asr_main(args)
        table = pd.read_csv(model_out, sep="\t")
        row = table.iloc[0]
        assert row["model"] == "ER"
        assert bool(row["rate_estimated"]) is False
        assert row["states"] == "red,blue"
        assert row["rate"] == pytest.approx(0.2)
        assert row["rate_red_to_blue"] == pytest.approx(0.2)
        assert row["rate_blue_to_red"] == pytest.approx(0.2)

    def test_sym_and_ard_models_estimate_rates(self, tmp_nwk, tmp_path):
        infile = tmp_nwk("((A:1,B:1):1,(C:1,D:1):1);", "tree.nwk")
        trait = _write_trait(
            tmp_path,
            [
                {"leaf_name": "A", "state": "red"},
                {"leaf_name": "B", "state": "blue"},
                {"leaf_name": "C", "state": "green"},
                {"leaf_name": "D", "state": "red"},
            ],
        )
        for model in ["SYM", "ARD"]:
            outfile = tmp_path / "{}.tsv".format(model.lower())
            model_out = tmp_path / "{}_model.tsv".format(model.lower())
            args = make_args(
                infile=infile,
                outfile=str(outfile),
                trait=trait,
                state_column="state",
                states="red,blue,green",
                missing_values=None,
                model=model,
                rate=None,
                rate_bounds="1e-4,10",
                root_prior="equal",
                target="intnode",
                output="probabilities",
                model_out=str(model_out),
            )
            asr_main(args)
            model_table = pd.read_csv(model_out, sep="\t")
            assert model_table.loc[0, "model"] == model
            assert bool(model_table.loc[0, "rate_estimated"]) is True
            rate_columns = [
                column
                for column in model_table.columns
                if column.startswith("rate_")
                and column != "rate_bounds"
                and column != "rate_estimated"
            ]
            assert all(
                float(model_table.loc[0, column]) >= 0.0 for column in rate_columns
            )
            out_table = pd.read_csv(outfile, sep="\t")
            probability_columns = [
                column for column in out_table.columns if column.startswith("p_")
            ]
            for _, row in out_table.iterrows():
                assert (
                    abs(sum(float(row[column]) for column in probability_columns) - 1.0)
                    < 1e-6
                )

    def test_ambiguous_states_use_multi_hot_tip_likelihood(self, tmp_nwk, tmp_path):
        infile = tmp_nwk("((A:1,B:1):1,C:2);", "tree.nwk")
        trait = _write_trait(
            tmp_path,
            [
                {"leaf_name": "A", "state": "red|blue"},
                {"leaf_name": "B", "state": "red"},
                {"leaf_name": "C", "state": "blue"},
            ],
        )
        outfile = tmp_path / "ambiguous.tsv"
        args = make_args(
            infile=infile,
            outfile=str(outfile),
            trait=trait,
            state_column="state",
            states="red,blue",
            missing_values=None,
            model="ER",
            rate=0.2,
            rate_bounds=None,
            root_prior="empirical",
            target="tip",
            output="probabilities",
            ambiguous_separator="|",
        )
        asr_main(args)
        table = pd.read_csv(outfile, sep="\t")
        a_row = table.loc[table["name"] == "A"].iloc[0]
        assert a_row["observed_state"] == "red|blue"
        assert a_row["p_red"] > 0.0
        assert a_row["p_blue"] > 0.0

    def test_custom_ambiguous_separator_preserves_structured_states(
        self, tmp_nwk, tmp_path
    ):
        infile = tmp_nwk("(A:1,B:1);", "tree.nwk")
        trait = _write_trait(
            tmp_path,
            [
                {"leaf_name": "A", "state": "red|warm;blue"},
                {"leaf_name": "B", "state": "red|warm"},
            ],
        )
        outfile = tmp_path / "custom_ambiguous.tsv"
        args = make_args(
            infile=infile,
            outfile=str(outfile),
            trait=trait,
            state_column="state",
            states="red|warm,blue",
            missing_values=None,
            model="ER",
            rate=0.2,
            rate_bounds=None,
            root_prior="equal",
            target="leaf",
            output="probabilities",
            ambiguous_separator=";",
        )

        asr_main(args)

        table = pd.read_csv(outfile, sep="\t")
        a_row = table.loc[table["name"] == "A"].iloc[0]
        warm_column = "p_{}".format(asr._safe_column_state("red|warm"))
        assert a_row["observed_state"] == "red|warm;blue"
        assert a_row[warm_column] > 0.0
        assert a_row["p_blue"] > 0.0

    def test_state_identifiers_are_injective_across_output_artifacts(self):
        states = ["a-b", "a b", "state_612d62", "A"]
        state_ids = [asr._safe_column_state(state) for state in states]
        assert len(state_ids) == len(set(state_ids))

        tree = Tree("(A:1,B:1);", parser=1)
        posterior_by_node = {
            node: pd.Series([0.1, 0.2, 0.3, 0.4]).to_numpy() for node in tree.traverse()
        }
        observed_state_by_leaf = {"A": "a-b", "B": "a b"}
        table = asr._build_output_table(
            tree=tree,
            states=states,
            observed_state_by_leaf=observed_state_by_leaf,
            posterior_by_node=posterior_by_node,
            targets={"all"},
            output_mode="probabilities",
        )
        probability_columns = ["p_{}".format(state_id) for state_id in state_ids]
        assert all(column in table.columns for column in probability_columns)
        assert len(table.columns) == len(set(table.columns))

        asr._annotate_tree(tree, states, posterior_by_node, observed_state_by_leaf)
        for node in tree.traverse():
            assert all(
                "asr_p_{}".format(state_id) in node.props for state_id in state_ids
            )

    def test_state_identifiers_preserve_legacy_safe_underscores(self):
        assert asr._safe_column_state("dark_blue") == "dark_blue"
        assert asr._safe_column_state("café_2") == "café_2"

    def test_optimizer_non_success_is_rejected(self, monkeypatch):
        tree = Tree("(A:1,B:1);", parser=1)
        states = ["red", "blue"]
        likelihood_by_leaf = {
            "A": pd.Series([1.0, 0.0]).to_numpy(),
            "B": pd.Series([0.0, 1.0]).to_numpy(),
        }

        class FailedResult:
            success = False
            fun = 1.0
            x = pd.Series([0.0]).to_numpy()
            message = "iteration limit reached"

        monkeypatch.setattr(asr, "minimize", lambda *args, **kwargs: FailedResult())
        with pytest.raises(ValueError, match="iteration limit reached"):
            asr._fit_rate_matrix(
                tree=tree,
                model="ER",
                states=states,
                likelihood_by_leaf=likelihood_by_leaf,
                root_prior=pd.Series([0.5, 0.5]).to_numpy(),
            )

    def test_tree_out_writes_nhx_annotations(self, tmp_nwk, tmp_path):
        infile = tmp_nwk("((A:1,B:1):1,C:2);", "tree.nwk")
        trait = _write_trait(
            tmp_path,
            [
                {"leaf_name": "A", "state": "red"},
                {"leaf_name": "B", "state": "red"},
                {"leaf_name": "C", "state": "blue"},
            ],
        )
        tree_out = tmp_path / "annotated.nwk"
        args = make_args(
            infile=infile,
            outfile=str(tmp_path / "asr.tsv"),
            trait=trait,
            state_column="state",
            states="red,blue",
            missing_values=None,
            model="ER",
            rate=0.2,
            rate_bounds=None,
            root_prior="equal",
            target="intnode",
            output="map",
            tree_out=str(tree_out),
            tree_outformat="auto",
            tree_annotation="all",
        )
        asr_main(args)
        annotated = tree_out.read_text()
        assert "&&NHX" in annotated
        assert "asr_state=" in annotated
        assert "asr_probability=" in annotated
        assert "asr_p_red=" in annotated

    def test_stochastic_map_out_is_reproducible_with_seed(self, tmp_nwk, tmp_path):
        infile = tmp_nwk("((A:1,B:1):1,C:2);", "tree.nwk")
        trait = _write_trait(
            tmp_path,
            [
                {"leaf_name": "A", "state": "red"},
                {"leaf_name": "B", "state": "red"},
                {"leaf_name": "C", "state": "blue"},
            ],
        )
        first_out = tmp_path / "map1.tsv"
        second_out = tmp_path / "map2.tsv"
        common = dict(
            infile=infile,
            trait=trait,
            state_column="state",
            states="red,blue",
            missing_values=None,
            model="ER",
            rate=0.5,
            rate_bounds=None,
            root_prior="equal",
            target="intnode",
            output="map",
            n_sim=20,
            seed=7,
        )
        asr_main(
            make_args(
                outfile=str(tmp_path / "asr1.tsv"),
                stochastic_map_out=str(first_out),
                **common,
            )
        )
        asr_main(
            make_args(
                outfile=str(tmp_path / "asr2.tsv"),
                stochastic_map_out=str(second_out),
                **common,
            )
        )
        first = pd.read_csv(first_out, sep="\t")
        second = pd.read_csv(second_out, sep="\t")
        pd.testing.assert_frame_equal(first, second)
        assert {
            "from_state",
            "to_state",
            "mean_count",
            "posterior_frequency",
            "num_simulations",
        }.issubset(first.columns)
        assert set(first["num_simulations"]) == {20}

    def test_stochastic_map_reuses_uniformization_context_for_equal_branch_lengths(
        self, tmp_nwk, tmp_path, monkeypatch
    ):
        infile = tmp_nwk("((A:1,B:1):1,C:1);", "tree.nwk")
        trait = _write_trait(
            tmp_path,
            [
                {"leaf_name": "A", "state": "red"},
                {"leaf_name": "B", "state": "red"},
                {"leaf_name": "C", "state": "blue"},
            ],
        )
        call_lengths = list()
        original_builder = asr._build_uniformization_context

        def counted_builder(rate_matrix, branch_length):
            call_lengths.append(float(branch_length))
            return original_builder(rate_matrix, branch_length)

        monkeypatch.setattr(asr, "_build_uniformization_context", counted_builder)
        args = make_args(
            infile=infile,
            outfile=str(tmp_path / "asr.tsv"),
            trait=trait,
            state_column="state",
            states="red,blue",
            missing_values=None,
            model="ER",
            rate=0.5,
            rate_bounds=None,
            root_prior="equal",
            target="intnode",
            output="map",
            stochastic_map_out=str(tmp_path / "smap.tsv"),
            n_sim=5,
            seed=7,
            threads=1,
        )
        asr_main(args)
        assert call_lengths == [1.0]

    def test_threaded_stochastic_map_matches_single_thread_with_seed(
        self,
        tmp_nwk,
        tmp_path,
        monkeypatch,
    ):
        executor_calls = dict()

        class InlineExecutor:
            def __init__(self, **kwargs):
                executor_calls["kwargs"] = kwargs

            def __enter__(self):
                return self

            def __exit__(self, exc_type, exc, traceback):
                return False

            def map(self, function, payloads):
                payloads = list(payloads)
                executor_calls["chunk_sizes"] = [
                    len(seed_sequences) for _, seed_sequences in payloads
                ]
                return [function(payload) for payload in payloads]

        monkeypatch.setattr(asr, "ProcessPoolExecutor", InlineExecutor)
        monkeypatch.setattr(asr, "_get_process_pool_context", lambda: None)
        infile = tmp_nwk("((A:1,B:1):1,C:2);", "tree.nwk")
        trait = _write_trait(
            tmp_path,
            [
                {"leaf_name": "A", "state": "red"},
                {"leaf_name": "B", "state": "red"},
                {"leaf_name": "C", "state": "blue"},
            ],
        )
        single_out = tmp_path / "single.tsv"
        threaded_out = tmp_path / "threaded.tsv"
        common = dict(
            infile=infile,
            trait=trait,
            state_column="state",
            states="red,blue",
            missing_values=None,
            model="ER",
            rate=0.5,
            rate_bounds=None,
            root_prior="equal",
            target="intnode",
            output="map",
            n_sim=8,
            seed=11,
        )
        asr_main(
            make_args(
                outfile=str(tmp_path / "asr_single.tsv"),
                stochastic_map_out=str(single_out),
                threads=1,
                **common,
            )
        )
        asr_main(
            make_args(
                outfile=str(tmp_path / "asr_threaded.tsv"),
                stochastic_map_out=str(threaded_out),
                threads=4,
                **common,
            )
        )
        single = pd.read_csv(single_out, sep="\t")
        threaded = pd.read_csv(threaded_out, sep="\t")
        pd.testing.assert_frame_equal(single, threaded)
        assert executor_calls == {
            "kwargs": {"max_workers": 4},
            "chunk_sizes": [2, 2, 2, 2],
        }

    def test_stochastic_map_rejects_non_positive_threads(self, tmp_nwk, tmp_path):
        infile = tmp_nwk("(A:1,B:1);", "tree.nwk")
        trait = _write_trait(
            tmp_path,
            [
                {"leaf_name": "A", "state": "red"},
                {"leaf_name": "B", "state": "blue"},
            ],
        )
        args = make_args(
            infile=infile,
            outfile=str(tmp_path / "asr.tsv"),
            trait=trait,
            state_column="state",
            states="red,blue",
            missing_values=None,
            model="ER",
            rate=0.5,
            rate_bounds=None,
            root_prior="equal",
            target="all",
            output="map",
            stochastic_map_out=str(tmp_path / "smap.tsv"),
            n_sim=5,
            seed=1,
            threads=0,
        )
        with pytest.raises(ValueError, match="--threads"):
            asr_main(args)

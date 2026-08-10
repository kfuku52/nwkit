import math

import pandas as pd
import pytest

from nwkit.consensus import consensus_main
from nwkit.util import read_tree
from tests.helpers import make_args
from tests.helpers import write_tree_collection as _write_tree_collection


class TestConsensusMain:
    def test_threaded_consensus_matches_single_thread(self, tmp_path):
        infile = _write_tree_collection(
            tmp_path,
            [
                "((A:1,B:2):1,(C:3,D:4):5);",
                "((A:5,B:6):3,(C:7,D:8):9);",
                "((A:1,C:1):1,(B:1,D:1):1);",
            ],
            name="threaded.nwk",
        )
        single_out = str(tmp_path / "single.nwk")
        threaded_out = str(tmp_path / "threaded.nwk")
        common = dict(
            infile=infile,
            min_freq=0.5,
            reference=None,
            reference_format="auto",
            support_scale="percent",
            method="greedy",
            branch_length="mean",
            weight_tsv=None,
        )

        consensus_main(make_args(outfile=single_out, threads=1, **common))
        consensus_main(make_args(outfile=threaded_out, threads=2, **common))

        assert (
            read_tree(
                single_out, format="auto", quoted_node_names=True, quiet=True
            ).write()
            == read_tree(
                threaded_out,
                format="auto",
                quoted_node_names=True,
                quiet=True,
            ).write()
        )

    def test_strict_consensus_removes_conflicting_clades(self, tmp_path):
        infile = _write_tree_collection(
            tmp_path,
            [
                "((A:1,B:1):1,(C:1,D:1):1);",
                "((A:1,C:1):1,(B:1,D:1):1);",
            ],
            name="strict.nwk",
        )
        outfile = str(tmp_path / "strict_consensus.nwk")
        args = make_args(
            infile=infile,
            outfile=outfile,
            min_freq=0.5,
            reference=None,
            reference_format="auto",
            support_scale="percent",
            method="strict",
            branch_length="none",
            weight_tsv=None,
        )
        consensus_main(args)
        tree = read_tree(outfile, format="auto", quoted_node_names=True, quiet=True)
        assert {frozenset(child.leaf_names()) for child in tree.get_children()} == {
            frozenset({"A"}),
            frozenset({"B"}),
            frozenset({"C"}),
            frozenset({"D"}),
        }

    def test_builds_majority_consensus_tree(self, tmp_path):
        infile = _write_tree_collection(
            tmp_path,
            [
                "((A:1,B:1):1,(C:1,D:1):1);",
                "((A:1,B:1):1,(C:1,D:1):1);",
                "((A:1,C:1):1,(B:1,D:1):1);",
            ],
        )
        outfile = str(tmp_path / "consensus.nwk")
        args = make_args(
            infile=infile,
            outfile=outfile,
            min_freq=0.5,
            reference=None,
            reference_format="auto",
            support_scale="percent",
            method="greedy",
            branch_length="none",
            weight_tsv=None,
        )
        consensus_main(args)
        tree = read_tree(outfile, format="auto", quoted_node_names=True, quiet=True)
        root_children = {frozenset(child.leaf_names()) for child in tree.get_children()}
        assert root_children == {frozenset({"A", "B"}), frozenset({"C", "D"})}
        assert abs(tree.common_ancestor(["A", "B"]).support - (200.0 / 3.0)) < 1e-4
        assert abs(tree.common_ancestor(["C", "D"]).support - (200.0 / 3.0)) < 1e-4

    def test_transfers_consensus_support_to_reference_tree(self, tmp_path, tmp_nwk):
        infile = _write_tree_collection(
            tmp_path,
            [
                "((A:1,B:1):1,(C:1,D:1):1);",
                "((A:1,B:1):1,(C:1,D:1):1);",
                "((A:1,C:1):1,(B:1,D:1):1);",
            ],
        )
        reference = tmp_nwk("((A:1,C:1):1,(B:1,D:1):1);", "reference.nwk")
        outfile = str(tmp_path / "reference_supported.nwk")
        args = make_args(
            infile=infile,
            outfile=outfile,
            min_freq=0.5,
            reference=reference,
            reference_format="auto",
            support_scale="percent",
            method="greedy",
            branch_length="none",
            weight_tsv=None,
        )
        consensus_main(args)
        tree = read_tree(outfile, format="auto", quoted_node_names=True, quiet=True)
        assert abs(tree.common_ancestor(["A", "C"]).support - (100.0 / 3.0)) < 1e-4
        assert abs(tree.common_ancestor(["B", "D"]).support - (100.0 / 3.0)) < 1e-4

    def test_rejects_mismatched_leaf_sets(self, tmp_path):
        infile = _write_tree_collection(
            tmp_path,
            [
                "((A:1,B:1):1,(C:1,D:1):1);",
                "((A:1,B:1):1,(C:1,E:1):1);",
            ],
        )
        args = make_args(
            infile=infile,
            outfile="-",
            min_freq=0.5,
            reference=None,
            reference_format="auto",
            support_scale="percent",
            method="greedy",
            branch_length="none",
            weight_tsv=None,
        )
        with pytest.raises(ValueError, match="Leaf labels must be identical"):
            consensus_main(args)

    def test_weighted_consensus_prefers_heavier_tree(self, tmp_path):
        infile = _write_tree_collection(
            tmp_path,
            [
                "((A:1,B:1):1,(C:1,D:1):1);",
                "((A:1,C:1):1,(B:1,D:1):1);",
            ],
            name="weighted.nwk",
        )
        weight_tsv = tmp_path / "weights.tsv"
        pd.DataFrame({"weight": [3.0, 1.0]}).to_csv(weight_tsv, sep="\t", index=False)
        outfile = str(tmp_path / "weighted_consensus.nwk")
        args = make_args(
            infile=infile,
            outfile=outfile,
            min_freq=0.5,
            reference=None,
            reference_format="auto",
            support_scale="percent",
            method="greedy",
            branch_length="none",
            weight_tsv=str(weight_tsv),
        )
        consensus_main(args)
        tree = read_tree(outfile, format="auto", quoted_node_names=True, quiet=True)
        assert {frozenset(child.leaf_names()) for child in tree.get_children()} == {
            frozenset({"A", "B"}),
            frozenset({"C", "D"}),
        }
        assert abs(tree.common_ancestor(["A", "B"]).support - 75.0) < 1e-6

    def test_consensus_can_average_branch_lengths(self, tmp_path):
        infile = _write_tree_collection(
            tmp_path,
            [
                "((A:1,B:2):1,(C:3,D:4):5);",
                "((A:5,B:6):3,(C:7,D:8):9);",
            ],
            name="branch_lengths.nwk",
        )
        outfile = str(tmp_path / "branch_length_consensus.nwk")
        args = make_args(
            infile=infile,
            outfile=outfile,
            min_freq=0.5,
            reference=None,
            reference_format="auto",
            support_scale="percent",
            method="majority",
            branch_length="mean",
            weight_tsv=None,
        )
        consensus_main(args)
        tree = read_tree(outfile, format="auto", quoted_node_names=True, quiet=True)
        assert abs(tree.common_ancestor(["A", "B"]).dist - 2.0) < 1e-6
        assert abs(next(tree.search_nodes(name="A")).dist - 3.0) < 1e-6

    def test_extreme_finite_weights_are_normalized_before_aggregation(
        self,
        tmp_path,
    ):
        infile = _write_tree_collection(
            tmp_path,
            [
                "((A:1,B:2):1,(C:3,D:4):5);",
                "((A:5,B:6):3,(C:7,D:8):9);",
            ],
            name="extreme_weights.nwk",
        )
        weight_tsv = tmp_path / "extreme_weights.tsv"
        pd.DataFrame({"weight": [1e308, 1e308]}).to_csv(
            weight_tsv,
            sep="\t",
            index=False,
        )
        outfile = str(tmp_path / "extreme_weight_consensus.nwk")

        consensus_main(
            make_args(
                infile=infile,
                outfile=outfile,
                min_freq=0.5,
                reference=None,
                reference_format="auto",
                support_scale="proportion",
                method="strict",
                branch_length="mean",
                weight_tsv=str(weight_tsv),
            )
        )

        tree = read_tree(
            outfile,
            format="auto",
            quoted_node_names=True,
            quiet=True,
        )
        assert tree.common_ancestor(["A", "B"]).support == pytest.approx(1.0)
        assert tree.common_ancestor(["A", "B"]).dist == pytest.approx(2.0)
        assert next(tree.search_nodes(name="A")).dist == pytest.approx(3.0)

    @pytest.mark.parametrize("comparison", ["rooted", "unrooted"])
    @pytest.mark.parametrize("threads", [1, 2])
    def test_mean_branch_lengths_avoid_overflow(
        self,
        tmp_path,
        comparison,
        threads,
    ):
        infile = _write_tree_collection(
            tmp_path,
            [
                "(A:1e308,B:1e308,(C:1e308,D:1e308):1e308);",
                "(A:1e308,B:1e308,(C:1e308,D:1e308):1e308);",
            ],
            name="large_branch_lengths_{}_{}.nwk".format(
                comparison,
                threads,
            ),
        )
        outfile = str(
            tmp_path
            / "large_branch_consensus_{}_{}.nwk".format(
                comparison,
                threads,
            )
        )

        consensus_main(
            make_args(
                infile=infile,
                outfile=outfile,
                min_freq=0.5,
                reference=None,
                reference_format="auto",
                support_scale="proportion",
                method="strict",
                comparison=comparison,
                branch_length="mean",
                weight_tsv=None,
                threads=threads,
            )
        )

        tree = read_tree(
            outfile,
            format="auto",
            quoted_node_names=True,
            quiet=True,
        )
        assert tree.common_ancestor(["C", "D"]).dist == pytest.approx(1e308)
        assert all(
            math.isfinite(float(node.dist))
            for node in tree.traverse()
            if node.dist is not None
        )

    def test_singleton_clades_are_counted_once_and_lengths_are_collapsed(
        self, tmp_path
    ):
        infile = _write_tree_collection(
            tmp_path,
            [
                "((((A:1,B:1):2):3,C:1):1,D:1);",
                "(((A:1,B:1):4,C:1):1,D:1);",
            ],
            name="singleton_clades.nwk",
        )
        outfile = str(tmp_path / "singleton_consensus.nwk")
        args = make_args(
            infile=infile,
            outfile=outfile,
            min_freq=0.5,
            reference=None,
            reference_format="auto",
            support_scale="percent",
            method="strict",
            branch_length="mean",
            weight_tsv=None,
        )

        consensus_main(args)

        tree = read_tree(outfile, format="auto", quoted_node_names=True, quiet=True)
        ab = tree.common_ancestor(["A", "B"])
        assert ab.support == pytest.approx(100.0)
        assert ab.dist == pytest.approx(4.5)

    def test_consensus_accepts_its_own_root_polytomy_output(self, tmp_path):
        first_input = _write_tree_collection(
            tmp_path,
            [
                "((A:1,B:1):1,(C:1,D:1):1);",
                "((A:1,C:1):1,(B:1,D:1):1);",
            ],
            name="conflicting.nwk",
        )
        first_output = str(tmp_path / "polytomy.nwk")
        common = dict(
            min_freq=0.5,
            reference=None,
            reference_format="auto",
            support_scale="percent",
            method="strict",
            branch_length="none",
            weight_tsv=None,
        )
        consensus_main(make_args(infile=first_input, outfile=first_output, **common))

        second_output = str(tmp_path / "polytomy_roundtrip.nwk")
        consensus_main(make_args(infile=first_output, outfile=second_output, **common))

        tree = read_tree(
            second_output, format="auto", quoted_node_names=True, quiet=True
        )
        assert {frozenset(child.leaf_names()) for child in tree.get_children()} == {
            frozenset({"A"}),
            frozenset({"B"}),
            frozenset({"C"}),
            frozenset({"D"}),
        }

    def test_unrooted_consensus_is_independent_of_newick_root_position(
        self,
        tmp_path,
    ):
        infile = _write_tree_collection(
            tmp_path,
            [
                "(A,B,(C,(D,E)));",
                "(D,E,(C,(A,B)));",
            ],
            name="rerooted-unrooted.nwk",
        )
        outfile = str(tmp_path / "unrooted-consensus.nwk")
        args = make_args(
            infile=infile,
            outfile=outfile,
            min_freq=0.5,
            reference=None,
            reference_format="auto",
            support_scale="proportion",
            method="greedy",
            comparison="unrooted",
            branch_length="none",
            weight_tsv=None,
        )

        consensus_main(args)

        tree = read_tree(outfile, format=0, quoted_node_names=True, quiet=True)
        assert tree.common_ancestor(["D", "E"]).support == pytest.approx(1.0)
        assert tree.common_ancestor(["C", "D", "E"]).support == pytest.approx(1.0)

    def test_weight_tsv_rejects_non_integer_tree_ids(self, tmp_path):
        infile = _write_tree_collection(
            tmp_path,
            [
                "((A:1,B:1):1,(C:1,D:1):1);",
                "((A:1,C:1):1,(B:1,D:1):1);",
            ],
            name="weighted_invalid.nwk",
        )
        weight_tsv = tmp_path / "weights.tsv"
        pd.DataFrame({"tree_id": [1.5, 2.0], "weight": [3.0, 1.0]}).to_csv(
            weight_tsv, sep="\t", index=False
        )
        args = make_args(
            infile=infile,
            outfile="-",
            min_freq=0.5,
            reference=None,
            reference_format="auto",
            support_scale="percent",
            method="greedy",
            branch_length="none",
            weight_tsv=str(weight_tsv),
        )
        with pytest.raises(ValueError, match="tree_id"):
            consensus_main(args)

    def test_weight_tsv_rejects_nan_weights(self, tmp_path):
        infile = _write_tree_collection(
            tmp_path,
            [
                "((A:1,B:1):1,(C:1,D:1):1);",
                "((A:1,C:1):1,(B:1,D:1):1);",
            ],
            name="weighted_nan.nwk",
        )
        weight_tsv = tmp_path / "weights.tsv"
        pd.DataFrame({"weight": [float("nan"), 1.0]}).to_csv(
            weight_tsv, sep="\t", index=False
        )
        args = make_args(
            infile=infile,
            outfile="-",
            min_freq=0.5,
            reference=None,
            reference_format="auto",
            support_scale="percent",
            method="greedy",
            branch_length="none",
            weight_tsv=str(weight_tsv),
        )
        with pytest.raises(ValueError, match="weight"):
            consensus_main(args)

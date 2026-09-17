import gzip
import io

import pandas as pd
import pytest

from nwkit.nwk2table import nwk2table_main
from nwkit.table2nwk import table2nwk_main
from nwkit.util import read_tree
from tests.helpers import make_args


class TestTable2NwkMain:
    @pytest.mark.parametrize("mode", ["stdin", "gzip", "plain"])
    def test_input_containers_preserve_literal_node_names(
        self, tmp_path, monkeypatch, mode
    ):
        text = "\ufeffbranch_id\tparent\tname\n0\t-1\tR\n1\t0\t001\n2\t0\tNA\n"
        source = tmp_path / ("input.tsv.gz" if mode == "gzip" else "input.tsv")
        if mode == "gzip":
            with gzip.open(source, "wt", encoding="utf-8") as handle:
                handle.write(text)
        elif mode == "stdin":
            monkeypatch.setattr("sys.stdin", io.StringIO(text))
        else:
            source.write_text(text, encoding="utf-8")
        output = tmp_path / "out.nwk"
        table2nwk_main(
            make_args(
                infile="-" if mode == "stdin" else str(source),
                outfile=str(output),
            )
        )
        tree = read_tree(str(output), "1", True, quiet=True)
        assert list(tree.leaf_names()) == ["001", "NA"]

    @pytest.mark.parametrize("column", ["branch_id", "parent", "name", "dist"])
    @pytest.mark.parametrize("prefix", ["", "\ufeff", "\n \n"])
    def test_rejects_ambiguous_duplicate_headers_before_output(
        self, tmp_path, column, prefix
    ):
        table_path = tmp_path / "ambiguous.tsv"
        table_path.write_text(
            f"{prefix}branch_id\tparent\tname\tdist\t{column}\n"
            "0\t-1\troot\t0\t0\n1\t0\t001\t1\t99\n",
            encoding="utf-8",
        )
        outfile = tmp_path / "output.nwk"
        outfile.write_text("original tree")
        with pytest.raises(ValueError, match="duplicat.*column"):
            table2nwk_main(make_args(infile=str(table_path), outfile=str(outfile)))
        assert outfile.read_text() == "original tree"

    def test_roundtrip_support_tree(self, tmp_nwk, tmp_path):
        infile = tmp_nwk("((A:1,B:2)80:3,C:4);", "tree.nwk")
        table_path = tmp_path / "tree.tsv"
        outfile = tmp_path / "roundtrip.nwk"
        nwk2table_main(
            make_args(infile=infile, outfile=str(table_path), age=False, sister=True)
        )
        table2nwk_main(
            make_args(infile=str(table_path), outfile=str(outfile), outformat="auto")
        )
        tree = read_tree(
            str(outfile), format="auto", quoted_node_names=True, quiet=True
        )
        assert set(tree.leaf_names()) == {"A", "B", "C"}
        assert abs(next(tree.search_nodes(name="A")).dist - 1.0) < 1e-6
        assert abs(next(tree.search_nodes(name="B")).dist - 2.0) < 1e-6
        assert abs(tree.common_ancestor(["A", "B"]).support - 80.0) < 1e-6

    def test_roundtrip_named_internal_nodes(self, tmp_nwk, tmp_path):
        infile = tmp_nwk("((A:1,B:2)AB:3,C:4)root;", "tree.nwk")
        table_path = tmp_path / "tree.tsv"
        outfile = tmp_path / "roundtrip_named.nwk"
        nwk2table_main(
            make_args(infile=infile, outfile=str(table_path), age=False, sister=True)
        )
        table2nwk_main(
            make_args(infile=str(table_path), outfile=str(outfile), outformat="auto")
        )
        tree = read_tree(
            str(outfile), format="auto", quoted_node_names=True, quiet=True
        )
        assert tree.name == "root"
        assert tree.common_ancestor(["A", "B"]).name == "AB"
        assert abs(next(tree.search_nodes(name="A")).dist - 1.0) < 1e-6

    def test_rejects_multiple_roots(self, tmp_path):
        table_path = tmp_path / "invalid.tsv"
        pd.DataFrame(
            [
                {"branch_id": 0, "parent": -1, "name": "root1"},
                {"branch_id": 1, "parent": -1, "name": "root2"},
            ]
        ).to_csv(table_path, sep="\t", index=False)
        args = make_args(infile=str(table_path), outfile="-", outformat="auto")
        with pytest.raises(ValueError, match="exactly one root"):
            table2nwk_main(args)

    @pytest.mark.parametrize(
        "column,value",
        [
            ("dist", "inf"),
            ("dist", "-Infinity"),
            ("support", "NAN"),
            ("support", "1e9999"),
        ],
    )
    def test_rejects_non_finite_numeric_values(self, tmp_path, column, value):
        table_path = tmp_path / "invalid.tsv"
        rows = [
            {"branch_id": 0, "parent": -1, "name": "root", column: ""},
            {"branch_id": 1, "parent": 0, "name": "A", column: value},
        ]
        pd.DataFrame(rows).to_csv(table_path, sep="\t", index=False)

        with pytest.raises(
            ValueError, match="Column '{}' must contain finite".format(column)
        ):
            table2nwk_main(
                make_args(
                    infile=str(table_path),
                    outfile=str(tmp_path / "output.nwk"),
                    outformat="auto",
                )
            )

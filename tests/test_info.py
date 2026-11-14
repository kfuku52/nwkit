import io

import pytest

from nwkit.cli import main
from nwkit.info import info_main
from tests.helpers import make_args


@pytest.mark.parametrize("source", ["stdin", "inline", "file"])
def test_cli_reports_actual_input_source(source, tmp_path, monkeypatch, capsys):
    monkeypatch.chdir(tmp_path)
    text = "(A:1,B:1);"
    if source == "stdin":
        # '-' always means stdin, even when a file with that name exists.
        (tmp_path / "-").write_text("(C:2,D:2,E:2);")
        monkeypatch.setattr("sys.stdin", io.StringIO(text))
        infile, expected = "-", "Tree input: stdin"
    elif source == "inline":
        infile, expected = text, "Tree input: inline text"
    else:
        (tmp_path / "tree.nwk").write_text(text)
        infile = "tree.nwk"
        expected = f"Tree file PATH: {(tmp_path / infile).resolve()}"
    main(["info", "-i", infile])
    output = capsys.readouterr().out
    assert output.splitlines()[0] == expected
    assert "Number of leaves: 2\n" in output
    assert "Tree length: 2.0\n" in output


class TestInfoMain:
    def test_basic_info(self, tmp_nwk, capsys):
        path = tmp_nwk("((A:1,B:1):1,(C:1,D:1):1);")
        args = make_args(infile=path)
        info_main(args)
        captured = capsys.readouterr()
        assert "Number of leaves: 4" in captured.out
        assert "Number of nodes: 7" in captured.out

    def test_singleton_count(self, tmp_nwk, capsys):
        path = tmp_nwk("(((A:1,B:1):1):1,(C:1,D:1):1);")
        args = make_args(infile=path)
        info_main(args)
        captured = capsys.readouterr()
        assert "Number of singleton nodes: 1" in captured.out

    def test_multifurcation_count(self, tmp_nwk, capsys):
        path = tmp_nwk("(A:1,B:1,C:1,D:1);")
        args = make_args(infile=path)
        info_main(args)
        captured = capsys.readouterr()
        assert "Number of multifurcation nodes: 1" in captured.out

    def test_zero_branch_length(self, tmp_nwk, capsys):
        path = tmp_nwk("((A:0,B:1):1,(C:1,D:0):1);")
        args = make_args(infile=path)
        info_main(args)
        captured = capsys.readouterr()
        assert "Number of nodes with zero branch length: 2" in captured.out

    def test_species_names(self, tmp_nwk, capsys):
        nwk = "((Homo_sapiens_G1:1,Mus_musculus_G1:1):1,Danio_rerio_G1:1);"
        path = tmp_nwk(nwk)
        args = make_args(infile=path)
        info_main(args)
        captured = capsys.readouterr()
        assert "Number of species" in captured.out
        assert "Homo_sapiens" in captured.out
        assert "Mus_musculus" in captured.out
        assert "Danio_rerio" in captured.out

    def test_tree_length(self, tmp_nwk, capsys):
        path = tmp_nwk("((A:1,B:2):3,(C:4,D:5):6);")
        args = make_args(infile=path)
        info_main(args)
        captured = capsys.readouterr()
        # Tree length = 1+2+3+4+5+6 = 21
        assert "Tree length: 21" in captured.out

    def test_negative_branch_length(self, tmp_nwk, capsys):
        path = tmp_nwk("((A:-1,B:1):1,(C:1,D:1):1);")
        args = make_args(infile=path)
        info_main(args)
        captured = capsys.readouterr()
        assert "Number of nodes with negative branch length: 1" in captured.out

    def test_unnamed_leaf_not_counted_as_empty_species(self, tmp_nwk, capsys):
        path = tmp_nwk("(:1,B:1);")
        args = make_args(infile=path)
        info_main(args)
        captured = capsys.readouterr()
        assert (
            "Number of species in the leaf name convention of GENUS_SPECIES_GENEID: 1"
            in captured.out
        )
        assert "Species names: B" in captured.out

    def test_custom_species_regex(self, tmp_nwk, capsys):
        path = tmp_nwk("((Homo.sapiens|A:1,Mus.musculus|A:1):1,Homo.sapiens|B:1);")
        args = make_args(infile=path, species_regex=r"^([A-Za-z]+)\.([A-Za-z]+)\|")
        info_main(args)
        captured = capsys.readouterr()
        assert "Number of species parsed by --species-regex: 2" in captured.out
        assert "Species names: Homo_sapiens, Mus_musculus" in captured.out

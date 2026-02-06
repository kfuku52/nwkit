import os
import pytest
from ete4 import Tree

from nwkit.info import info_main
from nwkit.util import read_tree
from tests.helpers import make_args, DATA_DIR


class TestInfoMain:
    def test_basic_info(self, tmp_nwk, capsys):
        path = tmp_nwk('((A:1,B:1):1,(C:1,D:1):1);')
        args = make_args(infile=path)
        info_main(args)
        captured = capsys.readouterr()
        assert 'Number of leaves: 4' in captured.out
        assert 'Number of nodes: 7' in captured.out

    def test_singleton_count(self, tmp_nwk, capsys):
        path = tmp_nwk('(((A:1,B:1):1):1,(C:1,D:1):1);')
        args = make_args(infile=path)
        info_main(args)
        captured = capsys.readouterr()
        assert 'Number of singleton nodes: 1' in captured.out

    def test_multifurcation_count(self, tmp_nwk, capsys):
        path = tmp_nwk('(A:1,B:1,C:1,D:1);')
        args = make_args(infile=path)
        info_main(args)
        captured = capsys.readouterr()
        assert 'Number of multifurcation nodes: 1' in captured.out

    def test_zero_branch_length(self, tmp_nwk, capsys):
        path = tmp_nwk('((A:0,B:1):1,(C:1,D:0):1);')
        args = make_args(infile=path)
        info_main(args)
        captured = capsys.readouterr()
        assert 'Number of nodes with zero branch length: 2' in captured.out

    def test_species_names(self, tmp_nwk, capsys):
        nwk = '((Homo_sapiens_G1:1,Mus_musculus_G1:1):1,Danio_rerio_G1:1);'
        path = tmp_nwk(nwk)
        args = make_args(infile=path)
        info_main(args)
        captured = capsys.readouterr()
        assert 'Number of species' in captured.out
        assert 'Homo_sapiens' in captured.out
        assert 'Mus_musculus' in captured.out
        assert 'Danio_rerio' in captured.out

    def test_tree_length(self, tmp_nwk, capsys):
        path = tmp_nwk('((A:1,B:2):3,(C:4,D:5):6);')
        args = make_args(infile=path)
        info_main(args)
        captured = capsys.readouterr()
        # Tree length = 1+2+3+4+5+6 = 21
        assert 'Tree length: 21' in captured.out

    def test_with_data_file(self, capsys):
        infile = os.path.join(DATA_DIR, 'info1', 'input.nwk')
        if not os.path.exists(infile):
            pytest.skip('Test data not found')
        args = make_args(infile=infile)
        info_main(args)
        captured = capsys.readouterr()
        assert 'Number of leaves:' in captured.out

    def test_wiki_exact_output(self, capsys):
        """Wiki example: 29-leaf tree with Danio_rerio, Gadus_morhua, Xenopus_tropicalis.

        Expected output:
        Number of leaves: 29
        Number of nodes: 57
        Number of singleton nodes: 0
        Number of multifurcation nodes: 0
        Number of nodes with zero branch length: 0
        Number of species: 3
        """
        infile = os.path.join(DATA_DIR, 'info1', 'input.nwk')
        if not os.path.exists(infile):
            pytest.skip('Test data not found')
        args = make_args(infile=infile)
        info_main(args)
        captured = capsys.readouterr()
        assert 'Number of leaves: 29' in captured.out
        assert 'Number of nodes: 57' in captured.out
        assert 'Number of singleton nodes: 0' in captured.out
        assert 'Number of multifurcation nodes: 0' in captured.out
        assert 'Number of nodes with zero branch length: 0' in captured.out
        assert 'Number of species in the leaf name convention of GENUS_SPECIES_GENEID: 3' in captured.out
        assert 'Danio_rerio' in captured.out
        assert 'Gadus_morhua' in captured.out
        assert 'Xenopus_tropicalis' in captured.out

    def test_negative_branch_length(self, tmp_nwk, capsys):
        path = tmp_nwk('((A:-1,B:1):1,(C:1,D:1):1);')
        args = make_args(infile=path)
        info_main(args)
        captured = capsys.readouterr()
        assert 'Number of nodes with negative branch length: 1' in captured.out

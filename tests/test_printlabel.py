import os
import pytest
from ete4 import Tree

from nwkit.printlabel import printlabel_main
from nwkit.util import read_tree
from tests.helpers import make_args, DATA_DIR


class TestPrintlabelMain:
    def test_print_all_labels(self, tmp_nwk, capsys):
        path = tmp_nwk('((A:1,B:1)AB:1,(C:1,D:1)CD:1)root;')
        args = make_args(
            infile=path, format='1',
            pattern='.*', target='all', sister=False,
        )
        printlabel_main(args)
        captured = capsys.readouterr()
        assert 'A' in captured.out
        assert 'B' in captured.out
        assert 'AB' in captured.out

    def test_print_leaf_labels(self, tmp_nwk, capsys):
        path = tmp_nwk('((A:1,B:1)AB:1,(C:1,D:1)CD:1)root;')
        args = make_args(
            infile=path, format='1',
            pattern='.*', target='leaf', sister=False,
        )
        printlabel_main(args)
        captured = capsys.readouterr()
        assert 'A' in captured.out
        assert 'B' in captured.out
        assert 'AB' not in captured.out

    def test_print_intnode_labels(self, tmp_nwk, capsys):
        path = tmp_nwk('((A:1,B:1)AB:1,(C:1,D:1)CD:1)root;')
        args = make_args(
            infile=path, format='1',
            pattern='.*', target='intnode', sister=False,
        )
        printlabel_main(args)
        captured = capsys.readouterr()
        assert 'AB' in captured.out
        assert 'CD' in captured.out

    def test_pattern_filter(self, tmp_nwk, capsys):
        path = tmp_nwk('((A1:1,A2:1):1,(B1:1,B2:1):1);')
        args = make_args(
            infile=path,
            pattern='A.*', target='leaf', sister=False,
        )
        printlabel_main(args)
        captured = capsys.readouterr()
        assert 'A1' in captured.out
        assert 'A2' in captured.out
        assert 'B1' not in captured.out

    def test_sister_output(self, tmp_nwk, capsys):
        path = tmp_nwk('((A:1,B:1):1,(C:1,D:1):1);')
        args = make_args(
            infile=path,
            pattern='A', target='leaf', sister=True,
        )
        printlabel_main(args)
        captured = capsys.readouterr()
        assert 'B' in captured.out

    def test_no_match(self, tmp_nwk, capsys):
        path = tmp_nwk('((A:1,B:1):1,(C:1,D:1):1);')
        args = make_args(
            infile=path,
            pattern='Z.*', target='leaf', sister=False,
        )
        printlabel_main(args)
        captured = capsys.readouterr()
        assert captured.out.strip() == ''

    def test_exact_leaf_label_set(self, tmp_nwk, capsys):
        """Verify exact set of leaf labels printed with pattern filter."""
        path = tmp_nwk('((A1:1,A2:1):1,(B1:1,B2:1):1);')
        args = make_args(
            infile=path,
            pattern='A.*', target='leaf', sister=False,
        )
        printlabel_main(args)
        captured = capsys.readouterr()
        printed = set(captured.out.strip().split('\n'))
        assert printed == {'A1', 'A2'}

    def test_sister_exact_labels(self, tmp_nwk, capsys):
        """Verify exact sister labels for pattern-matched leaves."""
        path = tmp_nwk('((A:1,B:1):1,(C:1,D:1):1);')
        args = make_args(
            infile=path,
            pattern='A', target='leaf', sister=True,
        )
        printlabel_main(args)
        captured = capsys.readouterr()
        # Sister of A is B
        assert 'B' in captured.out.strip()
        # Should not contain C or D
        assert 'C' not in captured.out
        assert 'D' not in captured.out

    def test_intnode_exact_labels(self, tmp_nwk, capsys):
        """Verify exact set of internal node names printed."""
        path = tmp_nwk('((A:1,B:1)AB:1,(C:1,D:1)CD:1)root;')
        args = make_args(
            infile=path, format='1',
            pattern='.*', target='intnode', sister=False,
        )
        printlabel_main(args)
        captured = capsys.readouterr()
        printed = set(captured.out.strip().split('\n'))
        assert printed == {'AB', 'CD', 'root'}

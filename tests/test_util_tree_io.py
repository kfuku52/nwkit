import os
import sys
import pytest
import io
from ete4 import Tree

from nwkit import util as util_mod
from nwkit.util import (
    _compile_node_placeholder_pattern,
    iter_newick_stream,
    split_newick_stream,
    read_tree,
    write_tree,
)
from tests.helpers import make_args


class TestReadTree:
    def test_read_from_file(self, tmp_nwk):
        path = tmp_nwk('((A:1,B:1):1,(C:1,D:1):1);')
        tree = read_tree(path, format='auto', quoted_node_names=True, quiet=True)
        assert set(tree.leaf_names()) == {'A', 'B', 'C', 'D'}

    def test_read_with_explicit_format(self, tmp_nwk):
        path = tmp_nwk('((A:1,B:1):1,(C:1,D:1):1);')
        tree = read_tree(path, format='1', quoted_node_names=True, quiet=True)
        assert len(list(tree.leaves())) == 4

    def test_streaming_newick_parser_handles_quote_boundaries(self):
        text = "('A''quoted':1,B:1);\n(C:1,D:1);\n"
        assert list(iter_newick_stream(io.StringIO(text), chunk_size=1)) == [
            "('A''quoted':1,B:1);",
            '(C:1,D:1);',
        ]

    def test_newick_collection_parsers_handle_double_quoted_semicolons(self):
        text = '("A;""quoted":1,B:1);\n(C:1,D:1);\n'
        expected = [
            '("A;""quoted":1,B:1);',
            '(C:1,D:1);',
        ]
        assert split_newick_stream(text) == expected
        assert list(iter_newick_stream(io.StringIO(text), chunk_size=1)) == expected

    def test_newick_collection_parsers_ignore_semicolons_in_nested_comments(self):
        text = '(A[outer;[inner;comment]]:1,B:1);\n(C:1,D:1);\n'
        expected = [
            '(A[outer;[inner;comment]]:1,B:1);',
            '(C:1,D:1);',
        ]
        assert split_newick_stream(text) == expected
        assert list(iter_newick_stream(io.StringIO(text), chunk_size=1)) == expected

    @pytest.mark.parametrize('streaming', [False, True])
    def test_newick_collection_parsers_reject_unterminated_comments(
        self,
        streaming,
    ):
        text = '(A[unterminated;comment'
        with pytest.raises(ValueError, match='unterminated Newick comment'):
            if streaming:
                list(iter_newick_stream(io.StringIO(text), chunk_size=1))
            else:
                split_newick_stream(text)

    def test_quoted_name_detection_ignores_quotes_inside_comments(self):
        text = "(A:1[comment 'text' \"double\"],B:1);"
        tree = read_tree(
            text,
            format='1',
            quoted_node_names=False,
            quiet=True,
        )
        assert set(tree.leaf_names()) == {'A', 'B'}
        inspection = util_mod.inspect_tree_text(
            text,
            format='1',
            quoted_node_names=False,
        )
        assert inspection['parse_ok'] is True
        assert inspection['has_quoted_node_names'] is False
        assert inspection['has_quoted_internal_node_names'] is False

    def test_read_tree_with_internal_names(self, tmp_nwk):
        path = tmp_nwk('((A:1,B:1)AB:1,(C:1,D:1)CD:1)root;')
        tree = read_tree(path, format='auto', quoted_node_names=True, quiet=True)
        leaf_names = set(tree.leaf_names())
        assert leaf_names == {'A', 'B', 'C', 'D'}

    def test_read_tree_auto_format_detection(self, tmp_nwk):
        # Format 0: support values
        path = tmp_nwk('((A:1,B:1)90:1,(C:1,D:1)85:1);')
        tree = read_tree(path, format='auto', quoted_node_names=True, quiet=True)
        assert set(tree.leaf_names()) == {'A', 'B', 'C', 'D'}
        assert abs(tree.children[0].support - 90.0) < 1e-9
        assert str(tree.children[0].name or '') == ''

    def test_read_tree_auto_prefers_quoted_numeric_internal_names(self, tmp_nwk):
        path = tmp_nwk("((A:1,B:1)'42':1,C:1)'99':1;")
        tree = read_tree(path, format='auto', quoted_node_names=True, quiet=True)
        assert tree.name == '99'
        child_names = {str(child.name or '') for child in tree.children}
        assert '42' in child_names

    def test_read_tree_auto_warns_on_ambiguous_numeric_internal_labels(self, tmp_nwk, capsys):
        path = tmp_nwk('((A:1,B:1)42:1,C:1)99:1;')
        tree = read_tree(path, format='auto', quoted_node_names=True, quiet=False)
        captured = capsys.readouterr()
        assert 'Ambiguous tree format' in captured.err
        assert tree.name in ('', None)

    def test_read_tree_preserves_subunit_support_values_when_mixed_with_percent_support(self, tmp_nwk):
        path = tmp_nwk('((A:1,B:1)0.95:1,(C:1,D:1)80:1);')
        tree = read_tree(path, format='0', quoted_node_names=True, quiet=True)
        assert abs(tree.common_ancestor(['A', 'B']).support - 0.95) < 1e-9
        assert abs(tree.common_ancestor(['C', 'D']).support - 80.0) < 1e-9

    def test_read_tree_preserves_explicit_one_percent_support(self, tmp_nwk):
        path = tmp_nwk('(((A:1,B:1)1:1,C:1)99:1,D:1);')
        tree = read_tree(path, format='0', quoted_node_names=True, quiet=True)
        assert tree.common_ancestor(['A', 'B']).support == pytest.approx(1.0)

    def test_read_tree_marks_missing_support_on_proportion_scale(self, tmp_nwk):
        path = tmp_nwk('(((A:1,B:1)0.8:1,C:1):1,D:1);')
        tree = read_tree(path, format='0', quoted_node_names=True, quiet=True)
        assert tree.common_ancestor(['A', 'B']).support == pytest.approx(0.8)
        assert tree.common_ancestor(['A', 'B', 'C']).support == pytest.approx(-999999.0)

    def test_read_tree_auto_strict_raises_on_ambiguous_numeric_internal_labels(self, tmp_nwk):
        path = tmp_nwk('((A:1,B:1)42:1,C:1)99:1;')
        with pytest.raises(ValueError, match='Ambiguous tree format'):
            read_tree(path, format='auto-strict', quoted_node_names=True, quiet=True)

    def test_read_tree_rejects_quoted_names_when_flag_disabled(self, tmp_nwk):
        path = tmp_nwk("('A,B':1,C:2);")
        with pytest.raises(ValueError, match='--quoted-node-names yes'):
            read_tree(path, format='auto', quoted_node_names=False, quiet=True)

    def test_read_tree_rejects_double_quoted_names_when_flag_disabled(self, tmp_nwk):
        path = tmp_nwk('("A,B":1,C:2);')
        with pytest.raises(ValueError, match='--quoted-node-names yes'):
            read_tree(path, format='auto', quoted_node_names=False, quiet=True)

    def test_read_tree_invalid_raises(self, tmp_nwk):
        path = tmp_nwk('not_a_tree')
        with pytest.raises(Exception, match='Failed to parse'):
            read_tree(path, format='auto', quoted_node_names=True, quiet=True)

    def test_read_tree_rejects_non_finite_branch_lengths(self, tmp_nwk):
        path = tmp_nwk('(A:nan,B:1);')
        with pytest.raises(ValueError, match='branch lengths must be finite'):
            read_tree(path, format='auto', quoted_node_names=True, quiet=True)

    def test_read_tree_rejects_non_finite_root_support(self, tmp_nwk):
        path = tmp_nwk('(A:1,B:1)nan;')
        with pytest.raises(ValueError, match='support values must be finite'):
            read_tree(path, format=0, quoted_node_names=True, quiet=True)

    def test_read_from_stdin_multiline(self, monkeypatch):
        monkeypatch.setattr(sys, 'stdin', io.StringIO('((A:1,B:1):1,\n(C:1,D:1):1);\n'))
        tree = read_tree('-', format='auto', quoted_node_names=True, quiet=True)
        assert set(tree.leaf_names()) == {'A', 'B', 'C', 'D'}


class TestWriteTree:
    def test_placeholder_regex_size_is_independent_of_node_count(self):
        prefix = 'NWKITNODE0123456789abcdef_'
        pattern = _compile_node_placeholder_pattern(prefix)

        assert pattern.fullmatch(prefix + '0000000042')
        assert not pattern.fullmatch(prefix + '0000000042extra')
        assert len(pattern.pattern) < 100

    def test_write_to_file(self, tmp_nwk, tmp_outfile):
        path = tmp_nwk('((A:1,B:1):1,(C:1,D:1):1);')
        tree = read_tree(path, format='1', quoted_node_names=True, quiet=True)
        args = make_args(outfile=tmp_outfile)
        write_tree(tree, args, format='1', quiet=True)
        assert os.path.exists(tmp_outfile)
        tree2 = read_tree(tmp_outfile, format='1', quoted_node_names=True, quiet=True)
        assert set(tree2.leaf_names()) == {'A', 'B', 'C', 'D'}

    def test_write_to_stdout(self, tmp_nwk, capsys):
        path = tmp_nwk('((A:1,B:1):1,(C:1,D:1):1);')
        tree = read_tree(path, format='1', quoted_node_names=True, quiet=True)
        args = make_args(outfile='-')
        write_tree(tree, args, format='1', quiet=True)
        captured = capsys.readouterr()
        assert 'A' in captured.out
        assert 'B' in captured.out

    def test_write_preserves_node_names(self, tmp_nwk, tmp_outfile):
        path = tmp_nwk('((A:1,B:1)AB:1,(C:1,D:1)CD:1)root;')
        tree = read_tree(path, format='1', quoted_node_names=True, quiet=True)
        args = make_args(outfile=tmp_outfile)
        write_tree(tree, args, format='1', quiet=True)
        with open(tmp_outfile) as f:
            content = f.read()
        assert 'AB' in content
        assert 'CD' in content

    def test_write_auto_no_subcommand_in_argv(self, tmp_nwk, tmp_outfile, monkeypatch):
        path = tmp_nwk('((A:1,B:1):1,(C:1,D:1):1);')
        tree = read_tree(path, format='auto', quoted_node_names=True, quiet=True)
        args = make_args(outfile=tmp_outfile)
        monkeypatch.setattr(sys, 'argv', ['nwkit'])
        write_tree(tree, args, format='auto', quiet=True)
        tree2 = read_tree(tmp_outfile, format='auto', quoted_node_names=True, quiet=True)
        assert set(tree2.leaf_names()) == {'A', 'B', 'C', 'D'}

    def test_write_auto_without_infile_format_global(self, tmp_outfile, monkeypatch):
        import nwkit.util as util_mod
        tree = Tree('((A:1,B:1):1,(C:1,D:1):1);', parser=1)
        args = make_args(outfile=tmp_outfile)
        monkeypatch.delattr(util_mod, 'INFILE_FORMAT', raising=False)
        monkeypatch.setattr(sys, 'argv', ['nwkit'])
        write_tree(tree, args, format='auto', quiet=True)
        tree2 = read_tree(tmp_outfile, format='auto', quoted_node_names=True, quiet=True)
        assert set(tree2.leaf_names()) == {'A', 'B', 'C', 'D'}

    def test_write_auto_respects_explicit_input_format(self, tmp_nwk, tmp_outfile, monkeypatch):
        path = tmp_nwk('((A:1,B:1)AB:1,(C:1,D:1)CD:1)ROOT;')
        tree = read_tree(path, format='1', quoted_node_names=True, quiet=True)
        args = make_args(outfile=tmp_outfile)
        monkeypatch.setattr(sys, 'argv', ['nwkit'])
        write_tree(tree, args, format='auto', quiet=True)
        with open(tmp_outfile) as f:
            out = f.read()
        assert 'AB' in out
        assert 'CD' in out

    def test_write_auto_uses_tree_specific_format_after_multiple_reads(self, tmp_nwk, tmp_outfile, monkeypatch):
        path1 = tmp_nwk('((A:1,B:1)AB:1,(C:1,D:1)CD:1)ROOT;', 'tree1.nwk')
        path2 = tmp_nwk('((A:1,B:1)90:1,(C:1,D:1)80:1);', 'tree2.nwk')
        tree1 = read_tree(path1, format='1', quoted_node_names=True, quiet=True)
        _ = read_tree(path2, format='auto', quoted_node_names=True, quiet=True)
        args = make_args(outfile=tmp_outfile)
        monkeypatch.setattr(sys, 'argv', ['nwkit'])
        write_tree(tree1, args, format='auto', quiet=True)
        with open(tmp_outfile) as f:
            out = f.read()
        assert 'AB' in out
        assert 'CD' in out

    def test_write_preserves_quoted_leaf_names(self, tmp_nwk, tmp_outfile):
        path = tmp_nwk("('A,B':1,C:1);")
        tree = read_tree(path, format='1', quoted_node_names=True, quiet=True)
        args = make_args(outfile=tmp_outfile)
        write_tree(tree, args, format='1', quiet=True)
        with open(tmp_outfile) as f:
            out = f.read()
        assert "'A,B'" in out
        tree2 = read_tree(tmp_outfile, format='1', quoted_node_names=True, quiet=True)
        assert set(tree2.leaf_names()) == {'A,B', 'C'}

    def test_write_auto_quotes_embedded_double_quote_for_round_trip(
        self,
        tmp_outfile,
    ):
        tree = Tree('(A:1,B:1);', parser=1)
        tree['A'].name = 'a"b'
        write_tree(
            tree,
            make_args(outfile=tmp_outfile),
            format='1',
            quiet=True,
            name_quote='auto',
        )
        with open(tmp_outfile) as handle:
            output = handle.read()
        assert "'a\"b'" in output
        round_tripped = read_tree(
            tmp_outfile,
            format='1',
            quoted_node_names=True,
            quiet=True,
        )
        assert set(round_tripped.leaf_names()) == {'a"b', 'B'}

    def test_write_preserves_numeric_internal_names_with_quotes(self, tmp_nwk, tmp_outfile):
        path = tmp_nwk("((A:1,B:1)'42':1,C:1)'99':1;")
        tree = read_tree(path, format='auto', quoted_node_names=True, quiet=True)
        args = make_args(outfile=tmp_outfile)
        write_tree(tree, args, format='1', quiet=True)
        with open(tmp_outfile) as f:
            out = f.read()
        assert "'42'" in out
        assert "'99'" in out
        tree2 = read_tree(tmp_outfile, format='auto', quoted_node_names=True, quiet=True)
        assert tree2.name == '99'
        child_names = {str(child.name or '') for child in tree2.children}
        assert '42' in child_names

    def test_write_does_not_replace_placeholder_like_property_values(self, tmp_outfile):
        tree = Tree('(A:1,B:1);', parser=1)
        marker = 'NODENAME_PLACEHOLDER0000000000'
        tree['A'].props['tag'] = marker
        args = make_args(outfile=tmp_outfile)
        write_tree(tree, args, format='1', quiet=True, props=['tag'])
        output = read_tree(tmp_outfile, format='1', quoted_node_names=True, quiet=True)
        assert output['A'].props['tag'] == marker

    def test_write_leaves_unknown_placeholder_like_values_unchanged(self, tmp_outfile):
        tree = Tree('(A:1,B:1);', parser=1)
        marker = 'NODENAME_PLACEHOLDER9999999999'
        tree['A'].props['tag'] = marker
        args = make_args(outfile=tmp_outfile)
        write_tree(tree, args, format='1', quiet=True, props=['tag'])
        with open(tmp_outfile) as handle:
            assert marker in handle.read()

    def test_write_restores_tree_when_name_serialization_fails(
        self,
        tmp_nwk,
        tmp_outfile,
        monkeypatch,
    ):
        tree = read_tree(
            tmp_nwk('((A:1,B:1)AB:1,C:1)ROOT;'),
            format='1',
            quoted_node_names=True,
            quiet=True,
        )
        original_values = [
            (node, node.name, node.support)
            for node in tree.traverse()
        ]
        original_serializer = util_mod._serialize_newick_node_name
        calls = 0

        def fail_on_third_name(*args, **kwargs):
            nonlocal calls
            calls += 1
            if calls == 3:
                raise ValueError('injected name serialization failure')
            return original_serializer(*args, **kwargs)

        monkeypatch.setattr(
            util_mod,
            '_serialize_newick_node_name',
            fail_on_third_name,
        )

        with pytest.raises(ValueError, match='injected name serialization failure'):
            write_tree(tree, make_args(outfile=tmp_outfile), format='1', quiet=True)

        assert calls == 3
        assert [
            (node, node.name, node.support)
            for node in tree.traverse()
        ] == original_values
        assert not os.path.exists(tmp_outfile)

    @pytest.mark.parametrize('attribute,value,error_match', [
        ('dist', float('inf'), 'branch lengths must be finite'),
        ('support', float('nan'), 'support values must be finite'),
    ])
    def test_write_rejects_non_finite_tree_values_before_opening_output(
        self,
        tmp_outfile,
        attribute,
        value,
        error_match,
    ):
        tree = Tree('((A:1,B:1):1,C:1);', parser=1)
        target_node = tree['A'] if attribute == 'dist' else tree.common_ancestor(['A', 'B'])
        setattr(target_node, attribute, value)
        with open(tmp_outfile, 'w') as handle:
            handle.write('existing output')

        with pytest.raises(ValueError, match=error_match):
            write_tree(tree, make_args(outfile=tmp_outfile), format='1', quiet=True)

        with open(tmp_outfile) as handle:
            assert handle.read() == 'existing output'

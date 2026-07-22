import os
import sys
import pytest
import io
import pickle
import sqlite3
import tarfile
from contextlib import contextmanager
from ete4 import Tree

from nwkit.util import (
    _break_stale_lock_if_needed,
    _build_ete_taxonomy_database,
    acquire_exclusive_lock,
    extract_species_label,
    extract_taxonomy_query,
    get_ete_ncbitaxa,
    get_monophyletic_species_groups,
    iter_newick_stream,
    resolve_download_dir,
    resolve_ete_data_dir,
    read_tree,
    write_tree,
    remove_singleton,
    label2sciname,
    read_item_per_line_file,
    annotate_scientific_names,
    annotate_duplication_confidence_scores,
    get_subtree_leaf_name_sets,
    get_subtree_leaf_bitmasks,
    validate_unique_named_leaves,
    validate_distinct_output_paths,
    is_all_leaf_names_identical,
    get_target_nodes,
    is_rooted,
)
from tests.helpers import make_args


def _write_valid_taxonomy_db(path):
    connection = sqlite3.connect(path)
    connection.executescript(
        'CREATE TABLE stats (version INT PRIMARY KEY);'
        'CREATE TABLE species (taxid INT PRIMARY KEY, parent INT, spname TEXT, common TEXT, rank TEXT, track TEXT);'
        'CREATE TABLE synonym (taxid INT, spname TEXT);'
        'CREATE TABLE merged (taxid_old INT, taxid_new INT);'
        'INSERT INTO stats VALUES (2);'
        "INSERT INTO species VALUES (1, 1, 'root', '', 'no rank', '1');"
    )
    connection.commit()
    connection.close()
    with open(str(path) + '.traverse.pkl', 'wb') as handle:
        pickle.dump([1, 1], handle)


def _write_valid_taxdump(path):
    source_dir = path.parent / '{}-members'.format(path.name)
    source_dir.mkdir()
    for name in ('nodes.dmp', 'names.dmp', 'merged.dmp'):
        (source_dir / name).write_text('1\t|\t1\t|\n')
    with tarfile.open(path, 'w:gz') as archive:
        for name in ('nodes.dmp', 'names.dmp', 'merged.dmp'):
            archive.add(source_dir / name, arcname=name)


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

    def test_read_tree_auto_strict_raises_on_ambiguous_numeric_internal_labels(self, tmp_nwk):
        path = tmp_nwk('((A:1,B:1)42:1,C:1)99:1;')
        with pytest.raises(ValueError, match='Ambiguous tree format'):
            read_tree(path, format='auto-strict', quoted_node_names=True, quiet=True)

    def test_read_tree_rejects_quoted_names_when_flag_disabled(self, tmp_nwk):
        path = tmp_nwk("('A,B':1,C:2);")
        with pytest.raises(ValueError, match='--quoted-node-names yes'):
            read_tree(path, format='auto', quoted_node_names=False, quiet=True)

    def test_read_tree_invalid_raises(self, tmp_nwk):
        path = tmp_nwk('not_a_tree')
        with pytest.raises(Exception, match='Failed to parse'):
            read_tree(path, format='auto', quoted_node_names=True, quiet=True)

    def test_read_from_stdin_multiline(self, monkeypatch):
        monkeypatch.setattr(sys, 'stdin', io.StringIO('((A:1,B:1):1,\n(C:1,D:1):1);\n'))
        tree = read_tree('-', format='auto', quoted_node_names=True, quiet=True)
        assert set(tree.leaf_names()) == {'A', 'B', 'C', 'D'}


class TestWriteTree:
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


class TestDownloadDirHelpers:
    def test_resolve_download_dir_returns_none_for_auto(self, tmp_path):
        outfile = tmp_path / 'out' / 'tree.nwk'
        args = make_args(outfile=str(outfile), download_dir='auto')
        resolved = resolve_download_dir(args)
        assert resolved is None

    def test_resolve_download_dir_uses_outfile_parent_for_inferred(self, tmp_path):
        outfile = tmp_path / 'out' / 'tree.nwk'
        args = make_args(outfile=str(outfile), download_dir='inferred')
        resolved = resolve_download_dir(args)
        assert resolved == os.path.join(os.path.realpath(outfile.parent), 'downloads')

    def test_resolve_download_dir_uses_cwd_for_stdout_when_inferred(self, monkeypatch, tmp_path):
        monkeypatch.chdir(tmp_path)
        args = make_args(outfile='-', download_dir='inferred')
        resolved = resolve_download_dir(args)
        assert resolved == os.path.join(os.path.realpath(tmp_path), 'downloads')

    def test_resolve_download_dir_respects_explicit_path(self, tmp_path):
        explicit_dir = tmp_path / 'shared-cache'
        args = make_args(outfile='-', download_dir=str(explicit_dir))
        resolved = resolve_download_dir(args)
        assert resolved == os.path.realpath(explicit_dir)

    def test_resolve_ete_data_dir_appends_ete4(self, tmp_path):
        explicit_dir = tmp_path / 'shared-cache'
        args = make_args(download_dir=str(explicit_dir))
        assert resolve_ete_data_dir(args) == os.path.join(os.path.realpath(explicit_dir), 'ete4')

    def test_resolve_ete_data_dir_returns_none_for_auto(self):
        args = make_args(download_dir='auto')
        assert resolve_ete_data_dir(args) is None

    def test_acquire_exclusive_lock_creates_and_removes_file(self, tmp_path):
        lock_path = tmp_path / '.lock'
        with acquire_exclusive_lock(str(lock_path), timeout_seconds=1):
            assert lock_path.exists()
        assert not lock_path.exists()

    def test_acquire_exclusive_lock_breaks_stale_lock(self, tmp_path):
        lock_path = tmp_path / '.lock'
        lock_path.write_text('999999999\n')
        with acquire_exclusive_lock(str(lock_path), timeout_seconds=1):
            assert lock_path.exists()
        assert not lock_path.exists()

    def test_new_empty_lock_receives_metadata_grace_period(self, tmp_path):
        lock_path = tmp_path / '.lock'
        lock_path.write_bytes(b'')
        assert _break_stale_lock_if_needed(str(lock_path)) is False
        assert lock_path.exists()

    def test_distinct_output_paths_rejects_realpath_aliases(self, tmp_path):
        output = tmp_path / 'result.tsv'
        with pytest.raises(ValueError, match='Output paths must be distinct'):
            validate_distinct_output_paths([
                ('--outfile', str(output)),
                ('--report', str(tmp_path / '.' / 'result.tsv')),
            ])

    def test_get_ete_ncbitaxa_uses_download_dir_and_lock(self, monkeypatch, tmp_path):
        explicit_dir = tmp_path / 'shared-cache'
        args = make_args(download_dir=str(explicit_dir))
        calls = dict()

        @contextmanager
        def fake_lock(lock_path, lock_label='Lock', poll_seconds=1, timeout_seconds=3600):
            calls['lock_path'] = lock_path
            calls['lock_label'] = lock_label
            yield

        def fake_download_ete_taxdump(taxdump_file):
            calls['downloaded_taxdump_file'] = taxdump_file
            _write_valid_taxdump(tmp_path / 'shared-cache' / 'ete4' / 'taxdump.tar.gz')
            return True

        def fake_build(dbfile, taxdump_file):
            calls['built_dbfile'] = dbfile
            calls['built_taxdump_file'] = taxdump_file
            _write_valid_taxonomy_db(dbfile)

        class FakeNCBI:
            def __init__(self, dbfile=None, taxdump_file=None, memory=False, update=True):
                calls['dbfile'] = dbfile
                calls['taxdump_file'] = taxdump_file
                calls['taxdump_exists_at_init'] = taxdump_file is not None and os.path.exists(taxdump_file)
                calls['memory'] = memory
                calls['update'] = update
                self.db = None

        monkeypatch.setattr('nwkit.util.acquire_exclusive_lock', fake_lock)
        monkeypatch.setattr('nwkit.util._find_existing_ete_taxonomy_assets', lambda exclude_dbfile=None: None)
        monkeypatch.setattr('nwkit.util._download_ete_taxdump', fake_download_ete_taxdump)
        monkeypatch.setattr('nwkit.util._build_ete_taxonomy_database', fake_build)
        monkeypatch.setattr('nwkit.util.ete4.NCBITaxa', FakeNCBI)
        ncbi = get_ete_ncbitaxa(args=args)
        assert isinstance(ncbi, FakeNCBI)
        expected_ete_dir = os.path.join(os.path.realpath(explicit_dir), 'ete4')
        assert calls['lock_path'] == os.path.join(expected_ete_dir, '.ete4_taxonomy.lock')
        assert calls['lock_label'] == 'ETE4 taxonomy DB'
        assert calls['dbfile'] == os.path.join(expected_ete_dir, 'taxa.sqlite')
        assert calls['taxdump_file'] is None
        assert calls['downloaded_taxdump_file'] == os.path.join(expected_ete_dir, 'taxdump.tar.gz')
        assert calls['built_dbfile'] == os.path.join(expected_ete_dir, 'taxa.sqlite')
        assert calls['built_taxdump_file'] == os.path.join(expected_ete_dir, 'taxdump.tar.gz')

    def test_taxonomy_build_isolated_from_calling_working_directory(self, monkeypatch, tmp_path):
        caller_dir = tmp_path / 'caller'
        target_dir = tmp_path / 'cache'
        caller_dir.mkdir()
        target_dir.mkdir()
        sentinel = caller_dir / 'taxa.tab'
        sentinel.write_text('keep me')
        taxdump = target_dir / 'taxdump.tar.gz'
        taxdump.write_bytes(b'placeholder')
        observed = {}

        def fake_run(command, cwd, check):
            observed['cwd'] = cwd
            observed['command'] = command
            dbfile = command[-2]
            with open(dbfile, 'wb') as handle:
                handle.write(b'database')
            with open(dbfile + '.traverse.pkl', 'wb') as handle:
                pickle.dump([1], handle)

        monkeypatch.chdir(caller_dir)
        monkeypatch.setattr('nwkit.util.subprocess.run', fake_run)
        monkeypatch.setattr('nwkit.util._validate_ete_taxonomy_db', lambda *args, **kwargs: True)
        _build_ete_taxonomy_database(str(target_dir / 'taxa.sqlite'), str(taxdump))
        assert os.path.realpath(observed['cwd']) != os.path.realpath(caller_dir)
        assert sentinel.read_text() == 'keep me'
        assert (target_dir / 'taxa.sqlite').read_bytes() == b'database'

    def test_get_ete_ncbitaxa_auto_uses_ete_default_location(self, monkeypatch, tmp_path):
        calls = dict()

        default_db = tmp_path / 'ete-default' / 'taxa.sqlite'
        default_db.parent.mkdir()
        _write_valid_taxonomy_db(default_db)

        @contextmanager
        def fake_lock(lock_path, **kwargs):
            calls['lock_path'] = lock_path
            yield

        class FakeNCBI:
            def __init__(self, *args, **kwargs):
                calls['args'] = args
                calls['kwargs'] = kwargs
                self.db = None

        monkeypatch.setattr('ete4.ncbi_taxonomy.ncbiquery.DEFAULT_TAXADB', str(default_db))
        monkeypatch.setattr('nwkit.util.acquire_exclusive_lock', fake_lock)
        monkeypatch.setattr('nwkit.util.ete4.NCBITaxa', FakeNCBI)
        ncbi = get_ete_ncbitaxa(args=make_args(download_dir='auto'))
        assert isinstance(ncbi, FakeNCBI)
        assert calls['args'] == ()
        assert calls['kwargs'] == {'dbfile': str(default_db), 'update': False}
        assert calls['lock_path'] == str(default_db.parent / '.ete4_taxonomy.lock')

    def test_get_ete_ncbitaxa_checks_existing_taxdump_before_build(self, monkeypatch, tmp_path):
        explicit_dir = tmp_path / 'shared-cache'
        expected_ete_dir = os.path.join(os.path.realpath(explicit_dir), 'ete4')
        os.makedirs(expected_ete_dir, exist_ok=True)
        _write_valid_taxdump(tmp_path / 'shared-cache' / 'ete4' / 'taxdump.tar.gz')
        args = make_args(download_dir=str(explicit_dir))
        calls = dict(download_count=0)

        @contextmanager
        def fake_lock(lock_path, lock_label='Lock', poll_seconds=1, timeout_seconds=3600):
            calls['lock_path'] = lock_path
            yield

        def fake_download_ete_taxdump(taxdump_file):
            calls['download_count'] += 1
            return False

        def fake_build(dbfile, taxdump_file):
            calls['build_count'] = calls.get('build_count', 0) + 1
            _write_valid_taxonomy_db(dbfile)

        class FakeNCBI:
            def __init__(self, dbfile=None, taxdump_file=None, memory=False, update=True):
                calls['dbfile'] = dbfile
                calls['taxdump_file'] = taxdump_file
                self.db = None

        monkeypatch.setattr('nwkit.util.acquire_exclusive_lock', fake_lock)
        monkeypatch.setattr('nwkit.util._find_existing_ete_taxonomy_assets', lambda exclude_dbfile=None: None)
        monkeypatch.setattr('nwkit.util._download_ete_taxdump', fake_download_ete_taxdump)
        monkeypatch.setattr('nwkit.util._build_ete_taxonomy_database', fake_build)
        monkeypatch.setattr('nwkit.util.ete4.NCBITaxa', FakeNCBI)
        ncbi = get_ete_ncbitaxa(args=args)
        assert isinstance(ncbi, FakeNCBI)
        assert calls['lock_path'] == os.path.join(expected_ete_dir, '.ete4_taxonomy.lock')
        assert calls['dbfile'] == os.path.join(expected_ete_dir, 'taxa.sqlite')
        assert calls['taxdump_file'] is None
        assert calls['download_count'] == 1
        assert calls['build_count'] == 1

    def test_get_ete_ncbitaxa_reuses_existing_db_without_taxdump_update(self, monkeypatch, tmp_path):
        explicit_dir = tmp_path / 'shared-cache'
        expected_ete_dir = os.path.join(os.path.realpath(explicit_dir), 'ete4')
        os.makedirs(expected_ete_dir, exist_ok=True)
        existing_dbfile = os.path.join(expected_ete_dir, 'taxa.sqlite')
        _write_valid_taxonomy_db(existing_dbfile)
        args = make_args(download_dir=str(explicit_dir))
        calls = dict(download_count=0)

        @contextmanager
        def fake_lock(lock_path, lock_label='Lock', poll_seconds=1, timeout_seconds=3600):
            calls['lock_path'] = lock_path
            yield

        def fake_download_ete_taxdump(taxdump_file):
            calls['download_count'] += 1

        class FakeNCBI:
            def __init__(self, dbfile=None, taxdump_file=None, memory=False, update=True):
                calls['dbfile'] = dbfile
                calls['taxdump_file'] = taxdump_file
                calls['update'] = update
                self.db = None

        monkeypatch.setattr('nwkit.util.acquire_exclusive_lock', fake_lock)
        monkeypatch.setattr('nwkit.util._find_existing_ete_taxonomy_assets', lambda exclude_dbfile=None: None)
        monkeypatch.setattr('nwkit.util._download_ete_taxdump', fake_download_ete_taxdump)
        monkeypatch.setattr('nwkit.util.ete4.NCBITaxa', FakeNCBI)
        ncbi = get_ete_ncbitaxa(args=args)
        assert isinstance(ncbi, FakeNCBI)
        assert calls['lock_path'] == os.path.join(expected_ete_dir, '.ete4_taxonomy.lock')
        assert calls['dbfile'] == existing_dbfile
        assert calls['taxdump_file'] is None
        assert calls['update'] is False
        assert calls['download_count'] == 0

    def test_get_ete_ncbitaxa_seeds_local_cache_from_existing_user_db(self, monkeypatch, tmp_path):
        explicit_dir = tmp_path / 'shared-cache'
        expected_ete_dir = os.path.join(os.path.realpath(explicit_dir), 'ete4')
        source_dir = tmp_path / 'user-cache'
        source_dir.mkdir()
        source_dbfile = source_dir / 'taxa.sqlite'
        source_traverse = source_dir / 'taxa.sqlite.traverse.pkl'
        source_taxdump = source_dir / 'taxdump.tar.gz'
        _write_valid_taxonomy_db(source_dbfile)
        _write_valid_taxdump(source_taxdump)
        args = make_args(download_dir=str(explicit_dir))
        calls = dict(download_count=0)

        @contextmanager
        def fake_lock(lock_path, lock_label='Lock', poll_seconds=1, timeout_seconds=3600):
            calls['lock_path'] = lock_path
            yield

        def fake_download_ete_taxdump(taxdump_file):
            calls['download_count'] += 1

        class FakeNCBI:
            def __init__(self, dbfile=None, taxdump_file=None, memory=False, update=True):
                calls['dbfile'] = dbfile
                calls['taxdump_file'] = taxdump_file
                calls['update'] = update
                self.db = None

        monkeypatch.setattr('nwkit.util.acquire_exclusive_lock', fake_lock)
        monkeypatch.setattr(
            'nwkit.util._find_existing_ete_taxonomy_assets',
            lambda exclude_dbfile=None: {
                'dbfile': str(source_dbfile),
                'traverse_file': str(source_traverse),
                'taxdump_file': str(source_taxdump),
            },
        )
        monkeypatch.setattr('nwkit.util._download_ete_taxdump', fake_download_ete_taxdump)
        monkeypatch.setattr('nwkit.util.ete4.NCBITaxa', FakeNCBI)
        ncbi = get_ete_ncbitaxa(args=args)
        assert isinstance(ncbi, FakeNCBI)
        assert calls['lock_path'] == os.path.join(expected_ete_dir, '.ete4_taxonomy.lock')
        assert calls['dbfile'] == os.path.join(expected_ete_dir, 'taxa.sqlite')
        assert calls['taxdump_file'] is None
        assert calls['update'] is False
        assert calls['download_count'] == 0
        assert os.path.getsize(os.path.join(expected_ete_dir, 'taxa.sqlite')) > 0
        assert os.path.getsize(os.path.join(expected_ete_dir, 'taxa.sqlite.traverse.pkl')) > 0
        assert os.path.getsize(os.path.join(expected_ete_dir, 'taxdump.tar.gz')) > 0

    def test_get_ete_ncbitaxa_rebuilds_corrupt_database(self, monkeypatch, tmp_path):
        explicit_dir = tmp_path / 'shared-cache'
        ete_dir = explicit_dir / 'ete4'
        ete_dir.mkdir(parents=True)
        dbfile = ete_dir / 'taxa.sqlite'
        dbfile.write_bytes(b'not a sqlite database')
        calls = {'download': 0, 'build': 0}

        def fake_download(path):
            calls['download'] += 1
            _write_valid_taxdump(ete_dir / 'taxdump.tar.gz')
            return True

        def fake_build(dbfile, taxdump_file):
            calls['build'] += 1
            os.remove(dbfile)
            _write_valid_taxonomy_db(dbfile)

        monkeypatch.setattr('nwkit.util._find_existing_ete_taxonomy_assets', lambda exclude_dbfile=None: None)
        monkeypatch.setattr('nwkit.util._download_ete_taxdump', fake_download)
        monkeypatch.setattr('nwkit.util._build_ete_taxonomy_database', fake_build)
        ncbi = get_ete_ncbitaxa(args=make_args(download_dir=str(explicit_dir)))
        try:
            assert calls == {'download': 1, 'build': 1}
            assert ncbi.get_taxid_translator([1])[1] == 'root'
        finally:
            ncbi.db.close()


class TestSpeciesGrouping:
    def test_monophyletic_duplicate_species_are_grouped(self):
        tree = Tree('(((Homo_sapiens_gene1:1,Homo_sapiens_gene2:1):1,Pan_troglodytes_gene1:1):1,Mus_musculus_gene1:1);', parser=1)
        leaf_name_to_sci_name, species_to_leaf_names = get_monophyletic_species_groups(tree)
        assert leaf_name_to_sci_name['Homo_sapiens_gene1'] == 'Homo_sapiens'
        assert leaf_name_to_sci_name['Homo_sapiens_gene2'] == 'Homo_sapiens'
        assert species_to_leaf_names['Homo_sapiens'] == ['Homo_sapiens_gene1', 'Homo_sapiens_gene2']

    def test_split_duplicate_species_raise(self):
        tree = Tree('((Homo_sapiens_gene1:1,Pan_troglodytes_gene1:1):1,(Homo_sapiens_gene2:1,Mus_musculus_gene1:1):1);', parser=1)
        with pytest.raises(ValueError, match='not monophyletic'):
            get_monophyletic_species_groups(tree)


class TestRemoveSingleton:
    def test_remove_singleton_node(self):
        # Create tree with a singleton: ((A,B)) -> should become (A,B)
        tree = Tree('(((A:1,B:1):1):1,(C:1,D:1):1);', parser=1)
        num_nodes_before = len(list(tree.traverse()))
        tree = remove_singleton(tree, verbose=False, preserve_branch_length=True)
        num_nodes_after = len(list(tree.traverse()))
        assert num_nodes_after < num_nodes_before
        assert set(tree.leaf_names()) == {'A', 'B', 'C', 'D'}

    def test_no_singleton(self, simple_tree):
        num_nodes_before = len(list(simple_tree.traverse()))
        tree = remove_singleton(simple_tree, verbose=False)
        num_nodes_after = len(list(tree.traverse()))
        assert num_nodes_before == num_nodes_after

    def test_preserve_branch_length(self):
        tree = Tree('(((A:1,B:1):2):3,C:6);', parser=1)
        tree = remove_singleton(tree, verbose=False, preserve_branch_length=True)
        # After removing singleton, branch lengths should be preserved (summed)
        assert set(tree.leaf_names()) == {'A', 'B', 'C'}

    def test_remove_singleton_root_wrapper(self):
        tree = Tree('((A:1,B:1):1);', parser=1)
        tree = remove_singleton(tree, verbose=False, preserve_branch_length=True)
        assert set(tree.leaf_names()) == {'A', 'B'}
        assert len(tree.get_children()) == 2
        assert abs(tree.get_distance('A', 'B') - 2.0) < 1e-6
        for node in tree.traverse():
            if not node.is_leaf:
                assert len(node.get_children()) != 1


class TestLabel2Sciname:
    def test_single_string(self):
        result = label2sciname('Homo_sapiens_GENE1')
        assert result == 'Homo_sapiens'

    def test_list_input(self):
        result = label2sciname(['Homo_sapiens_GENE1', 'Mus_musculus_GENE2'])
        assert result == ['Homo_sapiens', 'Mus_musculus']

    def test_no_species_info(self):
        result = label2sciname('SingleWord')
        assert result is None

    def test_custom_delimiter(self):
        result = label2sciname('Homo-sapiens-GENE1', in_delim='-')
        assert result == 'Homo_sapiens'

    def test_custom_out_delimiter(self):
        result = label2sciname('Homo_sapiens_GENE1', out_delim=' ')
        assert result == 'Homo sapiens'

    def test_empty_list(self):
        result = label2sciname([])
        assert result == []

    def test_none_scalar(self):
        result = label2sciname(None)
        assert result is None

    def test_list_with_none(self):
        result = label2sciname(['Homo_sapiens_GENE1', None, 'Mus_musculus_GENE2'])
        assert result == ['Homo_sapiens', None, 'Mus_musculus']


class TestExtractSpeciesLabel:
    def test_default_regex_parses_with_suffix(self):
        assert extract_species_label('Homo_sapiens_GENE1') == 'Homo_sapiens'

    def test_default_regex_parses_exact_binomial(self):
        assert extract_species_label('Homo_sapiens') == 'Homo_sapiens'

    def test_custom_regex_uses_capture_groups(self):
        assert extract_species_label(
            'Homo.sapiens|GENE1',
            species_regex=r'^([A-Za-z]+)\.([A-Za-z]+)\|',
        ) == 'Homo_sapiens'

    def test_custom_regex_allows_space_output(self):
        assert extract_species_label(
            'Homo.sapiens|GENE1',
            species_regex=r'^([A-Za-z]+)\.([A-Za-z]+)\|',
            out_delim=' ',
        ) == 'Homo sapiens'

    def test_taxonomic_keeps_species_label_qualifier(self):
        args = make_args(species_parser='taxonomic')
        assert extract_species_label('Dictyostelium_cf_discoideum', args=args) == 'Dictyostelium_cf_discoideum'
        assert extract_species_label('Amoeba_sp_JDSRuffled', args=args) == 'Amoeba_sp_JDSRuffled'
        assert extract_species_label('Solanum_lycopersicum_cultivar_Heinz1706_gene1', args=args) == 'Solanum_lycopersicum_cultivar_Heinz1706'
        assert extract_species_label('Escherichia_coli_serovar_O157_gene1', args=args) == 'Escherichia_coli_serovar_O157'

    def test_taxonomic_taxonomy_query_falls_back(self):
        args = make_args(species_parser='taxonomic')
        assert extract_taxonomy_query('Dictyostelium_cf_discoideum', args=args) == 'Dictyostelium discoideum'
        assert extract_taxonomy_query('Amoeba_sp_JDSRuffled', args=args) == 'Amoeba'
        assert extract_taxonomy_query('Solanum_lycopersicum_cultivar_Heinz1706_gene1', args=args) == 'Solanum lycopersicum'
        assert extract_taxonomy_query('Escherichia_coli_serovar_O157_gene1', args=args) == 'Escherichia coli'

    def test_species_map_tsv_overrides_species_label_and_taxonomy_query(self, tmp_path):
        map_path = tmp_path / 'species_map.tsv'
        map_path.write_text(
            'leaf_name\tspecies_label\ttaxonomy_query\n'
            'Sample42\tMapped_species\tMapped species\n'
        )
        args = make_args(species_map_tsv=str(map_path))
        assert extract_species_label('Sample42', args=args) == 'Mapped_species'
        assert extract_taxonomy_query('Sample42', args=args) == 'Mapped species'


class TestReadItemPerLineFile:
    def test_basic(self, tmp_path):
        f = tmp_path / 'items.txt'
        f.write_text('apple\nbanana\ncherry\n')
        result = read_item_per_line_file(str(f))
        assert result == ['apple', 'banana', 'cherry']

    def test_empty_lines_stripped(self, tmp_path):
        f = tmp_path / 'items.txt'
        f.write_text('apple\n\nbanana\n\n')
        result = read_item_per_line_file(str(f))
        assert result == ['apple', 'banana']

    def test_crlf_and_whitespace_are_normalized(self, tmp_path):
        f = tmp_path / 'items.txt'
        f.write_bytes(b' apple \r\nbanana\r\n\r\n  cherry  \r\n')
        result = read_item_per_line_file(str(f))
        assert result == ['apple', 'banana', 'cherry']


class TestAnnotateScientificNames:
    def test_annotate(self, species_tree):
        tree = annotate_scientific_names(species_tree)
        sci_names = [leaf.props.get('sci_name') for leaf in tree.leaves()]
        assert 'Homo_sapiens' in sci_names
        assert 'Mus_musculus' in sci_names
        assert 'Danio_rerio' in sci_names

    def test_annotate_with_unnamed_leaf(self):
        tree = Tree('((:1,Homo_sapiens_gene1:1):1,Mus_musculus_gene1:1);', parser=1)
        tree = annotate_scientific_names(tree)
        unnamed = [leaf for leaf in tree.leaves() if not leaf.name][0]
        assert unnamed.props.get('sci_name') is None


class TestAnnotateDuplicationConfidenceScores:
    def test_annotate(self, species_tree):
        tree = annotate_scientific_names(species_tree)
        tree = annotate_duplication_confidence_scores(tree)
        for node in tree.traverse():
            if not node.is_leaf:
                assert 'dup_conf_score' in node.props
                assert 0 <= node.props.get('dup_conf_score') <= 1

    def test_no_duplication(self):
        # Tree where each species appears once in each child clade
        nwk = '((Homo_sapiens_G1:1,Mus_musculus_G1:1):1,(Danio_rerio_G1:1,Xenopus_laevis_G1:1):1);'
        tree = Tree(nwk, parser=1)
        tree = annotate_scientific_names(tree)
        tree = annotate_duplication_confidence_scores(tree)
        # Root node has no species overlap -> dup_conf_score = 0
        assert tree.props.get('dup_conf_score') == 0.0

    def test_full_duplication(self):
        # Tree where both children have the same species
        nwk = '((Homo_sapiens_G1:1,Mus_musculus_G1:1):1,(Homo_sapiens_G2:1,Mus_musculus_G2:1):1);'
        tree = Tree(nwk, parser=1)
        tree = annotate_scientific_names(tree)
        tree = annotate_duplication_confidence_scores(tree)
        assert tree.props.get('dup_conf_score') == 1.0

    def test_polytomy_internal_node_gets_zero_dup_conf_score(self):
        nwk = '(Homo_sapiens_G1:1,Homo_sapiens_G2:1,Mus_musculus_G1:1);'
        tree = Tree(nwk, parser=1)
        tree = annotate_scientific_names(tree)
        tree = annotate_duplication_confidence_scores(tree)
        assert tree.props.get('dup_conf_score') == 0.0


class TestGetSubtreeLeafNameSets:
    def test_sets_are_correct(self):
        tree = Tree('((A:1,B:1):1,(C:1,D:1):1);', parser=1)
        sets = get_subtree_leaf_name_sets(tree)
        assert sets[tree] == {'A', 'B', 'C', 'D'}
        assert sets[tree.common_ancestor(['A', 'B'])] == {'A', 'B'}
        assert sets[tree.common_ancestor(['C', 'D'])] == {'C', 'D'}
        for leaf in tree.leaves():
            assert sets[leaf] == {leaf.name}


class TestGetSubtreeLeafBitmasks:
    def test_bitmasks_are_correct(self):
        tree = Tree('((A:1,B:1):1,(C:1,D:1):1);', parser=1)
        leaf_name_to_bit = {'A': 0, 'B': 1, 'C': 2, 'D': 3}
        masks = get_subtree_leaf_bitmasks(tree, leaf_name_to_bit)
        assert masks[tree] == 0b1111
        assert masks[tree.common_ancestor(['A', 'B'])] == 0b0011
        assert masks[tree.common_ancestor(['C', 'D'])] == 0b1100
        for leaf in tree.leaves():
            assert masks[leaf] == (1 << leaf_name_to_bit[leaf.name])

    def test_missing_leaf_in_mapping_raises_clear_error(self):
        tree = Tree('((A:1,B:1):1,(C:1,D:1):1);', parser=1)
        leaf_name_to_bit = {'A': 0, 'B': 1, 'C': 2}
        with pytest.raises(ValueError, match='not found in reference mapping'):
            get_subtree_leaf_bitmasks(tree, leaf_name_to_bit)


class TestValidateUniqueNamedLeaves:
    def test_duplicate_leaf_names_raise(self):
        tree = Tree('((A:1,A:1):1,B:1);', parser=1)
        with pytest.raises(ValueError, match='Duplicated leaf labels'):
            validate_unique_named_leaves(tree, '--infile', " for 'transfer'")

    def test_empty_leaf_names_raise(self):
        tree = Tree('(A:1,:1,B:1);', parser=1)
        with pytest.raises(ValueError, match='Empty leaf labels'):
            validate_unique_named_leaves(tree, '--infile', " for 'transfer'")


class TestIsAllLeafNamesIdentical:
    def test_identical(self):
        t1 = Tree('((A,B),(C,D));', parser=1)
        t2 = Tree('((A,C),(B,D));', parser=1)
        assert is_all_leaf_names_identical(t1, t2) is True

    def test_not_identical(self):
        t1 = Tree('((A,B),(C,D));', parser=1)
        t2 = Tree('((A,B),(C,E));', parser=1)
        assert is_all_leaf_names_identical(t1, t2) is False

    def test_verbose_mode_handles_none_leaf_names(self):
        t1 = Tree('((A,B),(C,D));', parser=1)
        t2 = Tree('((A,B),(C,D));', parser=1)
        first_leaf = next(iter(t2.leaves()))
        first_leaf.name = None
        assert is_all_leaf_names_identical(t1, t2, verbose=True) is False


class TestGetTargetNodes:
    def test_all(self, simple_tree):
        nodes = get_target_nodes(simple_tree, 'all')
        assert len(nodes) == len(list(simple_tree.traverse()))

    def test_root(self, simple_tree):
        nodes = get_target_nodes(simple_tree, 'root')
        assert len(nodes) == 1
        assert nodes[0].is_root

    def test_leaf(self, simple_tree):
        nodes = get_target_nodes(simple_tree, 'leaf')
        assert all(n.is_leaf for n in nodes)
        assert len(nodes) == 4

    def test_intnode(self, simple_tree):
        nodes = get_target_nodes(simple_tree, 'intnode')
        assert all(not n.is_leaf for n in nodes)

    def test_invalid_target_raises(self, simple_tree):
        with pytest.raises(ValueError, match='Unknown target'):
            get_target_nodes(simple_tree, 'unknown_target')


class TestIsRooted:
    def test_rooted(self, simple_tree):
        assert is_rooted(simple_tree) is True

    def test_unrooted(self, unrooted_tree):
        assert is_rooted(unrooted_tree) is False

    def test_single_leaf_tree_is_rooted(self):
        tree = Tree('A;', parser=1)
        assert is_rooted(tree) is True

    def test_polytomy_root_tree_is_unrooted(self):
        tree = Tree('(A:1,B:1,C:1,D:1);', parser=1)
        assert is_rooted(tree) is False

import csv
import errno
import hashlib
import math
import os
import pickle
import re
import secrets
import shutil
import sqlite3
import subprocess
import sys
import tarfile
import tempfile
import time
from itertools import islice
from urllib.parse import quote
from collections import Counter, defaultdict
from collections.abc import Set
from contextlib import contextmanager
from io import StringIO
import ete4
import pandas as pd
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
from ete4 import Tree
from nwkit.fasta import parse_fasta, write_fasta
from nwkit.conventions import DEFAULT_TABLE_MISSING_VALUES
from nwkit.species_parser import (
    extract_parsed_species,
    get_species_parser,
)

NODENAME_PLACEHOLDER_PATTERN = re.compile(r'NODENAME_PLACEHOLDER\d{10}')
QUOTED_NODE_NAME_PATTERN = re.compile(r"'(?:[^']|'')*'")
NUMERIC_NODE_NAME_PATTERN = re.compile(r'[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?')
UNQUOTED_NODE_NAME_PATTERN = re.compile(r"^[^\s\(\)\[\]':;,]+$")
TREE_FORMAT_PROP = '_nwkit_parser_format'
MISSING_SUPPORT_VALUE = -999999.0
DOWNLOAD_LOCK_POLL_SECONDS = 1
DOWNLOAD_LOCK_TIMEOUT_SECONDS = 3600
LOCK_METADATA_GRACE_SECONDS = 2
ETE_TAXDUMP_URL = 'https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz'
ETE_TAXDUMP_MAX_BYTES = 512 * 1024 * 1024
ETE_TAXONOMY_DEFAULT_MAX_AGE_DAYS = 30.0
COMMON_ETE_CACHE_DIRS = (
    os.path.join(os.path.expanduser('~'), '.local', 'share', 'ete'),
    os.path.join(os.path.expanduser('~'), '.etetoolkit'),
)
def read_input_text(infile):
    if infile == '-':
        return sys.stdin.read()
    if os.path.isfile(infile):
        with open(infile) as handle:
            return handle.read()
    return str(infile)


def _read_raw_tsv_column_values(text, column_name):
    reader = csv.reader(StringIO(text), delimiter='\t')
    try:
        header = next(reader)
    except StopIteration:
        return None
    if column_name not in header:
        return None
    column_index = header.index(column_name)
    values = list()
    for row in reader:
        if (len(row) == 0) or all(cell == '' for cell in row):
            continue
        values.append(row[column_index] if column_index < len(row) else '')
    return values


def read_tsv_preserving_leaf_name(path):
    text = read_input_text(path)
    dataframe = pd.read_csv(StringIO(text), sep='\t', keep_default_na=False)
    if 'leaf_name' not in dataframe.columns:
        return dataframe
    raw_leaf_names = _read_raw_tsv_column_values(text, 'leaf_name')
    if raw_leaf_names is None:
        return dataframe
    if len(raw_leaf_names) != len(dataframe.index):
        raise ValueError(
            "Failed to preserve the 'leaf_name' column while parsing: {}".format(path)
        )
    dataframe = dataframe.copy()
    dataframe['leaf_name'] = raw_leaf_names
    return dataframe


def parse_table_missing_values(value=None):
    if value is None:
        return set(DEFAULT_TABLE_MISSING_VALUES)
    return {item.strip() for item in str(value).split(',')}


def is_missing_table_value(value, missing_values=None):
    if pd.isna(value):
        return True
    markers = parse_table_missing_values(missing_values) if not isinstance(missing_values, set) else missing_values
    return str(value).strip() in markers


def read_tip_table(path, option_name='--trait', tree_leaf_names=None, required_columns=(),
                   unmatched='warn', missing_values=None):
    dataframe = read_tsv_preserving_leaf_name(path)
    if 'leaf_name' not in dataframe.columns:
        raise ValueError("Column 'leaf_name' is required in '{}'.".format(option_name))
    missing_columns = [column for column in required_columns if column not in dataframe.columns]
    if missing_columns:
        raise ValueError(
            "Missing required column(s) in '{}': {}".format(option_name, ', '.join(missing_columns))
        )
    dataframe = dataframe.copy()
    dataframe['leaf_name'] = [str(leaf_name) for leaf_name in dataframe['leaf_name'].tolist()]
    if any(leaf_name.strip() == '' for leaf_name in dataframe['leaf_name'].tolist()):
        raise ValueError("Column 'leaf_name' in '{}' must not contain empty values.".format(option_name))
    duplicated = dataframe.loc[
        dataframe['leaf_name'].duplicated(keep=False), 'leaf_name'
    ].unique().tolist()
    if duplicated:
        raise ValueError(
            "Duplicated 'leaf_name' entries in '{}': {}".format(
                option_name,
                ', '.join(sorted(str(name) for name in duplicated)),
            )
        )
    table_only = list()
    tree_only = list()
    if tree_leaf_names is not None:
        tree_leaf_names = [str(name) for name in tree_leaf_names]
        tree_leaf_set = set(tree_leaf_names)
        table_leaf_set = set(dataframe['leaf_name'])
        table_only = sorted(table_leaf_set - tree_leaf_set)
        tree_only = sorted(tree_leaf_set - table_leaf_set)
        if unmatched not in ('warn', 'error', 'ignore'):
            raise ValueError("Unsupported '--unmatched' policy: {}".format(unmatched))
        if unmatched == 'error' and (table_only or tree_only):
            raise ValueError(
                "{} and tree tips differ (table-only={}; tree-only={}).".format(
                    option_name,
                    ','.join(table_only),
                    ','.join(tree_only),
                )
            )
        if unmatched == 'warn':
            if table_only:
                sys.stderr.write("Rows in {} not found in tree: {}\n".format(option_name, ' '.join(table_only)))
            if tree_only:
                sys.stderr.write("Tree tips not found in {}: {}\n".format(option_name, ' '.join(tree_only)))
    markers = parse_table_missing_values(missing_values)
    for column in dataframe.columns:
        if column == 'leaf_name':
            continue
        dataframe[column] = [
            pd.NA if is_missing_table_value(value, markers) else value
            for value in dataframe[column].tolist()
        ]
    return dataframe, table_only, tree_only


def count_set_bits(value):
    return int(value).bit_count()

def warn_cleanup_failure(resource_label, exc):
    sys.stderr.write('Warning: failed to clean up {}: {}\n'.format(resource_label, exc))


def validate_distinct_output_paths(outputs):
    """Reject multiple output roles that resolve to the same filesystem path."""
    by_path = defaultdict(list)
    for option_name, path in outputs:
        if path in (None, '', '-'):
            continue
        by_path[os.path.realpath(os.fspath(path))].append(str(option_name))
    collisions = [
        (path, option_names)
        for path, option_names in by_path.items()
        if len(option_names) > 1
    ]
    if collisions:
        details = '; '.join(
            '{} -> {}'.format(', '.join(option_names), path)
            for path, option_names in sorted(collisions)
        )
        raise ValueError('Output paths must be distinct: {}.'.format(details))

def resolve_download_dir(args=None):
    raw_dir = getattr(args, 'download_dir', 'auto') if args is not None else 'auto'
    if raw_dir is None:
        return None
    normalized = str(raw_dir).strip()
    if normalized.lower() in ['', 'auto']:
        return None
    if normalized.lower() == 'inferred':
        outfile = getattr(args, 'outfile', '-') if args is not None else '-'
        if outfile not in ['', None, '-']:
            inferred_base = os.path.dirname(os.path.realpath(outfile))
        else:
            raw_out_dir = getattr(args, 'out_dir', None) if args is not None else None
            inferred_base = os.path.realpath(raw_out_dir if raw_out_dir not in ['', None] else os.getcwd())
        return os.path.join(inferred_base, 'downloads')
    return os.path.realpath(normalized)

def resolve_ete_data_dir(args=None):
    download_dir = resolve_download_dir(args)
    if download_dir is None:
        return None
    return os.path.join(download_dir, 'ete4')

def _find_existing_ete_taxonomy_assets(exclude_dbfile=None):
    excluded_realpath = os.path.realpath(exclude_dbfile) if exclude_dbfile else None
    for cache_dir in COMMON_ETE_CACHE_DIRS:
        dbfile = os.path.join(cache_dir, 'taxa.sqlite')
        if not os.path.isfile(dbfile):
            continue
        if not _validate_ete_taxonomy_db(dbfile, require_traverse=True, full_check=True):
            continue
        if excluded_realpath is not None and os.path.realpath(dbfile) == excluded_realpath:
            continue
        assets = {'dbfile': dbfile, 'validated': True}
        traverse_file = dbfile + '.traverse.pkl'
        taxdump_file = os.path.join(cache_dir, 'taxdump.tar.gz')
        if os.path.isfile(traverse_file):
            assets['traverse_file'] = traverse_file
        if os.path.isfile(taxdump_file):
            assets['taxdump_file'] = taxdump_file
        return assets
    return None

def _seed_ete_taxonomy_assets(target_dbfile, source_assets):
    target_dir = os.path.dirname(target_dbfile)
    os.makedirs(target_dir, exist_ok=True)
    source_traverse_file = source_assets.get('traverse_file')
    if source_traverse_file is not None:
        _atomic_copy_file(source_traverse_file, target_dbfile + '.traverse.pkl')
    source_taxdump_file = source_assets.get('taxdump_file')
    if source_taxdump_file is not None:
        if _validate_ete_taxdump(source_taxdump_file):
            _atomic_copy_file(source_taxdump_file, os.path.join(target_dir, 'taxdump.tar.gz'))
    _atomic_copy_sqlite(source_assets['dbfile'], target_dbfile)

def _contains_quoted_node_names(newick_text):
    return QUOTED_NODE_NAME_PATTERN.search(str(newick_text)) is not None

def _contains_quoted_internal_node_names(newick_text):
    return re.search(r"\)\s*'(?:[^']|'')*'", str(newick_text)) is not None

def _is_numeric_node_name(name):
    if name in (None, ''):
        return False
    return NUMERIC_NODE_NAME_PATTERN.fullmatch(str(name).strip()) is not None

def _should_quote_node_name(name, is_internal=False):
    if name in (None, ''):
        return False
    name_text = str(name)
    if UNQUOTED_NODE_NAME_PATTERN.fullmatch(name_text) is None:
        return True
    if is_internal and _is_numeric_node_name(name_text):
        return True
    return False

def _serialize_newick_node_name(name, is_internal=False):
    name_text = str(name)
    if _should_quote_node_name(name_text, is_internal=is_internal):
        return "'{}'".format(name_text.replace("'", "''"))
    return name_text

def _is_missing_support_value(support):
    if support is None:
        return False
    try:
        return abs(float(support) - MISSING_SUPPORT_VALUE) < 10 ** -9
    except (TypeError, ValueError):
        return False

def _named_internal_nodes(tree):
    return [
        node for node in tree.traverse()
        if (not node.is_leaf) and (str(node.name or '').strip() != '')
    ]

def _auto_format_ambiguity_message():
    return (
        'Ambiguous tree format: unquoted numeric internal labels can be interpreted as '
        'support values (format 0) or internal node names (format 1). '
        'Use --format 0, --format 1, or quote internal node names to disambiguate.'
    )

def _read_tree_auto(infile):
    quoted_internal_names = _contains_quoted_internal_node_names(infile)
    parsed = dict()
    for candidate_format in [1, 0]:
        try:
            parsed[candidate_format] = Tree(infile, parser=candidate_format)
        except Exception:
            parsed[candidate_format] = None
    parser1_tree = parsed[1]
    parser0_tree = parsed[0]
    ambiguity_message = None
    if parser1_tree is not None:
        parser1_named_internal_nodes = _named_internal_nodes(parser1_tree)
        if parser1_named_internal_nodes and (
            quoted_internal_names or
            any(not _is_numeric_node_name(node.name) for node in parser1_named_internal_nodes)
        ):
            return 1, parser1_tree, None
        if (
            parser0_tree is not None and
            parser1_named_internal_nodes and
            all(_is_numeric_node_name(node.name) for node in parser1_named_internal_nodes)
        ):
            ambiguity_message = _auto_format_ambiguity_message()
    if parser0_tree is not None:
        return 0, parser0_tree, ambiguity_message
    if parser1_tree is not None:
        return 1, parser1_tree, None
    for candidate_format in [2, 3, 4, 5, 6, 7, 8, 9, 100]:
        try:
            return candidate_format, Tree(infile, parser=candidate_format), None
        except Exception:
            pass
    raise Exception('Failed to parse the input tree.')

def _assert_lock_path_is_regular_file(lock_path, lock_label='Lock'):
    if not os.path.lexists(lock_path):
        return
    if os.path.islink(lock_path) or (not os.path.isfile(lock_path)):
        raise IsADirectoryError('{} path exists but is not a file: {}'.format(lock_label, lock_path))

def _try_create_lock_file(lock_path):
    _assert_lock_path_is_regular_file(lock_path)
    try:
        fd = os.open(lock_path, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
    except FileExistsError:
        return None
    token = secrets.token_hex(16)
    try:
        os.write(fd, '{} {}\n'.format(os.getpid(), token).encode('ascii'))
        os.fsync(fd)
    finally:
        os.close(fd)
    stat_result = os.stat(lock_path)
    return token, stat_result.st_dev, stat_result.st_ino

def _read_lock_owner_pid(lock_path):
    try:
        with open(lock_path) as lock_handle:
            first_line = lock_handle.readline().strip()
    except OSError:
        return None
    if first_line == '':
        return None
    try:
        pid = int(first_line.split()[0])
    except ValueError:
        return None
    if pid <= 0:
        return None
    return pid

def _is_process_alive(pid):
    try:
        os.kill(int(pid), 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    except OSError as exc:
        if getattr(exc, 'errno', None) == errno.ESRCH:
            return False
        return True
    return True

def _break_stale_lock_if_needed(lock_path, lock_label='Lock'):
    if not os.path.lexists(lock_path):
        return False
    _assert_lock_path_is_regular_file(lock_path, lock_label=lock_label)
    try:
        stat_before = os.stat(lock_path)
    except FileNotFoundError:
        return False
    owner_pid = _read_lock_owner_pid(lock_path)
    if owner_pid is None:
        age_seconds = max(0.0, time.time() - stat_before.st_mtime)
        if age_seconds < LOCK_METADATA_GRACE_SECONDS:
            return False
        stale_reason = 'missing/invalid owner PID'
    elif _is_process_alive(owner_pid):
        return False
    else:
        stale_reason = 'owner PID {} is not running'.format(owner_pid)
    try:
        stat_now = os.stat(lock_path)
    except FileNotFoundError:
        return False
    if (
        (stat_before.st_ino != stat_now.st_ino)
        or (stat_before.st_size != stat_now.st_size)
        or (stat_before.st_mtime_ns != stat_now.st_mtime_ns)
    ):
        return False
    try:
        os.remove(lock_path)
    except FileNotFoundError:
        return False
    sys.stderr.write('Removed stale {} lock: {} ({})\n'.format(lock_label, lock_path, stale_reason))
    return True


def _remove_owned_lock(lock_path, owner_record, lock_label='Lock'):
    owner_token, owner_device, owner_inode = owner_record
    try:
        stat_before = os.stat(lock_path)
        with open(lock_path) as lock_handle:
            fields = lock_handle.readline().strip().split()
    except FileNotFoundError:
        return
    _assert_lock_path_is_regular_file(lock_path, lock_label=lock_label)
    if (
        stat_before.st_dev != owner_device
        or stat_before.st_ino != owner_inode
        or len(fields) < 2
        or fields[1] != owner_token
    ):
        sys.stderr.write(
            'Warning: {} ownership changed before release; leaving lock in place: {}\n'.format(
                lock_label, lock_path
            )
        )
        return
    os.remove(lock_path)

@contextmanager
def acquire_exclusive_lock(
    lock_path,
    lock_label='Lock',
    poll_seconds=DOWNLOAD_LOCK_POLL_SECONDS,
    timeout_seconds=DOWNLOAD_LOCK_TIMEOUT_SECONDS,
):
    poll_seconds = int(poll_seconds)
    timeout_seconds = int(timeout_seconds)
    if poll_seconds <= 0:
        raise ValueError('poll_seconds must be > 0.')
    if timeout_seconds <= 0:
        raise ValueError('timeout_seconds must be > 0.')
    lock_path = os.path.realpath(lock_path)
    lock_dir = os.path.dirname(lock_path)
    if lock_dir != '':
        if os.path.exists(lock_dir) and (not os.path.isdir(lock_dir)):
            raise NotADirectoryError('Lock parent path exists but is not a directory: {}'.format(lock_dir))
        os.makedirs(lock_dir, exist_ok=True)
    wait_start = time.time()
    has_reported_wait = False
    while True:
        owner_record = _try_create_lock_file(lock_path)
        if owner_record is not None:
            try:
                yield
            finally:
                if os.path.lexists(lock_path):
                    _remove_owned_lock(
                        lock_path,
                        owner_record=owner_record,
                        lock_label=lock_label,
                    )
            return
        _assert_lock_path_is_regular_file(lock_path, lock_label=lock_label)
        if _break_stale_lock_if_needed(lock_path=lock_path, lock_label=lock_label):
            continue
        elapsed = time.time() - wait_start
        if not has_reported_wait:
            sys.stderr.write('Another process holds {}. Waiting for lock release: {}\n'.format(lock_label, lock_path))
            has_reported_wait = True
        if elapsed > timeout_seconds:
            raise TimeoutError(
                'Timed out after {:,} sec waiting for {} lock: {}'.format(
                    timeout_seconds,
                    lock_label,
                    lock_path,
                )
            )
        time.sleep(poll_seconds)

def _validate_ete_taxdump(taxdump_file):
    if not os.path.isfile(taxdump_file):
        return False
    try:
        if os.path.getsize(taxdump_file) <= 0 or os.path.getsize(taxdump_file) > ETE_TAXDUMP_MAX_BYTES:
            return False
        with tarfile.open(taxdump_file, mode='r:gz') as archive:
            members = {member.name: member for member in archive.getmembers()}
            for required_name in ('nodes.dmp', 'names.dmp', 'merged.dmp'):
                member = members.get(required_name)
                if member is None or not member.isfile() or member.size <= 0:
                    return False
    except (OSError, tarfile.TarError):
        return False
    return True


def _validate_ete_taxonomy_db(dbfile, require_traverse=False, full_check=False):
    if not os.path.isfile(dbfile) or os.path.getsize(dbfile) <= 0:
        return False
    try:
        uri = 'file:{}?mode=ro'.format(quote(os.path.abspath(dbfile), safe='/'))
        connection = sqlite3.connect(uri, uri=True)
        try:
            if full_check and connection.execute('PRAGMA quick_check').fetchone() != ('ok',):
                return False
            table_names = {
                row[0]
                for row in connection.execute("SELECT name FROM sqlite_master WHERE type='table'")
            }
            if not {'stats', 'species', 'synonym', 'merged'}.issubset(table_names):
                return False
            if connection.execute('SELECT version FROM stats').fetchone() is None:
                return False
            if connection.execute('SELECT 1 FROM species LIMIT 1').fetchone() is None:
                return False
        finally:
            connection.close()
    except (OSError, sqlite3.DatabaseError):
        return False
    if require_traverse:
        traverse_file = dbfile + '.traverse.pkl'
        try:
            with open(traverse_file, 'rb') as handle:
                traversal = pickle.load(handle)
            if not isinstance(traversal, list) or not traversal:
                return False
        except (OSError, EOFError, pickle.PickleError, ValueError, TypeError):
            return False
    return True


def _atomic_copy_file(source, target):
    target_dir = os.path.dirname(target)
    os.makedirs(target_dir, exist_ok=True)
    fd, temporary = tempfile.mkstemp(prefix='.{}.'.format(os.path.basename(target)), dir=target_dir)
    os.close(fd)
    try:
        shutil.copy2(source, temporary)
        os.replace(temporary, target)
    finally:
        if os.path.exists(temporary):
            os.remove(temporary)


def _atomic_copy_sqlite(source, target):
    target_dir = os.path.dirname(target)
    os.makedirs(target_dir, exist_ok=True)
    fd, temporary = tempfile.mkstemp(prefix='.{}.'.format(os.path.basename(target)), dir=target_dir)
    os.close(fd)
    source_uri = 'file:{}?mode=ro'.format(quote(os.path.abspath(source), safe='/'))
    try:
        source_connection = sqlite3.connect(source_uri, uri=True)
        target_connection = sqlite3.connect(temporary)
        try:
            source_connection.backup(target_connection)
        finally:
            target_connection.close()
            source_connection.close()
        os.replace(temporary, target)
    finally:
        if os.path.exists(temporary):
            os.remove(temporary)


def _download_ete_taxdump(taxdump_file):
    target_dir = os.path.dirname(taxdump_file)
    os.makedirs(target_dir, exist_ok=True)
    session = requests.Session()
    retry = Retry(
        total=3,
        connect=3,
        read=3,
        status=3,
        backoff_factor=0.5,
        status_forcelist=(429, 500, 502, 503, 504),
        allowed_methods=frozenset(('GET',)),
    )
    adapter = HTTPAdapter(max_retries=retry)
    session.mount('https://', adapter)
    session.mount('http://', adapter)
    try:
        with session.get(ETE_TAXDUMP_URL + '.md5', timeout=(10, 30)) as checksum_response:
            checksum_response.raise_for_status()
            checksum_text = checksum_response.text
        if len(checksum_text.encode('utf-8')) > 4096:
            raise ValueError('Unexpected NCBI taxonomy checksum response.')
        checksum_tokens = checksum_text.split()
        if not checksum_tokens:
            raise ValueError('Unexpected NCBI taxonomy checksum response.')
        expected_md5 = checksum_tokens[0].lower()
        if re.fullmatch(r'[0-9a-f]{32}', expected_md5) is None:
            raise ValueError('Unexpected NCBI taxonomy checksum response.')
        if _validate_ete_taxdump(taxdump_file):
            local_md5 = hashlib.md5()
            with open(taxdump_file, 'rb') as handle:
                for chunk in iter(lambda: handle.read(1024 * 1024), b''):
                    local_md5.update(chunk)
            if local_md5.hexdigest() == expected_md5:
                return False
        fd, temporary = tempfile.mkstemp(prefix='.taxdump.', suffix='.tar.gz', dir=target_dir)
        os.close(fd)
        try:
            downloaded = 0
            digest = hashlib.md5()
            with session.get(ETE_TAXDUMP_URL, stream=True, timeout=(10, 120)) as response:
                response.raise_for_status()
                content_length = response.headers.get('Content-Length')
                if content_length is not None and int(content_length) > ETE_TAXDUMP_MAX_BYTES:
                    raise ValueError('NCBI taxonomy archive exceeds the download size limit.')
                with open(temporary, 'wb') as handle:
                    for chunk in response.iter_content(chunk_size=1024 * 1024):
                        if not chunk:
                            continue
                        downloaded += len(chunk)
                        if downloaded > ETE_TAXDUMP_MAX_BYTES:
                            raise ValueError('NCBI taxonomy archive exceeds the download size limit.')
                        digest.update(chunk)
                        handle.write(chunk)
                    handle.flush()
                    os.fsync(handle.fileno())
            if digest.hexdigest() != expected_md5:
                raise ValueError('NCBI taxonomy archive checksum verification failed.')
            if not _validate_ete_taxdump(temporary):
                raise ValueError('Downloaded NCBI taxonomy archive is invalid.')
            os.replace(temporary, taxdump_file)
        finally:
            if os.path.exists(temporary):
                os.remove(temporary)
    finally:
        session.close()
    return True


def _build_ete_taxonomy_database(dbfile, taxdump_file):
    target_dir = os.path.dirname(dbfile)
    with tempfile.TemporaryDirectory(prefix='.taxonomy-build-', dir=target_dir) as work_dir:
        temporary_db = os.path.join(work_dir, 'taxa.sqlite')
        script = (
            'import sys; '
            'from ete4.ncbi_taxonomy.ncbiquery import update_db; '
            'update_db(sys.argv[1], sys.argv[2])'
        )
        try:
            subprocess.run(
                [sys.executable, '-c', script, temporary_db, os.path.realpath(taxdump_file)],
                cwd=work_dir,
                check=True,
            )
        except subprocess.CalledProcessError as exc:
            raise RuntimeError('Failed to build the ETE4 taxonomy database.') from exc
        if not _validate_ete_taxonomy_db(temporary_db, require_traverse=True, full_check=True):
            raise ValueError('ETE4 produced an invalid taxonomy database.')
        os.replace(temporary_db + '.traverse.pkl', dbfile + '.traverse.pkl')
        os.replace(temporary_db, dbfile)


def _taxonomy_cache_max_age_seconds(args):
    value = getattr(args, 'taxonomy_cache_max_age_days', ETE_TAXONOMY_DEFAULT_MAX_AGE_DAYS)
    try:
        days = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError("'--taxonomy-cache-max-age-days' must be numeric.") from exc
    if not math.isfinite(days) or days < 0:
        raise ValueError("'--taxonomy-cache-max-age-days' must be finite and >= 0.")
    return days * 24 * 60 * 60


def _taxonomy_cache_is_stale(dbfile, args):
    if bool(getattr(args, 'refresh_taxonomy_cache', False)):
        return True
    return time.time() - os.path.getmtime(dbfile) > _taxonomy_cache_max_age_seconds(args)

def get_ete_ncbitaxa(args=None):
    ete_data_dir = resolve_ete_data_dir(args)
    if ete_data_dir is None:
        from ete4.ncbi_taxonomy.ncbiquery import DEFAULT_TAXADB

        dbfile = os.path.realpath(DEFAULT_TAXADB)
        ete_data_dir = os.path.dirname(dbfile)
    else:
        dbfile = os.path.join(ete_data_dir, 'taxa.sqlite')
    os.makedirs(ete_data_dir, exist_ok=True)
    taxdump_file = os.path.join(ete_data_dir, 'taxdump.tar.gz')
    lock_path = os.path.join(ete_data_dir, '.ete4_taxonomy.lock')
    with acquire_exclusive_lock(lock_path=lock_path, lock_label='ETE4 taxonomy DB'):
        os.makedirs(ete_data_dir, exist_ok=True)
        db_valid = _validate_ete_taxonomy_db(dbfile, require_traverse=True)
        if db_valid and not _taxonomy_cache_is_stale(dbfile, args):
            return ete4.NCBITaxa(dbfile=dbfile, update=False)
        existing_assets = _find_existing_ete_taxonomy_assets(exclude_dbfile=dbfile)
        if (
            not db_valid
            and existing_assets is not None
            and (
                existing_assets.get('validated')
                or _validate_ete_taxonomy_db(
                    existing_assets['dbfile'],
                    require_traverse=True,
                    full_check=True,
                )
            )
            and not _taxonomy_cache_is_stale(existing_assets['dbfile'], args)
        ):
            _seed_ete_taxonomy_assets(target_dbfile=dbfile, source_assets=existing_assets)
            return ete4.NCBITaxa(dbfile=dbfile, update=False)
        archive_changed = _download_ete_taxdump(taxdump_file)
        if db_valid and not archive_changed:
            os.utime(dbfile, None)
            return ete4.NCBITaxa(dbfile=dbfile, update=False)
        _build_ete_taxonomy_database(dbfile=dbfile, taxdump_file=taxdump_file)
        return ete4.NCBITaxa(dbfile=dbfile, update=False)

def read_tree(infile, format, quoted_node_names, quiet=False):
    global INFILE_FORMAT
    infile = read_input_text(infile).strip()
    if infile == '':
        raise Exception('Failed to parse the input tree.')
    if (not quoted_node_names) and _contains_quoted_node_names(infile):
        raise ValueError('Quoted node names were found in the input tree. Re-run with --quoted-node-names yes.')
    if format in ('auto', 'auto-strict'):
        format_original = format
        format, tree, ambiguity_message = _read_tree_auto(infile)
        if ambiguity_message is not None:
            if format_original == 'auto-strict':
                raise ValueError(ambiguity_message)
            if not quiet:
                sys.stderr.write('Warning: {}\n'.format(ambiguity_message))
        INFILE_FORMAT = format
    else:
        format_original = format
        format = int(format)
        tree = Tree(infile, parser=format)
    INFILE_FORMAT = int(format)
    # Keep per-node format metadata so writing a subtree keeps the original parser.
    for node in tree.traverse():
        node.props[TREE_FORMAT_PROP] = INFILE_FORMAT
    if format==0: # flexible with support values
        # Single pass over nodes to detect max support and collect candidates.
        # If max support is > 1.0, treat non-root None/1.0 as missing support.
        max_support = None
        missing_support_candidates = list()
        for node in tree.traverse():
            support = node.support
            if support is not None:
                if (max_support is None) or (support > max_support):
                    max_support = support
            if not node.is_root:
                if support is None or (abs(float(support) - 1.0) < 10**-9): # 1.0 is default for missing support values
                    missing_support_candidates.append(node)
        if (max_support is not None) and (max_support > 1.0):
            for node in missing_support_candidates:
                node.support = -999999
    if format==1: # flexible with internal node names
        for node in tree.traverse():
            if (
                (not node.is_root) and
                ('support' not in node.props or node.support is None)
            ):  # Don't set support on root (breaks ete4 set_outgroup)
                node.support = -999999
    if not quiet:
        num_leaves = len(list(tree.leaves()))
        txt = 'Number of leaves in input tree = {:,}, Input tree format = {}\n'
        sys.stderr.write(txt.format(num_leaves, format))
    return tree

def write_tree(tree, args, format, quiet=False, props=None):
    if format=='auto':
        format = tree.props.get(TREE_FORMAT_PROP, globals().get('INFILE_FORMAT', 0))
        format = int(format)
    else:
        format = int(format)
    if not quiet:
        num_leaves = len(list(tree.leaves()))
        txt = 'Number of leaves in output tree = {:,}, Output tree format = {}\n'
        sys.stderr.write(txt.format(num_leaves, format))
    node_name_dict = dict()
    original_node_names = list()
    original_support_values = list()
    i = 0
    for node in tree.traverse():
        original_node_names.append((node, node.name))
        if (node.name is not None) and (str(node.name) != ''):
            placeholder_name = 'NODENAME_PLACEHOLDER'+str(i).zfill(10)
            node_name_dict[placeholder_name] = _serialize_newick_node_name(
                node.name,
                is_internal=(not node.is_leaf),
            )
            node.name = placeholder_name
            i += 1
        if _is_missing_support_value(node.support):
            original_support_values.append((node, node.support))
            node.support = None
    try:
        if props is None:
            props = getattr(args, 'output_properties', None)
        if props is not None:
            props = sorted(set(str(prop) for prop in props if prop != TREE_FORMAT_PROP))
        write_kwargs = {
            'parser': format,
            'format_root_node': True,
        }
        if props:
            write_kwargs['props'] = props
        tree_str = tree.write(**write_kwargs)
    finally:
        for node, node_name in original_node_names:
            node.name = node_name
        for node, support in original_support_values:
            node.support = support
    if tree_str.endswith(':;'):
        tree_str = tree_str[:-2]+';'
    if node_name_dict:
        tree_str = NODENAME_PLACEHOLDER_PATTERN.sub(lambda m: str(node_name_dict[m.group(0)]), tree_str)
    if args.outfile=='-':
        print(tree_str)
    else:
        with open(args.outfile, mode='w') as f:
            f.write(tree_str)

def split_newick_stream(newick_text):
    trees = list()
    buffer = list()
    in_quote = False
    text = str(newick_text)
    i = 0
    while i < len(text):
        char = text[i]
        buffer.append(char)
        if char == "'":
            if in_quote and (i + 1 < len(text)) and (text[i + 1] == "'"):
                buffer.append(text[i + 1])
                i += 1
            else:
                in_quote = not in_quote
        elif (char == ';') and (not in_quote):
            tree_text = ''.join(buffer).strip()
            if tree_text != '':
                trees.append(tree_text)
            buffer = list()
        i += 1
    if ''.join(buffer).strip() != '':
        raise ValueError('Input tree collection ended before a terminal semicolon.')
    return trees


def iter_newick_stream(handle, chunk_size=1024 * 1024):
    """Yield semicolon-terminated Newick records without loading the stream at once."""
    buffer = list()
    in_quote = False
    possible_closing_quote = False
    while True:
        chunk = handle.read(chunk_size)
        if chunk == '':
            break
        for char in chunk:
            if possible_closing_quote:
                if char == "'":
                    buffer.append(char)
                    possible_closing_quote = False
                    continue
                in_quote = False
                possible_closing_quote = False
            buffer.append(char)
            if char == "'":
                if in_quote:
                    possible_closing_quote = True
                else:
                    in_quote = True
            elif char == ';' and not in_quote:
                tree_text = ''.join(buffer).strip()
                if tree_text:
                    yield tree_text
                buffer = list()
    if possible_closing_quote:
        in_quote = False
    if in_quote or ''.join(buffer).strip():
        raise ValueError('Input tree collection ended before a terminal semicolon.')


def iter_tree_strings(infile):
    if infile == '-':
        yield from iter_newick_stream(sys.stdin)
        return
    if os.path.isfile(infile):
        with open(infile) as handle:
            yield from iter_newick_stream(handle)
        return
    yield from split_newick_stream(str(infile))

def read_trees(infile, format, quoted_node_names, quiet=False):
    tree_strings = read_tree_strings(infile)
    if len(tree_strings) == 0:
        raise Exception('Failed to parse the input trees.')
    trees = [read_tree(tree_string, format, quoted_node_names, quiet=True) for tree_string in tree_strings]
    if not quiet:
        sys.stderr.write('Number of input trees = {:,}\n'.format(len(trees)))
    return trees

def read_tree_strings(infile):
    tree_text = read_input_text(infile)
    return split_newick_stream(tree_text)

def inspect_tree_text(newick_text, format='auto', quoted_node_names=True):
    text = str(newick_text).strip()
    if text == '':
        return {
            'parse_ok': False,
            'parse_error': 'Failed to parse the input tree.',
            'input_format': '',
            'format_ambiguous': False,
            'has_quoted_node_names': False,
            'has_quoted_internal_node_names': False,
        }
    has_quoted_node_names = _contains_quoted_node_names(text)
    has_quoted_internal_node_names = _contains_quoted_internal_node_names(text)
    if (not quoted_node_names) and has_quoted_node_names:
        return {
            'parse_ok': False,
            'parse_error': 'Quoted node names were found in the input tree. Re-run with --quoted-node-names yes.',
            'input_format': '',
            'format_ambiguous': False,
            'has_quoted_node_names': has_quoted_node_names,
            'has_quoted_internal_node_names': has_quoted_internal_node_names,
        }
    try:
        if format in ('auto', 'auto-strict'):
            inferred_format, _, ambiguity_message = _read_tree_auto(text)
            input_format = inferred_format
            format_ambiguous = ambiguity_message is not None
            if format == 'auto-strict' and format_ambiguous:
                return {
                    'parse_ok': False,
                    'parse_error': ambiguity_message,
                    'input_format': '',
                    'format_ambiguous': True,
                    'has_quoted_node_names': has_quoted_node_names,
                    'has_quoted_internal_node_names': has_quoted_internal_node_names,
                }
        else:
            Tree(text, parser=int(format))
            input_format = int(format)
            format_ambiguous = False
    except Exception as exc:
        return {
            'parse_ok': False,
            'parse_error': str(exc),
            'input_format': '',
            'format_ambiguous': False,
            'has_quoted_node_names': has_quoted_node_names,
            'has_quoted_internal_node_names': has_quoted_internal_node_names,
        }
    return {
        'parse_ok': True,
        'parse_error': '',
        'input_format': input_format,
        'format_ambiguous': format_ambiguous,
        'has_quoted_node_names': has_quoted_node_names,
        'has_quoted_internal_node_names': has_quoted_internal_node_names,
    }

def assign_branch_ids(tree, strategy='levelorder'):
    node_to_branch_id = dict()
    for branch_id, node in enumerate(tree.traverse(strategy=strategy)):
        node_to_branch_id[node] = branch_id
    return node_to_branch_id


def get_node_class(node):
    if node.is_root:
        return 'root'
    if node.is_leaf:
        return 'leaf'
    return 'intnode'

def support_is_missing(support):
    if support is None:
        return True
    try:
        support_value = float(support)
    except (TypeError, ValueError):
        return False
    if math.isnan(support_value):
        return True
    return abs(support_value - MISSING_SUPPORT_VALUE) < 10 ** -9


def get_tree_property_names(tree):
    reserved = {TREE_FORMAT_PROP, 'name', 'dist', 'support'}
    return {
        str(prop)
        for node in tree.traverse()
        for prop in node.props
        if str(prop) not in reserved
    }

def compute_node_ages(tree, tolerance=10 ** -9):
    age_by_node = dict()
    for node in tree.traverse(strategy='postorder'):
        if node.is_leaf:
            age_by_node[node] = 0.0
            continue
        child_ages = list()
        for child in node.get_children():
            child_dist = 0.0 if (child.dist is None) else float(child.dist)
            child_ages.append(age_by_node[child] + child_dist)
        if len(child_ages) == 0:
            age_by_node[node] = 0.0
            continue
        if max(child_ages) - min(child_ages) > tolerance:
            raise ValueError('Tree must be ultrametric.')
        age_by_node[node] = child_ages[0]
    return age_by_node

def read_seqs(seqfile, seqformat, quiet):
    if str(seqformat).lower() != 'fasta':
        raise ValueError("Unsupported sequence format '{}'. Only 'fasta' is supported.".format(seqformat))
    if seqfile=='-':
        records = parse_fasta(sys.stdin)
    else:
        with open(seqfile, newline='') as fh:
            records = parse_fasta(fh)
    if not quiet:
        sys.stderr.write('Number of input sequences: {:,}\n'.format(len(records)))
    return records

def write_seqs(records, outfile, seqformat='fasta', quiet=False):
    if str(seqformat).lower() != 'fasta':
        raise ValueError("Unsupported sequence format '{}'. Only 'fasta' is supported.".format(seqformat))
    if not quiet:
        sys.stderr.write('Number of output sequences: {:,}\n'.format(len(records)))
    if outfile=='-':
        write_fasta(records, sys.stdout, normalize_newlines=True)
    else:
        with open(outfile, 'w', newline='') as fh:
            write_fasta(records, fh)

def remove_singleton(tree, verbose=False, preserve_branch_length=True):
    for node in tree.traverse():
        if node.is_leaf:
            continue
        num_children = len(node.get_children())
        if (num_children>1):
            continue
        if verbose:
            sys.stderr.write('Deleting a singleton node: {}\n'.format(node.name))
        node.delete(prevent_nondicotomic=False, preserve_branch_length=preserve_branch_length)
    # ete4 does not always collapse a singleton root in-place.
    # Explicitly promote the root child while more than one tip exists.
    while (len(tree.get_children()) == 1) and (len(list(tree.leaves())) > 1):
        child = tree.get_children()[0]
        if child.is_leaf:
            break
        if verbose:
            sys.stderr.write('Deleting a singleton node: {}\n'.format(tree.name))
        root_dist = tree.dist
        child_props = {k: v for k, v in child.props.items() if k != 'dist'}
        child_name = child.name
        tree.remove_child(child)
        for grandchild in list(child.get_children()):
            tree.add_child(grandchild)
        tree.dist = root_dist
        tree.name = child_name
        tree.props.clear()
        tree.props.update(child_props)
    return tree

def label2sciname(labels, in_delim='_', out_delim='_', args=None, species_parser=None, species_regex=None, species_map_tsv=None):
    if labels is None:
        return None
    is_str_input = isinstance(labels, str)
    if is_str_input:
        labels = [labels,]
    scinames = list()
    use_species_parser = any([
        args is not None,
        species_parser is not None,
        species_regex is not None,
        species_map_tsv is not None,
        in_delim == '_',
    ])
    for label in labels:
        if label is None:
            sciname = None
        elif use_species_parser:
            parsed_species = extract_parsed_species(
                label,
                args=args,
                species_parser=species_parser,
                species_regex=species_regex,
                species_map_tsv=species_map_tsv,
            )
            sciname = parsed_species.species_label
            if (sciname is not None) and (out_delim != '_'):
                sciname = sciname.replace('_', out_delim)
        else:
            label_str = str(label)
            splitted = label_str.split(in_delim)
            if len(splitted)>=2:
                sciname = splitted[0]+out_delim+splitted[1]
            else:
                sciname = None
        scinames.append(sciname)
    if is_str_input:
        scinames = scinames[0]
    return scinames

def extract_species_label(label, species_regex=None, out_delim='_', args=None, species_parser=None, species_map_tsv=None):
    parsed_species = extract_parsed_species(
        label,
        args=args,
        species_parser=species_parser,
        species_regex=species_regex,
        species_map_tsv=species_map_tsv,
    )
    species_label = parsed_species.species_label
    if (species_label is None) or (out_delim == '_'):
        return species_label
    return species_label.replace('_', out_delim)

def extract_taxonomy_query(label, out_delim=' ', args=None, species_parser=None, species_regex=None, species_map_tsv=None):
    parsed_species = extract_parsed_species(
        label,
        args=args,
        species_parser=species_parser,
        species_regex=species_regex,
        species_map_tsv=species_map_tsv,
    )
    taxonomy_query = parsed_species.taxonomy_query
    if taxonomy_query is None:
        return None
    if out_delim == ' ':
        return taxonomy_query
    return taxonomy_query.replace(' ', out_delim)

def get_species_group_records(tree, option_name='--infile', context='', args=None, species_parser=None, species_regex=None, species_map_tsv=None):
    leaf_name_to_species_label = dict()
    species_label_to_leaf_names = defaultdict(list)
    species_label_to_taxonomy_query = dict()
    unresolved_leaf_names = list()
    parser = get_species_parser(
        args=args,
        species_parser=species_parser,
        species_regex=species_regex,
        species_map_tsv=species_map_tsv,
    )
    for leaf in tree.leaves():
        parsed_species = parser.parse(leaf.name)
        species_label = parsed_species.species_label
        if species_label is None:
            unresolved_leaf_names.append(str(leaf.name))
            continue
        leaf_name_to_species_label[leaf.name] = species_label
        species_label_to_leaf_names[species_label].append(leaf.name)
        taxonomy_query = parsed_species.taxonomy_query
        previous_query = species_label_to_taxonomy_query.get(species_label)
        if previous_query is None:
            species_label_to_taxonomy_query[species_label] = taxonomy_query
        elif (taxonomy_query is not None) and (previous_query != taxonomy_query):
            raise ValueError(
                "Parsed species label maps to multiple taxonomy queries in '{}'{}: {}".format(
                    option_name,
                    context,
                    species_label,
                )
            )
    if unresolved_leaf_names:
        raise ValueError(
            "Leaf labels could not be parsed as species labels in '{}'{}: {}".format(
                option_name,
                context,
                ', '.join(sorted(unresolved_leaf_names)),
            )
        )
    for species_label, leaf_names in species_label_to_leaf_names.items():
        if len(leaf_names) <= 1:
            continue
        mrca = tree.common_ancestor(leaf_names)
        if set(mrca.leaf_names()) != set(leaf_names):
            raise ValueError(
                "Leaf labels for the same species are not monophyletic in '{}'{}: {}".format(
                    option_name,
                    context,
                    species_label,
                )
            )
    for species_label in species_label_to_leaf_names.keys():
        if species_label_to_taxonomy_query.get(species_label) is None:
            species_label_to_taxonomy_query[species_label] = species_label.replace('_', ' ')
    return (
        leaf_name_to_species_label,
        dict(species_label_to_leaf_names),
        species_label_to_taxonomy_query,
    )

def get_monophyletic_species_groups(tree, option_name='--infile', context='', args=None, species_parser=None, species_regex=None, species_map_tsv=None):
    leaf_name_to_species_label, species_label_to_leaf_names, _ = get_species_group_records(
        tree,
        option_name=option_name,
        context=context,
        args=args,
        species_parser=species_parser,
        species_regex=species_regex,
        species_map_tsv=species_map_tsv,
    )
    return leaf_name_to_species_label, species_label_to_leaf_names

class _LeafIntervalSet(Set):
    """Read-only set view over a contiguous DFS leaf interval."""
    __slots__ = ('_leaf_names', '_name_to_index', '_start', '_end')

    def __init__(self, leaf_names, name_to_index, start, end):
        self._leaf_names = leaf_names
        self._name_to_index = name_to_index
        self._start = start
        self._end = end

    def __contains__(self, value):
        index = self._name_to_index.get(value)
        return index is not None and self._start <= index < self._end

    def __iter__(self):
        return islice(self._leaf_names, self._start, self._end)

    def __len__(self):
        return self._end - self._start

    def __and__(self, other):
        return frozenset(value for value in self if value in other)

    def __rand__(self, other):
        return self.__and__(other)

    def __sub__(self, other):
        return frozenset(value for value in self if value not in other)

    def __repr__(self):
        return repr(set(self))


def get_subtree_leaf_name_sets(tree):
    leaf_nodes = list(tree.leaves())
    leaf_names = tuple(node.name for node in leaf_nodes)
    if len(set(leaf_names)) != len(leaf_names):
        subtree_leaf_name_sets = dict()
        for node in tree.traverse(strategy='postorder'):
            if node.is_leaf:
                subtree_leaf_name_sets[node] = {node.name}
            else:
                leaf_set = set()
                for child in node.get_children():
                    leaf_set.update(subtree_leaf_name_sets[child])
                subtree_leaf_name_sets[node] = leaf_set
        return subtree_leaf_name_sets
    name_to_index = {name: index for index, name in enumerate(leaf_names)}
    intervals = dict()
    for node in tree.traverse(strategy='postorder'):
        if node.is_leaf:
            index = name_to_index[node.name]
            intervals[node] = (index, index + 1)
        else:
            child_intervals = [intervals[child] for child in node.get_children()]
            intervals[node] = (
                min(start for start, _ in child_intervals),
                max(end for _, end in child_intervals),
            )
    return {
        node: _LeafIntervalSet(leaf_names, name_to_index, start, end)
        for node, (start, end) in intervals.items()
    }

def get_subtree_leaf_bitmasks(tree, leaf_name_to_bit):
    subtree_leaf_bitmasks = dict()
    for node in tree.traverse(strategy='postorder'):
        if node.is_leaf:
            if node.name not in leaf_name_to_bit:
                raise ValueError('Leaf label not found in reference mapping: {}'.format(node.name))
            subtree_leaf_bitmasks[node] = 1 << leaf_name_to_bit[node.name]
            continue
        bitmask = 0
        for child in node.get_children():
            bitmask |= subtree_leaf_bitmasks[child]
        subtree_leaf_bitmasks[node] = bitmask
    return subtree_leaf_bitmasks

def validate_unique_named_leaves(tree, option_name, context=''):
    leaf_names = list(tree.leaf_names())
    if len(leaf_names) != len(set(leaf_names)):
        raise ValueError("Duplicated leaf labels are not supported in '{}'{}.".format(option_name, context))
    has_missing = any((name is None) or (str(name) == '') for name in leaf_names)
    if has_missing:
        raise ValueError("Empty leaf labels are not supported in '{}'{}.".format(option_name, context))

def read_item_per_line_file(file):
    if file == '-':
        out = sys.stdin.read().splitlines()
    else:
        with open(file, 'r') as f:
            out = f.read().splitlines()
    out = [o.strip() for o in out if o.strip() != '']
    return out

def annotate_scientific_names(tree, species_regex=None, args=None, species_parser=None, species_map_tsv=None):
    for node in tree.leaves():
        parsed_species = extract_parsed_species(
            node.name,
            args=args,
            species_parser=species_parser,
            species_regex=species_regex,
            species_map_tsv=species_map_tsv,
        )
        node.add_props(
            species_label=parsed_species.species_label,
            taxonomy_query=parsed_species.taxonomy_query,
            sci_name=parsed_species.species_label,
        )
    return tree

def get_subtree_sci_name_sets(tree):
    subtree_sci_name_sets = dict()
    for node in tree.traverse(strategy='postorder'):
        if node.is_leaf:
            species_label = node.props.get('species_label', node.props.get('sci_name'))
            subtree_sci_name_sets[node] = {species_label}
            continue
        sci_names = set()
        for child in node.get_children():
            sci_names.update(subtree_sci_name_sets[child])
        subtree_sci_name_sets[node] = sci_names
    return subtree_sci_name_sets

def annotate_duplication_confidence_scores(tree, subtree_sci_name_sets=None):
    if subtree_sci_name_sets is None:
        subtree_sci_name_sets = get_subtree_sci_name_sets(tree)
    for node in tree.traverse():
        if node.is_leaf:
            continue
        children = node.get_children()
        if len(children) != 2:
            node.add_props(dup_conf_score=0.0)
            continue
        sp_child1 = subtree_sci_name_sets[children[0]]
        sp_child2 = subtree_sci_name_sets[children[1]]
        num_union = len(sp_child1.union(sp_child2))
        num_intersection = len(sp_child1.intersection(sp_child2))
        if num_union == 0:
            node.add_props(dup_conf_score=0.0)
        else:
            node.add_props(dup_conf_score=num_intersection / num_union)
    return tree

def is_all_leaf_names_identical(tree1, tree2, verbose=False):
    leaf_names1 = list(tree1.leaf_names())
    leaf_names2 = list(tree2.leaf_names())
    counts1 = Counter(leaf_names1)
    counts2 = Counter(leaf_names2)
    is_all_leaf_names_identical = (counts1 == counts2)
    if verbose:
        if not is_all_leaf_names_identical:
            all_names = set(counts1.keys()).union(set(counts2.keys()))
            unmatched_names = [
                str(name) for name in sorted(all_names, key=lambda x: str(x))
                if counts1.get(name, 0) != counts2.get(name, 0)
            ]
            sys.stderr.write('Unmatched leaf labels: {}\n'.format(' '.join(unmatched_names)))
    return is_all_leaf_names_identical

def get_target_nodes(tree, target):
    if (target=='all'):
        nodes = list(tree.traverse())
    elif (target=='root'):
        nodes = [tree]
    elif (target=='leaf'):
        nodes = list(tree.leaves())
    elif (target=='intnode'):
        nodes = [ node for node in tree.traverse() if not node.is_leaf ]
    else:
        raise ValueError('Unknown target: {}'.format(target))
    return nodes

def is_rooted(tree):
    num_subroot_nodes = len(tree.get_children())
    if num_subroot_nodes <= 2:
        return True
    return False

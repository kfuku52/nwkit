import hashlib
import io
import json
import os
import sys
import time
from datetime import datetime, timezone

from ete4 import Tree

from nwkit import __version__
from nwkit.conventions import get_stdin_input_options
from nwkit.util import inspect_tree_text, is_rooted, split_newick_stream


OUTPUT_ARGUMENTS = frozenset((
    'outfile',
    'report',
    'tree_out',
    'model_out',
    'stochastic_map_out',
    'output_table',
    'seqout',
    'out_dir',
    'manifest_out',
    'attribution_out',
    'group_table_prefix',
    'audit',
))


class _TeeTextWriter:
    def __init__(self, stream, capture=False):
        self.stream = stream
        self.capture = capture
        self.hasher = hashlib.sha256()
        self.byte_count = 0
        self.parts = list()

    def write(self, text):
        text = str(text)
        encoded = text.encode('utf-8')
        self.hasher.update(encoded)
        self.byte_count += len(encoded)
        if self.capture:
            self.parts.append(text)
        return self.stream.write(text)

    def flush(self):
        return self.stream.flush()

    def __getattr__(self, name):
        return getattr(self.stream, name)

    @property
    def sha256(self):
        return self.hasher.hexdigest()

    @property
    def text(self):
        return ''.join(self.parts)


def _sha256_bytes(data):
    return hashlib.sha256(data).hexdigest()


def _sha256_file(path):
    hasher = hashlib.sha256()
    size = 0
    with open(path, 'rb') as handle:
        while True:
            chunk = handle.read(1024 * 1024)
            if not chunk:
                break
            size += len(chunk)
            hasher.update(chunk)
    return hasher.hexdigest(), size


def _sha256_directory(path):
    hasher = hashlib.sha256()
    total_size = 0
    file_count = 0
    for root, dirnames, filenames in os.walk(path):
        dirnames.sort()
        for filename in sorted(filenames):
            filepath = os.path.join(root, filename)
            relative = os.path.relpath(filepath, path).replace(os.sep, '/')
            digest, size = _sha256_file(filepath)
            hasher.update(relative.encode('utf-8'))
            hasher.update(b'\0')
            hasher.update(digest.encode('ascii'))
            hasher.update(b'\n')
            total_size += size
            file_count += 1
    return hasher.hexdigest(), total_size, file_count


def _json_value(value):
    if isinstance(value, (str, int, float, bool, type(None))):
        return value
    if isinstance(value, (list, tuple)):
        return [_json_value(item) for item in value]
    if isinstance(value, dict):
        return {str(key): _json_value(item) for key, item in value.items()}
    return str(value)


def _argument_dict(args):
    return {
        key: _json_value(value)
        for key, value in vars(args).items()
        if key != 'handler'
    }


def _path_candidates_from_value(value):
    if isinstance(value, str):
        candidates = [value]
        if '@' in value:
            candidates.append(value.rsplit('@', 1)[1])
        return candidates
    if isinstance(value, (list, tuple)):
        candidates = list()
        for item in value:
            candidates.extend(_path_candidates_from_value(item))
        return candidates
    return []


def _input_file_records(args):
    records = list()
    seen = set()
    for argument, value in vars(args).items():
        if argument in OUTPUT_ARGUMENTS or argument == 'handler':
            continue
        if argument in ('manifest', 'attribution') and getattr(args, 'command', None) == 'image':
            continue
        for candidate in _path_candidates_from_value(value):
            if candidate in ('', '-') or not os.path.isfile(candidate):
                continue
            realpath = os.path.realpath(candidate)
            if realpath in seen:
                continue
            seen.add(realpath)
            digest, size = _sha256_file(realpath)
            records.append({
                'argument': argument,
                'path': realpath,
                'sha256': digest,
                'bytes': size,
            })
    return sorted(records, key=lambda record: (record['argument'], record['path']))


def _output_file_records(args):
    records = list()
    seen = set()
    output_arguments = set(OUTPUT_ARGUMENTS)
    if getattr(args, 'command', None) == 'image':
        output_arguments.update(('manifest', 'attribution'))
    for argument in output_arguments:
        if argument == 'audit':
            continue
        value = getattr(args, argument, None)
        if (
            argument == 'group_table_prefix'
            and value in (None, '')
            and getattr(args, 'command', None) == 'skim'
            and bool(getattr(args, 'output_groupfile', False))
        ):
            outfile = getattr(args, 'outfile', None)
            if outfile not in (None, '', '-'):
                value = outfile.removesuffix('.nwk')
        candidates = _path_candidates_from_value(value)
        if argument == 'group_table_prefix':
            candidates = [
                '{}.{}'.format(candidate, suffix)
                for candidate in candidates
                for suffix in ('all.tsv', 'sampled.tsv')
            ]
        for candidate in candidates:
            if candidate in ('', '-'):
                continue
            realpath = os.path.realpath(candidate)
            if realpath in seen:
                continue
            seen.add(realpath)
            if os.path.isfile(realpath):
                digest, size = _sha256_file(realpath)
                records.append({
                    'argument': argument,
                    'path': realpath,
                    'type': 'file',
                    'sha256': digest,
                    'bytes': size,
                })
            elif os.path.isdir(realpath):
                digest, size, file_count = _sha256_directory(realpath)
                records.append({
                    'argument': argument,
                    'path': realpath,
                    'type': 'directory',
                    'sha256': digest,
                    'bytes': size,
                    'file_count': file_count,
                })
    return sorted(records, key=lambda record: (record['argument'], record['path']))


def _primary_input_text(args, stdin_text):
    infile = getattr(args, 'infile', None)
    if infile == '-':
        return stdin_text
    if isinstance(infile, str) and os.path.isfile(infile):
        try:
            with open(infile) as handle:
                return handle.read()
        except UnicodeDecodeError:
            return None
    if isinstance(infile, str):
        return infile
    return None


def _input_summary(text, args):
    if text in (None, ''):
        return {}
    try:
        tree_strings = split_newick_stream(text)
    except ValueError:
        first_line = text.splitlines()[0] if text.splitlines() else ''
        if '\t' in first_line:
            return {
                'kind': 'table',
                'columns': first_line.split('\t'),
                'row_count': max(0, len(text.splitlines()) - 1),
            }
        return {'kind': 'text'}
    if not tree_strings:
        return {'kind': 'text'}
    inspection = inspect_tree_text(
        tree_strings[0],
        format=getattr(args, 'format', 'auto'),
        quoted_node_names=getattr(args, 'quoted_node_names', True),
    )
    summary = {
        'kind': 'newick',
        'tree_count': len(tree_strings),
        'parse_ok': inspection['parse_ok'],
        'input_format': inspection['input_format'],
        'format_ambiguous': inspection['format_ambiguous'],
    }
    if inspection['parse_ok']:
        try:
            tree = Tree(tree_strings[0], parser=int(inspection['input_format']))
            summary.update({
                'first_tree_tip_count': len(list(tree.leaves())),
                'first_tree_rooted': is_rooted(tree),
            })
        except Exception:
            pass
    return summary


def _seed_arguments(args):
    return {
        key: _json_value(value)
        for key, value in vars(args).items()
        if 'seed' in key.lower() and value is not None
    }


def _external_context(args):
    keys = (
        'download_dir',
        'taxonomy_source',
        'rank',
        'species_parser',
        'species_regex',
        'species_map_tsv',
    )
    return {
        key: _json_value(getattr(args, key))
        for key in keys
        if hasattr(args, key) and getattr(args, key) is not None
    }


def _write_audit(path, record):
    parent = os.path.dirname(os.path.realpath(path))
    if parent:
        os.makedirs(parent, exist_ok=True)
    with open(path, 'a') as handle:
        handle.write(json.dumps(record, sort_keys=True, ensure_ascii=False) + '\n')


def run_with_audit(args, argv, handler):
    audit_path = getattr(args, 'audit', None)
    if audit_path in (None, ''):
        return handler(args)
    if audit_path == '-':
        raise ValueError("'--audit' requires a file path; stdout is reserved for primary output.")
    original_stdin = sys.stdin
    original_stdout = sys.stdout
    original_stderr = sys.stderr
    stdin_text = None
    stdin_options = get_stdin_input_options(args)
    stdin_argument = stdin_options[0][0] if stdin_options else None
    if stdin_argument is not None:
        stdin_text = original_stdin.read()
        sys.stdin = io.StringIO(stdin_text)
    stdout_tee = _TeeTextWriter(original_stdout, capture=False)
    stderr_tee = _TeeTextWriter(original_stderr, capture=True)
    sys.stdout = stdout_tee
    sys.stderr = stderr_tee
    start_time = time.monotonic()
    started_at = datetime.now(timezone.utc).isoformat()
    status = 'ok'
    error = None
    return_value = None
    try:
        return_value = handler(args)
    except (Exception, SystemExit) as exc:
        status = 'error'
        error = {
            'type': type(exc).__name__,
            'message': str(exc),
        }
        raise
    finally:
        duration = time.monotonic() - start_time
        sys.stdin = original_stdin
        sys.stdout = original_stdout
        sys.stderr = original_stderr
        stderr_lines = stderr_tee.text.splitlines()
        record = {
            'schema': 'nwkit-audit-v1',
            'started_at_utc': started_at,
            'duration_seconds': round(duration, 6),
            'status': status,
            'error': error,
            'nwkit_version': __version__,
            'command': argv[0] if argv else '',
            'argv': list(argv),
            'arguments': _argument_dict(args),
            'random_seeds': _seed_arguments(args),
            'external_context': _external_context(args),
            'inputs': _input_file_records(args),
            'primary_input': _input_summary(_primary_input_text(args, stdin_text), args),
            'stdin': None if stdin_text is None else {
                'argument': stdin_argument,
                'sha256': _sha256_bytes(stdin_text.encode('utf-8')),
                'bytes': len(stdin_text.encode('utf-8')),
            },
            'outputs': _output_file_records(args),
            'stdout': {
                'sha256': stdout_tee.sha256,
                'bytes': stdout_tee.byte_count,
            },
            'warnings': [line for line in stderr_lines if 'warning' in line.lower()],
            'messages': stderr_lines[:500],
        }
        try:
            _write_audit(audit_path, record)
        except Exception as audit_exc:
            original_stderr.write('Warning: failed to write audit record: {}\n'.format(audit_exc))
    return return_value

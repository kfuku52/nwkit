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
from nwkit.util import (
    acquire_exclusive_lock,
    inspect_tree_text,
    is_rooted,
    split_newick_stream,
    validate_distinct_output_paths,
)


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
    def __init__(self, stream, capture=False, max_lines=500, max_line_chars=8192):
        self.stream = stream
        self.capture = capture
        self.max_lines = max_lines
        self.max_line_chars = max_line_chars
        self.hasher = hashlib.sha256()
        self.byte_count = 0
        self.lines = list()
        self.warning_lines = list()
        self._pending_line = ''
        self._pending_is_warning = False

    def _capture_complete_line(self, line, is_warning=False):
        line = line[:self.max_line_chars]
        if len(self.lines) < self.max_lines:
            self.lines.append(line)
        if is_warning and len(self.warning_lines) < self.max_lines:
            self.warning_lines.append(line)

    def _capture_text(self, text):
        segments = text.split('\n')
        for index, segment in enumerate(segments):
            if index == 0:
                self._pending_is_warning = (
                    self._pending_is_warning
                    or 'warning' in (self._pending_line + segment).lower()
                )
                self._pending_line = (self._pending_line + segment)[:self.max_line_chars]
            else:
                self._capture_complete_line(
                    self._pending_line,
                    is_warning=self._pending_is_warning,
                )
                self._pending_line = segment[:self.max_line_chars]
                self._pending_is_warning = 'warning' in segment.lower()

    def write(self, text):
        text = str(text)
        encoded = text.encode('utf-8')
        self.hasher.update(encoded)
        self.byte_count += len(encoded)
        if self.capture:
            self._capture_text(text)
        return self.stream.write(text)

    def flush(self):
        return self.stream.flush()

    def __getattr__(self, name):
        return getattr(self.stream, name)

    @property
    def sha256(self):
        return self.hasher.hexdigest()

    @property
    def captured_lines(self):
        lines = list(self.lines)
        if self._pending_line and len(lines) < self.max_lines:
            lines.append(self._pending_line)
        return lines

    @property
    def captured_warning_lines(self):
        lines = list(self.warning_lines)
        if (
            self._pending_line
            and self._pending_is_warning
            and len(lines) < self.max_lines
        ):
            lines.append(self._pending_line)
        return lines


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
        if key != 'handler' and not key.startswith('_nwkit_')
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
    records_by_path = dict()

    def add_candidate(argument, candidate):
        if candidate in ('', '-') or not os.path.isfile(candidate):
            return
        realpath = os.path.realpath(candidate)
        record = records_by_path.get(realpath)
        if record is None:
            digest, size = _sha256_file(realpath)
            record = {
                'path': realpath,
                'sha256': digest,
                'bytes': size,
                'arguments': set(),
            }
            records_by_path[realpath] = record
        record['arguments'].add(argument)

    for argument, value in vars(args).items():
        if (
            argument in OUTPUT_ARGUMENTS
            or argument == 'handler'
            or argument.startswith('_nwkit_')
        ):
            continue
        if argument in ('manifest', 'attribution') and getattr(args, 'command', None) == 'image':
            continue
        for candidate in _path_candidates_from_value(value):
            add_candidate(argument, candidate)
    if getattr(args, 'command', None) == 'compose':
        manifest_path = getattr(args, 'manifest', None)
        if manifest_path not in (None, '') and os.path.isfile(manifest_path):
            try:
                with open(manifest_path) as handle:
                    manifest = json.load(handle)
            except (OSError, ValueError, TypeError):
                manifest = None
            if isinstance(manifest, dict):
                base_dir = os.path.dirname(os.path.realpath(manifest_path))
                for key in ('root', 'name', 'support', 'length'):
                    value = manifest.get(key)
                    if value not in (None, '', '-'):
                        candidate = value if os.path.isabs(str(value)) else os.path.join(base_dir, str(value))
                        add_candidate('manifest:{}'.format(key), candidate)
                properties = manifest.get('properties', [])
                if isinstance(properties, list):
                    for index, entry in enumerate(properties):
                        if not isinstance(entry, dict) or entry.get('path') in (None, '', '-'):
                            continue
                        value = entry['path']
                        candidate = value if os.path.isabs(str(value)) else os.path.join(base_dir, str(value))
                        add_candidate('manifest:properties[{}]'.format(index), candidate)
    if getattr(args, 'command', None) == 'draw':
        for path in getattr(args, '_nwkit_tip_image_paths', ()):
            add_candidate('tip_image_manifest:asset', path)
    records = list()
    for record in records_by_path.values():
        arguments = sorted(record.pop('arguments'))
        record['argument'] = arguments[0]
        record['arguments'] = arguments
        records.append(record)
    return sorted(records, key=lambda record: (record['argument'], record['path']))


def _output_file_records(args):
    records_by_path = dict()
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
            record = records_by_path.get(realpath)
            if record is not None:
                record['arguments'].add(argument)
                continue
            if os.path.isfile(realpath):
                digest, size = _sha256_file(realpath)
                records_by_path[realpath] = {
                    'path': realpath,
                    'type': 'file',
                    'sha256': digest,
                    'bytes': size,
                    'arguments': {argument},
                }
            elif os.path.isdir(realpath):
                digest, size, file_count = _sha256_directory(realpath)
                records_by_path[realpath] = {
                    'path': realpath,
                    'type': 'directory',
                    'sha256': digest,
                    'bytes': size,
                    'file_count': file_count,
                    'arguments': {argument},
                }
    records = list()
    for record in records_by_path.values():
        arguments = sorted(record.pop('arguments'))
        record['argument'] = arguments[0]
        record['arguments'] = arguments
        records.append(record)
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
    payload = (json.dumps(record, sort_keys=True, ensure_ascii=False) + '\n').encode('utf-8')
    lock_path = os.path.realpath(path) + '.lock'
    with acquire_exclusive_lock(lock_path, lock_label='audit log'):
        fd = os.open(path, os.O_APPEND | os.O_CREAT | os.O_WRONLY, 0o666)
        try:
            view = memoryview(payload)
            while view:
                written = os.write(fd, view)
                view = view[written:]
            os.fsync(fd)
        finally:
            os.close(fd)


def run_with_audit(args, argv, handler):
    audit_path = getattr(args, 'audit', None)
    if audit_path in (None, ''):
        return handler(args)
    if audit_path == '-':
        raise ValueError("'--audit' requires a file path; stdout is reserved for primary output.")
    audit_collision_candidates = [('--audit', audit_path)]
    for argument in OUTPUT_ARGUMENTS:
        if argument == 'audit':
            continue
        value = getattr(args, argument, None)
        if argument == 'group_table_prefix':
            if (
                value in (None, '')
                and getattr(args, 'command', None) == 'skim'
                and bool(getattr(args, 'output_groupfile', False))
            ):
                outfile = getattr(args, 'outfile', None)
                if outfile not in (None, '', '-'):
                    value = outfile.removesuffix('.nwk')
            if value not in (None, ''):
                audit_collision_candidates.extend((
                    ('--group-table-prefix .all.tsv', '{}.all.tsv'.format(value)),
                    ('--group-table-prefix .sampled.tsv', '{}.sampled.tsv'.format(value)),
                ))
        else:
            for candidate in _path_candidates_from_value(value):
                audit_collision_candidates.append(('--{}'.format(argument.replace('_', '-')), candidate))
    validate_distinct_output_paths(audit_collision_candidates)
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
        stderr_lines = stderr_tee.captured_lines
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
            'warnings': stderr_tee.captured_warning_lines,
            'messages': stderr_lines,
        }
        try:
            _write_audit(audit_path, record)
        except Exception as audit_exc:
            original_stderr.write('Warning: failed to write audit record: {}\n'.format(audit_exc))
    return return_value

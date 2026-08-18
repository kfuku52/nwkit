import csv
import hashlib
import json
import os
import stat
import sys
import tempfile
import time
from datetime import datetime, timezone
from typing import Any

from ete4 import Tree

from nwkit import __version__
from nwkit.conventions import get_stdin_input_options
from nwkit.util import (
    acquire_exclusive_lock,
    inspect_tree_text,
    is_rooted,
    normalized_missing_path_key,
    split_newick_stream,
    validate_distinct_output_paths,
)

OUTPUT_ARGUMENTS = frozenset(
    (
        "outfile",
        "report",
        "tree_out",
        "model_out",
        "model_comparison_out",
        "gene_contrasts_out",
        "random_effects_out",
        "reconciliation_out",
        "response_sampling_covariance_out",
        "response_tip_summary_out",
        "predictor_sampling_covariance_out",
        "predictor_tip_summary_out",
        "sampling_covariance_out",
        "species_contrasts_out",
        "stochastic_map_out",
        "tip_summary_out",
        "output_table",
        "seqout",
        "out_dir",
        "manifest_out",
        "attribution_out",
        "group_table_prefix",
        "audit",
    )
)

INPUT_PATH_ARGUMENTS = frozenset(
    (
        "infile",
        "infile2",
        "data",
        "expression",
        "evolution_covariance",
        "gene_tree",
        "gene_tree_ensemble",
        "length_source",
        "manifest",
        "mcmctree_posterior",
        "name_source",
        "name_tsv",
        "property_source",
        "predictor_contrasts",
        "predictor_sampling_covariance",
        "reconciliation",
        "reconciliation_tree",
        "reference",
        "response_sampling_covariance",
        "root_source",
        "seqin",
        "species_list",
        "species_map_tsv",
        "species_name_tsv",
        "species_tree",
        "species_traits",
        "support_source",
        "table",
        "taxid_tsv",
        "tip_image_manifest",
        "trait",
        "tree",
        "weight_tsv",
        "densitree_trees",
    )
)

STDIN_SPOOL_MEMORY_BYTES = 8 * 1024 * 1024
INPUT_SUMMARY_MAX_CHARS = 4 * 1024 * 1024
STREAM_CHUNK_CHARS = 1024 * 1024


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
        self._pending_line = ""
        self._pending_is_warning = False

    def _capture_complete_line(self, line, is_warning=False):
        line = line[: self.max_line_chars]
        if len(self.lines) < self.max_lines:
            self.lines.append(line)
        if is_warning and len(self.warning_lines) < self.max_lines:
            self.warning_lines.append(line)

    def _capture_text(self, text):
        segments = text.split("\n")
        for index, segment in enumerate(segments):
            if index == 0:
                self._pending_is_warning = (
                    self._pending_is_warning
                    or "warning" in (self._pending_line + segment).lower()
                )
                self._pending_line = (self._pending_line + segment)[
                    : self.max_line_chars
                ]
            else:
                self._capture_complete_line(
                    self._pending_line,
                    is_warning=self._pending_is_warning,
                )
                self._pending_line = segment[: self.max_line_chars]
                self._pending_is_warning = "warning" in segment.lower()

    def write(self, text):
        text = str(text)
        encoded = text.encode("utf-8")
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
    flags = os.O_RDONLY
    if hasattr(os, "O_NOFOLLOW"):
        flags |= os.O_NOFOLLOW
    if hasattr(os, "O_NONBLOCK"):
        flags |= os.O_NONBLOCK
    opened_fd = os.open(path, flags)
    fd: int | None = opened_fd
    try:
        file_stat = os.fstat(opened_fd)
        if not stat.S_ISREG(file_stat.st_mode):
            raise ValueError(
                "Audit hashing only supports regular files: '{}'.".format(path)
            )
        handle = os.fdopen(opened_fd, "rb")
        fd = None
        with handle:
            while True:
                chunk = handle.read(1024 * 1024)
                if not chunk:
                    break
                size += len(chunk)
                hasher.update(chunk)
    finally:
        if fd is not None:
            os.close(fd)
    return hasher.hexdigest(), size


def _sha256_directory(path):
    hasher = hashlib.sha256()
    total_size = 0
    file_count = 0
    for root, dirnames, filenames in os.walk(path):
        dirnames.sort()
        for dirname in dirnames:
            directory_path = os.path.join(root, dirname)
            directory_stat = os.lstat(directory_path)
            if stat.S_ISLNK(directory_stat.st_mode):
                raise ValueError(
                    "Audit directory hashing does not follow symbolic links: '{}'.".format(
                        directory_path
                    )
                )
            if not stat.S_ISDIR(directory_stat.st_mode):
                raise ValueError(
                    "Audit directory hashing encountered a special entry: '{}'.".format(
                        directory_path
                    )
                )
        for filename in sorted(filenames):
            filepath = os.path.join(root, filename)
            file_stat = os.lstat(filepath)
            if not stat.S_ISREG(file_stat.st_mode):
                raise ValueError(
                    "Audit directory hashing only supports regular files: '{}'.".format(
                        filepath
                    )
                )
            relative = os.path.relpath(filepath, path).replace(os.sep, "/")
            digest, size = _sha256_file(filepath)
            hasher.update(relative.encode("utf-8"))
            hasher.update(b"\0")
            hasher.update(digest.encode("ascii"))
            hasher.update(b"\n")
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
        if key != "handler" and not key.startswith("_nwkit_")
    }


def _path_candidates_from_value(value):
    if isinstance(value, str):
        return [value]
    if isinstance(value, (list, tuple)):
        candidates: list[str] = []
        for item in value:
            candidates.extend(_path_candidates_from_value(item))
        return candidates
    return []


def _declared_input_path_candidates(args):
    candidates: list[tuple[str, str]] = []
    for argument in sorted(INPUT_PATH_ARGUMENTS):
        if not hasattr(args, argument):
            continue
        if (
            argument in ("manifest", "attribution")
            and getattr(args, "command", None) == "image"
        ):
            continue
        values = _path_candidates_from_value(getattr(args, argument))
        if argument == "property_source":
            values = [value.rsplit("@", 1)[1] for value in values if "@" in value]
        candidates.extend((argument, value) for value in values)
    return candidates


def _input_path_candidates(args):
    candidates = _declared_input_path_candidates(args)
    if getattr(args, "command", None) == "compose":
        manifest_path = getattr(args, "manifest", None)
        if manifest_path not in (None, "") and os.path.isfile(manifest_path):
            try:
                with open(manifest_path) as handle:
                    manifest = json.load(handle)
            except (OSError, ValueError, TypeError):
                manifest = None
            if isinstance(manifest, dict):
                base_dir = os.path.dirname(os.path.realpath(manifest_path))
                for key in ("root", "name", "support", "length"):
                    value = manifest.get(key)
                    if value not in (None, "", "-"):
                        candidate = (
                            value
                            if os.path.isabs(str(value))
                            else os.path.join(base_dir, str(value))
                        )
                        candidates.append(("manifest:{}".format(key), candidate))
                properties = manifest.get("properties", [])
                if isinstance(properties, list):
                    for index, entry in enumerate(properties):
                        if not isinstance(entry, dict) or entry.get("path") in (
                            None,
                            "",
                            "-",
                        ):
                            continue
                        value = entry["path"]
                        candidate = (
                            value
                            if os.path.isabs(str(value))
                            else os.path.join(base_dir, str(value))
                        )
                        candidates.append(
                            (
                                "manifest:properties[{}]".format(index),
                                candidate,
                            )
                        )
    if getattr(args, "command", None) == "draw":
        manifest_path = getattr(args, "tip_image_manifest", None)
        if manifest_path not in (None, "", "-") and os.path.isfile(manifest_path):
            tip_image_root = getattr(args, "tip_image_root", None)
            if tip_image_root in (None, ""):
                tip_image_root = os.path.dirname(os.path.realpath(manifest_path))
            try:
                with open(manifest_path, newline="") as handle:
                    for row in csv.DictReader(handle, delimiter="\t"):
                        value = row.get("local_path")
                        if value in (None, "", "-"):
                            continue
                        candidate = (
                            value
                            if os.path.isabs(value)
                            else os.path.join(tip_image_root, value)
                        )
                        candidates.append(("tip_image_manifest:asset", candidate))
            except (OSError, UnicodeError, csv.Error):
                pass
        for path in getattr(args, "_nwkit_tip_image_paths", ()):
            candidates.append(("tip_image_manifest:asset", path))
    return candidates


def _input_file_records(args, input_candidates=None):
    records_by_path: dict[str, dict[str, Any]] = {}
    if input_candidates is None:
        input_candidates = _input_path_candidates(args)
    for argument, candidate in input_candidates:
        if candidate in ("", "-"):
            continue
        try:
            candidate_exists = os.path.lexists(candidate)
        except (OSError, ValueError) as exc:
            raise ValueError(
                "Audit input path is invalid for '{}': {!r}.".format(
                    argument,
                    candidate,
                )
            ) from exc
        if not candidate_exists:
            continue
        try:
            candidate_stat = os.stat(candidate)
        except OSError as exc:
            raise ValueError(
                "Audit input path is unavailable for '{}': '{}'.".format(
                    argument,
                    candidate,
                )
            ) from exc
        if not stat.S_ISREG(candidate_stat.st_mode):
            raise ValueError(
                "Audit input paths must be regular files; use '-' for "
                "standard input ({}): '{}'.".format(argument, candidate)
            )
        realpath = os.path.realpath(candidate)
        record = records_by_path.get(realpath)
        if record is None:
            digest, size = _sha256_file(realpath)
            record = {
                "path": realpath,
                "sha256": digest,
                "bytes": size,
                "arguments": set(),
            }
            records_by_path[realpath] = record
        record["arguments"].add(argument)
    records = list()
    for record in records_by_path.values():
        arguments = sorted(record.pop("arguments"))
        record["argument"] = arguments[0]
        record["arguments"] = arguments
        records.append(record)
    return sorted(records, key=lambda record: (record["argument"], record["path"]))


def _output_file_records(args):
    records_by_path: dict[str, dict[str, Any]] = {}
    output_arguments = set(OUTPUT_ARGUMENTS)
    if getattr(args, "command", None) == "image":
        output_arguments.update(("manifest", "attribution"))
    for argument in output_arguments:
        if argument == "audit":
            continue
        value = getattr(args, argument, None)
        if (
            argument == "group_table_prefix"
            and value in (None, "")
            and getattr(args, "command", None) == "skim"
            and bool(getattr(args, "output_groupfile", False))
        ):
            outfile = getattr(args, "outfile", None)
            if outfile not in (None, "", "-"):
                value = outfile.removesuffix(".nwk")
        candidates = _path_candidates_from_value(value)
        if argument == "group_table_prefix":
            candidates = [
                "{}.{}".format(candidate, suffix)
                for candidate in candidates
                for suffix in ("all.tsv", "sampled.tsv")
            ]
        for candidate in candidates:
            if candidate in ("", "-"):
                continue
            realpath = os.path.realpath(candidate)
            record = records_by_path.get(realpath)
            if record is not None:
                record["arguments"].add(argument)
                continue
            if os.path.isfile(realpath):
                digest, size = _sha256_file(realpath)
                records_by_path[realpath] = {
                    "path": realpath,
                    "type": "file",
                    "sha256": digest,
                    "bytes": size,
                    "arguments": {argument},
                }
            elif os.path.isdir(realpath):
                digest, size, file_count = _sha256_directory(realpath)
                records_by_path[realpath] = {
                    "path": realpath,
                    "type": "directory",
                    "sha256": digest,
                    "bytes": size,
                    "file_count": file_count,
                    "arguments": {argument},
                }
    records = list()
    for record in records_by_path.values():
        arguments = sorted(record.pop("arguments"))
        record["argument"] = arguments[0]
        record["arguments"] = arguments
        records.append(record)
    return sorted(records, key=lambda record: (record["argument"], record["path"]))


def _planned_output_records(args):
    records_by_path: dict[str, dict[str, Any]] = {}
    for argument in OUTPUT_ARGUMENTS:
        if argument == "audit":
            continue
        value = getattr(args, argument, None)
        candidates = _path_candidates_from_value(value)
        if argument == "group_table_prefix":
            candidates = [
                "{}.{}".format(candidate, suffix)
                for candidate in candidates
                for suffix in ("all.tsv", "sampled.tsv")
            ]
        for candidate in candidates:
            if candidate in ("", "-"):
                continue
            realpath = os.path.realpath(candidate)
            record = records_by_path.setdefault(
                realpath,
                {"path": realpath, "arguments": set()},
            )
            record["arguments"].add(argument)
    records = []
    for record in records_by_path.values():
        arguments = sorted(record.pop("arguments"))
        record["argument"] = arguments[0]
        record["arguments"] = arguments
        records.append(record)
    return sorted(records, key=lambda record: (record["argument"], record["path"]))


def _primary_input_text(args, stdin_text):
    infile = getattr(args, "infile", None)
    gene_tree = getattr(args, "gene_tree", None)
    gene_tree_ensemble = getattr(args, "gene_tree_ensemble", None)
    ordinary_tree = getattr(args, "tree", None)
    primary_input = (
        ordinary_tree
        if getattr(args, "command", None) == "regress"
        and ordinary_tree not in (None, "")
        else (
            gene_tree
            if getattr(args, "command", None) == "regress"
            and gene_tree not in (None, "")
            else (
                gene_tree_ensemble
                if getattr(args, "command", None) == "regress"
                and gene_tree_ensemble not in (None, "")
                else infile
            )
        )
    )
    if primary_input == "-":
        return stdin_text
    if isinstance(primary_input, str) and os.path.isfile(primary_input):
        try:
            with open(primary_input, errors="replace") as handle:
                return handle.read(INPUT_SUMMARY_MAX_CHARS + 1)
        except UnicodeDecodeError:
            return None
    if isinstance(primary_input, str):
        return primary_input[: INPUT_SUMMARY_MAX_CHARS + 1]
    return None


def _input_summary(text, args):
    if text in (None, ""):
        return {}
    truncated = len(text) > INPUT_SUMMARY_MAX_CHARS
    if truncated:
        text = text[:INPUT_SUMMARY_MAX_CHARS]
    try:
        tree_strings = split_newick_stream(text)
    except ValueError:
        first_line = text.splitlines()[0] if text.splitlines() else ""
        if "\t" in first_line:
            summary = {
                "kind": "table",
                "columns": first_line.split("\t"),
                "row_count": max(0, len(text.splitlines()) - 1),
            }
            if truncated:
                summary["truncated"] = True
            return summary
        summary = {"kind": "text"}
        if truncated:
            summary["truncated"] = True
        return summary
    if not tree_strings:
        summary = {"kind": "text"}
        if truncated:
            summary["truncated"] = True
        return summary
    input_format = getattr(args, "format", "auto")
    if getattr(args, "command", None) == "regress":
        if getattr(args, "tree", None) not in (None, ""):
            input_format = getattr(args, "tree_format", "auto") or "auto"
        elif getattr(args, "gene_tree", None) not in (None, ""):
            input_format = getattr(args, "gene_tree_format", "auto") or "auto"
        elif getattr(args, "gene_tree_ensemble", None) not in (None, ""):
            input_format = getattr(args, "gene_tree_format", "auto") or "auto"
    inspection = inspect_tree_text(
        tree_strings[0],
        format=input_format,
        quoted_node_names=getattr(args, "quoted_node_names", True),
    )
    summary = {
        "kind": "newick",
        "tree_count": len(tree_strings),
        "parse_ok": inspection["parse_ok"],
        "input_format": inspection["input_format"],
        "format_ambiguous": inspection["format_ambiguous"],
    }
    if truncated:
        summary["truncated"] = True
    if inspection["parse_ok"]:
        try:
            tree = Tree(tree_strings[0], parser=int(inspection["input_format"]))
            summary.update(
                {
                    "first_tree_tip_count": len(list(tree.leaves())),
                    "first_tree_rooted": is_rooted(tree),
                }
            )
        except Exception:
            pass
    return summary


def _seed_arguments(args):
    return {
        key: _json_value(value)
        for key, value in vars(args).items()
        if "seed" in key.lower() and value is not None
    }


def _external_context(args):
    keys = (
        "download_dir",
        "taxonomy_source",
        "rank",
        "species_parser",
        "species_regex",
        "species_map_tsv",
    )
    return {
        key: _json_value(getattr(args, key))
        for key in keys
        if hasattr(args, key) and getattr(args, key) is not None
    }


def _open_secure_audit(path):
    path = os.fspath(path)
    if os.path.lexists(path):
        path_stat = os.lstat(path)
        if stat.S_ISLNK(path_stat.st_mode):
            raise ValueError(
                "Audit path must not be a symbolic link: '{}'.".format(path)
            )
        if not stat.S_ISREG(path_stat.st_mode):
            raise ValueError("Audit path must be a regular file: '{}'.".format(path))
    flags = os.O_APPEND | os.O_CREAT | os.O_WRONLY
    if hasattr(os, "O_NOFOLLOW"):
        flags |= os.O_NOFOLLOW
    fd = os.open(path, flags, 0o600)
    try:
        file_stat = os.fstat(fd)
        if not stat.S_ISREG(file_stat.st_mode):
            raise ValueError("Audit path must be a regular file: '{}'.".format(path))
        if hasattr(os, "geteuid") and file_stat.st_uid != os.geteuid():
            raise PermissionError(
                "Audit file must be owned by the current user: '{}'.".format(path)
            )
        if file_stat.st_nlink != 1:
            raise ValueError(
                "Audit file must not have multiple hard links: '{}'.".format(path)
            )
        if hasattr(os, "fchmod"):
            os.fchmod(fd, 0o600)
    except Exception:
        os.close(fd)
        raise
    return fd


def _prepare_audit(path):
    parent = os.path.dirname(os.path.abspath(os.fspath(path)))
    if parent:
        os.makedirs(parent, exist_ok=True)
    fd = _open_secure_audit(path)
    os.close(fd)


def _write_audit(path, record):
    parent = os.path.dirname(os.path.abspath(os.fspath(path)))
    if parent:
        os.makedirs(parent, exist_ok=True)
    payload = (json.dumps(record, sort_keys=True, ensure_ascii=False) + "\n").encode(
        "utf-8"
    )
    lock_path = os.path.realpath(path) + ".lock"
    with acquire_exclusive_lock(lock_path, lock_label="audit log"):
        fd = _open_secure_audit(path)
        try:
            view = memoryview(payload)
            while view:
                written = os.write(fd, view)
                view = view[written:]
            os.fsync(fd)
        finally:
            os.close(fd)


def _image_output_collision_candidates(args):
    if getattr(args, "command", None) != "image":
        return []
    out_dir = getattr(args, "out_dir", None)
    if out_dir in (None, ""):
        return []
    out_dir = os.path.realpath(out_dir)
    manifest_out = getattr(args, "manifest_out", None) or getattr(
        args, "manifest", None
    )
    attribution_out = getattr(args, "attribution_out", None) or getattr(
        args, "attribution", None
    )
    candidates = [
        ("--out-dir/unmatched.tsv", os.path.join(out_dir, "unmatched.tsv")),
    ]
    if manifest_out in (None, ""):
        candidates.append(("--manifest-out", os.path.join(out_dir, "manifest.tsv")))
    if attribution_out in (None, ""):
        candidates.append(
            ("--attribution-out", os.path.join(out_dir, "ATTRIBUTION.md"))
        )
    return candidates


def _audit_collision_candidates(args, audit_path):
    candidates = [("--audit", audit_path)]
    for argument in OUTPUT_ARGUMENTS:
        if argument == "audit":
            continue
        value = getattr(args, argument, None)
        if argument == "group_table_prefix":
            if (
                value in (None, "")
                and getattr(args, "command", None) == "skim"
                and bool(getattr(args, "output_groupfile", False))
            ):
                outfile = getattr(args, "outfile", None)
                if outfile not in (None, "", "-"):
                    value = outfile.removesuffix(".nwk")
            if value not in (None, ""):
                candidates.extend(
                    (
                        ("--group-table-prefix .all.tsv", "{}.all.tsv".format(value)),
                        (
                            "--group-table-prefix .sampled.tsv",
                            "{}.sampled.tsv".format(value),
                        ),
                    )
                )
        else:
            for candidate in _path_candidates_from_value(value):
                candidates.append(
                    ("--{}".format(argument.replace("_", "-")), candidate)
                )
    out_prefix = getattr(args, "out_prefix", None)
    if out_prefix not in (None, ""):
        from nwkit.conventions import (
            regression_bundle_lock_path,
            regression_bundle_paths,
        )

        candidates.extend(
            ("--out-prefix {}".format(argument), path)
            for argument, path in regression_bundle_paths(out_prefix).items()
        )
        candidates.append(
            ("--out-prefix transaction lock", regression_bundle_lock_path(out_prefix))
        )
    candidates.extend(_image_output_collision_candidates(args))
    return candidates


def _validate_audit_input_collision(audit_path, input_candidates):
    audit_realpath = os.path.realpath(audit_path)
    collision_arguments = set()
    for argument, candidate in input_candidates:
        if candidate in ("", "-"):
            continue
        candidate_realpath = os.path.realpath(candidate)
        same_path = normalized_missing_path_key(
            candidate_realpath
        ) == normalized_missing_path_key(audit_realpath)
        same_file = False
        if not same_path and os.path.exists(audit_path) and os.path.exists(candidate):
            try:
                same_file = os.path.samefile(audit_path, candidate)
            except OSError:
                same_file = False
        if same_path or same_file:
            collision_arguments.add(argument)
    if collision_arguments:
        raise ValueError(
            "Audit path must be distinct from input paths ({}): '{}'.".format(
                ", ".join(sorted(collision_arguments)),
                audit_realpath,
            )
        )


def _validate_audit_directory_collision(args, audit_path):
    audit_realpath = os.path.realpath(audit_path)
    out_dir = getattr(args, "out_dir", None)
    if out_dir in (None, "", "-"):
        return
    output_directory = os.path.realpath(out_dir)
    try:
        nested = (
            os.path.commonpath((audit_realpath, output_directory)) == output_directory
        )
    except ValueError:
        nested = False
    if nested:
        raise ValueError(
            "Audit path must not be inside output directory '--out-dir': '{}'.".format(
                output_directory
            )
        )


def _spool_stdin(stream):
    spool = tempfile.SpooledTemporaryFile(
        max_size=STDIN_SPOOL_MEMORY_BYTES,
        mode="w+t",
        encoding="utf-8",
        newline="",
    )
    hasher = hashlib.sha256()
    byte_count = 0
    summary_chunks = []
    summary_chars = 0
    try:
        while True:
            chunk = stream.read(STREAM_CHUNK_CHARS)
            if not chunk:
                break
            if not isinstance(chunk, str):
                raise TypeError("Standard input must be a text stream.")
            spool.write(chunk)
            encoded = chunk.encode("utf-8")
            hasher.update(encoded)
            byte_count += len(encoded)
            if summary_chars <= INPUT_SUMMARY_MAX_CHARS:
                remaining = INPUT_SUMMARY_MAX_CHARS + 1 - summary_chars
                if remaining > 0:
                    captured = chunk[:remaining]
                    summary_chunks.append(captured)
                    summary_chars += len(captured)
    except Exception:
        spool.close()
        raise
    spool.seek(0)
    return spool, "".join(summary_chunks), hasher.hexdigest(), byte_count


def _runtime_input_file_records(args):
    if getattr(args, "command", None) != "draw":
        return []
    records = []
    for path in getattr(args, "_nwkit_tip_image_paths", ()):
        if path in ("", "-") or not os.path.isfile(path):
            continue
        realpath = os.path.realpath(path)
        digest, size = _sha256_file(realpath)
        records.append(
            {
                "path": realpath,
                "sha256": digest,
                "bytes": size,
                "argument": "tip_image_manifest:asset",
                "arguments": ["tip_image_manifest:asset"],
            }
        )
    return sorted(records, key=lambda record: (record["argument"], record["path"]))


def _merge_input_records(snapshot_records, runtime_records):
    merged = {record["path"]: dict(record) for record in snapshot_records}
    for record in runtime_records:
        existing = merged.get(record["path"])
        if existing is None:
            merged[record["path"]] = dict(record)
            continue
        arguments = sorted(set(existing["arguments"]) | set(record["arguments"]))
        existing["argument"] = arguments[0]
        existing["arguments"] = arguments
    return sorted(
        merged.values(), key=lambda record: (record["argument"], record["path"])
    )


def run_with_audit(args, argv, handler):
    audit_path = getattr(args, "audit", None)
    if audit_path in (None, ""):
        return handler(args)
    if audit_path == "-":
        raise ValueError(
            "'--audit' requires a file path; stdout is reserved for primary output."
        )
    input_candidates = _input_path_candidates(args)
    input_records = _input_file_records(args, input_candidates=input_candidates)
    validate_distinct_output_paths(_audit_collision_candidates(args, audit_path))
    _validate_audit_input_collision(audit_path, input_candidates)
    _validate_audit_directory_collision(args, audit_path)
    _prepare_audit(audit_path)
    original_stdin = sys.stdin
    original_stdout = sys.stdout
    original_stderr = sys.stderr
    stdin_text = None
    stdin_sha256 = None
    stdin_bytes = None
    stdin_spool = None
    stdin_options = get_stdin_input_options(args)
    stdin_argument = stdin_options[0][0] if stdin_options else None
    if stdin_argument is not None:
        stdin_spool, stdin_text, stdin_sha256, stdin_bytes = _spool_stdin(
            original_stdin
        )
    try:
        primary_input = _input_summary(_primary_input_text(args, stdin_text), args)
    except BaseException:
        if stdin_spool is not None:
            stdin_spool.close()
        raise
    stdout_tee = _TeeTextWriter(original_stdout, capture=False)
    stderr_tee = _TeeTextWriter(original_stderr, capture=True)
    if stdin_spool is not None:
        sys.stdin = stdin_spool
    sys.stdout = stdout_tee
    sys.stderr = stderr_tee
    start_time = time.monotonic()
    started_at = datetime.now(timezone.utc).isoformat()
    status = "ok"
    error = None
    handler_exception = None
    return_value = None
    try:
        try:
            return_value = handler(args)
        except BaseException as exc:
            handler_exception = exc
            status = "error"
            error = {
                "type": type(exc).__name__,
                "message": str(exc),
            }
            raise
        finally:
            duration = time.monotonic() - start_time
            sys.stdin = original_stdin
            sys.stdout = original_stdout
            sys.stderr = original_stderr
            try:
                stderr_lines = stderr_tee.captured_lines
                runtime_records = _runtime_input_file_records(args)
                record = {
                    "schema": "nwkit-audit-v1",
                    "started_at_utc": started_at,
                    "duration_seconds": round(duration, 6),
                    "status": status,
                    "error": error,
                    "nwkit_version": __version__,
                    "command": argv[0] if argv else "",
                    "argv": list(argv),
                    "arguments": _argument_dict(args),
                    "random_seeds": _seed_arguments(args),
                    "external_context": _external_context(args),
                    "inputs": _merge_input_records(input_records, runtime_records),
                    "primary_input": primary_input,
                    "stdin": None
                    if stdin_text is None
                    else {
                        "argument": stdin_argument,
                        "sha256": stdin_sha256,
                        "bytes": stdin_bytes,
                    },
                    "planned_outputs": _planned_output_records(args),
                    "outputs": _output_file_records(args),
                    "stdout": {
                        "sha256": stdout_tee.sha256,
                        "bytes": stdout_tee.byte_count,
                    },
                    "warnings": stderr_tee.captured_warning_lines,
                    "messages": stderr_lines,
                }
                _write_audit(audit_path, record)
            except Exception as audit_exc:
                if handler_exception is None:
                    raise
                original_stderr.write(
                    "Warning: failed to finalize the audit record while "
                    "preserving the original command error: {}: {}\n".format(
                        type(audit_exc).__name__,
                        audit_exc,
                    )
                )
    finally:
        if stdin_spool is not None:
            stdin_spool.close()
    return return_value

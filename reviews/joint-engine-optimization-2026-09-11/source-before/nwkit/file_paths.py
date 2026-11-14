"""Filesystem identity checks shared by all input and output commands."""

import os
import unicodedata
from collections import defaultdict
from typing import Any


def _filesystem_is_case_insensitive(path):
    """Infer case behavior from the nearest existing ancestor without writes."""
    current = os.path.realpath(os.path.abspath(os.fspath(path)))
    while not os.path.exists(current):
        parent = os.path.dirname(current)
        if parent == current:
            return os.path.normcase("A") == os.path.normcase("a")
        current = parent
    swapped = current.swapcase()
    if swapped == current:
        return os.path.normcase("A") == os.path.normcase("a")
    try:
        return os.path.samefile(current, swapped)
    except OSError:
        return False


def normalized_missing_path_key(path):
    normalized = os.path.normcase(os.path.normpath(os.path.realpath(os.fspath(path))))
    if _filesystem_is_case_insensitive(normalized):
        normalized = unicodedata.normalize("NFC", normalized).casefold()
    return normalized


def validate_distinct_output_paths(outputs):
    """Reject output roles that resolve to the same path or filesystem object."""
    by_identity = defaultdict(list)
    for option_name, path in outputs:
        if path in (None, "", "-"):
            continue
        realpath = os.path.realpath(os.fspath(path))
        try:
            path_stat = os.stat(realpath)
        except OSError:
            identity: tuple[Any, ...] = (
                "path",
                normalized_missing_path_key(realpath),
            )
        else:
            identity = ("inode", path_stat.st_dev, path_stat.st_ino)
        by_identity[identity].append((str(option_name), realpath))
    collisions = [records for records in by_identity.values() if len(records) > 1]
    if collisions:
        details = list()
        for records in collisions:
            option_names = sorted(record[0] for record in records)
            paths = sorted({record[1] for record in records})
            details.append(
                "{} -> {}".format(
                    ", ".join(option_names),
                    " = ".join(paths),
                )
            )
        raise ValueError(
            "Output paths must be distinct: {}.".format("; ".join(sorted(details)))
        )


def paths_identify_same_file(first_path, second_path):
    """Return whether two existing or prospective paths identify one file."""
    if normalized_missing_path_key(first_path) == normalized_missing_path_key(
        second_path
    ):
        return True
    if not os.path.exists(first_path) or not os.path.exists(second_path):
        return False
    try:
        return os.path.samefile(first_path, second_path)
    except OSError:
        return False


def validate_outputs_do_not_replace_inputs(inputs, outputs, *, label="Output"):
    """Reject output paths that resolve to any protected input path."""
    for input_name, input_path in inputs:
        if input_path in (None, "", "-"):
            continue
        for output_name, output_path in outputs:
            if output_path in (None, "", "-"):
                continue
            if paths_identify_same_file(input_path, output_path):
                raise ValueError(
                    "{} '{}' must not overwrite input {} '{}'.".format(
                        label,
                        output_name,
                        input_name,
                        os.path.realpath(os.fspath(output_path)),
                    )
                )

"""Strict literal-keyed TSV input for gene-count and Ks family analyses."""

import csv
import os
from io import StringIO

from nwkit.util import read_input_text


def read_gene_family_rows(path, required):
    """Validate the literal TSV and return its header and positional data rows."""
    source = (
        open(path, encoding="utf-8-sig")
        if path != "-" and os.path.isfile(path)
        else StringIO(read_input_text(path))
    )
    with source:
        rows = list(csv.reader(source, delimiter="\t"))
    if not rows or not rows[0] or any(not key.strip() for key in rows[0]):
        raise ValueError(f"{path}: TSV needs nonempty named columns.")
    header = rows.pop(0)
    header_names = set(header)
    if len(header_names) != len(header):
        raise ValueError(f"{path}: duplicated TSV headers.")
    if any(key not in header_names for key in required):
        raise ValueError(f"{path}: required columns: {', '.join(required)}.")
    for line, row in enumerate(rows, 2):
        if len(row) != len(header):
            raise ValueError(f"{path}: row {line} has incorrect field count.")
    if not rows:
        raise ValueError(f"{path}: TSV needs at least one data row.")
    return header, rows


def read_gene_family_table(path, required):
    header, rows = read_gene_family_rows(path, required)
    return header, [dict(zip(header, row, strict=True)) for row in rows]

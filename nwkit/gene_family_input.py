"""Strict literal-keyed TSV input for gene-count and Ks family analyses."""

import csv
from io import StringIO

from nwkit.util import read_input_text


def read_gene_family_table(path, required):
    reader = csv.reader(StringIO(read_input_text(path)), delimiter="\t")
    rows = list(reader)
    if not rows or not rows[0] or any(not key.strip() for key in rows[0]):
        raise ValueError(f"{path}: TSV needs nonempty named columns.")
    header = rows[0]
    if len(set(header)) != len(header):
        raise ValueError(f"{path}: duplicated TSV headers.")
    if any(key not in header for key in required):
        raise ValueError(f"{path}: required columns: {', '.join(required)}.")
    result = []
    for line, row in enumerate(rows[1:], 2):
        if len(row) != len(header):
            raise ValueError(f"{path}: row {line} has incorrect field count.")
        result.append(dict(zip(header, row, strict=True)))
    if not result:
        raise ValueError(f"{path}: TSV needs at least one data row.")
    return header, result

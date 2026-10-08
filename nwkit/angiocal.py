"""Versioned AngioCal inputs; fossil ages remain in their original Ma units."""

import csv
import hashlib
import math
import os
from dataclasses import dataclass
from decimal import Decimal, InvalidOperation
from importlib.resources import files
from pathlib import Path

import requests

from nwkit.output_transaction import output_transaction
from nwkit.util import COMMON_ETE_CACHE_DIRS, resolve_download_dir, resolve_ete_data_dir

ANGIOCAL_VERSION = "v1.0"
ANGIOCAL_URL = (
    "https://raw.githubusercontent.com/eflowerproject/angiocal/"
    "263e31ce22cbd9df6c7e8644406e859f43f59d70/Data2b_CalibrationList.xls"
)
ANGIOCAL_SHA256 = "7e6e0b8ba8dfc8b58812a9c38cd2741cb3ca9f534bea096cf4b041279c3f1ee0"
# The normalized records belong to these exact source bytes, independently of
# the selected download's verification policy.
ANGIOCAL_NORMALIZED_SOURCE_SHA256 = ANGIOCAL_SHA256
ANGIOCAL_NORMALIZED_SHA256 = (
    "4f121728b4959c41f9c2fda6117f5ab311177c0b176e6380ad7b853ed86cb5ba"
)
ANGIOCAL_MAX_BYTES = 5 * 1024 * 1024
SUMMARY_COLUMNS = {
    "NFos": "fossil_id",
    "Fossil taxon": "fossil_taxon",
    "Minimum age": "minimum_age_ma",
    "Node calibrated": "node_calibrated",
    "Crown or stem?": "placement",
    "Clade calibrated (= node calibrated if by crown)": "clade",
    "Safe minimum age (Ma)": "safe_minimum_age",
    "Age quality score": "age_quality_score",
    "Node assignment score": "node_assignment_score",
    "Reconciliation score": "reconciliation_score",
    "Full reference (fossil relationships)": "relationship_reference",
    "Full reference (age)": "age_reference",
}
REQUIRED_TSV_COLUMNS = {
    "fossil_id",
    "fossil_taxon",
    "minimum_age_ma",
    "placement",
    "clade",
}


@dataclass(frozen=True)
class FossilCalibration:
    fossil_id: str
    fossil_taxon: str
    minimum_age_ma: float
    placement: str
    clade: str
    source_row: int
    node_calibrated: str = ""
    safe_minimum_age: str = ""
    age_quality_score: str = ""
    node_assignment_score: str = ""
    reconciliation_score: str = ""
    relationship_reference: str = ""
    age_reference: str = ""


@dataclass(frozen=True)
class AngioCalDataset:
    records: tuple[FossilCalibration, ...]
    path: str
    source: str
    sha256: str
    version: str = ANGIOCAL_VERSION


def positive_number(value, label):
    try:
        number = float(value)
    except (ValueError, TypeError, OverflowError) as exc:
        raise ValueError(f"{label} must be numeric.") from exc
    if not math.isfinite(number) or number <= 0:
        raise ValueError(f"{label} must be finite and positive.")
    return number


def fossil_id_text(value):
    # Keep the existing finite numeric range, but do not round identifiers or
    # accept fractional decimal text through a binary floating-point conversion.
    positive_number(value, "AngioCal fossil ID")
    try:
        number = Decimal(str(value).strip())
    except InvalidOperation as exc:
        raise ValueError("AngioCal fossil IDs must be positive integers.") from exc
    if number != number.to_integral_value():
        raise ValueError("AngioCal fossil IDs must be positive integers.")
    return str(int(number))


def _record(row, source_row):
    fields = {
        key: str(row[key]).strip() if row.get(key) is not None else ""
        for key in FossilCalibration.__dataclass_fields__
        if key not in {"source_row", "minimum_age_ma"}
    }
    fields["fossil_id"] = fossil_id_text(row.get("fossil_id"))
    if fields["placement"] not in {"crown", "stem"}:
        raise ValueError(f"Unknown AngioCal placement at row {source_row}.")
    node_words = fields["node_calibrated"].lower().split()
    if (
        node_words
        and node_words[0] in {"crown", "stem"}
        and node_words[0] != fields["placement"]
    ):
        raise ValueError(
            f"Conflicting AngioCal crown/stem placement at row {source_row}."
        )
    if not fields["clade"] or not fields["fossil_taxon"]:
        raise ValueError(f"Empty AngioCal clade or fossil taxon at row {source_row}.")
    age = positive_number(row.get("minimum_age_ma"), "AngioCal minimum age")
    return FossilCalibration(**fields, minimum_age_ma=age, source_row=source_row)


def _validate_headers(headers, required):
    if any(not name for name in headers) or len(set(headers)) != len(headers):
        raise ValueError("AngioCal headers must be nonempty and unique.")
    missing = required - set(headers)
    if missing:
        raise ValueError("Missing AngioCal columns: " + ", ".join(sorted(missing)))


def _read_xls(data):
    try:
        import xlrd
    except ModuleNotFoundError as exc:
        if exc.name != "xlrd":
            raise
        raise ValueError(
            "Custom AngioCal XLS files require the optional XLS reader: "
            "pip install 'nwkit[xls]'. Alternatively, use a normalized TSV. "
            "The pinned official v1.0 workbook does not require xlrd."
        ) from exc
    try:
        workbook = xlrd.open_workbook(file_contents=data, on_demand=True)
    except xlrd.XLRDError as exc:
        raise ValueError("Invalid AngioCal XLS workbook.") from exc
    try:
        if workbook.nsheets != 1:
            raise ValueError("AngioCal v1.0 requires exactly one worksheet.")
        sheet = workbook.sheet_by_index(0)
        if sheet.nrows < 3:
            raise ValueError("AngioCal worksheet has no calibration records.")
        headers = [str(value).strip() for value in sheet.row_values(1)]
        _validate_headers(headers, set(SUMMARY_COLUMNS))
        mapped_columns = [headers.index(name) for name in SUMMARY_COLUMNS]
        numeric_columns = [headers.index(name) for name in ("NFos", "Minimum age")]
        for index in range(2, sheet.nrows):
            values = sheet.row_values(index)
            if not any(value != "" for value in values):
                continue
            if any(
                sheet.cell_type(index, column) == xlrd.XL_CELL_ERROR
                for column in mapped_columns
            ):
                raise ValueError(f"Excel error in AngioCal worksheet row {index + 1}.")
            if any(
                sheet.cell_type(index, column)
                not in {xlrd.XL_CELL_NUMBER, xlrd.XL_CELL_TEXT}
                for column in numeric_columns
            ):
                raise ValueError(
                    f"AngioCal IDs and ages require numeric or decimal-text Excel cells at row {index + 1}; "
                    "dates and booleans are not fossil values."
                )
            row = dict(zip(headers, values, strict=True))
            yield _record(
                {dest: row[source] for source, dest in SUMMARY_COLUMNS.items()},
                index + 1,
            )
    finally:
        workbook.release_resources()


def _read_tsv(data):
    from io import StringIO

    reader = csv.DictReader(StringIO(data.decode("utf-8-sig")), delimiter="\t")
    _validate_headers(reader.fieldnames or [], REQUIRED_TSV_COLUMNS)
    for row in reader:
        if None in row or any(value is None for value in row.values()):
            raise ValueError(f"Malformed AngioCal TSV row {reader.line_num}.")
        yield _record(row, reader.line_num)


def _read_official_records():
    from io import StringIO

    data = files("nwkit").joinpath("data_angiocal").joinpath("v1.0.tsv").read_bytes()
    if hashlib.sha256(data).hexdigest() != ANGIOCAL_NORMALIZED_SHA256:
        raise ValueError("Bundled AngioCal v1.0 checksum verification failed.")
    reader = csv.DictReader(StringIO(data.decode("utf-8")), delimiter="\t")
    _validate_headers(
        reader.fieldnames or [], set(FossilCalibration.__dataclass_fields__)
    )
    for row in reader:
        yield _record(row, int(row["source_row"]))


def read_angiocal(path, *, official=False):
    declared_path = os.fspath(path)
    path = os.path.realpath(path)
    with open(path, "rb") as handle:
        data = handle.read(ANGIOCAL_MAX_BYTES + 1)
    if len(data) > ANGIOCAL_MAX_BYTES:
        raise ValueError("AngioCal input exceeds the size limit.")
    digest = hashlib.sha256(data).hexdigest()
    if official and digest != ANGIOCAL_SHA256:
        raise ValueError("AngioCal v1.0 checksum verification failed.")
    # The original BIFF workbook uses an OLE container. Detect it before text
    # decoding, including downloads or symlink targets without an extension.
    is_xls = data.startswith(
        b"\xd0\xcf\x11\xe0\xa1\xb1\x1a\xe1"
    ) or declared_path.lower().endswith(".xls")
    if is_xls and digest == ANGIOCAL_NORMALIZED_SOURCE_SHA256:
        records = tuple(_read_official_records())
    else:
        records = tuple(_read_xls(data) if is_xls else _read_tsv(data))
    if not records:
        raise ValueError("AngioCal input has no calibration records.")
    ids = [record.fossil_id for record in records]
    if len(set(ids)) != len(ids):
        raise ValueError("Duplicate AngioCal fossil IDs.")
    return AngioCalDataset(records, path, ANGIOCAL_URL if official else path, digest)


def _download_angiocal():
    try:
        with requests.get(ANGIOCAL_URL, stream=True, timeout=30) as response:
            response.raise_for_status()
            chunks = []
            size = 0
            for chunk in response.iter_content(chunk_size=64 * 1024):
                size += len(chunk)
                if size > ANGIOCAL_MAX_BYTES:
                    raise ValueError("AngioCal download exceeds the size limit.")
                chunks.append(chunk)
    except requests.RequestException as exc:
        raise RuntimeError("Failed to download AngioCal v1.0.") from exc
    data = b"".join(chunks)
    if hashlib.sha256(data).hexdigest() != ANGIOCAL_SHA256:
        raise ValueError("AngioCal v1.0 checksum verification failed.")
    return data


def angiocal_input_path(args):
    """Resolve the source before opening audit logs, outputs, or the network."""
    local_path = getattr(args, "angiocal_file", None)
    if local_path:
        return Path(local_path)
    base = resolve_download_dir(args)
    if base is None:
        base = os.path.join(os.path.expanduser("~"), ".cache", "nwkit")
    return Path(base) / "angiocal" / ANGIOCAL_VERSION / "Data2b_CalibrationList.xls"


def angiocal_protected_paths(args):
    paths = [("AngioCal source/cache", angiocal_input_path(args))]
    if getattr(args, "angiocal_taxonomy", True):
        from ete4.ncbi_taxonomy.ncbiquery import DEFAULT_TAXADB

        configured = resolve_ete_data_dir(args) or os.path.dirname(DEFAULT_TAXADB)
        for directory in set(COMMON_ETE_CACHE_DIRS) | {configured}:
            for filename in (
                "taxa.sqlite",
                "taxa.sqlite.traverse.pkl",
                "taxdump.tar.gz",
                ".ete4_taxonomy.lock",
            ):
                paths.append(("NCBI taxonomy cache", Path(directory) / filename))
    return paths


def load_angiocal(args):
    if getattr(args, "angiocal", "no") != ANGIOCAL_VERSION:
        raise ValueError("Supported AngioCal version: v1.0.")
    local_path = getattr(args, "angiocal_file", None)
    if local_path:
        if local_path == "-":
            raise ValueError("'--angiocal-file' requires an XLS or TSV file path.")
        return read_angiocal(local_path)
    path = angiocal_input_path(args)
    if not path.exists():
        data = _download_angiocal()
        with output_transaction([path], create_parents=True) as staged:
            staged.write_bytes(path, data)
    return read_angiocal(path, official=True)

"""Input and restart metadata for the image acquisition workflow."""

import csv
import hashlib
import json
import math
import threading
import time
from pathlib import Path
from typing import Any

from nwkit.image_metadata import normalize_species_name
from nwkit.output_transaction import output_transaction
from nwkit.taxonomic_labels import (
    TAXONOMIC_INFRASPECIFIC_RANKS,
    base_species_label,
    extract_species_label,
    taxonomic_token_key,
)

IMAGE_MANIFEST_FIELDS = (
    "leaf_name",
    "species_name",
    "provider",
    "provider_record_id",
    "matched_name",
    "matched_rank",
    "is_primary",
    "is_vector",
    "selection_reason",
    "license_code",
    "license_url",
    "attribution",
    "source_page_url",
    "media_url",
    "local_path",
    "score",
    "status",
)
UNMATCHED_FIELDS = ("leaf_name", "species_name", "reason", "details")
RANK_NAMES = {"subsp": "subspecies", "var": "variety", "forma": "form"}


def declared_rank(name):
    tokens = (normalize_species_name(name) or "").split()
    if len(tokens) >= 4:
        rank = taxonomic_token_key(tokens[2])
        if rank in TAXONOMIC_INFRASPECIFIC_RANKS:
            return RANK_NAMES.get(rank, rank)
        if tokens[-1].lower() in {"group", "subgroup"}:
            return tokens[-1].lower()
    if len(tokens) == 1:
        return "genus"
    if len(tokens) >= 3 and tokens[2].startswith(("'", "‘", "(")):
        return "cultivar"
    return "species"


def infraspecific_label(leaf):
    label = extract_species_label(str(leaf).replace(" ", "_"))
    if declared_rank(label) not in {"species", "genus"}:
        return normalize_species_name(label)
    full_name = normalize_species_name(leaf)
    if full_name and "(" in full_name and declared_rank(full_name) == "cultivar":
        return full_name
    return None


def parent_species(name):
    label = (
        "_".join(name.split()[:2])
        if declared_rank(name) == "cultivar"
        else name.replace(" ", "_")
    )
    return normalize_species_name(base_species_label(label))


def read_verified_rows(path):
    required = {"leaf_name", "source_page_url", "attribution"}
    with open(path, encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = reader.fieldnames or []
        if (
            len(set(fields)) != len(fields)
            or not all(field.strip() for field in fields)
            or not required.issubset(fields)
            or not {"local_path", "media_url"}.intersection(fields)
            or not {"license_code", "license_url"}.intersection(fields)
        ):
            raise ValueError(
                "--verified-images-tsv requires leaf_name, source_page_url, "
                "attribution, a media_url or local_path, and a license_code or license_url."
            )
        for row in reader:
            if None in row or None in row.values():
                raise ValueError(
                    f"--verified-images-tsv row {reader.line_num} has the wrong field count."
                )
            yield {key: value.strip() for key, value in row.items()}


def unique_taxon_aliases(connection, name, limit=2):
    owners_sql = "SELECT taxid FROM species WHERE spname = ? UNION SELECT taxid FROM synonym WHERE spname = ?"
    owners = connection.execute(owners_sql, (name, name)).fetchall()
    if len(owners) != 1:
        return []
    taxid = owners[0][0]
    names = connection.execute(
        """SELECT spname, 0 FROM species WHERE taxid = ?
        UNION ALL SELECT spname, 1 FROM synonym WHERE taxid = ?
        ORDER BY 2, 1""",
        (taxid, taxid),
    )
    aliases: list[str] = []
    seen = {name.casefold()}
    for alias, _ in names:
        if alias.casefold() in seen or len(alias.split()) < 2:
            continue
        if connection.execute(owners_sql, (alias, alias)).fetchall() == [(taxid,)]:
            aliases.append(alias)
            seen.add(alias.casefold())
        if len(aliases) == limit:
            break
    return aliases


def parse_provider_intervals(values, sources):
    intervals = {"inaturalist": 1.0}
    for value in values or ():
        provider, _, raw = value.partition("=")
        try:
            seconds = float(raw)
        except ValueError:
            seconds = -1
        if provider not in sources or not math.isfinite(seconds) or seconds < 0:
            raise ValueError(
                "--provider-interval requires SOURCE=SECONDS with finite, nonnegative seconds."
            )
        intervals[provider] = max(intervals.get(provider, 0), seconds)
    return intervals


class RequestPacer:
    """Share minimum API request spacing across this run's lookup threads."""

    def __init__(self, intervals):
        self.intervals = intervals
        self.locks = {
            source: threading.Lock() for source, delay in intervals.items() if delay > 0
        }
        self.next_request: dict[str, float] = {}

    def wait(self, source):
        if source in self.locks:
            with self.locks[source]:
                delay = self.next_request.get(source, 0) - time.monotonic()
                if delay > 0:
                    time.sleep(delay)
                self.next_request[source] = time.monotonic() + self.intervals[source]


def sha256_file(path):
    digest = hashlib.sha256()
    with open(path, "rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _fingerprint(value):
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, ensure_ascii=False).encode("utf-8")
    ).hexdigest()


class ImageCheckpoints:
    """Checkpoint validated output per query; failed or partial results stay retryable."""

    def __init__(self, out_dir, images_dir, args, sources):
        self.out_dir = Path(out_dir)
        self.images_dir = Path(images_dir).resolve()
        self.directory = self.out_dir / ".nwkit-cache" / "image-resume"
        self.refresh = bool(getattr(args, "refresh_cache", False))
        self.options = {
            option: getattr(args, option, None)
            for option in (
                "style",
                "fallback_rank",
                "license_max",
                "allow_nd",
                "max_per_species",
                "max_download_bytes",
                "output_format",
                "max_edge",
                "trim",
                "trim_shape",
                "canvas",
                "background",
                "synonyms",
                "wikimedia_depicts",
                "first_available",
                "download_dir",
                "taxonomy_cache_max_age_days",
                "refresh_taxonomy_cache",
            )
        }
        self.options["sources"] = list(sources)

    def _path(self, query, candidate):
        return self.directory / (
            _fingerprint([1, self.options, query, candidate]) + ".json"
        )

    def _image_path(self, row):
        path = (self.out_dir / row["local_path"]).resolve()
        if path.is_relative_to(self.images_dir) and path.is_file():
            return path
        raise ValueError("Checkpoint image is missing or outside the images directory.")

    def load(self, query, candidate, count):
        if self.refresh:
            return None
        try:
            data = json.loads(self._path(query, candidate).read_text(encoding="utf-8"))
            if not isinstance(data, list) or len(data) != count:
                return None
            rows = []
            for entry in data:
                row = entry["row"]
                if row.keys() != set(IMAGE_MANIFEST_FIELDS) or not all(
                    isinstance(value, str) for value in row.values()
                ):
                    return None
                if sha256_file(self._image_path(row)) != entry["sha256"]:
                    return None
                rows.append(row)
            return rows
        except (OSError, TypeError, ValueError, AttributeError, KeyError, RuntimeError):
            return None

    def save(self, query, candidate, rows, count):
        if len(rows) != count:
            return
        destination = self._path(query, candidate)
        entries = [
            {"row": row, "sha256": sha256_file(self._image_path(row))} for row in rows
        ]
        with output_transaction([destination], create_parents=True) as staged:
            staged.write_text(destination, lambda stream: json.dump(entries, stream))


def copy_selection_to_leaves(rows, leaves, *, resumed=False):
    selected: list[dict[str, Any]] = []
    for row in rows:
        for leaf in leaves:
            entry = dict(row, leaf_name=leaf)
            if resumed:
                entry["status"] = "resumed"
            selected.append(entry)
    return selected

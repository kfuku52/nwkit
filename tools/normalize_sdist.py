"""Normalize source-distribution archives for byte-for-byte reproducibility."""

import argparse
import copy
import gzip
import os
import tarfile
import tempfile
from pathlib import Path

MAX_GZIP_TIMESTAMP = (1 << 32) - 1


def parse_source_date_epoch(value: str) -> int:
    """Return a gzip-compatible non-negative SOURCE_DATE_EPOCH value."""
    try:
        epoch = int(value)
    except (TypeError, ValueError) as exc:
        raise ValueError("SOURCE_DATE_EPOCH must be an integer.") from exc
    if not 0 <= epoch <= MAX_GZIP_TIMESTAMP:
        raise ValueError(
            "SOURCE_DATE_EPOCH must be between 0 and {}.".format(MAX_GZIP_TIMESTAMP)
        )
    return epoch


def _normalized_mode(member: tarfile.TarInfo) -> int:
    if member.isdir():
        return 0o755
    if member.isfile():
        return 0o644
    if member.issym() or member.islnk():
        return 0o777
    raise ValueError(
        "Source distribution contains unsupported member type: '{}'".format(member.name)
    )


def _normalized_member(member: tarfile.TarInfo, epoch: int) -> tarfile.TarInfo:
    normalized = copy.copy(member)
    normalized.uid = 0
    normalized.gid = 0
    normalized.uname = ""
    normalized.gname = ""
    normalized.mtime = epoch
    normalized.mode = _normalized_mode(member)
    normalized.pax_headers = {}
    normalized.devmajor = 0
    normalized.devminor = 0
    return normalized


def normalize_sdist(path: Path, epoch: int) -> None:
    """Replace a ``.tar.gz`` sdist with a deterministic equivalent archive."""
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError("Source distribution does not exist: {}".format(path))

    temporary_path = None
    try:
        with tarfile.open(path, "r:gz") as source:
            members = sorted(source.getmembers(), key=lambda item: item.name)
            member_names = [member.name for member in members]
            if len(member_names) != len(set(member_names)):
                raise ValueError("Source distribution contains duplicate member names.")

            with tempfile.NamedTemporaryFile(
                dir=path.parent,
                prefix=".{}-".format(path.name),
                suffix=".tmp",
                delete=False,
            ) as raw_output:
                temporary_path = Path(raw_output.name)
                with gzip.GzipFile(
                    filename="",
                    mode="wb",
                    compresslevel=9,
                    fileobj=raw_output,
                    mtime=epoch,
                ) as compressed_output:
                    with tarfile.open(
                        fileobj=compressed_output,
                        mode="w",
                        format=tarfile.PAX_FORMAT,
                    ) as destination:
                        for member in members:
                            normalized = _normalized_member(member, epoch)
                            if member.isfile():
                                contents = source.extractfile(member)
                                if contents is None:
                                    raise RuntimeError(
                                        "Could not read sdist member: '{}'".format(
                                            member.name
                                        )
                                    )
                                with contents:
                                    destination.addfile(normalized, contents)
                            else:
                                destination.addfile(normalized)
        temporary_path.chmod(0o644)
        os.replace(temporary_path, path)
        temporary_path = None
    finally:
        if temporary_path is not None:
            temporary_path.unlink(missing_ok=True)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("sdist", type=Path)
    args = parser.parse_args()
    value = os.environ.get("SOURCE_DATE_EPOCH")
    if value is None:
        parser.error("SOURCE_DATE_EPOCH is required.")
    try:
        epoch = parse_source_date_epoch(value)
        normalize_sdist(args.sdist, epoch)
    except (FileNotFoundError, RuntimeError, ValueError) as exc:
        parser.error(str(exc))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

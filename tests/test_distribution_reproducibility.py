import gzip
import importlib.util
import os
import struct
import tarfile
from io import BytesIO
from pathlib import Path

import pytest


def _load_normalizer():
    path = Path(__file__).resolve().parents[1] / "tools" / "normalize_sdist.py"
    spec = importlib.util.spec_from_file_location("nwkit_normalize_sdist", path)
    if spec is None or spec.loader is None:
        raise RuntimeError("Could not load the sdist normalizer.")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


NORMALIZER = _load_normalizer()
normalize_sdist = NORMALIZER.normalize_sdist
parse_source_date_epoch = NORMALIZER.parse_source_date_epoch


def _write_archive(
    path: Path,
    *,
    reverse: bool,
    gzip_mtime: int,
    owner: tuple[int, int, str, str],
) -> None:
    members = [
        ("package-1.0", None, 0o700),
        ("package-1.0/data.txt", b"stable contents\n", 0o600),
        ("package-1.0/run", b"#!/bin/sh\n", 0o700),
    ]
    if reverse:
        members.reverse()
    with path.open("wb") as raw_output:
        with gzip.GzipFile(
            filename=path.name,
            mode="wb",
            fileobj=raw_output,
            mtime=gzip_mtime,
        ) as compressed_output:
            with tarfile.open(fileobj=compressed_output, mode="w") as archive:
                for index, (name, contents, mode) in enumerate(members):
                    member = tarfile.TarInfo(name)
                    member.type = (
                        tarfile.DIRTYPE if contents is None else tarfile.REGTYPE
                    )
                    member.size = 0 if contents is None else len(contents)
                    member.mode = mode
                    member.mtime = gzip_mtime + index
                    member.uid, member.gid, member.uname, member.gname = owner
                    archive.addfile(
                        member,
                        None if contents is None else BytesIO(contents),
                    )


def test_normalize_sdist_is_byte_reproducible_across_archive_metadata(tmp_path):
    first = tmp_path / "first.tar.gz"
    second = tmp_path / "second.tar.gz"
    _write_archive(
        first,
        reverse=False,
        gzip_mtime=123,
        owner=(501, 20, "developer", "staff"),
    )
    _write_archive(
        second,
        reverse=True,
        gzip_mtime=456,
        owner=(1001, 1001, "runner", "runner"),
    )

    normalize_sdist(first, 42)
    normalize_sdist(second, 42)

    assert first.read_bytes() == second.read_bytes()
    if os.name != "nt":
        assert first.stat().st_mode & 0o777 == 0o644
    assert struct.unpack("<I", first.read_bytes()[4:8])[0] == 42
    with tarfile.open(first, "r:gz") as archive:
        members = archive.getmembers()
        assert [member.name for member in members] == sorted(
            member.name for member in members
        )
        assert {member.uid for member in members} == {0}
        assert {member.gid for member in members} == {0}
        assert {member.uname for member in members} == {""}
        assert {member.gname for member in members} == {""}
        assert {member.mtime for member in members} == {42}
        assert archive.getmember("package-1.0").mode == 0o755
        assert archive.getmember("package-1.0/data.txt").mode == 0o644
        assert archive.getmember("package-1.0/run").mode == 0o644
        assert (
            archive.extractfile("package-1.0/data.txt").read() == b"stable contents\n"
        )


@pytest.mark.parametrize("value", ["not-an-integer", "-1", str(1 << 32)])
def test_parse_source_date_epoch_rejects_invalid_values(value):
    with pytest.raises(ValueError, match="SOURCE_DATE_EPOCH"):
        parse_source_date_epoch(value)

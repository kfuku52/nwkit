"""Exercise official records, custom XLS and TSV commands without xlrd."""

import builtins
import hashlib
import shutil
import struct
import subprocess
import sys
from pathlib import Path

import pytest

from nwkit import angiocal
from nwkit.cli import parser
from nwkit.mcmctree import mcmctree_main


def block_xlrd(monkeypatch):
    original_import = builtins.__import__

    def without_xlrd(name, *args, **kwargs):
        if name == "xlrd" or name.startswith("xlrd."):
            raise ModuleNotFoundError("No module named 'xlrd'", name="xlrd")
        return original_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", without_xlrd)


def test_official_records_keep_original_rows_and_metadata():
    records = tuple(angiocal._read_official_records())
    assert len(records) == 238
    assert [record.source_row for record in records] == list(range(3, 241))
    assert len({record.fossil_id for record in records}) == 238
    assert all(
        record.relationship_reference and record.age_reference for record in records
    )


def test_verified_official_xls_uses_records_without_xlrd(tmp_path, monkeypatch):
    fixture = Path(__file__).parent / "data" / "angiocal-v1.0-synthetic.xls"
    path = tmp_path / "official-workbook"
    path.write_bytes(fixture.read_bytes())
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    monkeypatch.setattr(angiocal, "ANGIOCAL_SHA256", digest)
    monkeypatch.setattr(angiocal, "ANGIOCAL_NORMALIZED_SOURCE_SHA256", digest)
    block_xlrd(monkeypatch)

    result = angiocal.read_angiocal(path, official=True)

    assert result.records == tuple(angiocal._read_official_records())
    assert result.path == str(path.resolve())
    assert result.source == angiocal.ANGIOCAL_URL
    assert result.sha256 == digest


def test_official_checksum_mismatch_cannot_use_bundled_records(tmp_path, monkeypatch):
    path = tmp_path / "modified.xls"
    path.write_bytes(b"not the pinned source workbook")
    block_xlrd(monkeypatch)
    with pytest.raises(ValueError, match="v1.0 checksum verification failed"):
        angiocal.read_angiocal(path, official=True)


def test_corrupt_bundled_records_are_rejected(tmp_path, monkeypatch):
    directory = tmp_path / "data_angiocal"
    directory.mkdir()
    (directory / "v1.0.tsv").write_text("corrupt records", encoding="utf-8")
    monkeypatch.setattr(angiocal, "files", lambda package: tmp_path)
    with pytest.raises(ValueError, match="Bundled AngioCal v1.0 checksum"):
        tuple(angiocal._read_official_records())


@pytest.mark.skipif(shutil.which("git") is None, reason="Git checkout requires Git")
def test_bundled_checksum_survives_windows_git_checkout(tmp_path, monkeypatch):
    repo = Path(__file__).resolve().parents[1]
    relative = "nwkit/data_angiocal/v1.0.tsv"
    subprocess.run(["git", "init", "-q", str(tmp_path)], check=True)
    shutil.copy2(repo / ".gitattributes", tmp_path / ".gitattributes")
    target = tmp_path / relative
    target.parent.mkdir(parents=True)
    shutil.copy2(repo / relative, target)
    command = ["git", "-C", str(tmp_path), "-c", "core.autocrlf=true"]
    subprocess.run(command + ["add", ".gitattributes", relative], check=True)
    target.unlink()
    subprocess.run(command + ["checkout-index", "--", relative], check=True)
    monkeypatch.setattr(angiocal, "files", lambda package: tmp_path / "nwkit")

    assert len(tuple(angiocal._read_official_records())) == 238


def test_custom_xls_uses_native_reader_without_xlrd(tmp_path, monkeypatch):
    fixture = Path(__file__).parent / "data" / "angiocal-v1.0-synthetic.xls"
    tree, report = tmp_path / "tree.tre", tmp_path / "report.tsv"
    tree.write_text("previous tree", encoding="utf-8")
    report.write_text("previous report", encoding="utf-8")
    args = parser.parse_args(
        [
            "mcmctree",
            "-i",
            "((A_a,B_b)Testaceae,C_c);",
            "--angiocal",
            "v1.0",
            "--angiocal-file",
            str(fixture),
            "--angiocal-taxonomy",
            "no",
            "--report",
            str(report),
            "-o",
            str(tree),
        ]
    )
    block_xlrd(monkeypatch)

    mcmctree_main(args)
    assert "L(12.5," in tree.read_text()
    assert "synthetic metadata" in report.read_text()


def test_fresh_cli_and_tsv_import_work_without_xlrd(tmp_path):
    path = tmp_path / "fossils.tsv"
    path.write_text(
        "fossil_id\tfossil_taxon\tminimum_age_ma\tplacement\tclade\n"
        "1\tSynthetic fossil\t12.5\tcrown\tTestaceae\n",
        encoding="utf-8",
    )
    script = """
import importlib.abc
import sys

class WithoutXlrd(importlib.abc.MetaPathFinder):
    def find_spec(self, fullname, path, target=None):
        if fullname == 'xlrd' or fullname.startswith('xlrd.'):
            raise ModuleNotFoundError("No module named 'xlrd'", name='xlrd')

sys.meta_path.insert(0, WithoutXlrd())
from nwkit import angiocal
from nwkit.cli import main
assert len(tuple(angiocal._read_official_records())) == 238
main(['mcmctree', '-i', '((A_a,B_b)Testaceae,C_c);',
      '--angiocal', 'v1.0', '--angiocal-file', sys.argv[1],
      '--angiocal-taxonomy', 'no', '-o', sys.argv[2]])
main(['mcmctree', '-i', '((A_a,B_b)Testaceae,C_c);',
      '--angiocal', 'v1.0', '--angiocal-file', sys.argv[3],
      '--angiocal-taxonomy', 'no', '-o', sys.argv[4]])
assert 'xlrd' not in sys.modules
"""
    result = subprocess.run(
        [
            sys.executable,
            "-c",
            script,
            str(path),
            str(tmp_path / "tree.tre"),
            str(Path(__file__).parent / "data" / "angiocal-v1.0-synthetic.xls"),
            str(tmp_path / "xls-tree.tre"),
        ],
        capture_output=True,
        text=True,
        timeout=30,
        check=False,
    )
    assert result.returncode == 0, result.stderr
    assert "L(12.5," in (tmp_path / "tree.tre").read_text()
    assert "L(12.5," in (tmp_path / "xls-tree.tre").read_text()


@pytest.mark.parametrize("corruption", ["truncated-sector", "overlapping-mini-stream"])
def test_corrupt_custom_xls_preserves_outputs_without_xlrd(
    tmp_path, monkeypatch, corruption
):
    fixture = Path(__file__).parent / "data" / "angiocal-v1.0-synthetic.xls"
    corrupt = tmp_path / "corrupt.xls"
    data = bytearray(fixture.read_bytes())
    if corruption == "truncated-sector":
        data = data[:-1]
        message = "truncated OLE sector"
    else:
        directory_start = struct.unpack_from("<I", data, 48)[0]
        directory_offset = (directory_start + 1) * 512
        workbook_start = struct.unpack_from("<I", data, directory_offset + 128 + 116)[0]
        workbook_size = struct.unpack_from("<Q", data, directory_offset + 128 + 120)[0]
        struct.pack_into(
            "<IQ", data, directory_offset + 116, workbook_start, workbook_size
        )
        message = "stream allocations overlap"
    corrupt.write_bytes(data)
    tree, report = tmp_path / "tree.tre", tmp_path / "report.tsv"
    tree.write_text("previous tree", encoding="utf-8")
    report.write_text("previous report", encoding="utf-8")
    args = parser.parse_args(
        [
            "mcmctree",
            "-i",
            "((A_a,B_b)Testaceae,C_c);",
            "--angiocal",
            "v1.0",
            "--angiocal-file",
            str(corrupt),
            "--angiocal-taxonomy",
            "no",
            "--report",
            str(report),
            "-o",
            str(tree),
        ]
    )
    block_xlrd(monkeypatch)
    with pytest.raises(ValueError, match=message):
        mcmctree_main(args)
    assert tree.read_text() == "previous tree"
    assert report.read_text() == "previous report"

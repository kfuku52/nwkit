#!/usr/bin/env python3
"""Build deterministic pre-submission archives for Zenodo and Dryad."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import hashlib
from pathlib import Path, PurePosixPath
import stat
import subprocess
import zipfile


PROJECT_ROOT = Path(__file__).resolve().parents[2]
PAPER_ROOT = PROJECT_ROOT / "paper"
FIXED_TIMESTAMP = (2026, 7, 16, 0, 0, 0)
VERSION = "0.27.0"
NWKIT_SOURCE_COMMIT = "f71dc345ac830a195992abadbf4b7a7394837991"


ZENODO_README = """# NWKIT manuscript code package

This archive contains the NWKIT 0.27.0 source snapshot and the scripts used to
generate the Systematic Biology manuscript analyses, result tables, figures,
supplement, and submission documents.

The `nwkit-source` directory is the exact Git tree at the commit recorded in
`NWKIT_SOURCE_COMMIT.txt`. The `manuscript-analysis` directory contains the
manuscript source, analysis scripts, environment specification, comparison
decisions, and reproduction instructions at `MANUSCRIPT_COMMIT.txt`. Data
inputs and generated result tables are in the separate Dryad package.

Create the environment and run the analyses in the order documented in
`manuscript-analysis/analysis/README.md`. All source and manuscript-analysis
code in this package is distributed under the NWKIT MIT License.
"""


DRYAD_README = """# NWKIT manuscript data package

This archive contains the inputs and generated results supporting the NWKIT
Systematic Biology manuscript. It includes checksum-pinned public PEPC inputs,
all tabular and tree outputs from the worked example, correctness and
predefined input-case results, benchmark measurements, generated manuscript figures,
and the human-readable supplement.

The source and redistribution information for the PEPC files is in
`data/README.md`; their upstream MIT license is in `data/CSUBST_LICENSE`.
Analysis code and the NWKIT source snapshot are provided in the separate Zenodo
code package. `MANIFEST.tsv` records the size and SHA-256 digest of every file
in this archive.
"""


@dataclass(frozen=True)
class ArchiveMember:
    archive_path: PurePosixPath
    data: bytes
    executable: bool = False


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def current_commit() -> str:
    return subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=PROJECT_ROOT,
        check=True,
        capture_output=True,
        text=True,
    ).stdout.strip()


def verify_git_state() -> None:
    paper_status = subprocess.run(
        ["git", "status", "--porcelain", "--untracked-files=all", "--", "paper"],
        cwd=PROJECT_ROOT,
        check=True,
        capture_output=True,
        text=True,
    ).stdout.strip()
    if paper_status:
        raise ValueError("Commit all paper/ changes before building deposit archives")

    changed = subprocess.run(
        ["git", "diff", "--name-only", NWKIT_SOURCE_COMMIT, "HEAD", "--"],
        cwd=PROJECT_ROOT,
        check=True,
        capture_output=True,
        text=True,
    ).stdout.splitlines()
    unexpected = [
        path for path in changed
        if path not in {".gitattributes", ".gitignore"} and not path.startswith("paper/")
    ]
    if unexpected:
        raise ValueError(
            "NWKIT source changed after the manuscript source commit: "
            + ", ".join(unexpected)
        )


def git_tree_members(commit: str, prefix: PurePosixPath) -> list[ArchiveMember]:
    output = subprocess.run(
        ["git", "ls-tree", "-r", "-z", commit],
        cwd=PROJECT_ROOT,
        check=True,
        capture_output=True,
    ).stdout
    members: list[ArchiveMember] = []
    for entry in output.split(b"\0"):
        if not entry:
            continue
        metadata, path_bytes = entry.split(b"\t", 1)
        mode, object_type, object_id = metadata.split()
        if object_type != b"blob":
            continue
        path_text = path_bytes.decode()
        data = subprocess.run(
            ["git", "cat-file", "blob", object_id.decode()],
            cwd=PROJECT_ROOT,
            check=True,
            capture_output=True,
        ).stdout
        members.append(
            ArchiveMember(
                prefix / path_text,
                data,
                executable=mode == b"100755",
            )
        )
    return members


def member_from_file(path: Path, archive_path: PurePosixPath) -> ArchiveMember:
    mode = path.stat().st_mode
    return ArchiveMember(archive_path, path.read_bytes(), bool(mode & stat.S_IXUSR))


def tree_members(root: Path, prefix: PurePosixPath) -> list[ArchiveMember]:
    members = []
    for path in sorted(item for item in root.rglob("*") if item.is_file()):
        if "__pycache__" in path.parts or path.suffix == ".pyc":
            continue
        members.append(member_from_file(path, prefix / path.relative_to(root).as_posix()))
    return members


def manifest(members: list[ArchiveMember]) -> bytes:
    rows = ["path\tbytes\tsha256"]
    for member in sorted(members, key=lambda item: str(item.archive_path)):
        rows.append(f"{member.archive_path}\t{len(member.data)}\t{sha256(member.data)}")
    return ("\n".join(rows) + "\n").encode()


def write_archive(path: Path, members: list[ArchiveMember]) -> None:
    unique_paths = [str(member.archive_path) for member in members]
    if len(unique_paths) != len(set(unique_paths)):
        raise ValueError(f"Duplicate paths requested for {path.name}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(path, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=9) as archive:
        for member in sorted(members, key=lambda item: str(item.archive_path)):
            info = zipfile.ZipInfo(str(member.archive_path), FIXED_TIMESTAMP)
            info.compress_type = zipfile.ZIP_DEFLATED
            permissions = 0o755 if member.executable else 0o644
            info.external_attr = (stat.S_IFREG | permissions) << 16
            archive.writestr(info, member.data)


def verify_archive(path: Path) -> None:
    with zipfile.ZipFile(path) as archive:
        bad_members = archive.testzip()
        if bad_members is not None:
            raise ValueError(f"CRC failure in {path.name}: {bad_members}")
        names = set(archive.namelist())
        if "MANIFEST.tsv" not in names:
            raise ValueError(f"Missing manifest in {path.name}")
        rows = archive.read("MANIFEST.tsv").decode().splitlines()
        if not rows or rows[0] != "path\tbytes\tsha256":
            raise ValueError(f"Invalid manifest header in {path.name}")
        manifest_names = set()
        for row in rows[1:]:
            name, size_text, expected_hash = row.split("\t")
            data = archive.read(name)
            if len(data) != int(size_text) or sha256(data) != expected_hash:
                raise ValueError(f"Manifest mismatch for {name} in {path.name}")
            manifest_names.add(name)
        if manifest_names != names - {"MANIFEST.tsv"}:
            raise ValueError(f"Manifest membership mismatch in {path.name}")


def zenodo_members(manuscript_commit: str) -> list[ArchiveMember]:
    members = [
        ArchiveMember(PurePosixPath("README.md"), ZENODO_README.encode()),
        ArchiveMember(
            PurePosixPath("NWKIT_SOURCE_COMMIT.txt"),
            f"{NWKIT_SOURCE_COMMIT}\n".encode(),
        ),
        ArchiveMember(
            PurePosixPath("MANUSCRIPT_COMMIT.txt"),
            f"{manuscript_commit}\n".encode(),
        ),
    ]
    members.extend(git_tree_members(NWKIT_SOURCE_COMMIT, PurePosixPath("nwkit-source")))
    for path in sorted((PAPER_ROOT / "analysis").iterdir()):
        if path.is_file() and path.suffix in {".py", ".tsv", ".yml", ".md"}:
            members.append(member_from_file(path, PurePosixPath("manuscript-analysis/analysis") / path.name))
    for name in (
        "README.md",
        "manuscript.md",
        "references.bib",
        "claims.tsv",
        "positioning.md",
        "published_uses.tsv",
        "table1.md",
        "supplement.md",
        "submission_readiness.md",
        "submission_metadata.md",
        "author_metadata_sources.md",
        "cover_letter.md",
        "deposit_metadata.md",
    ):
        path = PAPER_ROOT / name
        members.append(member_from_file(path, PurePosixPath("manuscript-analysis") / name))
    return members


def dryad_members() -> list[ArchiveMember]:
    members = [ArchiveMember(PurePosixPath("README.md"), DRYAD_README.encode())]
    members.extend(tree_members(PAPER_ROOT / "data", PurePosixPath("data")))
    members.extend(tree_members(PAPER_ROOT / "results", PurePosixPath("results")))
    members.extend(tree_members(PAPER_ROOT / "figures", PurePosixPath("figures")))
    members.append(member_from_file(PAPER_ROOT / "supplement.md", PurePosixPath("supplement.md")))
    members.append(
        member_from_file(PAPER_ROOT / "analysis" / "README.md", PurePosixPath("REPRODUCTION.md"))
    )
    members.append(member_from_file(PAPER_ROOT / "analysis" / "software_versions.tsv", PurePosixPath("provenance/software_versions.tsv")))
    members.append(member_from_file(PAPER_ROOT / "analysis" / "capability_matrix.tsv", PurePosixPath("provenance/capability_matrix.tsv")))
    return members


def add_manifest(members: list[ArchiveMember]) -> list[ArchiveMember]:
    return [*members, ArchiveMember(PurePosixPath("MANIFEST.tsv"), manifest(members))]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, default=PROJECT_ROOT / "output" / "archive")
    args = parser.parse_args()

    verify_git_state()
    manuscript_commit = current_commit()
    zenodo_path = args.output_dir / f"NWKIT_Zenodo_code_{VERSION}.zip"
    dryad_path = args.output_dir / f"NWKIT_Dryad_data_{VERSION}.zip"
    write_archive(zenodo_path, add_manifest(zenodo_members(manuscript_commit)))
    write_archive(dryad_path, add_manifest(dryad_members()))
    verify_archive(zenodo_path)
    verify_archive(dryad_path)

    checksums = [
        f"{sha256(zenodo_path.read_bytes())}  {zenodo_path.name}",
        f"{sha256(dryad_path.read_bytes())}  {dryad_path.name}",
    ]
    checksum_path = args.output_dir / "SHA256SUMS"
    checksum_path.write_text("\n".join(checksums) + "\n")
    print(checksum_path.read_text(), end="")


if __name__ == "__main__":
    main()

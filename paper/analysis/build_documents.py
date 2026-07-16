#!/usr/bin/env python3
"""Build and audit the Systematic Biology submission documents."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile
from urllib.request import urlopen
import zipfile


PROJECT_ROOT = Path(__file__).resolve().parents[2]
PAPER_ROOT = PROJECT_ROOT / "paper"
DEFAULT_OUTPUT_DIR = PROJECT_ROOT / "output" / "doc"
CSL_COMMIT = "439002f8dbcc44acd99e72c05fd6a6880c509d84"
CSL_SHA256 = "ad9701ec3f576f590661efcc5890201cb25bd9fb123dfa932b59fe2279f1bc95"
CSL_URL = (
    "https://raw.githubusercontent.com/citation-style-language/styles/"
    f"{CSL_COMMIT}/systematic-biology.csl"
)


def require_program(name: str, override: Path | None = None) -> str:
    if override is not None:
        if not override.is_file():
            raise FileNotFoundError(f"Program does not exist: {override}")
        return str(override)
    path = shutil.which(name)
    if path is None:
        raise FileNotFoundError(f"Required program is not on PATH: {name}")
    return path


def run(command: list[str], *, capture_output: bool = False) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        command,
        cwd=PROJECT_ROOT,
        check=True,
        capture_output=capture_output,
        text=True,
    )


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def fetch_csl(destination: Path, supplied_path: Path | None) -> None:
    if supplied_path is None:
        with urlopen(CSL_URL, timeout=60) as response:  # noqa: S310 - pinned HTTPS URL
            data = response.read()
    else:
        data = supplied_path.read_bytes()
    observed = hashlib.sha256(data).hexdigest()
    if observed != CSL_SHA256:
        raise ValueError(
            "Systematic Biology CSL checksum mismatch: "
            f"expected {CSL_SHA256}, observed {observed}"
        )
    destination.write_bytes(data)


def manuscript_word_count(pandoc: str, csl_path: Path) -> int:
    result = run([
        pandoc,
        str(PAPER_ROOT / "manuscript.md"),
        "--citeproc",
        f"--csl={csl_path}",
        f"--resource-path={PAPER_ROOT}",
        "--metadata",
        "link-citations=false",
        "--to=plain",
    ], capture_output=True)
    words: list[str] = []
    for line in result.stdout.splitlines():
        if line.strip() == "References":
            break
        words.extend(line.split())
    return len(words)


def pdf_pages(pdfinfo: str, path: Path) -> int:
    result = run([pdfinfo, str(path)], capture_output=True)
    match = re.search(r"^Pages:\s+(\d+)\s*$", result.stdout, flags=re.MULTILINE)
    if match is None:
        raise ValueError(f"Could not read page count from {path}")
    return int(match.group(1))


def docx_metrics(path: Path) -> dict[str, int]:
    with zipfile.ZipFile(path) as archive:
        bad_member = archive.testzip()
        if bad_member is not None:
            raise ValueError(f"CRC failure in {path}: {bad_member}")
        document_xml = archive.read("word/document.xml").decode()
        media_count = sum(name.startswith("word/media/") for name in archive.namelist())
    return {
        "line_number_sections": document_xml.count("w:lnNumType"),
        "embedded_media": media_count,
        "alt_text_attributes": document_xml.count("descr="),
    }


def build_documents(
    output_dir: Path,
    pandoc: str,
    soffice: str,
    pdfinfo: str,
    supplied_csl: Path | None,
) -> dict[str, object]:
    output_dir.mkdir(parents=True, exist_ok=True)
    main_docx = output_dir / "NWKIT_Systematic_Biology_draft.docx"
    main_pdf = output_dir / "NWKIT_Systematic_Biology_draft.pdf"
    supplement_docx = output_dir / "NWKIT_Systematic_Biology_supplement.docx"
    supplement_pdf = output_dir / "NWKIT_Systematic_Biology_supplement.pdf"

    with tempfile.TemporaryDirectory(prefix="nwkit-doc-build-") as temporary_directory:
        temporary_root = Path(temporary_directory)
        csl_path = temporary_root / "systematic-biology.csl"
        raw_main = temporary_root / "main-raw.docx"
        raw_supplement = temporary_root / "supplement-raw.docx"
        fetch_csl(csl_path, supplied_csl)

        run([
            pandoc,
            str(PAPER_ROOT / "manuscript.md"),
            "--citeproc",
            f"--csl={csl_path}",
            f"--resource-path={PAPER_ROOT}",
            "--metadata",
            "link-citations=false",
            f"--output={raw_main}",
        ])
        run([
            sys.executable,
            str(PAPER_ROOT / "analysis" / "format_docx.py"),
            str(raw_main),
            str(main_docx),
        ])
        run([
            pandoc,
            str(PAPER_ROOT / "supplement.md"),
            f"--resource-path={PAPER_ROOT}",
            f"--output={raw_supplement}",
        ])
        run([
            sys.executable,
            str(PAPER_ROOT / "analysis" / "format_docx.py"),
            "--supplement",
            str(raw_supplement),
            str(supplement_docx),
        ])

        for pdf_path in (main_pdf, supplement_pdf):
            pdf_path.unlink(missing_ok=True)
        run([
            soffice,
            "--headless",
            "--convert-to",
            "pdf",
            "--outdir",
            str(output_dir),
            str(main_docx),
        ])
        run([
            soffice,
            "--headless",
            "--convert-to",
            "pdf",
            "--outdir",
            str(output_dir),
            str(supplement_docx),
        ])

        word_count = manuscript_word_count(pandoc, csl_path)

    manuscript_text = (PAPER_ROOT / "manuscript.md").read_text()
    main_displays = len(re.findall(r"^!\[Figure \d+\]", manuscript_text, flags=re.MULTILINE))
    main_displays += len(re.findall(r"^\*\*Table \d+\.", manuscript_text, flags=re.MULTILINE))
    main_pages = pdf_pages(pdfinfo, main_pdf)
    supplement_pages = pdf_pages(pdfinfo, supplement_pdf)
    metrics = docx_metrics(main_docx)

    failures: list[str] = []
    if not 3000 <= word_count <= 4000:
        failures.append(f"main word count outside 3,000--4,000: {word_count}")
    if not 12 <= main_pages <= 14:
        failures.append(f"main page count outside 12--14: {main_pages}")
    if main_displays > 4:
        failures.append(f"more than four main displays: {main_displays}")
    if metrics["line_number_sections"] != 1:
        failures.append("main DOCX does not contain exactly one continuous line-number section")
    if metrics["embedded_media"] != 3:
        failures.append(f"main DOCX does not contain exactly three images: {metrics['embedded_media']}")
    if metrics["alt_text_attributes"] < 3:
        failures.append("main DOCX has fewer than three embedded alt-text attributes")
    if failures:
        raise ValueError("; ".join(failures))

    report: dict[str, object] = {
        "csl_commit": CSL_COMMIT,
        "csl_sha256": CSL_SHA256,
        "main_displays": main_displays,
        "main_docx_sha256": sha256_file(main_docx),
        "main_pages": main_pages,
        "main_pdf_sha256": sha256_file(main_pdf),
        "main_word_count_before_references": word_count,
        "supplement_docx_sha256": sha256_file(supplement_docx),
        "supplement_pages": supplement_pages,
        "supplement_pdf_sha256": sha256_file(supplement_pdf),
        **metrics,
    }
    (output_dir / "build_report.json").write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n"
    )
    return report


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--csl", type=Path, help="Use a local CSL file after checksum verification")
    parser.add_argument("--pandoc", type=Path)
    parser.add_argument("--soffice", type=Path)
    parser.add_argument("--pdfinfo", type=Path)
    args = parser.parse_args()

    report = build_documents(
        output_dir=args.output_dir,
        pandoc=require_program("pandoc", args.pandoc),
        soffice=require_program("soffice", args.soffice),
        pdfinfo=require_program("pdfinfo", args.pdfinfo),
        supplied_csl=args.csl,
    )
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

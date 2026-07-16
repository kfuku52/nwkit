# NWKIT manuscript workspace

This directory contains the working materials for a *Systematic Biology*
Software for Systematics and Evolution submission.

## Source of truth

`manuscript.md` is the canonical manuscript while the scientific argument,
analyses, and figures are changing rapidly. Citations use Pandoc-style keys
from `references.bib`, for example `[@Lemoine2021]`. Generated Word and PDF
files will be placed under `build/` and should not be edited during this phase.

Submission-oriented Word and PDF files are generated under `output/doc/` and
must not be edited directly while Markdown remains the source of truth. If
author or co-author revision later moves to Word with tracked changes, the Word
document should become the canonical source at a clearly recorded handoff.

## Working files

- `positioning.md`: target journal, central claim, scope, and claim boundaries.
- `claims.tsv`: auditable claim-to-evidence register.
- `published_uses.tsv`: candidate peer-reviewed NWKIT uses and example status.
- `outline.md`: section-level argument and word budget.
- `manuscript.md`: continuously revised manuscript draft.
- `references.bib`: verified bibliographic metadata.
- `analysis/`: scripts and instructions for comparisons, examples, and benchmarks.
- `data/`: data provenance and retrieval instructions.
- `figures/`: editable figure sources and rendered outputs.
- `submission_readiness.md`: completed checks and author-supplied items still
  required before ScholarOne submission.
- `submission_metadata.md`: fields and file assignments for ScholarOne.
- `cover_letter.md`: journal-specific cover-letter draft.
- `deposit_metadata.md`: ready-to-copy Zenodo and Dryad descriptions.

## Journal constraints

The target article type is Software for Systematics and Evolution. The current
journal instructions specify 3,000--4,000 words, 12--14 manuscript pages, and
no more than four figures and tables in total. The software must be open source,
well documented, and supported for at least two years after publication. Data
used in the manuscript must be deposited in Dryad and analysis scripts/code in
Zenodo before submission.

## Drafting rule

Claims marked `open`, `partial`, or `author-confirmation` in `claims.tsv` must
not be promoted to unqualified statements. Unresolved submission metadata is
recorded in `submission_readiness.md` rather than inserted into the manuscript.

## Building submission documents

Run the version-pinned document builder from the repository root:

```bash
python paper/analysis/build_documents.py
```

The script downloads the Systematic Biology CSL style at citation-style-language
styles commit `439002f8dbcc44acd99e72c05fd6a6880c509d84`, verifies its SHA-256
digest, builds both DOCX and PDF files, and writes `output/doc/build_report.json`.
It fails if the main document leaves the 3,000--4,000-word or 12--14-page range,
exceeds four main displays, or loses continuous line numbering, figures, or
embedded alternative text. Pandoc, Poppler (`pdfinfo`), and LibreOffice
(`soffice`) must be available on `PATH`; the first two are included in
`analysis/environment.yml`.

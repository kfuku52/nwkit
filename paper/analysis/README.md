# Reproducing the NWKIT manuscript analyses

Run all commands from the repository root. The manuscript environment was
created with Python 3.12, ETE 4.4.0, and Biopython 1.87. On macOS, the
conda-forge builds of ETE and Biopython avoid incorrect dynamic-library paths
that can be retained by locally built wheels.

```bash
conda env create -f paper/analysis/environment.yml
conda activate nwkit-paper
python -m pip install -e '.[test]'
python -m pip install dendropy==5.0.10 treeswift==1.1.45 phykit==2.3.0
```

Gotree 0.5.2 is built separately from its tagged source:

```bash
git clone --branch v0.5.2 https://github.com/evolbioinfo/gotree.git /tmp/gotree-v0.5.2
make -C /tmp/gotree-v0.5.2 build
```

Generate the reported outputs in this order:

```bash
python paper/analysis/inventory.py
python paper/analysis/pipeline_smoke.py
python paper/analysis/evaluate_correctness.py
python paper/analysis/worked_example.py
python paper/analysis/benchmark.py --repetitions 3 \
  --gotree /tmp/gotree-v0.5.2/gotree
python paper/analysis/make_figures.py
python paper/analysis/build_supplement.py
python paper/analysis/build_documents.py
python paper/analysis/build_archives.py
```

The scripts write raw results to `paper/results`, figures in PDF, SVG, and PNG
to `paper/figures`, pinned PEPC inputs to `paper/data/pepc`, and the generated
supplement to `paper/supplement.md`. `worked_example.py` verifies both downloaded
PEPC files against SHA-256 hashes before analysis. `evaluate_correctness.py`
exits nonzero if any predefined output or input-case check fails. Benchmark
failures are retained as rows with nonzero return codes rather than discarded.

`build_documents.py` pins and checksum-verifies the Systematic Biology CSL
style, builds the submission DOCX/PDF files, and audits their word count, page
count, main-display count, line numbering, embedded figures, and alternative
text. LibreOffice must be installed separately to provide `soffice`.

`build_archives.py` creates deterministic, checksum-manifested ZIP files under
`output/archive`: a Zenodo package containing the NWKIT source snapshot and
analysis code, and a Dryad package containing the pinned inputs, raw results,
figures, and supplement. The generated `SHA256SUMS` file is verified before
upload.

The capability comparison is a documented audit, not an automated inference
from command names. Its long-form source is `capability_matrix.tsv`; every main-
table classification should receive author review immediately before
submission because competitor software may change.

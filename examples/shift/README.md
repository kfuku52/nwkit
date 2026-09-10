# OU shift examples and validation evidence

`tree.nwk`, `traits.tsv` and `traits-with-se.tsv` are small command examples;
see [SHIFT.md](../../SHIFT.md) for executable examples and model conventions.

The `validation/` directory contains an 80-dataset paired BIC/pBIC simulation
pilot on eight-tip trees. `validation-16tip/` contains a 96-dataset pBIC extension
with weak and strong effects. See [SHIFT_VALIDATION.md](../../SHIFT_VALIDATION.md)
for results, denominators, caveats and reproduction commands.

Each evidence directory contains:

- `inputs.json`: trees, observations, seeds and full generating truth per case.
- `records.json` and `records.csv`: every attempted outer fit, including failures.
- `cell-summary.json`: condition-specific metrics and Monte Carlo intervals.
- `scenario-summary.json` and `table.md`: descriptive summaries across root/SE
  settings, keeping criterion, tip count and effect magnitude separate.
- `manifests.json`: generating grids, runtime metadata and execution-source hashes.
- `export.json` and `source-snapshot.json`: evidence-export source hashes and
  the actual source text indexed by SHA-256, including original simulation code
  and the exporter audit. Export checks the raw inputs, saved fits and summaries
  before publishing the complete bundle.

Canonical data are JSON. With pandas, read the CSV using
`pd.read_csv(path, keep_default_na=False)` to preserve the scenario name `null`.
Bootstrap frequencies use successful inner refits; failures and availability
have separate denominators. This pilot does not establish nominal error control.

The `joint-validation/` bundle contains the held-out small-tree joint-search
comparison described in [SHIFT_JOINT.md](../../SHIFT_JOINT.md). `records.json`
contains 14 outer fits; `candidates.json` preserves all 4,138 candidate attempts
including four failures. `models.json`, `inputs.json`, `manifest.json`,
`fit-manifests.json`, `source-snapshot.json` and `export.json` retain model
tables, inputs, execution sources and exporter provenance. pBIC comparisons
that fail representation equivalence have a null score improvement.

`pbic-correction/` stores the minimal before/after reproduction and backend
source hashes. `joint-validation-pbic-fixed/` replays the original joint grid
with the locally corrected backend; it is not independent replication. See
[SHIFT_PBIC.md](../../SHIFT_PBIC.md). Historical pBIC bundles describe the old
criterion and must not be treated as calibration of its corrected implementation.

The current default uses [calibrated selection](../../SHIFT_CALIBRATION.md).
`calibration-validation/` contains a separately seeded, frozen-protocol study of
that implementation. All earlier IC simulation commands explicitly request
`--selection ic` and retain their original interpretation.

Current alpha-envelope evidence is in `calibration-envelope/`,
`calibration-envelope-stress/` and `calibration-weak-null/`. Earlier calibration
directories retain the earlier plug-in method. See the
[review and next plan](../../reviews/shift-calibration-review.md).

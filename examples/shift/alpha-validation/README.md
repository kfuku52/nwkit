# Independent paired OU sensitivity evidence

Open [report.html](report.html) for the technical report. The full design,
definitions and reproduction commands are in [SHIFT_ALPHA.md](../../../SHIFT_ALPHA.md).
This bundle uses 520 new datasets and a locally corrected, unreleased kfl1ou
backend. It must not be confused with the earlier 12-dataset regression replay.

## Bundle contents

- `protocol.json`: frozen design; `jobs.jsonl`: completion status for every dataset.
- `inputs.jsonl.gz`: all generating truths, observations, and trees.
- `records.jsonl.gz`: all 4,160 selected-model records, mean effects, tip
  predictions, optima, and independent audits.
- `candidate-ledger.jsonl.gz`: every joint candidate attempt, criterion score, warning,
  and failure. Each fit provides both BIC and pBIC; do not double its sample size.
- `summary.json`: condition-level rates, Wilson intervals, failure bounds,
  mean RMSE, and alpha-boundary frequencies.
- `paired.json`: bound, criterion, and method comparisons on identical datasets.
- `audit.json`: completeness and independent numerical audit results.
- `null-profile-diagnostic.json`: post-hoc independent null-likelihood check,
  explicitly separate from the prespecified experiment.
- `backend-probe.tsv`, `installed-backend-sha256.json`, and
  `backend-correction-provenance.json`: exact backend identity and prior correction.
- `source/`: the driver, design, and generator snapshot saved before fitting.
- `analysis-source/`: final audit, aggregation, diagnostic, and report scripts.
- `artifact.json`, `report-data.json`: canonical report and reviewed plotted data.
- `report-delivery.json`: portable renderer verification receipt.
- `validation.json`: independent count/input checks and check-suite outcomes.
- `SHA256SUMS.json`: evidence file hashes, excluding the checksum file itself.

Raw RDS objects, full model tip/effect tables, and backend logs remain in the raw
run directory. The compressed evidence is intended for inspecting and reproducing
the study without committing thousands of RDS files. All selected point models
are auditable from the retained parameters, predictions, inputs and source, while a full
independent rerun requires the corrected R backend described in the main guide.

## Report and chart contract

The report audience is technical. Its section roles are title, technical summary,
metric definitions, quantitative findings, model/experimental design, numerical
validation, limitations, next steps and open questions. Definitions precede
findings to make denominators and effective-shift counting explicit. The null
profile check is a separate post-hoc validation section. Source metadata and
this file retain implementation details rather than adding an unrelated final
sources section to the reader-facing report.

The four grouped bar charts compare two alpha bounds within BIC/pBIC and
two-stage/joint methods. Root models use separate full-width panels. Repeated
bar charts are intentional: each asks the same bounded category-comparison
question, with false-positive or exact partition recovery as the response.
Each chart has eight plotted rows, paired sample context, denominator and interval
data, a percent scale, a 100% reference, and an adjacent interpretation. Two
series are identified by their legend and common left/right order. Native shared
renderer styling supplies the palette; no custom HTML chart runtime is used.
Tables serve exact lookup for the complete grid and paired discordance counts;
they are not intended as independent observations for statistical inference.

The HTML report uses the Data Analytics plugin's canonical portable reader.
The SQL recorded in its metadata is an actual read of reviewed Python aggregates
materialized in an in-memory SQLite table; the underlying fits and statistical
aggregation are performed by the saved Python/R scripts. The delivery receipt
records artifact validation, exact embedded-payload equality, chart/source
interaction, desktop and narrow-width rendering, and browser errors.

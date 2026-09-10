import argparse
import gzip
import hashlib
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "tools"))
from shift_known_error_envelope import known_error_envelope  # noqa: E402

from nwkit.shift_calibration import CalibratedSearch  # noqa: E402
from nwkit.util import read_tree  # noqa: E402

parser = argparse.ArgumentParser(
    description="Replay exploratory envelope probes into a new output file"
)
parser.add_argument("--output", type=Path, required=True)
args = parser.parse_args()
if args.output.exists():
    parser.error("Output already exists; choose a new path")

source = ROOT / "examples/shift/calibration-envelope-stress/records.jsonl.gz"
with gzip.open(source, "rt") as f:
    rows = [json.loads(x) for x in f]
output = []
for row in rows:
    if row["case"]["case_id"] not in (93, 105):
        continue
    tree = read_tree(row["tree"], "auto", True, quiet=True)
    search = CalibratedSearch(tree, convergence=False, variances=row["variances"])
    result = known_error_envelope(
        search, row["values"], seed=row["case"]["seed"] + 7000000000, max_evaluations=8
    )
    output.append(
        {
            "case_id": row["case"]["case_id"],
            "historical_plugin_p": row["fit"]["tests"][0]["p_value"],
            "envelope": result,
        }
    )
    print("case", row["case"]["case_id"], "done", flush=True)
result_text = (
    json.dumps(
        {
            "scope": "Exploratory probes selected because old plug-in tests rejected; not a calibration rate study. Budget 8.",
            "source_records_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
            "prototype_sha256": hashlib.sha256(
                (ROOT / "tools/shift_known_error_envelope.py").read_bytes()
            ).hexdigest(),
            "rows": output,
        },
        indent=2,
    )
    + "\n"
)

with args.output.open("x") as stream:
    stream.write(result_text)

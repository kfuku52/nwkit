"""Add fixed branch-Gaussian CLI and nucleotide RADTE measurements."""

import json
import sys
from pathlib import Path

from prepare import archive


def main():
    root = Path(sys.argv[1]).resolve()
    first = archive(root, "branch-cli-first", "97491cb")
    manifest = json.loads((root / "manifest.json").read_text())
    extra = []
    for n in [128, 512]:
        models = root / "inputs" / f"branch-models{n}.tsv"
        rows = ["branch_id\tmodel\tsigma2\talpha\ttheta\tjump_mean\tjump_variance"]
        for branch in range(1, 2 * n - 1):
            rows.append(
                f"{branch}\tOU\t1.2\t0.4\t0.2\t0.1\t0.05"
                if branch % 3 == 0
                else f"{branch}\tBM\t0.8\t\t\t\t"
            )
        models.write_text("\n".join(rows) + "\n")
        extra.append(
            dict(
                name=f"branch-gaussian-fixed-{n}",
                before=str(first),
                after=str(root / "after"),
                args=[
                    "asr",
                    "-i",
                    str(root / "inputs" / f"tree{n}.nwk"),
                    "--input-rooted",
                    "yes",
                    "--trait",
                    str(root / "inputs" / f"traits{n}.tsv"),
                    "--state-column",
                    "y",
                    "--model",
                    "BRANCH-GAUSSIAN",
                    "--branch-models",
                    str(models),
                    "--root-prior",
                    "gaussian",
                    "--root-mean",
                    "0",
                    "--root-variance",
                    "1",
                    "-o",
                    "result.tsv",
                ],
                tables=["result.tsv"],
            )
        )
    case = next(c for c in manifest["cases"] if c["name"] == "radte-branch-8")
    case["name"] = "radte-jc69-8"
    case["args"] += [
        "--alignment",
        str(root / "inputs/radte8/alignment.fasta"),
        "--substitution-model",
        "jc69",
    ]
    extra.append(case)
    manifest.update(cases=extra, runs=str(root / "additional-runs"))
    (root / "additional-manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n"
    )


if __name__ == "__main__":
    main()

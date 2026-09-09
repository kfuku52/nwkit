"""Time original R RADTE and native RADTE on the same bundled GeneRax input.

Models and output bundles differ. This is a workflow baseline, not an equivalent
estimator speedup comparison. R, ape, and the optional RADTE checkout are needed.
"""

import argparse
import json
import platform
import shutil
import sys
from pathlib import Path

from benchmark_radte import timed_run


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--radte-repo", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--repeats", type=int, default=3)
    args = parser.parse_args()
    repo, output = args.radte_repo.resolve(), args.outdir.resolve()
    if output.exists() or args.repeats < 1:
        parser.error("Use a new output directory and positive repeat count")
    rscript = shutil.which("Rscript")
    if rscript is None:
        parser.error("Rscript is not installed")
    data = repo / "data" / "example_generax_01"
    output.mkdir(parents=True)
    rows = []
    for repeat in range(-1, args.repeats):
        for method in ["R-RADTE", "native"]:
            prefix = output / f"{method}-{repeat}"
            if method == "R-RADTE":
                cmd = [
                    rscript,
                    str(repo / "radte.r"),
                    "--species_tree=" + str(data / "species_tree.nwk"),
                    "--generax_nhx=" + str(data / "gene_tree.nhx"),
                    "--species_node_bounds_tsv="
                    + str(data / "species_node_bounds.tsv"),
                    "--outdir=" + str(output),
                    "--prefix=" + prefix.name,
                    "--seed=1",
                    "--chronos_model=discrete",
                    "--chronos_lambda=1",
                    "--max_age=1000",
                    "--chronos_attempt_timeout_sec=5",
                    "--chronos_total_timeout_sec=15",
                ]
            else:
                cmd = [
                    sys.executable,
                    "-m",
                    "nwkit",
                    "radte",
                    "--generax-nhx",
                    str(data / "gene_tree.nhx"),
                    "--species-tree",
                    str(data / "species_tree.nwk"),
                    "--species-node-bounds-tsv",
                    str(data / "species_node_bounds.tsv"),
                    "--max-age",
                    "1000",
                    "--out-prefix",
                    str(prefix),
                ]
            row = timed_run(cmd, prefix)
            row.update(method=method, repeat=repeat)
            if repeat >= 0:
                rows.append(row)
            print(method, repeat, row["exit_code"], flush=True)
    (output / "results.json").write_text(
        json.dumps(dict(platform=platform.platform(), rows=rows), indent=2)
    )


if __name__ == "__main__":
    main()

"""MCMCtree multi-chain diagnostics using R posterior's reference estimators."""

import argparse
import math
import subprocess
from pathlib import Path


def diagnose_chains(
    paths, output, *, profile="chain", log=None, jobs=1, stopping_thresholds=None
):
    """Compute reference diagnostics, with an optional worker budget for chains."""
    if profile not in {"chain", "calibration"}:
        raise ValueError("Unknown posterior diagnostic profile")
    if not isinstance(jobs, int) or isinstance(jobs, bool) or not 1 <= jobs < 2**31:
        raise ValueError(
            "Posterior diagnostic jobs must be a positive integer below 2147483648"
        )
    if profile == "calibration" and jobs != 1:
        raise ValueError("Parallel diagnostics require the chain profile")
    if stopping_thresholds is not None:
        rhat, ess = stopping_thresholds
        if profile != "chain" or not all(
            math.isfinite(x) and x > 0 for x in (rhat, ess)
        ):
            raise ValueError(
                "Stopping thresholds must be positive and use the chain profile"
            )
    if len(paths) < 2:
        raise ValueError("At least two independent chains are required")
    name = (
        "mcmctree_diagnostics.R"
        if profile == "chain"
        else "mcmctree_calibration_diagnostics.R"
    )
    from nwkit.file_paths import (
        paths_identify_same_file,
        validate_outputs_do_not_replace_inputs,
    )
    from nwkit.output_transaction import output_transaction

    if any(
        paths_identify_same_file(left, right)
        for index, left in enumerate(paths)
        for right in paths[index + 1 :]
    ):
        raise ValueError("Independent chains must use distinct input files")
    validate_outputs_do_not_replace_inputs(
        [("chain", path) for path in paths], [("diagnostics", output)]
    )
    with output_transaction([output], create_parents=True) as staged:
        subprocess.run(
            [
                "Rscript",
                str(Path(__file__).with_name("data_r") / name),
                *(["--jobs", str(jobs)] if profile == "chain" else []),
                *(
                    ["--stop-at", *map(str, stopping_thresholds)]
                    if stopping_thresholds is not None
                    else []
                ),
                str(staged[output]),
                *map(str, paths),
            ],
            check=True,
            stdout=log,
            stderr=log,
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--profile", choices=("chain", "calibration"), default="chain")
    parser.add_argument(
        "--jobs", type=int, default=1, help="Worker budget for the chain profile"
    )
    parser.add_argument("output")
    parser.add_argument("chains", nargs="+")
    args = parser.parse_args()
    diagnose_chains(args.chains, args.output, profile=args.profile, jobs=args.jobs)


if __name__ == "__main__":
    main()

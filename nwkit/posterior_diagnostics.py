"""MCMCtree multi-chain diagnostics using R posterior's reference estimators."""

import argparse
import subprocess
from pathlib import Path


def diagnose_chains(paths, output, *, profile="chain"):
    if profile not in {"chain", "calibration"}:
        raise ValueError("Unknown posterior diagnostic profile")
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
                str(staged[output]),
                *map(str, paths),
            ],
            check=True,
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--profile", choices=("chain", "calibration"), default="chain")
    parser.add_argument("output")
    parser.add_argument("chains", nargs="+")
    args = parser.parse_args()
    diagnose_chains(args.chains, args.output, profile=args.profile)


if __name__ == "__main__":
    main()

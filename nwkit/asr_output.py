"""Publish reference and tree-ensemble ASR outputs only after successful fitting."""

import shutil
import sys
from contextlib import redirect_stdout
from copy import copy
from pathlib import Path
from tempfile import TemporaryFile

from nwkit.output_transaction import output_transaction


def run_ensemble_transaction(args, handler):
    from nwkit.asr import _validate_asr_output_paths

    _validate_asr_output_paths(args)
    paths = {
        name: path
        for name, path in vars(args).items()
        if (name == "outfile" or name.endswith("_out")) and path not in (None, "", "-")
    }
    staged_args = copy(args)
    staged_args._ensemble_output_staged = True
    if getattr(args, "figure_out", None):
        # The figure writer already uses this override for staged ASR outputs.
        staged_args._branch_figure_format = Path(args.figure_out).suffix.lower()[1:]
    # Keep large standard-output tables on disk until the entire fit succeeds.
    with TemporaryFile(mode="w+", encoding="utf-8") as captured:

        def publish_stdout():
            if not captured.tell():
                return
            captured.seek(0)
            shutil.copyfileobj(captured, sys.stdout)
            sys.stdout.flush()

        with output_transaction(paths.values(), after_install=publish_stdout) as staged:
            for name, path in paths.items():
                setattr(staged_args, name, staged[path])
            with redirect_stdout(captured):
                result = handler(staged_args)
    return result

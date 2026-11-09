"""Optional R/ggtree rendering backend preserving established plot styles."""

import subprocess
import sys
from pathlib import Path


def render_basic_tree(arguments):
    subprocess.run(
        [
            "Rscript",
            str(Path(__file__).with_name("data_r") / "basic_tree.R"),
            *arguments,
        ],
        check=True,
    )


if __name__ == "__main__":
    render_basic_tree(sys.argv[1:])

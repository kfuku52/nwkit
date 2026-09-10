"""Fit the fixed synthetic BM/OU layout and render conditional ASR histories.

Run from the repository root: python examples/branch_gaussian/fit/render.py
Observations, assignments and bounds are illustrative; diffusion parameters are fitted.
"""

from pathlib import Path

from nwkit.cli import main


def render():
    source = Path(__file__).resolve().parent
    output = source / "output"
    output.mkdir(exist_ok=True)
    common = [
        "asr",
        "--model",
        "BRANCH-GAUSSIAN",
        "-i",
        str(source / "tree.nwk"),
        "--input-rooted",
        "yes",
        "--branch-regimes",
        str(source / "branch_regimes.tsv"),
        "--regime-models",
        str(source / "regime_models.tsv"),
        "--root-prior",
        "gaussian",
        "--root-mean",
        "0",
        "--root-variance",
        "0.2",
        "--state-column",
        "x",
        "--figure-simulations",
        "6",
        "--figure-simulation-steps",
        "120",
        "--seed",
        "39",
        "--species-overlap-node-plot",
        "no",
        "--figure-width",
        "17",
        "--figure-height",
        "9",
    ]
    for extension in ("png", "svg", "pdf"):
        main(
            [
                *common,
                "--branch-fit",
                str(source / "fit.tsv"),
                "--trait",
                str(source / "traits.tsv"),
                "--standard-error-column",
                "se",
                "--figure-simulation-mode",
                "conditional",
                "--figure-tip-heatmap",
                "yes",
                "--figure-trait-tip-labels",
                "yes",
                "--figure-out",
                str(output / f"fitted.{extension}"),
                "-o",
                str(output / "asr.tsv"),
                "--model-out",
                str(output / "model.tsv"),
                "--process-out",
                str(output / "process.json"),
                "--branch-models-out",
                str(output / "branch-models.tsv"),
            ]
        )


if __name__ == "__main__":
    render()

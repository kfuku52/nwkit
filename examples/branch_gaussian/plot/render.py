"""Reproduce the synthetic mixed BM/OU/end-jump ASR and prior figures.

Run from the repository root: python examples/branch_gaussian/plot/render.py
Parameters and observed values are illustrative, not estimates from real data.
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
    for name, options in (
        (
            "posterior",
            [
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
            ],
        ),
        ("prior", ["--output", "prior-samples", "--prior-samples", "6"]),
    ):
        for extension in ("png", "svg", "pdf"):
            main(
                [
                    *common,
                    *options,
                    "--figure-out",
                    str(output / f"{name}.{extension}"),
                    "-o",
                    str(output / f"{name}.tsv"),
                    "--model-out",
                    str(output / f"{name}-model.tsv"),
                    "--process-out",
                    str(output / f"{name}-process.json"),
                    "--branch-models-out",
                    str(output / "branch-models.tsv"),
                ]
            )


if __name__ == "__main__":
    render()

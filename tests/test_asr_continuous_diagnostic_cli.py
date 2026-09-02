import pandas as pd
import pytest

from nwkit.cli import main


def run_bm_diagnostics(tmp_path, seed="12"):
    tmp_path.mkdir(parents=True, exist_ok=True)
    traits = tmp_path / "traits.tsv"
    output = tmp_path / "asr.tsv"
    samples = tmp_path / "samples.tsv"
    predictive = tmp_path / "predictive.tsv"
    bootstrap = tmp_path / "bootstrap.tsv"
    traits.write_text(
        "leaf_name\tvalue\tse\nA\t0\t0.1\nB\t1\t0.1\nC\t2\t0.2\nD\t4\t0.2\nE\tNA\tNA\n"
    )
    main(
        [
            "asr",
            "-i",
            "((A:1,B:1):1,(C:1,D:1):1,E:2)R;",
            "--input-rooted",
            "yes",
            "--trait",
            str(traits),
            "--state-column",
            "value",
            "--standard-error-column",
            "se",
            "--missing-values",
            "NA",
            "--posterior-samples-out",
            str(samples),
            "--posterior-samples",
            "12",
            "--posterior-predictive-out",
            str(predictive),
            "--posterior-predictive-simulations",
            "20",
            "--bootstrap-out",
            str(bootstrap),
            "--bootstrap-simulations",
            "3",
            "--seed",
            seed,
            "-o",
            str(output),
        ]
    )
    return output, samples, predictive, bootstrap


@pytest.mark.integration
def test_bm_cli_writes_reproducible_simulation_diagnostics(tmp_path):
    first = run_bm_diagnostics(tmp_path / "first")
    second = run_bm_diagnostics(tmp_path / "second")
    first_samples = pd.read_csv(first[1], sep="\t")
    second_samples = pd.read_csv(second[1], sep="\t")
    pd.testing.assert_frame_equal(first_samples, second_samples)
    assert len(first_samples) == 12 * 8
    assert set(first_samples.columns) == {
        "sample",
        "branch_id",
        "parent",
        "node_class",
        "name",
        "trait",
        "value",
    }
    predictive = pd.read_csv(first[2], sep="\t")
    assert set(predictive["statistic"]) == {
        "mean",
        "variance",
        "range",
        "max_abs_centered",
    }
    assert set(predictive["num_simulations"]) == {20}
    bootstrap = pd.read_csv(first[3], sep="\t")
    assert len(bootstrap) == 3
    assert set(bootstrap["fit_status"]) == {"ok"}
    assert bootstrap["sigma2"].notna().all()


def test_diagnostic_count_requires_its_output(tmp_path):
    traits = tmp_path / "traits.tsv"
    traits.write_text("leaf_name\tvalue\nA\t0\nB\t1\nC\t2\n")
    with pytest.raises(
        ValueError, match="--posterior-samples requires --posterior-samples-out"
    ):
        main(
            [
                "asr",
                "-i",
                "(A:1,B:1,C:1)R;",
                "--input-rooted",
                "yes",
                "--trait",
                str(traits),
                "--state-column",
                "value",
                "--posterior-samples",
                "10",
                "-o",
                str(tmp_path / "out.tsv"),
            ]
        )


def test_multivariate_diagnostics_are_rejected_explicitly(tmp_path):
    traits = tmp_path / "traits.tsv"
    traits.write_text("leaf_name\tx\ty\nA\t0\t1\nB\t1\t2\nC\t2\t3\n")
    with pytest.raises(ValueError, match="scalar continuous models"):
        main(
            [
                "asr",
                "-i",
                "(A:1,B:1,C:1)R;",
                "--input-rooted",
                "yes",
                "--trait",
                str(traits),
                "--state-column",
                "x,y",
                "--model",
                "MV-BM",
                "--posterior-samples-out",
                str(tmp_path / "samples.tsv"),
                "-o",
                str(tmp_path / "out.tsv"),
            ]
        )

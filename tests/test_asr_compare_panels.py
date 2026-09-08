"""Model panels reuse fits, preserve comparisons, and form one physical PDF page."""

import re

import numpy as np
import pandas as pd
import pytest

from nwkit.cli import main


def command(tmp_path, *extra):
    trait = tmp_path / "traits.tsv"
    trait.write_text(
        "leaf_name\tx\ty\nA\t0.1\t2.4\nB\t0.8\t1.7\nC\t1.2\t3.2\nD\t2.1\t2.1\nE\t-0.7\t0.4\nF\t-0.5\t2.0\n"
    )
    return [
        "asrcompare",
        "-i",
        "[&R]((A:1,B:1):1,(C:1,D:1):1,(E:1,F:1):1);",
        "--trait",
        str(trait),
        "--state-column",
        "x",
        "--models",
        "BM,OU",
        "--alpha",
        "0.8",
        "--sigma2",
        "0.5",
        "-o",
        str(tmp_path / "comparison.tsv"),
        *extra,
    ]


@pytest.mark.parametrize("tip_labels", ["no", "yes"])
def test_one_page_same_fit_values_and_shared_axes(tmp_path, monkeypatch, tip_labels):
    import nwkit.asr_compare as comparison
    import nwkit.asr_compare_panels as plotting

    main(command(tmp_path))
    before = pd.read_csv(tmp_path / "comparison.tsv", sep="\t")
    original_fit = comparison._fit_continuous
    original_build = plotting.build_comparison_panels
    calls = []

    def fit(context, candidate):
        calls.append(candidate.model_id)
        return original_fit(context, candidate)

    def build(context, table):
        fig = original_build(context, table)
        plots = [ax for ax in fig.axes if ax.axison]
        assert len(plots) == 6
        assert len({ax.get_ylim() for ax in plots}) == 1
        trait_plots = [ax for ax in plots if ax.get_xlabel() == "x"]
        assert len(trait_plots) == 4
        assert len({ax.get_xlim() for ax in trait_plots}) == 1
        if tip_labels == "yes":
            for ax in trait_plots:
                band = next(
                    child
                    for child in ax.child_axes
                    if child.get_label() == "trait-tip-labels"
                )
                assert band.get_xlim() == ax.get_xlim()
                labels = [text for text in band.texts if text.get_gid()]
                assert {text.get_text() for text in labels} == set("ABCDEF")
                if "simulation" in ax.get_title(loc="left"):
                    assert any(
                        "first history" in text.get_text() for text in band.texts
                    )
        # Cached marginals, rather than a second fit, supply every node.
        for posterior, _, _, _ in context.cache["figure_fits"].values():
            assert len(posterior) == len(list(context.tree.traverse()))
        return fig

    monkeypatch.setattr(comparison, "_fit_continuous", fit)
    monkeypatch.setattr(plotting, "build_comparison_panels", build)
    pdf = tmp_path / "panels.pdf"
    main(
        command(
            tmp_path,
            "--figure-layout",
            "panels",
            "--figure-out",
            str(pdf),
            "--figure-simulations",
            "2",
            "--figure-trait-tip-labels",
            tip_labels,
            "--figure-simulation-steps",
            "4",
            "--seed",
            "5",
        )
    )
    assert calls == ["BM", "OU[stationary]"]
    after = pd.read_csv(tmp_path / "comparison.tsv", sep="\t")
    times = [column for column in before if column.endswith("seconds")]
    pd.testing.assert_frame_equal(before.drop(columns=times), after.drop(columns=times))
    assert pdf.read_bytes().startswith(b"%PDF")
    assert len(re.findall(rb"/Type\s*/Page\b", pdf.read_bytes())) == 1


@pytest.mark.parametrize(
    "options,match",
    [
        (["--figure-layout", "panels"], "requires --figure-out"),
        (["--figure-simulations", "1"], "require --figure-layout"),
        (["--figure-width", "8"], "require --figure-layout"),
        (
            ["--figure-layout", "panels", "--figure-width", "-1"],
            "requires --figure-out",
        ),
    ],
)
def test_invalid_panel_controls_before_output(tmp_path, options, match):
    with pytest.raises(ValueError, match=match):
        main(command(tmp_path, *options))
    assert not (tmp_path / "comparison.tsv").exists()


def test_discrete_panels_rejected_before_fit(tmp_path):
    with pytest.raises(ValueError, match="continuous traits only"):
        main(
            command(
                tmp_path,
                "--trait-type",
                "discrete",
                "--models",
                "ER",
                "--figure-layout",
                "panels",
                "--figure-out",
                str(tmp_path / "a.pdf"),
            )
        )


def test_panel_failure_keeps_both_previous_outputs(tmp_path, monkeypatch):
    import nwkit.asr_compare_panels as plotting

    output = tmp_path / "comparison.tsv"
    pdf = tmp_path / "a.pdf"
    output.write_text("previous table")
    pdf.write_bytes(b"previous pdf")

    def fail(*args):
        raise RuntimeError("render failed")

    monkeypatch.setattr(plotting, "draw_comparison_panels", fail)
    with pytest.raises(RuntimeError, match="render failed"):
        main(command(tmp_path, "--figure-layout", "panels", "--figure-out", str(pdf)))
    assert output.read_text() == "previous table"
    assert pdf.read_bytes() == b"previous pdf"


def test_transformed_and_equivalent_candidates_remain_visible(tmp_path, monkeypatch):
    import nwkit.asr_compare_panels as plotting

    original = plotting.build_comparison_panels
    seen = []

    def build(context, table):
        # Also exercise an automatic-fit failure row without an expensive all-model fit.
        failed = dict(table.iloc[0])
        failed.update(
            model_id="failed-model",
            status="failed",
            message="Deliberate fit failure",
            comparison_group="",
            criterion_rank=pd.NA,
            criterion_value=np.nan,
        )
        equivalent = dict(
            failed,
            model_id="alias-model",
            status="equivalent",
            message="Statistically equivalent to BM",
        )
        fig = original(
            context,
            pd.concat([table, pd.DataFrame([failed, equivalent])], ignore_index=True),
        )
        text = "\n".join(t.get_text() for ax in fig.axes for t in ax.texts)
        assert "Branch simulation unavailable" in text
        assert "Deliberate fit failure" in text
        assert "equivalent" in text
        seen.append(True)
        return fig

    monkeypatch.setattr(plotting, "build_comparison_panels", build)
    main(
        command(
            tmp_path,
            "--models",
            "BM,OU,KAPPA",
            "--evolution-parameter",
            "0.5",
            "--alpha",
            "0.8",
            "--figure-layout",
            "panels",
            "--figure-out",
            str(tmp_path / "a.pdf"),
            "--figure-simulations",
            "1",
            "--figure-simulation-steps",
            "3",
        )
    )
    assert seen


def test_candidate_seed_is_stable_and_distinct():
    from nwkit.asr_compare_panels import _candidate_seed

    assert (
        np.random.default_rng(_candidate_seed(7, "BM")).normal()
        == np.random.default_rng(_candidate_seed(7, "BM")).normal()
    )
    assert (
        np.random.default_rng(_candidate_seed(7, "BM")).normal()
        != np.random.default_rng(_candidate_seed(7, "OU")).normal()
    )


def test_multivariate_comparison_retains_posterior_and_draws_each_trait(
    tmp_path, monkeypatch
):
    import nwkit.asr_compare_panels as plotting

    original = plotting.build_comparison_panels
    captured = []

    def build(context, table):
        figure = original(context, table)
        plots = [ax for ax in figure.axes if ax.axison]
        assert len(plots) == 5
        assert sum(ax.get_xlabel() == "x" for ax in plots) == 2
        assert sum(ax.get_xlabel() == "y" for ax in plots) == 2
        captured.append(True)
        return figure

    monkeypatch.setattr(plotting, "build_comparison_panels", build)
    args = command(
        tmp_path,
        "--models",
        "MV-OU-FULL",
        "--state-column",
        "x,y",
        "--attraction-matrix",
        "0.7,-0.2;0.3,1.1",
        "--diffusion-matrix",
        "1,0.3;0.3,1.2",
        "--figure-layout",
        "panels",
        "--figure-out",
        str(tmp_path / "multi.pdf"),
        "--figure-simulations",
        "1",
        "--figure-simulation-mode",
        "conditional",
        "--figure-simulation-steps",
        "3",
    )
    for flag in ("--alpha", "--sigma2"):
        index = args.index(flag)
        del args[index : index + 2]
    main(args)
    assert captured

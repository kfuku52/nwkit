"""Headless phylomorphospace and evolutionary-loading figures."""

from pathlib import Path

import numpy as np
from scipy.stats import chi2


def draw_pca(ancestors, loadings, eigenvalues, path, args):
    from matplotlib import rc_context
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.figure import Figure
    from matplotlib.patches import Ellipse

    from nwkit.asr_compare_figure import _font_family_for_text

    try:
        selected = [int(part) for part in args.figure_components.split(",")]
    except ValueError as exc:
        raise ValueError("--figure-components requires two component numbers.") from exc
    if (
        len(selected) != 2
        or len(set(selected)) != 2
        or any(i < 1 or i > len(eigenvalues) for i in selected)
    ):
        raise ValueError(
            "--figure-components requires two distinct available component numbers."
        )
    components = [f"PC{i}" for i in selected]
    first = ancestors[ancestors.component == components[0]].set_index("branch_id")
    second = (
        ancestors[ancestors.component == components[1]]
        .set_index("branch_id")
        .reindex(first.index)
    )
    family = _font_family_for_text(" ".join([*first.name, *loadings.trait]))
    settings = {"font.family": family} if family else {}
    with rc_context(settings):
        figure = Figure(
            figsize=(12, max(6, loadings.trait.nunique() * 0.23 + 2)),
            layout="constrained",
        )
        FigureCanvasAgg(figure)
        axes, heat = figure.subplots(1, 2, gridspec_kw={"width_ratios": [2.4, 1]})
        for identifier, row in first.iterrows():
            x, y = row["mean"], second.loc[identifier, "mean"]
            if row.parent_branch_id in first.index:
                parent = row.parent_branch_id
                axes.plot(
                    [first.loc[parent, "mean"], x],
                    [second.loc[parent, "mean"], y],
                    color="#8996a3",
                    linewidth=1,
                    zorder=1,
                )
            if row.node_class != "leaf":
                radius = chi2.ppf(args.ci_level, 2)
                ellipse = Ellipse(
                    (x, y),
                    2 * np.sqrt(radius * row.variance),
                    2 * np.sqrt(radius * second.loc[identifier, "variance"]),
                    facecolor="#0072B2",
                    edgecolor="none",
                    alpha=0.09,
                    zorder=0,
                )
                axes.add_patch(ellipse)
            elif args.figure_tip_labels == "yes":
                axes.annotate(
                    row["name"],
                    (x, y),
                    xytext=(4, 4),
                    textcoords="offset points",
                    fontsize=8,
                )
        tips = first.node_class == "leaf"
        axes.scatter(
            first.loc[tips, "mean"],
            second.loc[tips, "mean"],
            c="#D55E00",
            s=32,
            label="Observed tips",
            zorder=3,
        )
        axes.scatter(
            first.loc[~tips, "mean"],
            second.loc[~tips, "mean"],
            c="#0072B2",
            s=18,
            label="Conditional ancestors",
            zorder=2,
        )
        ratios = eigenvalues.set_index("component").explained_variance_ratio
        axes.set(
            xlabel=f"{components[0]} ({ratios[components[0]]:.1%})",
            ylabel=f"{components[1]} ({ratios[components[1]]:.1%})",
            title="Phylogenetic morphospace",
        )
        axes.set_aspect("equal", adjustable="datalim")
        axes.grid(alpha=0.15)
        axes.legend(loc="best", fontsize=8)
        axes.autoscale_view()
        traits = list(dict.fromkeys(loadings.trait))
        matrix = loadings.pivot(
            index="trait", columns="component", values="loading"
        ).loc[traits, components]
        plotted = heat.imshow(matrix, vmin=-1, vmax=1, cmap="RdBu_r", aspect="auto")
        heat.set(
            xticks=range(2),
            xticklabels=components,
            yticks=range(len(traits)),
            yticklabels=traits,
            title="Evolutionary loadings",
        )
        if len(traits) <= 15:
            for i in range(len(traits)):
                for j in range(2):
                    value = matrix.iloc[i, j]
                    heat.text(
                        j,
                        i,
                        f"{value:.2f}",
                        ha="center",
                        va="center",
                        color="white" if abs(value) > 0.6 else "black",
                        fontsize=9,
                    )
        figure.colorbar(plotted, ax=heat, shrink=0.65, label="Trait–PC correlation")
        figure.suptitle(f"Phylogenetic PCA · {args.model} · {args.mode}", fontsize=14)
        figure.supxlabel(
            f"Shading: conditional {args.ci_level:.0%} joint ellipses. Axes and evolutionary variances treated as fitted; lines join node means.",
            fontsize=8,
        )
        figure.savefig(path, format=Path(args.figure_out).suffix.lower()[1:], dpi=180)
        figure.clear()

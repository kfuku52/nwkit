"""Single-page PDF rendering for ASR model-comparison tables."""

import logging
import math

import numpy as np
import pandas as pd


def _criterion_label(criterion):
    return "AICc" if criterion == "aicc" else criterion.upper()


def _comparison_group_key(value):
    if value is None or bool(pd.isna(value)):
        return ""
    return str(value).strip()


def _comparison_set_label(group):
    labels = {
        "continuous": "Continuous",
        "discrete": "Discrete",
        "scalar": "single trait",
        "flat_root_integrated": "flat-root integrated likelihood",
        "stationary_root_ml": "stationary-root ML",
        "proper_root_ml": "proper-root ML",
        "discrete_ml": "CTMC ML",
        "equal": "equal root",
        "empirical": "empirical root",
        "stationary": "stationary root",
        "flat": "flat root",
        "fixed": "fixed root",
        "gaussian": "Gaussian root",
    }
    parts = []
    for token in str(group).split(":"):
        if token.endswith("char") and token[:-4].isdigit():
            parts.append(f"{token[:-4]} characters")
        elif token.endswith("d") and token[:-1].isdigit():
            parts.append(f"{token[:-1]} traits")
        else:
            parts.append(labels.get(token, token.replace("_", " ")))
    return " | ".join(parts)


def _status_summary(rows):
    counts = rows["status"].astype(str).value_counts(sort=False)
    return ", ".join(
        f"{status.replace('_', ' ')}: {count}" for status, count in counts.items()
    )


def _comparison_sections(table, criterion_label):
    working = table.copy()
    working["_source_order"] = np.arange(len(working))
    working["_plot_group"] = [
        _comparison_group_key(value) for value in working["comparison_group"]
    ]
    sections = []
    group_order = dict.fromkeys(
        group for group in working["_plot_group"].tolist() if group
    )
    for number, group in enumerate(group_order, start=1):
        rows = working[working["_plot_group"] == group].copy()
        rows["_rank_sort"] = pd.to_numeric(rows["criterion_rank"], errors="coerce")
        rows = rows.sort_values(
            ["_rank_sort", "_source_order"], kind="stable", na_position="last"
        )
        ranked_count = int(rows["criterion_rank"].notna().sum())
        unranked_count = len(rows) - ranked_count
        summary = f"{ranked_count} ranked by {criterion_label}"
        if unranked_count:
            summary += f", {unranked_count} not ranked"
        sections.append(
            (
                f"Comparison set {number}: {_comparison_set_label(group)} ({summary})",
                rows,
                False,
            )
        )
    unassigned = working[working["_plot_group"] == ""].sort_values(
        "_source_order", kind="stable"
    )
    if not unassigned.empty:
        sections.append(
            (
                "Not assigned to a comparison set "
                f"({len(unassigned)} candidates; {_status_summary(unassigned)})",
                unassigned,
                True,
            )
        )
    return sections


def _finite_number(value):
    try:
        number = float(value)
    except (TypeError, ValueError, OverflowError):
        return None
    return number if math.isfinite(number) else None


def _format_number(value):
    number = _finite_number(value)
    if number is None:
        return ""
    magnitude = abs(number)
    if magnitude >= 1e6 or (0.0 < magnitude < 1e-3):
        return f"{number:.3e}"
    return f"{number:.3f}"


def _format_count(value):
    number = _finite_number(value)
    return "" if number is None else str(int(number))


def _safe_text(value):
    if value is None:
        return ""
    try:
        if bool(pd.isna(value)):
            return ""
    except (TypeError, ValueError):
        pass
    return " ".join(str(value).split())


def _shorten(value, width):
    text = _safe_text(value)
    if len(text) <= width:
        return text
    if width <= 3:
        return "." * width
    return text[: width - 3].rstrip() + "..."


def _font_family_for_text(text):
    required = {ord(character) for character in str(text) if character.isprintable()}
    if not any(codepoint > 127 for codepoint in required):
        return None

    from matplotlib import font_manager, ft2font

    preferred = (
        "Noto Sans CJK JP",
        "Noto Sans JP",
        "Source Han Sans JP",
        "Arial Unicode MS",
        "Hiragino Sans",
        "Yu Gothic",
        "Meiryo",
        "Noto Sans CJK SC",
        "Noto Sans CJK TC",
    )
    priority = {name.casefold(): index for index, name in enumerate(preferred)}
    entries = sorted(
        font_manager.fontManager.ttflist,
        key=lambda entry: (
            priority.get(entry.name.casefold(), len(priority)),
            entry.name.casefold(),
            entry.fname,
        ),
    )
    seen = set()
    for entry in entries:
        if entry.fname in seen:
            continue
        seen.add(entry.fname)
        try:
            supported = set(ft2font.FT2Font(entry.fname).get_charmap())
        except (OSError, RuntimeError, ValueError):
            continue
        if required <= supported:
            return entry.name
    missing = ", ".join(f"U+{codepoint:04X}" for codepoint in sorted(required)[:12])
    raise ValueError(
        "The ASR comparison PDF contains characters for which no installed font "
        f"has complete coverage ({missing}). Install a Unicode/CJK font such as "
        "Noto Sans CJK and retry."
    )


def _figure_text(table, criterion):
    columns = (
        "model_id",
        "trait_type",
        "trait_columns",
        "status",
        "message",
        "comparison_group",
    )
    values = [
        str(value)
        for column in columns
        if column in table
        for value in table[column].dropna().tolist()
    ]
    values.extend(
        (
            "ASR model comparison",
            _criterion_label(criterion),
            "Comparison set Not assigned Rank Model Fit status Notes",
        )
    )
    return "\n".join(values)


def _fit_artist_to_width(artist, renderer, maximum_width):
    text = artist.get_text()
    if not text or artist.get_window_extent(renderer=renderer).width <= maximum_width:
        return
    suffix = "..."
    low, high = 0, len(text)
    while low < high:
        middle = (low + high + 1) // 2
        artist.set_text(text[:middle].rstrip() + suffix)
        if artist.get_window_extent(renderer=renderer).width <= maximum_width:
            low = middle
        else:
            high = middle - 1
    artist.set_text(text[:low].rstrip() + suffix if low else suffix)


def _fit_figure_text(figure, axis, title, bounded_axis_text):
    from matplotlib.backends.backend_agg import FigureCanvasAgg

    canvas = FigureCanvasAgg(figure)
    canvas.draw()
    renderer = canvas.get_renderer()
    _fit_artist_to_width(title, renderer, figure.bbox.width * 0.94)
    for artist, right_edge in bounded_axis_text:
        left = axis.transData.transform((artist.get_position()[0], 0.0))[0]
        right = axis.transData.transform((right_edge, 0.0))[0]
        _fit_artist_to_width(artist, renderer, max(1.0, right - left))


def _row_note(record, criterion, criterion_label):
    notes = []
    message = _safe_text(record.message)
    if message:
        notes.append(message)
    rank = _finite_number(record.criterion_rank)
    if rank is None and str(record.rankable).lower() in {"yes", "true", "1"}:
        criterion_value = _finite_number(getattr(record, criterion))
        comparable = _finite_number(record.num_comparable_models)
        if criterion_value is None:
            notes.append(f"{criterion_label} is unavailable for this model.")
        elif comparable is None or comparable < 2:
            notes.append(f"Fewer than two finite {criterion_label} values in this set.")
    return _shorten(" ".join(notes), 76)


def _draw_candidate_overview(pdf, table, criterion):
    from matplotlib.figure import Figure
    from matplotlib.patches import Rectangle

    criterion_label = _criterion_label(criterion)
    sections = _comparison_sections(table, criterion_label)
    display_row_count = sum(len(rows) + 1 for _label, rows, _plain in sections)
    height = max(4.8, 0.42 * display_row_count + 2.4)
    figure = Figure(figsize=(15.0, height))
    axis = figure.subplots()
    axis.set_xlim(0.0, 1.0)
    axis.set_ylim(0.0, display_row_count + 1.4)
    axis.axis("off")

    trait_types = ", ".join(dict.fromkeys(table["trait_type"].astype(str).tolist()))
    trait_columns = ", ".join(
        dict.fromkeys(table["trait_columns"].astype(str).tolist())
    )
    title = _shorten(
        f"ASR model comparison - {trait_types} ({trait_columns})",
        110,
    )
    title_artist = figure.suptitle(title, fontsize=15, fontweight="bold", y=0.98)

    status_summary = _shorten(_status_summary(table), 180)
    bounded_axis_text = []
    status_artist = axis.text(
        0.0,
        display_row_count + 1.05,
        status_summary,
        fontsize=8.5,
        color="#444444",
        va="bottom",
    )
    bounded_axis_text.append((status_artist, 0.985))
    headers = (
        (0.012, "Rank", "left"),
        (0.055, "Model", "left"),
        (0.17, "Fit status", "left"),
        (0.335, "n", "right"),
        (0.375, "k", "right"),
        (0.46, "logL", "right"),
        (0.54, criterion_label, "right"),
        (0.62, f"Delta {criterion_label}", "right"),
        (0.70, "Weight", "right"),
        (0.72, "Notes", "left"),
    )
    for x_position, label, alignment in headers:
        axis.text(
            x_position,
            display_row_count + 0.55,
            label,
            fontsize=9,
            fontweight="bold",
            ha=alignment,
            va="bottom",
        )

    status_colors = {
        "ok": "#009E73",
        "boundary": "#E69F00",
        "equivalent": "#4C78A8",
        "not_applicable": "#8A8A8A",
        "not_fitted": "#777777",
        "failed": "#D55E00",
        "nonconverged": "#D55E00",
        "nonregular": "#CC79A7",
        "no_likelihood": "#7B6FD0",
    }
    row_index = 0
    data_row_index = 0
    section_colors = ("#DCEAF7", "#E2F0E9", "#F1E5F0", "#F3EBDC")
    for section_index, (label, rows, plain) in enumerate(sections):
        section_y = display_row_count - row_index - 0.35
        axis.add_patch(
            Rectangle(
                (0.0, section_y - 0.37),
                1.0,
                0.74,
                facecolor=(
                    "#E5E5E5"
                    if plain
                    else section_colors[section_index % len(section_colors)]
                ),
                edgecolor="none",
                zorder=0,
            )
        )
        section_artist = axis.text(
            0.015,
            section_y,
            _shorten(label, 175),
            fontsize=8.5,
            fontweight="bold",
            color="#222222",
            va="center",
        )
        bounded_axis_text.append((section_artist, 0.985))
        row_index += 1
        for record in rows.itertuples(index=False):
            y_position = display_row_count - row_index - 0.35
            if data_row_index % 2:
                axis.add_patch(
                    Rectangle(
                        (0.0, y_position - 0.37),
                        1.0,
                        0.74,
                        facecolor="#F7F8F9",
                        edgecolor="none",
                        zorder=0,
                    )
                )
            rank = _finite_number(record.criterion_rank)
            is_best = str(record.is_best).lower() == "yes"
            axis.text(
                0.012,
                y_position,
                "-" if rank is None else f"#{int(rank)}",
                fontsize=7.5,
                fontweight="bold" if is_best else "normal",
                va="center",
            )
            model_artist = axis.text(
                0.055,
                y_position,
                _shorten(record.model_id, 20),
                fontsize=8.0,
                fontweight="bold" if is_best else "normal",
                va="center",
            )
            bounded_axis_text.append((model_artist, 0.165))
            status = str(record.status)
            status_color = status_colors.get(status, "#666666")
            status_text_color = "#111111" if status == "boundary" else "white"
            fit_status_artist = axis.text(
                0.17,
                y_position,
                _shorten(status.replace("_", " "), 17),
                fontsize=7.0,
                color=status_text_color,
                va="center",
                bbox={
                    "boxstyle": "round,pad=0.22",
                    "facecolor": status_color,
                    "edgecolor": status_color,
                },
            )
            bounded_axis_text.append((fit_status_artist, 0.325))
            numeric_values = (
                (0.335, _format_count(record.sample_size)),
                (0.375, _format_count(record.num_parameters)),
                (0.46, _format_number(record.log_likelihood)),
                (0.54, _format_number(getattr(record, criterion))),
                (0.62, _format_number(getattr(record, f"delta_{criterion}"))),
                (0.70, _format_number(getattr(record, f"{criterion}_weight"))),
            )
            for x_position, value in numeric_values:
                axis.text(
                    x_position,
                    y_position,
                    value,
                    fontsize=7.2,
                    ha="right",
                    va="center",
                )
            note_artist = axis.text(
                0.72,
                y_position,
                _row_note(record, criterion, criterion_label),
                fontsize=6.6,
                color="#555555",
                va="center",
            )
            bounded_axis_text.append((note_artist, 0.985))
            row_index += 1
            data_row_index += 1

    figure.text(
        0.5,
        0.018,
        "n = sample size; k = estimated parameter count. Ranks, deltas, and weights "
        "are calculated independently within each shaded comparison set.",
        ha="center",
        fontsize=8,
        color="#444444",
    )
    figure.subplots_adjust(left=0.035, right=0.985, top=0.88, bottom=0.065)
    _fit_figure_text(figure, axis, title_artist, bounded_axis_text)
    pdf.savefig(figure)


def draw_comparison_figure(table, path, criterion):
    import matplotlib as mpl
    from matplotlib.backends.backend_pdf import PdfPages

    font_loggers = [
        logging.getLogger(name) for name in ("fontTools", "fontTools.subset")
    ]
    previous_levels = [logger.level for logger in font_loggers]
    try:
        for logger in font_loggers:
            logger.setLevel(logging.WARNING)
        font_family = _font_family_for_text(_figure_text(table, criterion))
        settings = {"pdf.fonttype": 42, "ps.fonttype": 42}
        if font_family is not None:
            settings["font.family"] = font_family
        with mpl.rc_context(settings):
            with PdfPages(
                path,
                metadata={
                    "Title": "nwkit ASR model comparison",
                    "Subject": f"Grouped {_criterion_label(criterion)} comparison",
                    "Creator": "nwkit",
                    "Keywords": "ancestral state reconstruction model comparison",
                },
            ) as pdf:
                _draw_candidate_overview(pdf, table, criterion)
    finally:
        for logger, level in zip(font_loggers, previous_levels, strict=True):
            logger.setLevel(level)

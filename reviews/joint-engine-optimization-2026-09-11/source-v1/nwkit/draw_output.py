"""Drawing serialization and recoverable image/report output."""

import matplotlib.pyplot as plt

from nwkit import __version__
from nwkit.draw_quality import write_layout_report
from nwkit.output_transaction import output_transaction


def save_drawing(
    figure,
    outfile,
    image_format,
    transparent,
    layout_report,
    quality_report,
    *,
    transactional=True,
):
    metadata: dict[str, str | None] = {"Creator": f"NWKIT {__version__}"}
    if image_format == "svg":
        metadata["Date"] = None
    elif image_format == "pdf":
        metadata["CreationDate"] = None
        metadata["ModDate"] = None

    def write_image(path):
        figure.savefig(
            path,
            format=image_format,
            dpi=300,
            transparent=bool(transparent),
            metadata=metadata,
        )

    try:
        if not transactional:
            # A containing transaction (rootcompare) owns these staging files.
            write_image(outfile)
            write_layout_report(layout_report, quality_report)
            return
        paths = [outfile]
        if layout_report not in (None, "", "-"):
            paths.append(layout_report)
        write_stdout = (
            (lambda: write_layout_report("-", quality_report))
            if layout_report == "-"
            else None
        )
        with output_transaction(
            paths, create_parents=True, after_install=write_stdout
        ) as staged:
            write_image(staged[outfile])
            if layout_report not in (None, "", "-"):
                write_layout_report(staged[layout_report], quality_report)
    finally:
        plt.close(figure)

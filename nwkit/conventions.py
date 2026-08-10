from typing import Any

DEFAULT_TABLE_MISSING_VALUES = ("", "NA", "NaN", "nan", "?", "missing", "unknown")
DEFAULT_TABLE_MISSING_VALUES_CSV = ",".join(DEFAULT_TABLE_MISSING_VALUES)

STDIN_INPUT_DESTS = (
    "infile",
    "infile2",
    "seqin",
    "trait",
    "table",
    "taxid_tsv",
    "weight_tsv",
    "species_map_tsv",
    "species_name_tsv",
    "name_tsv",
    "species_list",
    "reference",
    "root_source",
    "name_source",
    "support_source",
    "length_source",
    "tip_image_manifest",
)


def get_stdin_input_options(args: Any) -> list[tuple[str, str]]:
    options = [
        (dest, "--{}".format(dest.replace("_", "-")))
        for dest in STDIN_INPUT_DESTS
        if getattr(args, dest, None) == "-"
    ]
    options.extend(
        ("property_source", "--property-source")
        for value in (getattr(args, "property_source", None) or [])
        if str(value).rsplit("@", 1)[-1] == "-"
    )
    return options

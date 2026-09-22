import io

import pytest

from nwkit.image import read_name_tsv as read_image_names
from nwkit.rename import read_name_tsv as read_rename_names
from nwkit.species_parser import ParsedSpecies, read_species_map_tsv

READERS = [
    (
        read_species_map_tsv,
        "leaf_name",
        "species_label",
        ParsedSpecies("Mapped_species"),
    ),
    (read_image_names, "leaf_name", "species_name", "Mapped species"),
    (read_rename_names, "old_name", "new_name", "Mapped_species"),
]


@pytest.mark.parametrize("reader,key,value,expected", READERS)
@pytest.mark.parametrize("source", ["file", "stdin"])
@pytest.mark.parametrize("malformation", ["duplicate", "empty", "extra", "short"])
def test_mapping_tables_reject_ambiguous_structure(
    reader, key, value, expected, source, malformation, tmp_path, monkeypatch
):
    header = [key, value, "note"]
    row = ["A", "Mapped_species", "annotation"]
    if malformation == "duplicate":
        header[2] = value
    elif malformation == "empty":
        header[2] = " "
    elif malformation == "extra":
        row.append("unexpected")
    else:
        row.pop()
    text = "\t".join(header) + "\n" + "\t".join(row) + "\n"
    stream = io.StringIO(text)
    if source == "stdin":
        monkeypatch.setattr("sys.stdin", stream)
        path = "-"
    else:
        path = tmp_path / "mapping.tsv"
        path.write_text(text, encoding="utf-8")
    with pytest.raises(ValueError, match="headers|number of fields"):
        reader(path)
    assert not stream.closed


@pytest.mark.parametrize("reader,key,value,expected", READERS)
@pytest.mark.parametrize("encoding", ["utf-8", "utf-8-sig"])
def test_mapping_files_accept_utf8_and_preserve_literal_keys(
    reader, key, value, expected, encoding, tmp_path
):
    path = tmp_path / "mapping.tsv"
    labels = ["001", "NA", " sample ", "caf\u00e9"]
    path.write_text(
        f"{key}\t{value}\tnote\n"
        + "".join(f"{label}\tMapped_species\t\n" for label in labels),
        encoding=encoding,
    )
    assert reader(path) == dict.fromkeys(labels, expected)

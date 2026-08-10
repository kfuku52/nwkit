from io import StringIO

from hypothesis import given
from hypothesis import strategies as st

from nwkit.fasta import FastaRecord, parse_fasta, write_fasta
from nwkit.image_metadata import sanitize_filename_component
from nwkit.util import iter_newick_stream, split_newick_stream

FASTA_IDENTIFIERS = st.text(
    alphabet="ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789_.-",
    min_size=1,
    max_size=24,
)
SEQUENCE_LINES = st.lists(
    st.text(alphabet="ACGTNacgtn-?", min_size=1, max_size=40),
    min_size=1,
    max_size=6,
)


@given(st.lists(st.tuples(FASTA_IDENTIFIERS, SEQUENCE_LINES), min_size=1, max_size=12))
def test_fasta_records_round_trip_without_information_loss(record_values):
    records = [
        FastaRecord(
            name=name,
            raw=">{} generated description\n{}\n".format(name, "\n".join(lines)),
        )
        for name, lines in record_values
    ]
    output = StringIO()

    assert write_fasta(records, output) == len(records)
    assert parse_fasta(StringIO(output.getvalue())) == records


QUOTED_NEWICK_LABELS = st.text(
    alphabet="abcXYZ 0123;,:()[]'\"",
    min_size=1,
    max_size=30,
)


@given(
    st.lists(QUOTED_NEWICK_LABELS, min_size=1, max_size=10),
    st.integers(min_value=1, max_value=16),
)
def test_newick_collection_parsers_agree_at_every_chunk_boundary(labels, chunk_size):
    trees = [
        "('{}':1,B{}:2);".format(label.replace("'", "''"), index)
        for index, label in enumerate(labels)
    ]
    text = "\n".join(trees) + "\n"

    assert split_newick_stream(text) == trees
    assert list(iter_newick_stream(StringIO(text), chunk_size=chunk_size)) == trees


@given(st.text(min_size=0, max_size=500))
def test_sanitized_filename_components_are_bounded_and_portable(value):
    component = sanitize_filename_component(value)

    assert component
    assert len(component) <= 64
    assert component == component.strip(".")
    assert all(character.isalnum() or character in "._-" for character in component)

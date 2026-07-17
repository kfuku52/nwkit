from io import BytesIO, StringIO, TextIOWrapper

import pytest

from nwkit.fasta import FastaRecord, parse_fasta, write_fasta
from nwkit.util import read_seqs, write_seqs


def test_parse_fasta_preserves_description_wrapping_and_final_newline():
    text = '>A description with spaces\nACGT\nTGCA\n>B\nNN--\n'
    records = parse_fasta(StringIO(text))
    assert [record.name for record in records] == ['A', 'B']
    assert records[0].raw == '>A description with spaces\nACGT\nTGCA\n'
    output = StringIO()
    assert write_fasta(records, output) == 2
    assert output.getvalue() == text


def test_write_fasta_separates_records_after_missing_final_newline():
    records = [
        FastaRecord(name='A', raw='>A\nACGT'),
        FastaRecord(name='B', raw='>B\nTGCA\n'),
    ]
    output = StringIO()
    write_fasta(records, output)
    assert output.getvalue() == '>A\nACGT\n>B\nTGCA\n'


def test_write_fasta_normalizes_newlines_before_platform_translation():
    raw_output = BytesIO()
    output = TextIOWrapper(raw_output, encoding='utf-8', newline='\r\n')
    records = [FastaRecord(name='A', raw='>A\r\nACGT\r\n')]
    write_fasta(records, output, normalize_newlines=True)
    output.flush()
    assert raw_output.getvalue() == b'>A\r\nACGT\r\n'


@pytest.mark.parametrize('text', ['>\nACGT\n', '>   \nACGT\n'])
def test_parse_fasta_rejects_empty_identifier(text):
    with pytest.raises(ValueError, match='does not contain a sequence identifier'):
        parse_fasta(StringIO(text))


def test_parse_fasta_rejects_content_before_first_header():
    with pytest.raises(ValueError, match="before the first '>'"):
        parse_fasta(StringIO('ACGT\n>A\nACGT\n'))


def test_parse_fasta_ignores_comment_preamble():
    text = '; legacy comment\n  # metadata\n\n>A description\nACGT\n'
    records = parse_fasta(StringIO(text))
    assert records == [FastaRecord(name='A', raw='>A description\nACGT\n')]


def test_read_and_write_seqs_reject_non_fasta_format(tmp_path):
    path = tmp_path / 'input.fasta'
    path.write_text('>A\nACGT\n')
    with pytest.raises(ValueError, match="Only 'fasta' is supported"):
        read_seqs(str(path), 'phylip', quiet=True)
    with pytest.raises(ValueError, match="Only 'fasta' is supported"):
        write_seqs([], str(tmp_path / 'output.phy'), 'phylip', quiet=True)

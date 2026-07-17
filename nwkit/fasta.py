from dataclasses import dataclass


@dataclass(frozen=True)
class FastaRecord:
    name: str
    raw: str


def parse_fasta(handle):
    records = list()
    current_name = None
    current_lines = list()
    for line_number, line in enumerate(handle, start=1):
        if line.startswith('>'):
            if current_name is not None:
                records.append(FastaRecord(name=current_name, raw=''.join(current_lines)))
            header = line[1:].strip()
            if not header:
                raise ValueError(
                    'FASTA header on line {} does not contain a sequence identifier.'.format(
                        line_number
                    )
                )
            current_name = header.split(maxsplit=1)[0]
            current_lines = [line]
        elif current_name is None:
            stripped = line.lstrip()
            if stripped.startswith((';', '#')):
                continue
            if stripped:
                raise ValueError(
                    "FASTA content before the first '>' header on line {}.".format(
                        line_number
                    )
                )
        else:
            current_lines.append(line)
    if current_name is not None:
        records.append(FastaRecord(name=current_name, raw=''.join(current_lines)))
    return records


def write_fasta(records, handle, normalize_newlines=False):
    previous_ended_with_newline = True
    count = 0
    for record in records:
        raw = record.raw
        if normalize_newlines:
            raw = raw.replace('\r\n', '\n').replace('\r', '\n')
        if not previous_ended_with_newline:
            handle.write('\n')
        handle.write(raw)
        previous_ended_with_newline = raw.endswith(('\n', '\r'))
        count += 1
    return count

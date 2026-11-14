"""Dependency-free BIFF cell-value reader for AngioCal workbooks.

Uses MS-XLS record layouts and the OpenOffice Excel file-format reference.
Formatting is read only to distinguish dates from numbers. Formula expressions,
macros, drawings and external references are never evaluated.
"""

import struct
from dataclasses import dataclass, field

from nwkit.xls_compound import require, uint16, uint32, workbook_stream
from nwkit.xls_strings import (
    Segments,
    codepage_encoding,
    is_date_format,
    legacy_string,
    shared_strings,
)

BOF_CODES = {0x0009, 0x0209, 0x0409, 0x0809}
CELL_CODES = {
    0x0002,
    0x0003,
    0x0004,
    0x0005,
    0x0006,
    0x0203,
    0x0204,
    0x0205,
    0x0206,
    0x0406,
    0x027E,
    0x00BD,
    0x00FD,
    0x00D6,
}


@dataclass(frozen=True)
class Record:
    code: int
    data: bytes
    offset: int


@dataclass(frozen=True)
class Cell:
    value: str | float | int
    kind: str


EMPTY_CELL = Cell("", "empty")


@dataclass
class Sheet:
    cells: dict[tuple[int, int], Cell] = field(default_factory=dict)
    nrows: int = 0
    ncols: int = 0
    max_rows: int = 65536

    def put(self, row, column, cell):
        require(
            0 <= row < self.max_rows and 0 <= column < 256, "cell is out of bounds."
        )
        require((row, column) not in self.cells, "duplicate cell records.")
        self.cells[row, column] = cell
        self.nrows = max(self.nrows, row + 1)
        self.ncols = max(self.ncols, column + 1)

    def cell(self, row, column):
        return self.cells.get((row, column), EMPTY_CELL)

    def row_values(self, row):
        return [self.cell(row, column).value for column in range(self.ncols)]


@dataclass
class Context:
    version: int
    encoding: str
    strings: list[str] = field(default_factory=list)
    xfs: list[int] = field(default_factory=list)
    formats: dict[int, str] = field(default_factory=dict)
    ixfe: int | None = None
    date_kinds: dict[int, str] = field(default_factory=dict)

    def format_id(self, data):
        if self.version == 2:
            if not self.xfs:
                return data[5] & 0x3F
            index = data[4] & 0x3F
            if index == 0x3F:
                require(self.ixfe is not None, "missing extended BIFF2 format index.")
                index = self.ixfe
        else:
            index = uint16(data, 4)
        if not self.xfs:
            require(index == 0, "missing cell format.")
            return 0
        require(
            index is not None and 0 <= index < len(self.xfs),
            "cell format index is out of bounds.",
        )
        return self.xfs[index]

    def number(self, value, data):
        format_id = self.format_id(data)
        kind = self.date_kinds.get(format_id)
        if kind is None:
            kind = "date" if is_date_format(format_id, self.formats) else "number"
            self.date_kinds[format_id] = kind
        return Cell(value, kind)


def records(data, start=0, stop=None):
    result: list[Record] = []
    position = start
    depth = 0
    substreams: list[int] = []
    stop = len(data) if stop is None else stop
    require(0 <= start < stop <= len(data), "invalid BIFF substream bounds.")
    while position < stop:
        require(position + 4 <= stop, "truncated BIFF record header.")
        code, length = struct.unpack_from("<HH", data, position)
        end = position + 4 + length
        require(end <= stop, "truncated BIFF record payload.")
        payload = data[position + 4 : end]
        if not result:
            require(code in BOF_CODES, "missing BIFF BOF record.")
        if code == 0x002F:
            raise ValueError("Encrypted AngioCal XLS workbooks are not supported.")
        if code in BOF_CODES:
            require(len(payload) >= 4, "truncated BIFF BOF.")
            kind = uint16(payload, 2)
            if depth:
                require(
                    (substreams[-1] == 0x0100 and depth == 1)
                    or (substreams[-1] in {0x0010, 0x0020, 0x0040} and kind == 0x0020),
                    "invalid nested BIFF substream.",
                )
            substreams.append(kind)
            depth += 1
        if depth == 1:
            result.append(Record(code, payload, position))
        if code == 0x000A:
            require(length == 0, "invalid BIFF EOF.")
            depth -= 1
            substreams.pop()
            if depth == 0:
                return result
        position = end
    raise ValueError("Invalid AngioCal XLS workbook: missing BIFF EOF record.")


def biff_version(record):
    if record.code == 0x0809:
        version = uint16(record.data)
        versions = {
            0x0600: 8,
            0x0500: 5,
            0x0400: 4,
            0x0300: 3,
            0x0200: 2,
            0x0000: 2,
            0x0007: 2,
        }
        require(version in versions, "unsupported BIFF version.")
        return versions[version]
    return {0x0009: 2, 0x0209: 3, 0x0409: 4}[record.code]


def continued_parts(items, index):
    parts = [items[index].data]
    index += 1
    while index < len(items) and items[index].code == 0x003C:
        parts.append(items[index].data)
        index += 1
    return parts, index


def _format_record(context, record, parts):
    if context.version == 8:
        reader = Segments(parts)
        key = uint16(reader.read(2))
        value = reader.unicode()
    else:
        key = uint16(record.data) if context.version >= 5 else len(context.formats)
        offset = 2 if context.version >= 4 and record.code != 0x001E else 0
        require(offset < len(record.data), "truncated number format.")
        value = legacy_string(record.data, offset, context.encoding, 1)
    require(key not in context.formats, "duplicate number format.")
    context.formats[key] = value


def _xf_record(context, data):
    if context.version >= 5:
        value = uint16(data, 2)
    else:
        offset = 2 if context.version == 2 else 1
        require(len(data) > offset, "truncated cell format.")
        value = data[offset] & (0x3F if context.version == 2 else 0xFF)
    context.xfs.append(value)


def make_context(items, version):
    codepages = [uint16(item.data) for item in items if item.code == 0x0042]
    encoding = (
        "latin1"
        if version == 8
        else codepage_encoding(codepages[-1] if codepages else None)
    )
    context = Context(version, encoding)
    index = 0
    seen_sst = False
    while index < len(items):
        item = items[index]
        parts, next_index = continued_parts(items, index)
        if item.code == 0x00FC:
            require(version == 8 and not seen_sst, "invalid shared-string table.")
            context.strings = shared_strings(parts)
            seen_sst = True
        elif item.code in {0x00E0, 0x0043, 0x0243, 0x0443}:
            _xf_record(context, item.data)
        elif item.code in {0x041E, 0x001E}:
            _format_record(context, item, parts)
        index = next_index
    return context


def rk_number(value):
    if value & 2:
        number = float(struct.unpack("<i", struct.pack("<I", value))[0] >> 2)
    else:
        number = struct.unpack("<d", struct.pack("<II", 0, value & 0xFFFFFFFC))[0]
    return number / 100 if value & 1 else number


def _boolean_error(data, offset):
    require(offset + 2 <= len(data), "truncated boolean/error cell.")
    value, error = data[offset : offset + 2]
    require(error in {0, 1} and (error or value in {0, 1}), "invalid boolean cell.")
    return Cell(value, "error" if error else "boolean")


def _formula(data, offset, context):
    require(offset + 8 <= len(data), "truncated formula result.")
    cached = data[offset : offset + 8]
    if cached[6:8] != b"\xff\xff":
        return context.number(struct.unpack("<d", cached)[0], data)
    kind = cached[0]
    require(kind in {0, 1, 2, 3}, "invalid formula result type.")
    if kind == 0:
        return None
    if kind == 3:
        return Cell("", "text")
    return _boolean_error(bytes([cached[2], int(kind == 2)]), 0)


def _text_cell(item, parts, context, offset):
    if item.code == 0x00FD:
        index = uint32(item.data, offset)
        require(index < len(context.strings), "shared-string index is out of bounds.")
        value = context.strings[index]
    elif context.version == 8:
        value = Segments([parts[0][offset:]] + parts[1:]).unicode()
    else:
        require(offset < len(item.data), "truncated text cell.")
        value = legacy_string(
            b"".join(parts), offset, context.encoding, 1 if context.version == 2 else 2
        )
    return Cell(value, "text")


def _single_cell(item, parts, context):
    data = item.data
    offset = 7 if context.version == 2 else 6
    require(len(data) >= offset, "truncated cell header.")
    code = item.code
    if code in {0x0004, 0x0204, 0x00D6, 0x00FD}:
        return _text_cell(item, parts, context, offset)
    if code in {0x0005, 0x0205}:
        return _boolean_error(data, offset)
    if code in {0x0006, 0x0206, 0x0406}:
        return _formula(data, offset, context)
    if code == 0x027E:
        return context.number(rk_number(uint32(data, offset)), data)
    if code == 0x0002:
        return context.number(float(uint16(data, offset)), data)
    require(offset + 8 <= len(data), "truncated numeric cell.")
    return context.number(struct.unpack_from("<d", data, offset)[0], data)


def _multiple_rk(sheet, data, context):
    require(
        len(data) >= 12 and (len(data) - 6) % 6 == 0, "invalid MULRK record length."
    )
    row, first, last = uint16(data), uint16(data, 2), uint16(data, len(data) - 2)
    count = (len(data) - 6) // 6
    require(last < 256 and last - first + 1 == count, "invalid MULRK columns.")
    for index in range(count):
        offset = 4 + index * 6
        prefix = struct.pack("<HH", row, first + index) + data[offset : offset + 2]
        sheet.put(
            row,
            first + index,
            context.number(rk_number(uint32(data, offset + 2)), prefix),
        )


def read_sheet(items, context):
    sheet = Sheet(max_rows=65536 if context.version == 8 else 16384)
    pending: tuple[int, int] | None = None
    index = 1
    while index < len(items):
        item = items[index]
        parts, next_index = continued_parts(items, index)
        if item.code in {0x0007, 0x0207}:
            if pending is None:
                raise ValueError(
                    "Invalid AngioCal XLS workbook: string result has no preceding formula."
                )
            value = (
                Segments(parts).unicode()
                if context.version == 8
                else legacy_string(
                    b"".join(parts),
                    0,
                    context.encoding,
                    1 if context.version == 2 else 2,
                )
            )
            sheet.put(*pending, Cell(value, "text"))
            pending = None
        elif item.code == 0x0044:
            context.ixfe = uint16(item.data)
        elif item.code in CELL_CODES:
            require(pending is None, "missing formula string result.")
            if item.code == 0x00BD:
                _multiple_rk(sheet, item.data, context)
            else:
                cell = _single_cell(item, parts, context)
                position = uint16(item.data), uint16(item.data, 2)
                if cell is None:
                    pending = position
                else:
                    sheet.put(*position, cell)
        index = next_index
    require(pending is None, "missing formula string result.")
    return sheet


def _legacy_workbook(stream, items):
    sheets = []
    headers = [item for item in items if item.code == 0x008F]
    offsets = [uint32(item.data) for item in items if item.code == 0x008E]
    require(
        len(offsets) <= 1
        and (not offsets or headers and offsets[0] == headers[0].offset),
        "invalid BIFF4 sheet offset.",
    )
    codepages = [item for item in items if item.code == 0x0042]
    encoding = make_context(codepages, 4).encoding
    names = [
        legacy_string(item.data, 0, encoding, 1)
        for item in items
        if item.code == 0x0085
    ]
    require(
        not names or len(names) == len(headers), "BIFF4 sheet count does not match."
    )
    for index, item in enumerate(headers):
        require(len(item.data) >= 5, "truncated BIFF4 sheet header.")
        name = legacy_string(item.data, 4, encoding, 1)
        require(not names or name == names[index], "BIFF4 sheet name does not match.")
        start = item.offset + 4 + len(item.data)
        size = uint32(item.data)
        require(0 < size <= len(stream) - start, "invalid BIFF4 sheet size.")
        sheet_records = records(stream, start)
        require(
            sheet_records[-1].offset + 4 == start + size
            and biff_version(sheet_records[0]) == 4,
            "BIFF4 sheet size or version does not match.",
        )
        if uint16(sheet_records[0].data, 2) == 0x0010:
            sheets.append(
                read_sheet(sheet_records, make_context(codepages + sheet_records, 4))
            )
    return tuple(sheets)


def read_xls(data):
    stream = workbook_stream(data)
    global_records = records(stream)
    version = biff_version(global_records[0])
    kind = uint16(global_records[0].data, 2)
    if kind == 0x0100 and version == 4:
        return _legacy_workbook(stream, global_records)
    require(kind in {0x0005, 0x0010}, "unsupported BIFF substream type.")
    context = make_context(global_records, version)
    if kind == 0x0010:
        return (read_sheet(global_records, context),)
    directory = {}
    global_end = global_records[-1].offset + 4
    for item in global_records:
        if item.code != 0x0085:
            continue
        require(len(item.data) >= 6, "truncated worksheet directory.")
        offset = uint32(item.data)
        require(
            global_end <= offset < len(stream) and offset not in directory,
            "invalid worksheet offset.",
        )
        directory[offset] = item.data[5]
    offsets = sorted(directory)
    sheets = {}
    for offset, stop in zip(offsets, offsets[1:] + [len(stream)], strict=True):
        if directory[offset] != 0:
            continue
        sheet_records = records(stream, offset, stop)
        require(
            biff_version(sheet_records[0]) == version
            and uint16(sheet_records[0].data, 2) == 0x0010,
            "invalid worksheet BOF.",
        )
        sheets[offset] = read_sheet(sheet_records, context)
    return tuple(sheets[offset] for offset in directory if offset in sheets)

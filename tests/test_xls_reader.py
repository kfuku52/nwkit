"""Small, independently constructed BIFF/CFB fixtures and malformed-input cases."""

import struct
import subprocess
import sys

import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

from nwkit.xls_compound import END, FAT, FREE, SIGNATURE, workbook_stream
from nwkit.xls_reader import read_xls
from nwkit.xls_strings import is_date_format


def record(code, data=b""):
    return struct.pack("<HH", code, len(data)) + data


def bof(kind=5, version=8):
    code = {2: 0x0009, 3: 0x0209, 4: 0x0409, 5: 0x0809, 8: 0x0809}[version]
    return record(
        code,
        struct.pack("<HHHHII", 0x0600 if version == 8 else 0x0500, kind, 0, 0, 0, 0),
    )


def unicode_text(value):
    data = value.encode("utf-16le")
    return struct.pack("<HB", len(data) // 2, 1) + data


def prefix(row=0, column=0, xf=0):
    return struct.pack("<HHH", row, column, xf)


def number(value, row=0, column=0, xf=0):
    return record(0x0203, prefix(row, column, xf) + struct.pack("<d", value))


def label(value, row=0, column=0):
    return record(0x0204, prefix(row, column) + unicode_text(value))


def workbook(cells=b"", strings=(), extra=b"", xfs=(0,), sheets=1):
    metadata = b"".join(
        record(0x00E0, struct.pack("<HH", 0, fmt) + bytes(16)) for fmt in xfs
    )
    if strings:
        metadata += record(
            0x00FC,
            struct.pack("<II", len(strings), len(strings))
            + b"".join(unicode_text(value) for value in strings),
        )
    metadata += extra
    bounds = [
        record(0x0085, struct.pack("<IBBB", 0, 0, 0, 1) + b"\0A") for _ in range(sheets)
    ]
    offset = len(bof()) + len(metadata) + sum(map(len, bounds)) + len(record(0x000A))
    worksheets = []
    for index in range(sheets):
        bounds[index] = record(0x0085, struct.pack("<IBBB", offset, 0, 0, 1) + b"\0A")
        worksheet = bof(0x10) + cells + record(0x000A)
        worksheets.append(worksheet)
        offset += len(worksheet)
    return bof() + metadata + b"".join(bounds) + record(0x000A) + b"".join(worksheets)


def directory_entry(name, kind, start=END, size=0, child=FREE):
    data = bytearray(128)
    encoded = (name + "\0").encode("utf-16le")
    data[: len(encoded)] = encoded
    struct.pack_into("<HBBIII", data, 64, len(encoded), kind, 1, FREE, FREE, child)
    struct.pack_into("<IQ", data, 116, start, size)
    return bytes(data)


def compound(stream, *, mini=False, fragmented=False, version=3, name="Workbook"):
    size = 512 if version == 3 else 4096
    header = bytearray(size)
    header[:8] = SIGNATURE
    struct.pack_into(
        "<HHHHH", header, 24, 0x003E, version, 0xFFFE, 9 if version == 3 else 12, 6
    )
    struct.pack_into(
        "<IIIIIIIII",
        header,
        40,
        0 if version == 3 else 1,
        1,
        0,
        0,
        4096,
        2 if mini else END,
        1 if mini else 0,
        END,
        0,
    )
    struct.pack_into("<109I", header, 76, 1, *([FREE] * 108))
    if mini:
        mini_count = (len(stream) + 63) // 64
        root_data = stream.ljust(mini_count * 64, b"\0")
        mini_fat = list(range(1, mini_count)) + [END]
        mini_fat.extend([FREE] * (size // 4 - len(mini_fat)))
        root_start = 3
        wb_start, wb_size = 0, len(stream)
    else:
        root_data = stream.ljust(max(4096, len(stream)), b"\0")
        root_start = 2
        wb_start, wb_size = root_start, len(root_data)
    count = (len(root_data) + size - 1) // size
    indices = [root_start + index * (2 if fragmented else 1) for index in range(count)]
    sectors = [bytes(size) for _ in range(indices[-1] + 1)]
    table = [FREE] * (size // 4)
    table[0], table[1] = END, FAT
    if mini:
        table[2] = END
        sectors[2] = struct.pack("<" + "I" * len(mini_fat), *mini_fat)
    for index, sector in enumerate(indices):
        table[sector] = indices[index + 1] if index + 1 < count else END
        sectors[sector] = root_data[index * size : (index + 1) * size].ljust(
            size, b"\0"
        )
    root = directory_entry(
        "Root Entry",
        5,
        root_start if mini else END,
        len(root_data) if mini else 0,
        child=1,
    )
    sectors[0] = (root + directory_entry(name, 2, wb_start, wb_size)).ljust(size, b"\0")
    sectors[1] = struct.pack("<" + "I" * len(table), *table)
    return bytes(header) + b"".join(sectors)


@pytest.mark.parametrize(
    "mini,fragmented,version",
    [
        (False, False, 3),
        (False, True, 3),
        (True, False, 3),
        (True, True, 3),
        (False, False, 4),
        (True, False, 4),
    ],
)
@pytest.mark.parametrize("name", ["Workbook", "Book"])
def test_compound_streams_and_sector_sizes(mini, fragmented, version, name):
    stream = workbook(
        label("Hervé 🌿") + number(12.5, column=1) + label("filler" * 100, row=1)
    )
    data = compound(
        stream, mini=mini, fragmented=fragmented, version=version, name=name
    )
    assert workbook_stream(data).startswith(stream)
    sheet = read_xls(data)[0]
    assert sheet.row_values(0) == ["Hervé 🌿", 12.5]


@pytest.mark.parametrize(
    "packed,expected",
    [
        (2, 0.0),
        (0xFFFFFFFE, -1.0),
        ((1234 << 2) | 3, 12.34),
        ((-1234 << 2) & 0xFFFFFFFF | 3, -12.34),
        (0x3FF80000, 1.5),
        (0x3FF80001, 0.015),
    ],
)
def test_rk_integer_double_sign_and_scale(packed, expected):
    sheet = read_xls(workbook(record(0x027E, prefix() + struct.pack("<I", packed))))[0]
    assert sheet.cell(0, 0).value == expected


def test_multiple_rk_and_sparse_rows():
    data = struct.pack("<HHHIHIH", 2, 1, 0, (25 << 2) | 2, 0, (1250 << 2) | 3, 2)
    sheet = read_xls(workbook(record(0x00BD, data)))[0]
    assert (sheet.nrows, sheet.ncols) == (3, 3)
    assert sheet.row_values(2) == ["", 25.0, 12.5]


def test_shared_strings_and_unicode_continuation_encoding_switch():
    # A compressed string continues as UTF-16, then a new string begins.
    first = struct.pack("<IIHB", 2, 2, 3, 0) + b"A"
    continuation = b"\1" + "éΩ".encode("utf-16le") + unicode_text("next")
    sst = record(0x00FC, first) + record(0x003C, continuation)
    cells = record(0x00FD, prefix() + struct.pack("<I", 0)) + record(
        0x00FD, prefix(column=1) + struct.pack("<I", 1)
    )
    sheet = read_xls(workbook(cells, extra=sst))[0]
    assert sheet.row_values(0) == ["AéΩ", "next"]


def test_surrogate_pair_across_string_continuation():
    data = "🌿".encode("utf-16le")
    sst = record(0x00FC, struct.pack("<IIHB", 1, 1, 2, 1) + data[:2]) + record(
        0x003C, b"\1" + data[2:]
    )
    assert (
        read_xls(workbook(record(0x00FD, prefix() + bytes(4)), extra=sst))[0]
        .cell(0, 0)
        .value
        == "🌿"
    )


def test_rich_text_and_extension_continuation_has_no_character_flag():
    first = struct.pack("<IIHBHI", 2, 2, 1, 0x0C, 1, 3) + b"A" + b"\0\0"
    continuation = b"\0\0xyz" + unicode_text("B")
    extra = record(0x00FC, first) + record(0x003C, continuation)
    cells = record(0x00FD, prefix() + bytes(4)) + record(
        0x00FD, prefix(column=1) + struct.pack("<I", 1)
    )
    assert read_xls(workbook(cells, extra=extra))[0].row_values(0) == ["A", "B"]


@pytest.mark.parametrize(
    "kind,value,expected", [(1, 1, "boolean"), (2, 7, "error"), (3, 0, "text")]
)
def test_formula_cached_non_numeric_results(kind, value, expected):
    cached = bytes([kind, 0, value, 0, 0, 0, 255, 255])
    sheet = read_xls(workbook(record(0x0006, prefix() + cached + bytes(8))))[0]
    assert sheet.cell(0, 0).kind == expected
    assert sheet.cell(0, 0).value == ("" if kind == 3 else value)


def test_formula_cached_numeric_and_string_results():
    numeric = record(0x0006, prefix() + struct.pack("<d", 12.5) + bytes(8))
    text = record(0x0006, prefix(column=1) + bytes(6) + b"\xff\xff" + bytes(8))
    text += record(0x0207, unicode_text("saved result"))
    assert read_xls(workbook(numeric + text))[0].row_values(0) == [12.5, "saved result"]


@pytest.mark.parametrize(
    "format_id,format_text,expected",
    [
        (14, None, True),
        (46, None, True),
        (164, "yyyy-mm-dd", True),
        (165, "[h]:mm:ss", True),
        (166, '0.0 "Ma"', False),
        (167, "0.0\\m", False),
        (168, "[Red][>=0]0.0", False),
        (169, "0.0E+00", False),
        (170, "[$-409]m/d/yyyy", True),
        (171, "0.0_ m", True),
        (172, "0.0_m", False),
    ],
)
def test_number_formats_preserve_dates_without_matching_literals(
    format_id, format_text, expected
):
    formats = {} if format_text is None else {format_id: format_text}
    assert is_date_format(format_id, formats) is expected
    extra = (
        b""
        if format_text is None
        else record(0x041E, struct.pack("<H", format_id) + unicode_text(format_text))
    )
    cell = read_xls(workbook(number(12.5), extra=extra, xfs=(format_id,)))[0].cell(0, 0)
    assert cell.kind == ("date" if expected else "number")


@pytest.mark.parametrize(
    "format_text,expected",
    [("[" * 30000 + "0.0", "number"), ("[" * 30000 + "yyyy", "date")],
    ids=["numeric-format", "date-format"],
)
def test_long_unclosed_number_formats_do_not_stall_workbook_reading(
    tmp_path, format_text, expected
):
    extra = record(0x041E, struct.pack("<H", 164) + unicode_text(format_text))
    path = tmp_path / "malformed-format.xls"
    path.write_bytes(
        workbook(
            b"".join(number(row + 0.5, row=row) for row in range(10000)),
            extra=extra,
            xfs=(164,),
        )
    )
    # Keep malformed-format nontermination confined to a disposable process.
    result = subprocess.run(
        [
            sys.executable,
            "-c",
            "from pathlib import Path; import sys; "
            "from nwkit.xls_reader import read_xls; "
            "sheet = read_xls(Path(sys.argv[1]).read_bytes())[0]; "
            "assert sheet.nrows == 10000; "
            "assert sheet.cell(9999, 0).value == 9999.5; "
            "assert all(cell.kind == sys.argv[2] for cell in sheet.cells.values())",
            str(path),
            expected,
        ],
        capture_output=True,
        text=True,
        timeout=30,
        check=False,
    )
    assert result.returncode == 0, result.stderr


@pytest.mark.parametrize("version", [2, 3, 4, 5])
def test_legacy_codepage_text_and_number_cells(version):
    attributes = bytes(3) if version == 2 else bytes(2)
    header = struct.pack("<HH", 0, 0) + attributes
    value = "Hervé".encode("cp1252")
    text = record(
        0x0004 if version == 2 else 0x0204,
        header
        + (bytes([len(value)]) if version == 2 else struct.pack("<H", len(value)))
        + value,
    )
    numeric_header = struct.pack("<HH", 0, 1) + attributes
    numeric = record(
        0x0003 if version == 2 else 0x0203, numeric_header + struct.pack("<d", 12.5)
    )
    data = (
        bof(0x10, version)
        + record(0x0042, struct.pack("<H", 1252))
        + text
        + numeric
        + record(0x000A)
    )
    assert read_xls(data)[0].row_values(0) == ["Hervé", 12.5]


def test_biff5_workbook_globals_and_custom_date_format():
    codepage = record(0x0042, struct.pack("<H", 1252))
    xf = record(0x00E0, struct.pack("<HH", 0, 164) + bytes(12))
    fmt = record(0x041E, struct.pack("<HB", 164, 10) + b"yyyy-mm-dd")
    bounds = record(0x0085, struct.pack("<IBBB", 0, 0, 0, 1) + b"A")
    offset = len(bof(version=5) + codepage + xf + fmt + bounds + record(0x000A))
    bounds = record(0x0085, struct.pack("<IBBB", offset, 0, 0, 1) + b"A")
    stream = bof(version=5) + codepage + xf + fmt + bounds + record(0x000A)
    stream += bof(0x10, 5) + number(12.5) + record(0x000A)
    cell = read_xls(compound(stream, mini=True, name="Book"))[0].cell(0, 0)
    assert cell.value == 12.5 and cell.kind == "date"


@pytest.mark.parametrize("version", [2, 3, 4, 5])
def test_legacy_formula_string_continuation(version):
    header = struct.pack("<HH", 0, 0) + (bytes(3) if version == 2 else bytes(2))
    formula = record(6, header + bytes(6) + b"\xff\xff" + bytes(8))
    length = b"\x05" if version == 2 else b"\x05\x00"
    cached = record(0x0007 if version == 2 else 0x0207, length + b"He")
    cached += record(0x003C, b"rv\xe9")
    stream = (
        bof(0x10, version)
        + record(0x0042, struct.pack("<H", 1252))
        + formula
        + cached
        + record(0x000A)
    )
    assert read_xls(stream)[0].cell(0, 0).value == "Hervé"


def test_biff4_workbook_with_embedded_sheet_substreams():
    sheets = []
    bounds = b""
    for name, value in [(b"A", 12.5), (b"B", 20.0)]:
        worksheet = bof(0x10, 4) + number(value) + record(0x000A)
        sheets.append(
            record(0x008F, struct.pack("<IB", len(worksheet), len(name)) + name)
            + worksheet
        )
        bounds += record(0x0085, bytes([len(name)]) + name)
    codepage = record(0x0042, struct.pack("<H", 1252))
    offset = len(bof(0x100, 4) + codepage + bounds) + 8
    stream = (
        bof(0x100, 4) + codepage + bounds + record(0x008E, struct.pack("<I", offset))
    )
    stream += b"".join(sheets) + record(0x000A)
    result = read_xls(stream)
    assert [sheet.cell(0, 0).value for sheet in result] == [12.5, 20.0]


def test_multiple_worksheets_are_preserved_for_caller_validation():
    assert len(read_xls(workbook(number(12.5), sheets=2))) == 2


@pytest.mark.parametrize("reverse", [False, True])
@pytest.mark.parametrize("second_kind", [0, 2])
def test_sheet_directory_cannot_alias_nested_substreams(reverse, second_kind):
    # Both directory entries point into one nested BOF/EOF sequence. Parsing
    # each entry independently would revisit the same bytes repeatedly.
    first = len(bof()) + 2 * 13 + 4
    entries = [(first, 0), (first + len(bof(0x10)), second_kind)]
    if reverse:
        entries.reverse()
    bounds = b"".join(
        record(0x0085, struct.pack("<IBBB", offset, 0, kind, 1) + b"\0A")
        for offset, kind in entries
    )
    stream = (
        bof()
        + bounds
        + record(0x000A)
        + bof(0x10)
        + bof(0x10)
        + number(12.5)
        + record(0x000A) * 2
    )
    with pytest.raises(ValueError, match="missing BIFF EOF|invalid nested BIFF"):
        read_xls(stream)


def test_physical_sheet_order_can_differ_from_directory_order():
    worksheets = [bof(0x10) + number(value) + record(0x000A) for value in (12.5, 20.0)]
    first = len(bof()) + 2 * 13 + 4
    bounds = b"".join(
        record(0x0085, struct.pack("<IBBB", offset, 0, 0, 1) + b"\0A")
        for offset in (first + len(worksheets[0]), first)
    )
    stream = bof() + bounds + record(0x000A) + b"".join(worksheets)
    assert [sheet.cell(0, 0).value for sheet in read_xls(stream)] == [20.0, 12.5]


def test_nested_chart_does_not_hide_following_worksheet_cells():
    chart = bof(0x20) + number(99.0) + record(0x000A)
    cells = number(12.5) + chart + number(20.0, column=1)
    assert read_xls(workbook(cells, sheets=2))[0].row_values(0) == [12.5, 20.0]


def test_nested_worksheet_cannot_silently_hide_cell_records():
    nested = bof(0x10) + number(99.0, row=2) + record(0x000A)
    with pytest.raises(ValueError, match="invalid nested BIFF substream"):
        read_xls(workbook(number(12.5) + nested))


def test_biff2_extended_format_index_and_integer_cell():
    xf = record(0x0043, bytes([0, 0, 14, 0]))
    attributes = bytes([63, 0, 0])
    cells = b"".join(
        record(
            0x0002, struct.pack("<HH", 0, column) + attributes + struct.pack("<H", 25)
        )
        for column in (0, 1)
    )
    stream = bof(0x10, 2) + xf + record(0x0044, bytes(2)) + cells + record(0x000A)
    sheet = read_xls(stream)[0]
    assert sheet.row_values(0) == [25.0, 25.0]
    assert all(sheet.cell(0, column).kind == "date" for column in (0, 1))
    with pytest.raises(ValueError, match="missing extended BIFF2 format index"):
        read_xls(bof(0x10, 2) + xf + cells + record(0x000A))


@pytest.mark.parametrize(
    "version,rows", [(2, 16384), (3, 16384), (4, 16384), (5, 16384), (8, 65536)]
)
def test_worksheet_row_and_column_bounds_depend_on_biff_version(version, rows):
    def worksheet(row, column):
        attributes = bytes(3) if version == 2 else bytes(2)
        cell = record(
            0x0003 if version == 2 else 0x0203,
            struct.pack("<HH", row, column) + attributes + struct.pack("<d", 12.5),
        )
        return bof(0x10, version) + cell + record(0x000A)

    sheet = read_xls(worksheet(rows - 1, 255))[0]
    assert (sheet.nrows, sheet.ncols) == (rows, 256)
    assert sheet.cell(rows - 1, 255).value == 12.5
    with pytest.raises(ValueError, match="cell is out of bounds"):
        read_xls(worksheet(0, 256))
    if version < 8:
        with pytest.raises(ValueError, match="cell is out of bounds"):
            read_xls(worksheet(rows, 0))


@pytest.mark.parametrize(
    "codepage,text,message",
    [
        (65535, b"", "unsupported legacy code page"),
        (932, b"\x81", "invalid legacy string"),
    ],
)
def test_invalid_legacy_codepages_and_text_are_rejected(codepage, text, message):
    stream = (
        bof(0x10, 3)
        + record(0x0042, struct.pack("<H", codepage))
        + record(0x0204, prefix() + struct.pack("<H", len(text)) + text)
        + record(0x000A)
    )
    with pytest.raises(ValueError, match=message):
        read_xls(stream)


@pytest.mark.parametrize("allocation", ["directory", "root-mini-stream", "mini-fat"])
def test_regular_workbook_cannot_overlap_structural_streams(allocation):
    stream = bof(0x10) + number(12.5) + record(0x000A)
    data = bytearray(compound(stream))
    if allocation == "directory":
        # A directory chain also claims all Workbook sectors. Their bytes happen
        # to look like unused directory entries, which must not conceal overlap.
        struct.pack_into("<I", data, 1024, 2)
    elif allocation == "root-mini-stream":
        struct.pack_into("<IQ", data, 512 + 116, 2, 4096)
    else:
        struct.pack_into("<II", data, 60, 2, 8)
    with pytest.raises(ValueError, match="stream allocations overlap"):
        read_xls(bytes(data))


@pytest.mark.parametrize(
    "codepage,encoding,value",
    [(10079, "mac_iceland", "Þór"), (10081, "mac_turkish", "ağaç")],
)
def test_legacy_mac_codepages(codepage, encoding, value):
    text = value.encode(encoding)
    stream = (
        bof(0x10, 5)
        + record(0x0042, struct.pack("<H", codepage))
        + record(0x0204, prefix() + struct.pack("<H", len(text)) + text)
        + record(0x000A)
    )
    assert read_xls(stream)[0].cell(0, 0).value == value


@pytest.mark.parametrize(
    "mutation",
    [
        "fat-cycle",
        "fat-out-of-bounds",
        "directory-cycle",
        "mini-cycle",
        "mini-out-of-bounds",
        "oversized-stream",
        "fat-marker",
        "bad-byte-order",
        "bad-sector-shift",
        "missing-workbook",
    ],
)
def test_invalid_compound_allocations_are_rejected(mutation):
    mini = mutation.startswith("mini-")
    data = bytearray(compound(workbook(label("x")), mini=mini))
    if mutation == "fat-cycle":
        struct.pack_into("<I", data, 1024 + 2 * 4, 2)
    elif mutation == "fat-out-of-bounds":
        struct.pack_into("<I", data, 1024 + 2 * 4, 99999)
    elif mutation == "directory-cycle":
        struct.pack_into("<I", data, 512 + 128 + 72, 1)
    elif mutation == "mini-cycle":
        struct.pack_into("<I", data, 1536, 0)
    elif mutation == "mini-out-of-bounds":
        struct.pack_into("<I", data, 512 + 128 + 116, 100000)
    elif mutation == "oversized-stream":
        struct.pack_into("<Q", data, 512 + 128 + 120, 0xFFFFFFFF)
    elif mutation == "fat-marker":
        struct.pack_into("<I", data, 1024 + 4, FREE)
    elif mutation == "bad-byte-order":
        struct.pack_into("<H", data, 28, 0xFEFF)
    elif mutation == "bad-sector-shift":
        struct.pack_into("<H", data, 30, 15)
    else:
        data[512 + 128 : 512 + 128 + 2] = "Z".encode("utf-16le")
    with pytest.raises(ValueError, match="Invalid AngioCal XLS"):
        read_xls(bytes(data))


@pytest.mark.parametrize(
    "data",
    [
        b"",
        b"PK\x03\x04not-xls",
        SIGNATURE,
        bof()[:-1],
        bof() + number(1),
        bof() + record(0x002F) + record(0x000A),
        workbook(record(0x00FD, prefix() + struct.pack("<I", 123))),
        workbook(number(1) + number(2)),
        workbook(record(0x00BD, b"bad")),
        workbook(record(0x0006, prefix() + bytes(6) + b"\xff\xff")),
        workbook(record(0x0207, unicode_text("orphan"))),
        workbook(extra=record(0x00FC, struct.pack("<II", 0xFFFFFFFF, 0xFFFFFFFF))),
    ],
)
def test_truncated_invalid_or_unsupported_biff_is_rejected(data):
    with pytest.raises(ValueError):
        read_xls(data)


@pytest.mark.parametrize("continuation", [b"\x02A", b"\x01A", b"\x01\x00\xd8"])
def test_malformed_string_continuations_are_rejected(continuation):
    sst = record(0x00FC, struct.pack("<IIHB", 1, 1, 2, 0) + b"A") + record(
        0x003C, continuation
    )
    with pytest.raises(ValueError):
        read_xls(workbook(extra=sst))


@given(st.binary(max_size=600))
@settings(max_examples=150, deadline=None)
def test_arbitrary_bytes_produce_only_a_validation_error(data):
    try:
        read_xls(data)
    except ValueError:
        pass


@pytest.mark.parametrize("container", [False, True])
@given(
    index=st.integers(min_value=0, max_value=6000),
    value=st.integers(min_value=0, max_value=255),
)
@settings(max_examples=250, deadline=None)
def test_mutated_valid_inputs_never_escape_validation(container, index, value):
    data = workbook(label("Fossil 🌿") + number(12.5, column=1))
    if container:
        data = compound(data, mini=True)
    mutated = bytearray(data)
    mutated[index % len(data)] = value
    try:
        read_xls(bytes(mutated))
    except ValueError:
        pass

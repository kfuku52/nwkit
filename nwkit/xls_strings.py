"""BIFF text and number-format decoding used by the AngioCal XLS reader.

References: [MS-XLS] XLUnicodeRichExtendedString/Continue, and OpenOffice's
Excel file-format documentation for pre-Unicode BIFF strings.
"""

import codecs
import re

from nwkit.xls_compound import require, uint16, uint32


class Segments:
    """Keep CONTINUE boundaries: character continuations have an extra flag."""

    def __init__(self, parts):
        self.parts = parts
        self.index = 0
        self.offset = 0
        self.remaining = sum(len(part) for part in parts)

    def _advance(self):
        while self.index < len(self.parts) and self.offset == len(
            self.parts[self.index]
        ):
            self.index += 1
            self.offset = 0
        require(self.index < len(self.parts), "truncated string data.")

    def read(self, size):
        require(0 <= size <= self.remaining, "truncated string data.")
        result = []
        while size:
            self._advance()
            part = self.parts[self.index]
            count = min(size, len(part) - self.offset)
            result.append(part[self.offset : self.offset + count])
            self.offset += count
            self.remaining -= count
            size -= count
        return b"".join(result)

    def characters(self, count, wide):
        result = []
        while count:
            if self.index >= len(self.parts) or self.offset == len(
                self.parts[self.index]
            ):
                self._advance()
                flag = self.read(1)[0]
                require(flag in {0, 1}, "invalid string continuation flag.")
                wide = bool(flag)
            part = self.parts[self.index]
            width = 2 if wide else 1
            available = (len(part) - self.offset) // width
            require(available > 0, "split UTF-16 code unit.")
            take = min(count, available)
            data = self.read(take * width)
            if not wide:
                data = data.decode("latin1").encode("utf-16le")
            result.append(data)
            count -= take
        try:
            return b"".join(result).decode("utf-16le")
        except UnicodeError as exc:
            raise ValueError(
                "Invalid AngioCal XLS workbook: invalid UTF-16 string."
            ) from exc

    def unicode(self, length_bytes=2):
        count_data = self.read(length_bytes)
        count = uint16(count_data) if length_bytes == 2 else count_data[0]
        flag = self.read(1)[0]
        require(flag & ~0x0D == 0, "invalid Unicode string flags.")
        runs = uint16(self.read(2)) if flag & 8 else 0
        extension = uint32(self.read(4)) if flag & 4 else 0
        value = self.characters(count, bool(flag & 1))
        self.read(4 * runs)
        self.read(extension)
        return value


def shared_strings(parts):
    reader = Segments(parts)
    total, unique = uint32(reader.read(4)), uint32(reader.read(4))
    require(
        unique <= total and unique <= reader.remaining // 3,
        "invalid shared-string count.",
    )
    strings = [reader.unicode() for _ in range(unique)]
    require(reader.remaining == 0, "extra shared-string data.")
    return strings


def codepage_encoding(codepage):
    special = {
        10000: "mac_roman",
        10006: "mac_greek",
        10007: "mac_cyrillic",
        10029: "mac_latin2",
        10079: "mac_iceland",
        10081: "mac_turkish",
        20127: "ascii",
        32768: "mac_roman",
        32769: "cp1252",
        65001: "utf-8",
    }
    encoding = (
        "latin1" if codepage is None else special.get(codepage, "cp" + str(codepage))
    )
    try:
        codecs.lookup(encoding)
    except LookupError as exc:
        raise ValueError(
            "Invalid AngioCal XLS workbook: unsupported legacy code page."
        ) from exc
    return encoding


def legacy_string(data, offset, encoding, length_bytes=2):
    require(
        offset >= 0 and offset + length_bytes <= len(data),
        "truncated legacy string length.",
    )
    count = uint16(data, offset) if length_bytes == 2 else data[offset]
    start = offset + length_bytes
    require(start + count <= len(data), "truncated legacy string.")
    try:
        return data[start : start + count].decode(encoding)
    except UnicodeError as exc:
        raise ValueError(
            "Invalid AngioCal XLS workbook: invalid legacy string."
        ) from exc


def is_date_format(format_id, formats):
    if format_id not in formats:
        return format_id in set(range(14, 23)) | set(range(27, 37)) | {
            45,
            46,
            47,
        } | set(range(50, 59)) | set(range(71, 82))
    text = formats[format_id]
    # Ignore literals, escaped/fill characters, colors, conditions and locale
    # directives. Bracketed elapsed-time tokens are date/time formats.
    text = re.sub(r'"[^\"]*"|\\.|_.|\*.', "", text)
    # Find disjoint bracket spans explicitly. A regex with an unbounded
    # negated class retries every opening bracket when no closing bracket
    # exists, making malformed number formats quadratic to classify.
    pieces = []
    start = 0
    while True:
        left = text.find("[", start)
        right = text.find("]", left + 1) if left >= 0 else -1
        if right < 0:
            pieces.append(text[start:])
            break
        pieces.append(text[start:left])
        if re.fullmatch(r"[hms]+", text[left + 1 : right], re.IGNORECASE):
            return True
        start = right + 1
    return re.search(r"[ymdhs]", "".join(pieces), re.IGNORECASE) is not None

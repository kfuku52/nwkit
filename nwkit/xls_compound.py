"""Read bounded workbook streams from MS-CFB containers, without dependencies.

Layout reference: Microsoft [MS-CFB], sections 2.2 through 2.6. This reader
only extracts root-level Workbook/Book streams; it does not execute OLE objects.
"""

import struct
from dataclasses import dataclass

SIGNATURE = b"\xd0\xcf\x11\xe0\xa1\xb1\x1a\xe1"
MAX_XLS_BYTES = 5 * 1024 * 1024
FREE = 0xFFFFFFFF
END = 0xFFFFFFFE
FAT = 0xFFFFFFFD
DIFAT = 0xFFFFFFFC


def require(condition, message):
    if not condition:
        raise ValueError("Invalid AngioCal XLS workbook: " + message)


def uint16(data, offset=0):
    require(offset >= 0 and offset + 2 <= len(data), "truncated integer.")
    return struct.unpack_from("<H", data, offset)[0]


def uint32(data, offset=0):
    require(offset >= 0 and offset + 4 <= len(data), "truncated integer.")
    return struct.unpack_from("<I", data, offset)[0]


def chain(start, table, limit):
    result: list[int] = []
    seen = set()
    current = start
    while current != END:
        require(0 <= current < len(table), "sector chain is out of bounds.")
        require(current not in seen, "cyclic sector chain.")
        require(len(result) < limit, "sector chain exceeds its size limit.")
        seen.add(current)
        result.append(current)
        current = table[current]
    return result


@dataclass(frozen=True)
class DirectoryEntry:
    name: str
    kind: int
    left: int
    right: int
    child: int
    start: int
    size: int


class CompoundFile:
    def __init__(self, data):
        require(len(data) <= MAX_XLS_BYTES, "input exceeds the size limit.")
        require(len(data) >= 512 and data.startswith(SIGNATURE), "invalid OLE header.")
        self.data = data
        self.version = uint16(data, 26)
        require(uint16(data, 28) == 0xFFFE, "unsupported OLE byte order.")
        shift = uint16(data, 30)
        require((self.version, shift) in {(3, 9), (4, 12)}, "unsupported OLE version.")
        require(uint16(data, 32) == 6, "invalid mini-sector size.")
        require(uint32(data, 56) == 4096, "invalid mini-stream cutoff.")
        self.sector_size = 1 << shift
        require(len(data) % self.sector_size == 0, "truncated OLE sector.")
        self.sector_count = len(data) // self.sector_size - 1
        self.reserved = set()
        self.used = set()
        self.fat = self._fat()
        directory = self._regular(uint32(data, 48))
        if self.version == 4:
            require(
                len(directory) // self.sector_size == uint32(data, 40),
                "directory sector count does not match.",
            )
        self.entries = self._directory(directory)
        require(self.entries and self.entries[0].kind == 5, "missing root directory.")
        self.mini_stream, self.mini_fat = self._mini_storage()

    def _sector(self, index):
        require(0 <= index < self.sector_count, "OLE sector is out of bounds.")
        offset = (index + 1) * self.sector_size
        return self.data[offset : offset + self.sector_size]

    def _fat(self):
        count = uint32(self.data, 44)
        require(0 < count <= self.sector_count, "invalid FAT sector count.")
        indices = [
            value
            for value in struct.unpack_from("<109I", self.data, 76)
            if value != FREE
        ]
        next_sector = uint32(self.data, 68)
        difat_count = uint32(self.data, 72)
        require(difat_count <= self.sector_count, "invalid DIFAT sector count.")
        difat_sectors = set()
        for _ in range(difat_count):
            require(next_sector not in difat_sectors, "cyclic DIFAT chain.")
            difat_sectors.add(next_sector)
            values = struct.unpack(
                "<" + "I" * (self.sector_size // 4), self._sector(next_sector)
            )
            indices.extend(value for value in values[:-1] if value != FREE)
            next_sector = values[-1]
        require(
            next_sector in {END, FREE} if not difat_count else next_sector == END,
            "invalid DIFAT terminator.",
        )
        require(
            len(indices) == count and len(set(indices)) == count,
            "invalid FAT sector list.",
        )
        require(not set(indices) & difat_sectors, "FAT and DIFAT sectors overlap.")
        table: list[int] = []
        for index in indices:
            table.extend(
                struct.unpack("<" + "I" * (self.sector_size // 4), self._sector(index))
            )
        require(len(table) >= self.sector_count, "FAT does not cover the file.")
        require(
            all(table[index] == FAT for index in indices), "invalid FAT sector marker."
        )
        require(
            all(table[index] == DIFAT for index in difat_sectors),
            "invalid DIFAT sector marker.",
        )
        self.reserved = set(indices) | difat_sectors
        return table[: self.sector_count]

    def _regular(self, start, size=None):
        if size is not None:
            require(0 <= size <= len(self.data), "stream exceeds the file size.")
            limit = (size + self.sector_size - 1) // self.sector_size
        else:
            limit = self.sector_count
        indices = chain(start, self.fat, limit)
        require(
            not self.reserved.intersection(indices),
            "stream overlaps allocation tables.",
        )
        if size is not None:
            require(
                len(indices) == limit, "stream sector count does not match its size."
            )
        require(not self.used.intersection(indices), "stream allocations overlap.")
        self.used.update(indices)
        result = b"".join(self._sector(index) for index in indices)
        return result if size is None else result[:size]

    def _directory(self, data):
        entries = []
        for offset in range(0, len(data), 128):
            entry = data[offset : offset + 128]
            kind = entry[66]
            if kind == 0:
                entries.append(DirectoryEntry("", 0, FREE, FREE, FREE, END, 0))
                continue
            require(kind in {1, 2, 5}, "invalid directory entry type.")
            length = uint16(entry, 64)
            require(
                2 <= length <= 64 and length % 2 == 0, "invalid directory name length."
            )
            require(
                entry[length - 2 : length] == b"\0\0", "unterminated directory name."
            )
            try:
                name = entry[: length - 2].decode("utf-16le")
            except UnicodeError as exc:
                raise ValueError(
                    "Invalid AngioCal XLS workbook: invalid directory name."
                ) from exc
            size = (
                uint32(entry, 120)
                if self.version == 3
                else struct.unpack_from("<Q", entry, 120)[0]
            )
            entries.append(
                DirectoryEntry(
                    name,
                    kind,
                    uint32(entry, 68),
                    uint32(entry, 72),
                    uint32(entry, 76),
                    uint32(entry, 116),
                    size,
                )
            )
        return entries

    def _root_streams(self):
        pending = [self.entries[0].child]
        seen = set()
        result = []
        while pending:
            index = pending.pop()
            if index == FREE:
                continue
            require(0 < index < len(self.entries), "directory child is out of bounds.")
            require(index not in seen, "cyclic directory tree.")
            seen.add(index)
            entry = self.entries[index]
            require(entry.kind in {1, 2}, "invalid root child.")
            pending.extend([entry.left, entry.right])
            if entry.kind == 2:
                result.append(entry)
        return result

    def _mini_storage(self):
        root = self.entries[0]
        count = uint32(self.data, 64)
        require(count <= self.sector_count, "invalid mini-FAT sector count.")
        allocation = (
            self._regular(uint32(self.data, 60), count * self.sector_size)
            if count
            else b""
        )
        table = struct.unpack("<" + "I" * (len(allocation) // 4), allocation)
        mini_stream = self._regular(root.start, root.size) if root.size else b""
        return mini_stream, table

    def _mini(self, entry):
        mini_stream, table = self.mini_stream, self.mini_fat
        require(entry.size <= len(mini_stream), "mini-stream exceeds the root stream.")
        limit = (entry.size + 63) // 64
        indices = chain(entry.start, table, limit)
        require(len(indices) == limit, "mini-stream sector count does not match.")
        require(
            all((index + 1) * 64 <= len(mini_stream) for index in indices),
            "mini-sector is out of bounds.",
        )
        return b"".join(
            mini_stream[index * 64 : (index + 1) * 64] for index in indices
        )[: entry.size]

    def workbook(self):
        matches = [
            entry
            for entry in self._root_streams()
            if entry.name.casefold() in {"workbook", "book"}
        ]
        require(len(matches) == 1, "expected one root Workbook or Book stream.")
        entry = matches[0]
        require(entry.size > 0, "empty workbook stream.")
        return (
            self._mini(entry)
            if entry.size < 4096
            else self._regular(entry.start, entry.size)
        )


def workbook_stream(data):
    require(len(data) <= MAX_XLS_BYTES, "input exceeds the size limit.")
    if data.startswith(SIGNATURE):
        return CompoundFile(data).workbook()
    return data

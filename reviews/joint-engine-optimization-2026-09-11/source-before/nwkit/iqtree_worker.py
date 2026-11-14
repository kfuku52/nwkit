"""Client for the optional external worker linked to an unmodified IQ-TREE library."""

import queue
import subprocess
import threading
import time
import weakref
from collections import deque

import numpy as np

PREFIX = "NWKIT_IQTREE_V1 "


def _read_lines(pipe, messages):
    try:
        for line in pipe:
            messages.put(line)
    finally:
        messages.put(None)
        pipe.close()


def _stop(process):
    if process.poll() is None:
        try:
            process.stdin.write("QUIT\n")
            process.stdin.flush()
            process.wait(timeout=2)
        except (OSError, subprocess.TimeoutExpired):
            process.kill()
            process.wait()
    try:
        process.stdin.close()
    except OSError:
        pass


class IQTreeWorker:
    """Fail closed on protocol errors; never silently restart/reload a session."""

    def __init__(self, command, *, timeout=600):
        self.timeout = timeout
        self.messages: queue.Queue[str | None] = queue.Queue()
        self.tail: deque[str] = deque(maxlen=20)
        self.process = subprocess.Popen(
            command,
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,
        )
        self._cleanup = weakref.finalize(self, _stop, self.process)
        self.reader = threading.Thread(
            target=_read_lines,
            args=(self.process.stdout, self.messages),
            daemon=True,
        )
        self.reader.start()
        try:
            ready = self._receive().split()
            if len(ready) != 2 or ready[0] != "READY":
                raise ValueError("Invalid IQ-TREE library worker handshake.")
            count = int(ready[1])
            if count < 3:
                raise ValueError("Invalid IQ-TREE library worker branch count.")
            self.splits, self.lengths = [], []
            for _ in range(count):
                fields = self._receive().split()
                if len(fields) < 3 or fields[0] != "EDGE":
                    raise ValueError("Invalid IQ-TREE library worker edge record.")
                self.lengths.append(float(fields[1]))
                self.splits.append(frozenset(fields[2:]))
            if self._receive() != "END":
                raise ValueError("Incomplete IQ-TREE library worker handshake.")
            if not np.isfinite(self.lengths).all() or np.any(
                np.array(self.lengths) <= 0
            ):
                raise ValueError("Invalid IQ-TREE library worker initial lengths.")
        except BaseException:
            self.close()
            raise

    def _receive(self):
        deadline = time.monotonic() + self.timeout
        while True:
            try:
                line = self.messages.get(timeout=max(0, deadline - time.monotonic()))
            except queue.Empty as exc:
                raise ValueError("IQ-TREE library worker timed out.") from exc
            if line is None:
                raise ValueError(
                    "IQ-TREE library worker exited before completing a request. "
                    + "".join(self.tail)[-6000:]
                )
            if line.startswith(PREFIX):
                message = line[len(PREFIX) :].strip()
                if message.startswith("ERROR "):
                    raise ValueError("IQ-TREE library worker: " + message[6:])
                return message
            if line.startswith("NWKIT_IQTREE_"):
                raise ValueError("Unsupported IQ-TREE library worker protocol version.")
            self.tail.append(line)

    def evaluate(self, lengths):
        if not self._cleanup.alive or self.process.poll() is not None:
            self.close()
            raise ValueError("IQ-TREE library worker is closed.")
        lengths = np.asarray(lengths, dtype=float)
        if (
            lengths.shape != (len(self.splits),)
            or not np.isfinite(lengths).all()
            or np.any(lengths <= 0)
        ):
            raise ValueError(
                "IQ-TREE library worker requires one positive length per edge."
            )
        try:
            assert self.process.stdin is not None
            self.process.stdin.write(
                "EVAL " + " ".join(format(x, ".17g") for x in lengths) + "\n"
            )
            self.process.stdin.flush()
            fields = self._receive().split()
            expected = 2 + 2 * len(lengths)
            if len(fields) != expected or fields[0] != "VALUE":
                raise ValueError("Invalid IQ-TREE library worker derivative response.")
            values = np.array([float(x) for x in fields[1:]])
            if not np.isfinite(values).all():
                raise ValueError("Nonfinite IQ-TREE library worker derivatives.")
        except (ValueError, OSError) as exc:
            self.close()
            raise ValueError(str(exc)) from exc
        diagonal = -values[1 + len(lengths) :]
        return -values[0], -values[1 : 1 + len(lengths)], diagonal

    def close(self):
        self._cleanup()
        self.reader.join(timeout=2)

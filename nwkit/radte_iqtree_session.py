"""Persistent IQ-TREE protocol client; each instance owns one loaded tree."""

import queue
import subprocess
import threading
import time
import weakref
from collections import deque

import numpy as np

PREFIX = "IQTREE_SESSION_V1 "


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


class IQTreeSession:
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
                raise ValueError("Invalid IQ-TREE session handshake.")
            count = int(ready[1])
            if count < 3:
                raise ValueError("Invalid IQ-TREE session branch count.")
            self.splits, self.lengths = [], []
            for _ in range(count):
                fields = self._receive().split()
                if len(fields) < 3 or fields[0] != "EDGE":
                    raise ValueError("Invalid IQ-TREE session edge record.")
                self.lengths.append(float(fields[1]))
                self.splits.append(frozenset(fields[2:]))
            if self._receive() != "END":
                raise ValueError("Incomplete IQ-TREE session handshake.")
            if not np.isfinite(self.lengths).all() or np.any(
                np.array(self.lengths) <= 0
            ):
                raise ValueError("Invalid IQ-TREE session initial lengths.")
        except BaseException:
            self.close()
            raise

    def _receive(self):
        deadline = time.monotonic() + self.timeout
        while True:
            try:
                line = self.messages.get(timeout=max(0, deadline - time.monotonic()))
            except queue.Empty as exc:
                raise ValueError("IQ-TREE session timed out.") from exc
            if line is None:
                raise ValueError(
                    "IQ-TREE session exited; requires an IQ-TREE build with "
                    "--likelihood-session. " + "".join(self.tail)[-6000:]
                )
            if line.startswith(PREFIX):
                message = line[len(PREFIX) :].strip()
                if message.startswith("ERROR "):
                    raise ValueError("IQ-TREE session: " + message[6:])
                return message
            if line.startswith("IQTREE_SESSION_"):
                raise ValueError("Unsupported IQ-TREE session protocol version.")
            self.tail.append(line)

    def evaluate(self, lengths, *, second_derivatives=True):
        if not self._cleanup.alive or self.process.poll() is not None:
            self.close()
            raise ValueError("IQ-TREE session is closed.")
        lengths = np.asarray(lengths, dtype=float)
        if (
            lengths.shape != (len(self.splits),)
            or not np.isfinite(lengths).all()
            or np.any(lengths <= 0)
        ):
            raise ValueError("IQ-TREE session requires one positive length per edge.")
        try:
            assert self.process.stdin is not None
            operation = "EVAL" if second_derivatives else "SCORE"
            self.process.stdin.write(
                operation + " " + " ".join(format(x, ".17g") for x in lengths) + "\n"
            )
            self.process.stdin.flush()
            fields = self._receive().split()
            expected = 2 + (2 if second_derivatives else 1) * len(lengths)
            response = "VALUE" if second_derivatives else "SCORE"
            if len(fields) != expected or fields[0] != response:
                raise ValueError("Invalid IQ-TREE session derivative response.")
            values = np.array([float(x) for x in fields[1:]])
            if not np.isfinite(values).all():
                raise ValueError("Nonfinite IQ-TREE session derivatives.")
        except (ValueError, OSError) as exc:
            self.close()
            raise ValueError(str(exc)) from exc
        diagonal = -values[1 + len(lengths) :] if second_derivatives else None
        return -values[0], -values[1 : 1 + len(lengths)], diagonal

    def statistics(self):
        """Return cumulative worker timings; no likelihood evaluation is performed."""
        names = (
            "requests",
            "score_requests",
            "spectral_attempts",
            "stable_attempts",
            "tip_cache_reuses",
            "likelihood_seconds",
            "derivative_seconds",
        )
        try:
            if not self._cleanup.alive or self.process.poll() is not None:
                raise ValueError("IQ-TREE session is closed.")
            assert self.process.stdin is not None
            self.process.stdin.write("STATS\n")
            self.process.stdin.flush()
            fields = self._receive().split()
            if len(fields) != len(names) + 1 or fields[0] != "STATS":
                raise ValueError("Invalid IQ-TREE session statistics response.")
            values = [float(x) for x in fields[1:]]
            if not np.isfinite(values).all() or any(x < 0 for x in values):
                raise ValueError("Invalid IQ-TREE session statistics values.")
            return dict(zip(names, values, strict=True))
        except (ValueError, OSError) as exc:
            self.close()
            raise ValueError(str(exc)) from exc

    def close(self):
        self._cleanup()
        self.reader.join(timeout=2)

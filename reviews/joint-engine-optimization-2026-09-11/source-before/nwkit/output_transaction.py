"""Stage and recover related file outputs without exposing partial writes.

Files are individually replaced atomically; a handled commit failure restores
the whole set. This is not a crash-atomic multi-file filesystem transaction.
Locks serialize NWKIT writers, while inode checks reject replaced staging files.
"""

import hashlib
import os
import secrets
import shutil
import stat
import tempfile
from contextlib import ExitStack, contextmanager
from dataclasses import dataclass

from nwkit.file_paths import normalized_missing_path_key, validate_distinct_output_paths


def regular_output_mode(path):
    try:
        info = os.lstat(path)
    except FileNotFoundError:
        return None
    if not stat.S_ISREG(info.st_mode):
        raise ValueError(f"Existing output target must be a regular file: '{path}'.")
    return stat.S_IMODE(info.st_mode)


def new_output_mode(directory):
    """Read the effective creation mode without changing the process umask."""
    for _ in range(100):
        path = os.path.join(directory, f".nwkit-mode-probe-{secrets.token_hex(16)}")
        try:
            descriptor = os.open(path, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o666)
        except FileExistsError:
            continue
        try:
            return stat.S_IMODE(os.fstat(descriptor).st_mode)
        finally:
            os.close(descriptor)
            os.remove(path)
    raise FileExistsError("Could not allocate an output-mode probe.")


def validate_output_targets(paths, *, follow_symlinks=True):
    """Validate every target and parent before staging or modifying any output."""
    paths = list(paths)
    validate_distinct_output_paths([(str(path), path) for path in paths])
    targets = {}
    for path in paths:
        target = os.path.abspath(os.fspath(path))
        if follow_symlinks:
            target = os.path.realpath(target)
        regular_output_mode(target)
        parent = os.path.dirname(target)
        while not os.path.exists(parent):
            parent = os.path.dirname(parent)
        if not os.path.isdir(parent):
            raise NotADirectoryError(f"Output parent is not a directory: '{parent}'.")
        targets[path] = target
    return targets


def output_lock_path(path):
    identity = hashlib.sha256(
        normalized_missing_path_key(path).encode("utf-8")
    ).hexdigest()
    return os.path.join(os.path.dirname(path), f".nwkit-output-{identity}.lock")


@dataclass(frozen=True)
class _Stage:
    target: str
    path: str
    device: int
    inode: int
    mode: int


def _new_stage(target):
    directory = os.path.dirname(target)
    mode = regular_output_mode(target)
    if mode is None:
        mode = new_output_mode(directory)
    descriptor, path = tempfile.mkstemp(
        prefix=f".{os.path.basename(target)}.stage.", dir=directory
    )
    try:
        info = os.fstat(descriptor)
        return _Stage(target, path, info.st_dev, info.st_ino, mode)
    finally:
        os.close(descriptor)


def _open_regular(path, *, writable=False):
    access = os.O_RDWR if writable else os.O_RDONLY
    flags = access | getattr(os, "O_NOFOLLOW", 0) | getattr(os, "O_NONBLOCK", 0)
    descriptor = os.open(path, flags)
    if not stat.S_ISREG(os.fstat(descriptor).st_mode):
        os.close(descriptor)
        raise ValueError(f"Existing output target must be a regular file: '{path}'.")
    return descriptor


def _open_stage(stage):
    descriptor = _open_regular(stage.path, writable=True)
    try:
        info = os.fstat(descriptor)
        if (info.st_dev, info.st_ino) != (stage.device, stage.inode):
            raise RuntimeError("An output staging file was replaced before commit.")
        return descriptor
    except BaseException:
        os.close(descriptor)
        raise


class _StagedOutputs(dict):
    def __init__(self, targets, stages):
        self._stages = dict(zip(targets, stages, strict=True))
        super().__init__((target, stage.path) for target, stage in self._stages.items())

    def write_text(self, target, writer):
        """Write through a verified descriptor without following a replaced path."""
        descriptor = _open_stage(self._stages[target])
        try:
            os.ftruncate(descriptor, 0)
            stream = os.fdopen(descriptor, "w", encoding="utf-8", newline="")
        except BaseException:
            os.close(descriptor)
            raise
        with stream:
            writer(stream)


def _finish_stage(stage):
    # Windows FlushFileBuffers requires a handle opened with write access.
    descriptor = _open_stage(stage)
    try:
        os.fsync(descriptor)
        if hasattr(os, "fchmod"):
            os.fchmod(descriptor, stage.mode)
        else:
            os.chmod(stage.path, stage.mode)
        path_info = os.lstat(stage.path)
        if not stat.S_ISREG(path_info.st_mode) or (
            path_info.st_dev,
            path_info.st_ino,
        ) != (stage.device, stage.inode):
            raise RuntimeError("An output staging file was replaced before commit.")
    finally:
        os.close(descriptor)


def _backup_output(target):
    descriptor, path = tempfile.mkstemp(
        prefix=f".{os.path.basename(target)}.backup.", dir=os.path.dirname(target)
    )
    try:
        with os.fdopen(descriptor, "wb") as backup:
            info = os.fstat(backup.fileno())
            with os.fdopen(_open_regular(target), "rb") as source:
                mode = stat.S_IMODE(os.fstat(source.fileno()).st_mode)
                shutil.copyfileobj(source, backup, length=1024 * 1024)
            backup.flush()
            os.fsync(backup.fileno())
            if hasattr(os, "fchmod"):
                os.fchmod(backup.fileno(), mode)
        path_info = os.lstat(path)
        if not stat.S_ISREG(path_info.st_mode) or (
            path_info.st_dev,
            path_info.st_ino,
        ) != (info.st_dev, info.st_ino):
            raise RuntimeError("An output backup file was replaced before commit.")
        if not hasattr(os, "fchmod"):
            os.chmod(path, mode)
        return path
    except BaseException:
        if os.path.lexists(path):
            os.remove(path)
        raise


def replace_output(source, target):
    os.replace(source, target)


@dataclass
class _Installation:
    target: str
    staged_path: str
    backup: str | None = None
    installed: bool = False


def _restore_outputs(installations):
    failures = []
    for entry in reversed(installations):
        try:
            if entry.installed:
                if entry.backup is not None:
                    replace_output(entry.backup, entry.target)
                elif os.path.lexists(entry.target):
                    os.remove(entry.target)
        except BaseException as exc:
            failures.append((entry.target, exc))
    if failures:
        raise RuntimeError(
            "Could not restore: " + ", ".join(path for path, _ in failures)
        ) from failures[0][1]


def commit_outputs(staged_outputs, *, after_install=None):
    """Install (target, staged-path) pairs; preserve backups if recovery fails."""
    installations = []
    cleanup_backups = False
    try:
        for target, staged_path in staged_outputs:
            entry = _Installation(target, staged_path)
            installations.append(entry)
            if os.path.lexists(target):
                entry.backup = _backup_output(target)
            replace_output(staged_path, target)
            entry.installed = True
        if after_install is not None:
            after_install()
        cleanup_backups = True
    except BaseException:
        # An interrupt can arrive just after replace, before flag assignment.
        for entry in installations:
            if not entry.installed and not os.path.lexists(entry.staged_path):
                entry.installed = True
        try:
            _restore_outputs(installations)
            cleanup_backups = True
        except BaseException as exc:
            raise RuntimeError(
                "Failed to restore outputs after a commit error; backup files were preserved."
            ) from exc
        raise
    finally:
        if cleanup_backups:
            for entry in installations:
                if entry.backup is not None and os.path.lexists(entry.backup):
                    os.remove(entry.backup)
        for _, staged_path in staged_outputs:
            if os.path.lexists(staged_path):
                os.remove(staged_path)


@contextmanager
def output_transaction(
    paths, *, follow_symlinks=True, create_parents=False, after_install=None
):
    """Yield a target-to-staging-path map and commit it on successful exit."""
    # Locks also serve taxonomy/download callers; keep their implementation
    # shared, and avoid an import cycle with util's path-helper reexports.
    from nwkit.util import acquire_exclusive_lock

    targets = validate_output_targets(paths, follow_symlinks=follow_symlinks)
    stages = []
    with ExitStack() as locks:
        try:
            for target in sorted(targets.values(), key=normalized_missing_path_key):
                directory = os.path.dirname(target)
                if create_parents:
                    os.makedirs(directory, exist_ok=True)
                elif not os.path.isdir(directory):
                    raise FileNotFoundError(
                        f"Output directory does not exist: '{directory}'."
                    )
                locks.enter_context(
                    acquire_exclusive_lock(
                        output_lock_path(target), lock_label="NWKIT output"
                    )
                )
            for target in targets.values():
                stages.append(_new_stage(target))
            yield _StagedOutputs(targets, stages)
            for stage in stages:
                _finish_stage(stage)
            staged_outputs = [(stage.target, stage.path) for stage in stages]
            if after_install is None:
                commit_outputs(staged_outputs)
            else:
                commit_outputs(staged_outputs, after_install=after_install)
        finally:
            for stage in stages:
                if os.path.lexists(stage.path):
                    os.remove(stage.path)

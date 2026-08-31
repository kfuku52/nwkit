import os
import stat
from pathlib import Path

import pytest

from nwkit import output_transaction as output_module
from nwkit.output_transaction import output_transaction


@pytest.mark.parametrize("fail_after_replace", [False, True])
def test_commit_failure_recovers_existing_and_new_outputs(
    monkeypatch, tmp_path, fail_after_replace
):
    existing, new = tmp_path / "old.txt", tmp_path / "new.txt"
    existing.write_text("original")
    original_replace = output_module.replace_output
    calls = 0

    def fail_second(source, target):
        nonlocal calls
        calls += 1
        if calls == 2:
            if fail_after_replace:
                original_replace(source, target)
            raise OSError("injected commit failure")
        original_replace(source, target)

    monkeypatch.setattr(output_module, "replace_output", fail_second)
    with pytest.raises(OSError, match="injected commit failure"):
        with output_transaction([existing, new]) as staged:
            Path(staged[existing]).write_text("replacement")
            Path(staged[new]).write_text("new output")
    assert existing.read_text() == "original"
    assert not new.exists()
    assert set(tmp_path.iterdir()) == {existing}


def test_restore_failure_keeps_recovery_backups(monkeypatch, tmp_path):
    targets = [tmp_path / "first.txt", tmp_path / "second.txt"]
    for target in targets:
        target.write_text("original")
    original_replace = output_module.replace_output
    calls = 0

    def fail_commit_and_restore(source, target):
        nonlocal calls
        calls += 1
        if calls >= 2:
            raise OSError("injected failure")
        original_replace(source, target)

    monkeypatch.setattr(output_module, "replace_output", fail_commit_and_restore)
    with pytest.raises(RuntimeError, match="backup files were preserved"):
        with output_transaction(targets) as staged:
            for path in staged.values():
                Path(path).write_text("replacement")
    backups = list(tmp_path.glob(".*.backup.*"))
    assert len(backups) == 2
    assert all(path.read_text() == "original" for path in backups)
    assert not list(tmp_path.glob(".*.stage.*"))
    assert not list(tmp_path.glob("*.lock"))


def test_recovery_continues_when_one_target_cannot_be_restored(monkeypatch, tmp_path):
    first, second = tmp_path / "first", tmp_path / "second"
    for path in (first, second):
        path.write_text("original")
    original_replace = output_module.replace_output

    def fail_one_restore(source, target):
        if ".backup." in str(source) and str(target) == str(second):
            raise PermissionError("second target cannot be restored")
        original_replace(source, target)

    def fail_after_install():
        raise OSError("commit callback failed")

    monkeypatch.setattr(output_module, "replace_output", fail_one_restore)
    with pytest.raises(RuntimeError, match="backup files were preserved"):
        with output_transaction(
            [first, second], after_install=fail_after_install
        ) as staged:
            for path in staged.values():
                Path(path).write_text("new")
    assert first.read_text() == "original"
    assert second.read_text() == "new"
    backups = list(tmp_path.glob(".*.backup.*"))
    assert len(backups) == 1
    assert backups[0].read_text() == "original"


def test_staging_failure_does_not_change_either_target(tmp_path):
    paths = [tmp_path / "first", tmp_path / "second"]
    for path in paths:
        path.write_text("original")
    with pytest.raises(OSError, match="write failure"):
        with output_transaction(paths) as staged:
            Path(staged[paths[0]]).write_text("partial")
            raise OSError("write failure")
    assert all(path.read_text() == "original" for path in paths)
    assert set(tmp_path.iterdir()) == set(paths)


def test_all_destinations_are_validated_before_staging(tmp_path):
    original = tmp_path / "original"
    original.write_text("original")
    with pytest.raises(ValueError, match="regular file"):
        with output_transaction([original, tmp_path]):
            pytest.fail("invalid target was accepted")
    assert original.read_text() == "original"
    assert set(tmp_path.iterdir()) == {original}


@pytest.mark.skipif(os.name == "nt", reason="POSIX permission bits")
def test_transaction_preserves_modes_and_obeys_umask(tmp_path):
    existing, new = tmp_path / "existing", tmp_path / "new"
    existing.write_text("old")
    existing.chmod(0o640)
    previous = os.umask(0o027)
    try:
        with output_transaction([existing, new]) as staged:
            for path in staged.values():
                Path(path).write_text("replacement")
    finally:
        os.umask(previous)
    assert stat.S_IMODE(existing.stat().st_mode) == 0o640
    assert stat.S_IMODE(new.stat().st_mode) == 0o640


@pytest.mark.skipif(os.name == "nt", reason="symlink creation may need privileges")
def test_transaction_preserves_output_symlink_and_detects_staging_alias(tmp_path):
    target, alias, victim = [tmp_path / name for name in ("target", "alias", "victim")]
    target.write_text("old")
    victim.write_text("untouched")
    alias.symlink_to(target)
    with output_transaction([alias]) as staged:
        Path(staged[alias]).write_text("updated")
    assert alias.is_symlink()
    assert target.read_text() == "updated"
    with pytest.raises((OSError, RuntimeError)):
        with output_transaction([target]) as staged:
            path = Path(staged[target])
            path.unlink()
            path.symlink_to(victim)
    assert victim.read_text() == "untouched"
    assert target.read_text() == "updated"


def test_stdout_failure_recovers_files(tmp_path):
    target = tmp_path / "output"
    target.write_text("original")

    def broken_stdout():
        raise BrokenPipeError("stdout closed")

    with pytest.raises(BrokenPipeError):
        with output_transaction([target], after_install=broken_stdout) as staged:
            Path(staged[target]).write_text("new")
    assert target.read_text() == "original"


def test_staging_flush_has_write_access_for_windows(monkeypatch, tmp_path):
    original_fsync = os.fsync

    def writable_fsync(descriptor):
        # FlushFileBuffers requires write access; exercise that contract on
        # POSIX too, without changing the staged bytes.
        os.write(descriptor, b"")
        original_fsync(descriptor)

    monkeypatch.setattr(os, "fsync", writable_fsync)
    destination = tmp_path / "result"
    with output_transaction([destination]) as staged:
        Path(staged[destination]).write_text("complete")
    assert destination.read_text() == "complete"


def test_regression_rejects_a_replaced_stage_before_writing_through_it(
    monkeypatch, tmp_path
):
    import pandas as pd

    from nwkit.regression_pipeline import _write_dataframes_transactionally

    first, second, victim = (tmp_path / name for name in ("first", "second", "victim"))
    for path in (first, second, victim):
        path.write_text("original")
    original_stage = output_module._new_stage

    def replace_stage(target):
        stage = original_stage(target)
        if os.path.basename(target) == second.name:
            os.remove(stage.path)
            os.link(victim, stage.path)
        return stage

    monkeypatch.setattr(output_module, "_new_stage", replace_stage)
    with pytest.raises(RuntimeError, match="staging file was replaced"):
        _write_dataframes_transactionally(
            [
                (str(first), pd.DataFrame({"value": [1]})),
                (str(second), pd.DataFrame({"value": [2]})),
            ]
        )
    assert all(path.read_text() == "original" for path in (first, second, victim))
    assert set(tmp_path.iterdir()) == {first, second, victim}

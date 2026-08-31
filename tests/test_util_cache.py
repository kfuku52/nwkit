import os
import pickle
from contextlib import contextmanager

import pytest
import requests

from nwkit import util as util_mod
from nwkit.util import (
    _break_stale_lock_if_needed,
    _build_ete_taxonomy_database,
    _validate_ete_taxonomy_db,
    acquire_exclusive_lock,
    get_ete_ncbitaxa,
    resolve_download_dir,
    resolve_ete_data_dir,
    validate_distinct_output_paths,
)
from tests.helpers import make_args
from tests.util_test_support import (
    write_valid_taxdump as _write_valid_taxdump,
)
from tests.util_test_support import (
    write_valid_taxonomy_db as _write_valid_taxonomy_db,
)


class TestDownloadDirHelpers:
    def test_resolve_download_dir_returns_none_for_auto(self, tmp_path):
        outfile = tmp_path / "out" / "tree.nwk"
        args = make_args(outfile=str(outfile), download_dir="auto")
        resolved = resolve_download_dir(args)
        assert resolved is None

    def test_resolve_download_dir_uses_outfile_parent_for_inferred(self, tmp_path):
        outfile = tmp_path / "out" / "tree.nwk"
        args = make_args(outfile=str(outfile), download_dir="inferred")
        resolved = resolve_download_dir(args)
        assert resolved == os.path.join(os.path.realpath(outfile.parent), "downloads")

    def test_resolve_download_dir_uses_cwd_for_stdout_when_inferred(
        self, monkeypatch, tmp_path
    ):
        monkeypatch.chdir(tmp_path)
        args = make_args(outfile="-", download_dir="inferred")
        resolved = resolve_download_dir(args)
        assert resolved == os.path.join(os.path.realpath(tmp_path), "downloads")

    def test_resolve_download_dir_respects_explicit_path(self, tmp_path):
        explicit_dir = tmp_path / "shared-cache"
        args = make_args(outfile="-", download_dir=str(explicit_dir))
        resolved = resolve_download_dir(args)
        assert resolved == os.path.realpath(explicit_dir)

    def test_resolve_ete_data_dir_appends_ete4(self, tmp_path):
        explicit_dir = tmp_path / "shared-cache"
        args = make_args(download_dir=str(explicit_dir))
        assert resolve_ete_data_dir(args) == os.path.join(
            os.path.realpath(explicit_dir), "ete4"
        )

    def test_resolve_ete_data_dir_returns_none_for_auto(self):
        args = make_args(download_dir="auto")
        assert resolve_ete_data_dir(args) is None

    def test_acquire_exclusive_lock_creates_and_removes_file(self, tmp_path):
        lock_path = tmp_path / ".lock"
        with acquire_exclusive_lock(str(lock_path), timeout_seconds=1):
            assert lock_path.exists()
        assert not lock_path.exists()

    def test_acquire_exclusive_lock_breaks_stale_lock(self, tmp_path):
        lock_path = tmp_path / ".lock"
        lock_path.write_text("999999999\n")
        with acquire_exclusive_lock(str(lock_path), timeout_seconds=1):
            assert lock_path.exists()
        assert not lock_path.exists()

    def test_acquire_exclusive_lock_rejects_final_symlink_without_deleting_target(
        self,
        tmp_path,
    ):
        victim = tmp_path / "victim.txt"
        lock_path = tmp_path / ".lock"
        victim.write_text("999999999\n")
        lock_path.symlink_to(victim)

        with pytest.raises(IsADirectoryError, match="not a file"):
            with acquire_exclusive_lock(str(lock_path), timeout_seconds=1):
                pass

        assert victim.read_text() == "999999999\n"
        assert lock_path.is_symlink()

    def test_new_empty_lock_receives_metadata_grace_period(self, tmp_path):
        lock_path = tmp_path / ".lock"
        lock_path.write_bytes(b"")
        assert _break_stale_lock_if_needed(str(lock_path)) is False
        assert lock_path.exists()

    def test_windows_liveness_check_never_calls_terminating_os_kill(
        self,
        monkeypatch,
    ):
        monkeypatch.setattr(util_mod.os, "name", "nt")
        monkeypatch.setattr(
            util_mod,
            "_is_windows_process_alive",
            lambda pid: pid == 123,
        )

        def unexpected_kill(pid, signal):
            raise AssertionError("os.kill must not be used for Windows liveness checks")

        monkeypatch.setattr(util_mod.os, "kill", unexpected_kill)

        assert util_mod._is_process_alive(123) is True
        assert util_mod._is_process_alive(124) is False

    def test_taxdump_checksum_response_is_streamed_with_a_hard_limit(
        self,
        monkeypatch,
        tmp_path,
    ):
        class OversizedResponse:
            status_code = 200
            headers = {}

            def raise_for_status(self):
                return None

            def iter_content(self, chunk_size):
                yield b"x" * (util_mod.ETE_CHECKSUM_MAX_BYTES + 1)

            def close(self):
                return None

        class FakeSession:
            def __init__(self):
                self.calls = []

            def mount(self, *args, **kwargs):
                return None

            def get(self, url, **kwargs):
                self.calls.append((url, kwargs))
                return OversizedResponse()

            def close(self):
                return None

        session = FakeSession()
        monkeypatch.setattr(requests, "Session", lambda: session)

        with pytest.raises(ValueError, match="checksum response exceeds"):
            util_mod._download_ete_taxdump(str(tmp_path / "taxdump.tar.gz"))

        assert len(session.calls) == 1
        _, request_kwargs = session.calls[0]
        assert request_kwargs["allow_redirects"] is False
        assert request_kwargs["stream"] is True
        assert request_kwargs["headers"]["Accept"] == "*/*"

    def test_taxdump_md5_is_explicitly_nonsecurity_for_fips(
        self,
        monkeypatch,
        tmp_path,
    ):
        existing_payload = b"old taxonomy archive"
        downloaded_payload = b"new taxonomy archive"
        real_md5 = util_mod.hashlib.md5
        expected_md5 = real_md5(downloaded_payload).hexdigest()
        taxdump_file = tmp_path / "taxdump.tar.gz"
        taxdump_file.write_bytes(existing_payload)

        class Response:
            status_code = 200

            def __init__(self, payload):
                self.payload = payload
                self.headers = {"Content-Length": str(len(payload))}

            def raise_for_status(self):
                return None

            def iter_content(self, chunk_size):
                yield self.payload

            def close(self):
                return None

        class FakeSession:
            def mount(self, *args, **kwargs):
                return None

            def get(self, url, **kwargs):
                if url.endswith(".md5"):
                    return Response(
                        (expected_md5 + "  taxdump.tar.gz\n").encode("ascii")
                    )
                return Response(downloaded_payload)

            def close(self):
                return None

        md5_security_flags = []

        def fips_md5(*args, **kwargs):
            used_for_security = kwargs.get("usedforsecurity")
            md5_security_flags.append(used_for_security)
            if used_for_security is not False:
                raise ValueError("MD5 is unavailable for security use in FIPS mode")
            return real_md5(*args, **kwargs)

        monkeypatch.setattr(requests, "Session", FakeSession)
        monkeypatch.setattr(util_mod, "_validate_ete_taxdump", lambda path: True)
        monkeypatch.setattr(util_mod.hashlib, "md5", fips_md5)

        assert util_mod._download_ete_taxdump(str(taxdump_file)) is True
        assert taxdump_file.read_bytes() == downloaded_payload
        assert md5_security_flags == [False, False]

    def test_taxdump_redirect_cannot_leave_expected_https_host(self):
        class RedirectResponse:
            status_code = 302
            headers = {"Location": "https://attacker.example/taxdump.tar.gz"}

            def close(self):
                return None

        class FakeSession:
            def get(self, url, **kwargs):
                return RedirectResponse()

        with pytest.raises(ValueError, match="expected HTTPS host"):
            util_mod._open_ete_download_response(
                FakeSession(),
                util_mod.ETE_TAXDUMP_URL,
                stream=True,
                timeout=(1, 1),
            )

    def test_distinct_output_paths_rejects_realpath_aliases(self, tmp_path):
        output = tmp_path / "result.tsv"
        with pytest.raises(ValueError, match="Output paths must be distinct"):
            validate_distinct_output_paths(
                [
                    ("--outfile", str(output)),
                    ("--report", str(tmp_path / "." / "result.tsv")),
                ]
            )

    def test_distinct_missing_outputs_respect_case_insensitive_filesystems(
        self,
        monkeypatch,
        tmp_path,
    ):
        monkeypatch.setattr(
            "nwkit.file_paths._filesystem_is_case_insensitive",
            lambda path: True,
        )

        with pytest.raises(ValueError, match="Output paths must be distinct"):
            validate_distinct_output_paths(
                [
                    ("--outfile", str(tmp_path / "Result.tsv")),
                    ("--report", str(tmp_path / "result.tsv")),
                ]
            )

    def test_distinct_output_paths_rejects_hardlink_aliases(self, tmp_path):
        output = tmp_path / "result.tsv"
        hardlink = tmp_path / "result-hardlink.tsv"
        symlink = tmp_path / "result-symlink.tsv"
        output.write_text("existing")
        os.link(output, hardlink)
        symlink.symlink_to(hardlink)

        with pytest.raises(ValueError, match="Output paths must be distinct"):
            validate_distinct_output_paths(
                [
                    ("--outfile", str(output)),
                    ("--report", str(symlink)),
                ]
            )

    def test_get_ete_ncbitaxa_uses_download_dir_and_lock(self, monkeypatch, tmp_path):
        explicit_dir = tmp_path / "shared-cache"
        args = make_args(download_dir=str(explicit_dir))
        calls = dict()

        @contextmanager
        def fake_lock(
            lock_path, lock_label="Lock", poll_seconds=1, timeout_seconds=3600
        ):
            calls["lock_path"] = lock_path
            calls["lock_label"] = lock_label
            yield

        def fake_download_ete_taxdump(taxdump_file):
            calls["downloaded_taxdump_file"] = taxdump_file
            _write_valid_taxdump(tmp_path / "shared-cache" / "ete4" / "taxdump.tar.gz")
            return True

        def fake_build(dbfile, taxdump_file):
            calls["built_dbfile"] = dbfile
            calls["built_taxdump_file"] = taxdump_file
            _write_valid_taxonomy_db(dbfile)

        class FakeNCBI:
            def __init__(
                self, dbfile=None, taxdump_file=None, memory=False, update=True
            ):
                calls["dbfile"] = dbfile
                calls["taxdump_file"] = taxdump_file
                calls["taxdump_exists_at_init"] = (
                    taxdump_file is not None and os.path.exists(taxdump_file)
                )
                calls["memory"] = memory
                calls["update"] = update
                self.db = None

        monkeypatch.setattr("nwkit.util.acquire_exclusive_lock", fake_lock)
        monkeypatch.setattr(
            "nwkit.util._find_existing_ete_taxonomy_assets",
            lambda exclude_dbfile=None: None,
        )
        monkeypatch.setattr(
            "nwkit.util._download_ete_taxdump", fake_download_ete_taxdump
        )
        monkeypatch.setattr("nwkit.util._build_ete_taxonomy_database", fake_build)
        monkeypatch.setattr("nwkit.util.ete4.NCBITaxa", FakeNCBI)
        ncbi = get_ete_ncbitaxa(args=args)
        assert isinstance(ncbi, FakeNCBI)
        expected_ete_dir = os.path.join(os.path.realpath(explicit_dir), "ete4")
        assert calls["lock_path"] == os.path.join(
            expected_ete_dir, ".ete4_taxonomy.lock"
        )
        assert calls["lock_label"] == "ETE4 taxonomy DB"
        assert calls["dbfile"] == os.path.join(expected_ete_dir, "taxa.sqlite")
        assert calls["taxdump_file"] is None
        assert calls["downloaded_taxdump_file"] == os.path.join(
            expected_ete_dir, "taxdump.tar.gz"
        )
        assert calls["built_dbfile"] == os.path.join(expected_ete_dir, "taxa.sqlite")
        assert calls["built_taxdump_file"] == os.path.join(
            expected_ete_dir, "taxdump.tar.gz"
        )

    def test_taxonomy_validation_does_not_execute_pickle_globals(self, tmp_path):
        dbfile = tmp_path / "taxa.sqlite"
        marker = tmp_path / "executed"
        _write_valid_taxonomy_db(dbfile)

        class Payload:
            def __reduce__(self):
                return os.system, ("touch {}".format(marker),)

        with open(str(dbfile) + ".traverse.pkl", "wb") as handle:
            pickle.dump(Payload(), handle)

        assert _validate_ete_taxonomy_db(str(dbfile), require_traverse=True) is False
        assert not marker.exists()

    def test_taxonomy_build_isolated_from_calling_working_directory(
        self, monkeypatch, tmp_path
    ):
        caller_dir = tmp_path / "caller"
        target_dir = tmp_path / "cache"
        caller_dir.mkdir()
        target_dir.mkdir()
        sentinel = caller_dir / "taxa.tab"
        sentinel.write_text("keep me")
        taxdump = target_dir / "taxdump.tar.gz"
        taxdump.write_bytes(b"placeholder")
        observed = {}

        def fake_run(command, cwd, check):
            observed["cwd"] = cwd
            observed["command"] = command
            dbfile = command[-2]
            with open(dbfile, "wb") as handle:
                handle.write(b"database")
            with open(dbfile + ".traverse.pkl", "wb") as handle:
                pickle.dump([1], handle)

        monkeypatch.chdir(caller_dir)
        monkeypatch.setattr("nwkit.util.subprocess.run", fake_run)
        monkeypatch.setattr(
            "nwkit.util._validate_ete_taxonomy_db", lambda *args, **kwargs: True
        )
        _build_ete_taxonomy_database(str(target_dir / "taxa.sqlite"), str(taxdump))
        assert os.path.realpath(observed["cwd"]) != os.path.realpath(caller_dir)
        assert sentinel.read_text() == "keep me"
        assert (target_dir / "taxa.sqlite").read_bytes() == b"database"

    def test_get_ete_ncbitaxa_auto_uses_ete_default_location(
        self, monkeypatch, tmp_path
    ):
        calls = dict()

        default_db = tmp_path / "ete-default" / "taxa.sqlite"
        default_db.parent.mkdir()
        _write_valid_taxonomy_db(default_db)

        @contextmanager
        def fake_lock(lock_path, **kwargs):
            calls["lock_path"] = lock_path
            yield

        class FakeNCBI:
            def __init__(self, *args, **kwargs):
                calls["args"] = args
                calls["kwargs"] = kwargs
                self.db = None

        monkeypatch.setattr(
            "ete4.ncbi_taxonomy.ncbiquery.DEFAULT_TAXADB", str(default_db)
        )
        monkeypatch.setattr("nwkit.util.acquire_exclusive_lock", fake_lock)
        monkeypatch.setattr("nwkit.util.ete4.NCBITaxa", FakeNCBI)
        ncbi = get_ete_ncbitaxa(args=make_args(download_dir="auto"))
        assert isinstance(ncbi, FakeNCBI)
        assert calls["args"] == ()
        assert calls["kwargs"] == {"dbfile": str(default_db), "update": False}
        assert calls["lock_path"] == str(default_db.parent / ".ete4_taxonomy.lock")

    def test_get_ete_ncbitaxa_checks_existing_taxdump_before_build(
        self, monkeypatch, tmp_path
    ):
        explicit_dir = tmp_path / "shared-cache"
        expected_ete_dir = os.path.join(os.path.realpath(explicit_dir), "ete4")
        os.makedirs(expected_ete_dir, exist_ok=True)
        _write_valid_taxdump(tmp_path / "shared-cache" / "ete4" / "taxdump.tar.gz")
        args = make_args(download_dir=str(explicit_dir))
        calls = dict(download_count=0)

        @contextmanager
        def fake_lock(
            lock_path, lock_label="Lock", poll_seconds=1, timeout_seconds=3600
        ):
            calls["lock_path"] = lock_path
            yield

        def fake_download_ete_taxdump(taxdump_file):
            calls["download_count"] += 1
            return False

        def fake_build(dbfile, taxdump_file):
            calls["build_count"] = calls.get("build_count", 0) + 1
            _write_valid_taxonomy_db(dbfile)

        class FakeNCBI:
            def __init__(
                self, dbfile=None, taxdump_file=None, memory=False, update=True
            ):
                calls["dbfile"] = dbfile
                calls["taxdump_file"] = taxdump_file
                self.db = None

        monkeypatch.setattr("nwkit.util.acquire_exclusive_lock", fake_lock)
        monkeypatch.setattr(
            "nwkit.util._find_existing_ete_taxonomy_assets",
            lambda exclude_dbfile=None: None,
        )
        monkeypatch.setattr(
            "nwkit.util._download_ete_taxdump", fake_download_ete_taxdump
        )
        monkeypatch.setattr("nwkit.util._build_ete_taxonomy_database", fake_build)
        monkeypatch.setattr("nwkit.util.ete4.NCBITaxa", FakeNCBI)
        ncbi = get_ete_ncbitaxa(args=args)
        assert isinstance(ncbi, FakeNCBI)
        assert calls["lock_path"] == os.path.join(
            expected_ete_dir, ".ete4_taxonomy.lock"
        )
        assert calls["dbfile"] == os.path.join(expected_ete_dir, "taxa.sqlite")
        assert calls["taxdump_file"] is None
        assert calls["download_count"] == 1
        assert calls["build_count"] == 1

    def test_get_ete_ncbitaxa_reuses_existing_db_without_taxdump_update(
        self, monkeypatch, tmp_path
    ):
        explicit_dir = tmp_path / "shared-cache"
        expected_ete_dir = os.path.join(os.path.realpath(explicit_dir), "ete4")
        os.makedirs(expected_ete_dir, exist_ok=True)
        existing_dbfile = os.path.join(expected_ete_dir, "taxa.sqlite")
        _write_valid_taxonomy_db(existing_dbfile)
        args = make_args(download_dir=str(explicit_dir))
        calls = dict(download_count=0)

        @contextmanager
        def fake_lock(
            lock_path, lock_label="Lock", poll_seconds=1, timeout_seconds=3600
        ):
            calls["lock_path"] = lock_path
            yield

        def fake_download_ete_taxdump(taxdump_file):
            calls["download_count"] += 1

        class FakeNCBI:
            def __init__(
                self, dbfile=None, taxdump_file=None, memory=False, update=True
            ):
                calls["dbfile"] = dbfile
                calls["taxdump_file"] = taxdump_file
                calls["update"] = update
                self.db = None

        monkeypatch.setattr("nwkit.util.acquire_exclusive_lock", fake_lock)
        monkeypatch.setattr(
            "nwkit.util._find_existing_ete_taxonomy_assets",
            lambda exclude_dbfile=None: None,
        )
        monkeypatch.setattr(
            "nwkit.util._download_ete_taxdump", fake_download_ete_taxdump
        )
        monkeypatch.setattr("nwkit.util.ete4.NCBITaxa", FakeNCBI)
        ncbi = get_ete_ncbitaxa(args=args)
        assert isinstance(ncbi, FakeNCBI)
        assert calls["lock_path"] == os.path.join(
            expected_ete_dir, ".ete4_taxonomy.lock"
        )
        assert calls["dbfile"] == existing_dbfile
        assert calls["taxdump_file"] is None
        assert calls["update"] is False
        assert calls["download_count"] == 0

    def test_get_ete_ncbitaxa_seeds_local_cache_from_existing_user_db(
        self, monkeypatch, tmp_path
    ):
        explicit_dir = tmp_path / "shared-cache"
        expected_ete_dir = os.path.join(os.path.realpath(explicit_dir), "ete4")
        source_dir = tmp_path / "user-cache"
        source_dir.mkdir()
        source_dbfile = source_dir / "taxa.sqlite"
        source_traverse = source_dir / "taxa.sqlite.traverse.pkl"
        source_taxdump = source_dir / "taxdump.tar.gz"
        _write_valid_taxonomy_db(source_dbfile)
        _write_valid_taxdump(source_taxdump)
        args = make_args(download_dir=str(explicit_dir))
        calls = dict(download_count=0)

        @contextmanager
        def fake_lock(
            lock_path, lock_label="Lock", poll_seconds=1, timeout_seconds=3600
        ):
            calls["lock_path"] = lock_path
            yield

        def fake_download_ete_taxdump(taxdump_file):
            calls["download_count"] += 1

        class FakeNCBI:
            def __init__(
                self, dbfile=None, taxdump_file=None, memory=False, update=True
            ):
                calls["dbfile"] = dbfile
                calls["taxdump_file"] = taxdump_file
                calls["update"] = update
                self.db = None

        monkeypatch.setattr("nwkit.util.acquire_exclusive_lock", fake_lock)
        monkeypatch.setattr(
            "nwkit.util._find_existing_ete_taxonomy_assets",
            lambda exclude_dbfile=None: {
                "dbfile": str(source_dbfile),
                "traverse_file": str(source_traverse),
                "taxdump_file": str(source_taxdump),
            },
        )
        monkeypatch.setattr(
            "nwkit.util._download_ete_taxdump", fake_download_ete_taxdump
        )
        monkeypatch.setattr("nwkit.util.ete4.NCBITaxa", FakeNCBI)
        ncbi = get_ete_ncbitaxa(args=args)
        assert isinstance(ncbi, FakeNCBI)
        assert calls["lock_path"] == os.path.join(
            expected_ete_dir, ".ete4_taxonomy.lock"
        )
        assert calls["dbfile"] == os.path.join(expected_ete_dir, "taxa.sqlite")
        assert calls["taxdump_file"] is None
        assert calls["update"] is False
        assert calls["download_count"] == 0
        assert os.path.getsize(os.path.join(expected_ete_dir, "taxa.sqlite")) > 0
        assert (
            os.path.getsize(os.path.join(expected_ete_dir, "taxa.sqlite.traverse.pkl"))
            > 0
        )
        assert os.path.getsize(os.path.join(expected_ete_dir, "taxdump.tar.gz")) > 0

    def test_get_ete_ncbitaxa_rebuilds_corrupt_database(self, monkeypatch, tmp_path):
        explicit_dir = tmp_path / "shared-cache"
        ete_dir = explicit_dir / "ete4"
        ete_dir.mkdir(parents=True)
        dbfile = ete_dir / "taxa.sqlite"
        dbfile.write_bytes(b"not a sqlite database")
        calls = {"download": 0, "build": 0}

        def fake_download(path):
            calls["download"] += 1
            _write_valid_taxdump(ete_dir / "taxdump.tar.gz")
            return True

        def fake_build(dbfile, taxdump_file):
            calls["build"] += 1
            os.remove(dbfile)
            _write_valid_taxonomy_db(dbfile)

        monkeypatch.setattr(
            "nwkit.util._find_existing_ete_taxonomy_assets",
            lambda exclude_dbfile=None: None,
        )
        monkeypatch.setattr("nwkit.util._download_ete_taxdump", fake_download)
        monkeypatch.setattr("nwkit.util._build_ete_taxonomy_database", fake_build)
        ncbi = get_ete_ncbitaxa(args=make_args(download_dir=str(explicit_dir)))
        try:
            assert calls == {"download": 1, "build": 1}
            assert ncbi.get_taxid_translator([1])[1] == "root"
        finally:
            ncbi.db.close()

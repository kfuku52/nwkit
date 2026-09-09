"""External library discovery and explicit setup never build during analysis."""

import json
import sys
from pathlib import Path

import pytest

from nwkit import iqtree_library


def capabilities():
    return dict(
        protocol=1,
        engine="iqtree3",
        library_version="3.1.4",
        library_sha256="1" * 64,
        adapter_sha256="2" * 64,
    )


def fake_worker(tmp_path, data):
    script = tmp_path / "worker"
    script.write_text(f"#!{sys.executable}\nprint({json.dumps(data)!r})\n")
    script.chmod(0o755)
    return str(script)


def test_auto_absent_is_cli_but_explicit_library_fails(monkeypatch):
    monkeypatch.delenv("NWKIT_IQTREE_WORKER", raising=False)
    monkeypatch.setattr(iqtree_library.shutil, "which", lambda _: None)
    assert iqtree_library.select_worker() is None
    with pytest.raises(ValueError, match="python -m nwkit.iqtree_library build"):
        iqtree_library.select_worker("library")
    with pytest.raises(ValueError, match="not found"):
        iqtree_library.select_worker("auto", "missing-worker")


def test_explicit_cli_does_not_probe_or_build(monkeypatch):
    monkeypatch.setenv("NWKIT_IQTREE_WORKER", "/missing/worker")
    monkeypatch.setattr(
        iqtree_library.subprocess,
        "run",
        lambda *a, **k: pytest.fail("must not execute"),
    )
    assert iqtree_library.select_worker("cli") is None
    with pytest.raises(ValueError, match="cannot be combined"):
        iqtree_library.select_worker("cli", "/missing/worker")


def test_external_worker_is_discovered_and_identified(tmp_path, monkeypatch):
    executable = fake_worker(tmp_path, capabilities())
    monkeypatch.setenv("NWKIT_IQTREE_WORKER", executable)
    found = iqtree_library.select_worker()
    assert found["executable"] == executable
    assert found["worker_sha256"] == iqtree_library.file_hash(executable)
    assert found["library_sha256"] == "1" * 64


@pytest.mark.parametrize(
    "key,value",
    [
        ("protocol", 2),
        ("engine", "iqtree2"),
        ("library_version", "2.4.0"),
        ("library_version", "unknown"),
        ("library_sha256", "missing"),
    ],
)
def test_broken_worker_never_falls_back_to_cli(tmp_path, monkeypatch, key, value):
    data = capabilities()
    data[key] = value
    monkeypatch.setenv("NWKIT_IQTREE_WORKER", fake_worker(tmp_path, data))
    with pytest.raises(ValueError, match="Cannot use IQ-TREE library worker"):
        iqtree_library.select_worker("auto")


def test_no_arbitrary_upper_version_bound(tmp_path):
    data = capabilities()
    data["library_version"] = "12.0.0"
    assert (
        iqtree_library.find_worker(fake_worker(tmp_path, data))["library_version"]
        == "12.0.0"
    )


def test_library_setup_requires_a_library_not_an_executable(tmp_path):
    (tmp_path / "iqtree3").write_bytes(b"not a library")
    with pytest.raises(ValueError, match="No libiqtree.a"):
        iqtree_library.build_worker(tmp_path, tmp_path / "install")


def test_adapter_build_reuses_matching_compile_settings(tmp_path, monkeypatch):
    upstream = tmp_path / "official source"
    (upstream / "tree").mkdir(parents=True)
    (upstream / "tree/phylotree.h").touch()
    build = tmp_path / "build"
    build.mkdir()
    (build / "CMakeCache.txt").write_text(
        "BUILD_LIB:BOOL=ON\nCMAKE_BUILD_TYPE:STRING=Release\n"
        f"CMAKE_HOME_DIRECTORY:INTERNAL={upstream}\n"
        "CMAKE_EXE_LINKER_FLAGS:STRING=-pthread\n"
    )
    original = str(upstream / "main/main.cpp")
    (build / "compile_commands.json").write_text(
        json.dumps(
            [
                dict(
                    directory=str(build),
                    file=original,
                    arguments=[
                        "c++",
                        "-DBUILD_LIB",
                        "-I",
                        str(upstream),
                        "-fPIC",
                        "-fopenmp",
                        "-c",
                        original,
                        "-o",
                        "main.o",
                    ],
                )
            ]
        )
    )
    monkeypatch.setenv("LDFLAGS", "-L'library path'")
    source, output, library = (
        tmp_path / name for name in ["worker.cpp", "worker", "libiqtree.a"]
    )
    command, cwd, found = iqtree_library._compile_command(
        build, source, output, library, {}
    )
    assert found == upstream and cwd == str(build)
    assert "-DBUILD_LIB" in command and "-fopenmp" in command and "-fPIC" in command
    assert original not in command and "main.o" not in command and "-c" not in command
    assert str(source) in command and str(library) in command
    assert "-Llibrary path" in command
    assert command[command.index("-o") + 1] == str(output)


def test_check_auto_reports_absence_without_installing(monkeypatch, capsys):
    monkeypatch.delenv("NWKIT_IQTREE_WORKER", raising=False)
    monkeypatch.setattr(iqtree_library.shutil, "which", lambda _: None)
    iqtree_library.main(["check", "--interface", "auto"])
    assert json.loads(capsys.readouterr().out) is None


def test_worker_failure_is_a_setup_error(tmp_path):
    script = tmp_path / "worker"
    script.write_text(f"#!{sys.executable}\nraise SystemExit(2)\n")
    script.chmod(0o755)
    with pytest.raises(ValueError, match="Cannot use"):
        iqtree_library.worker_capabilities(script)


def test_source_distribution_contains_no_iqtree_implementation():
    source = Path(iqtree_library.__file__).with_name("data_iqtree")
    assert [p.name for p in source.iterdir()] == ["worker.cpp"]
    assert "SPDX-License-Identifier: MIT" in (source / "worker.cpp").read_text()

"""Discover or explicitly build an external worker from a user's IQ-TREE library.

This module never downloads IQ-TREE. Analysis calls only discover installed
workers; compiling requires the separate ``build`` command below.
"""

import argparse
import hashlib
import json
import os
import re
import shlex
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

WORKER = "nwkit-iqtree-worker"
BUILD_HELP = (
    "Build and install the optional worker with "
    "python -m nwkit.iqtree_library build --build-dir IQTREE_BUILD --prefix PREFIX. "
    "See IQTREE_LIBRARY.md or the IQ-TREE-library wiki page."
)


def file_hash(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def worker_capabilities(executable):
    try:
        result = subprocess.run(
            [str(executable), "--capabilities"],
            capture_output=True,
            text=True,
            timeout=30,
            check=True,
        )
        data = json.loads(result.stdout)
        version = re.match(r"(\d+)\.", data.get("library_version", ""))
        if (
            data.get("protocol") != 1
            or data.get("engine") != "iqtree3"
            or version is None
            or int(version.group(1)) < 3
            or any(
                not re.fullmatch(r"[0-9a-f]{64}", data.get(key, ""))
                for key in ("library_sha256", "adapter_sha256")
            )
        ):
            raise ValueError("Unsupported worker protocol or library metadata.")
    except (
        OSError,
        subprocess.SubprocessError,
        ValueError,
        AttributeError,
        TypeError,
    ) as exc:
        raise ValueError(
            f"Cannot use IQ-TREE library worker {executable}: {exc}. {BUILD_HELP}"
        ) from exc
    return data


def find_worker(executable=None, *, required=False):
    requested = executable or os.environ.get("NWKIT_IQTREE_WORKER")
    resolved = shutil.which(requested or WORKER)
    if resolved is None:
        if requested or required:
            raise ValueError(
                f"IQ-TREE library worker not found: {requested or WORKER}. {BUILD_HELP}"
            )
        return None
    # An installed but broken worker is an error, not a reason to change engines.
    metadata = worker_capabilities(resolved)
    return dict(executable=resolved, **metadata, worker_sha256=file_hash(resolved))


def select_worker(interface="auto", executable=None):
    if interface not in {"auto", "cli", "library"}:
        raise ValueError("Unknown IQ-TREE interface: " + interface)
    if interface == "cli":
        if executable is not None:
            raise ValueError(
                "--iqtree-worker cannot be combined with --iqtree-interface cli."
            )
        return None
    return find_worker(executable, required=interface == "library")


def _cache(build):
    values = {}
    for line in (build / "CMakeCache.txt").read_text().splitlines():
        if not line.startswith(("#", "//")) and ":" in line and "=" in line:
            key, value = line.split("=", 1)
            values[key.split(":", 1)[0]] = value
    return values


def _compile_command(build, source, output, library, metadata):
    cache = _cache(build)
    if cache.get("BUILD_LIB", "").upper() not in {"ON", "TRUE", "1"}:
        raise ValueError("IQ-TREE must first be built with -DBUILD_LIB=ON.")
    entries = json.loads((build / "compile_commands.json").read_text())
    upstream = Path(cache["CMAKE_HOME_DIRECTORY"])
    entry = next(
        (
            item
            for item in entries
            if Path(item["file"]).resolve() == (upstream / "main/main.cpp").resolve()
        ),
        None,
    )
    if entry is None:
        raise ValueError("IQ-TREE compile_commands.json lacks main/main.cpp.")
    arguments = entry.get("arguments") or shlex.split(entry["command"])
    command = []
    skip = False
    for argument in arguments:
        if skip:
            skip = False
        elif argument in {"-o", "-MF", "-MT", "-MQ"}:
            skip = True
        elif argument not in {"-c", "-MD", "-MMD", entry["file"]}:
            command.append(argument)
    if not command or not (upstream / "tree/phylotree.h").is_file():
        raise ValueError(
            "Keep the IQ-TREE source/headers at the location recorded by CMake."
        )
    for key, value in metadata.items():
        command.append(f'-D{key}="{value}"')
    command += [str(source), str(library), "-o", str(output)]
    command += shlex.split(cache.get("CMAKE_EXE_LINKER_FLAGS", ""))
    command += shlex.split(
        cache.get(
            "CMAKE_EXE_LINKER_FLAGS_" + cache.get("CMAKE_BUILD_TYPE", "").upper(), ""
        )
    )
    # IQ-TREE folds its bundled zlib into libiqtree.a, but a system zlib
    # selected by CMake remains an external dependency of the static archive.
    # Use the same configuration-specific library and put it after the archive.
    if not re.search(r"nozlib|static", cache.get("IQTREE_FLAGS", "")):
        configurations = (
            ("DEBUG", "RELEASE")
            if cache.get("CMAKE_BUILD_TYPE", "").upper() == "DEBUG"
            else ("RELEASE", "DEBUG")
        )
        for configuration in configurations:
            zlib = cache.get("ZLIB_LIBRARY_" + configuration, "")
            if zlib and not zlib.endswith("-NOTFOUND"):
                command.extend(zlib.split(";"))
                break
    command += shlex.split(os.environ.get("LDFLAGS", ""))
    if sys.platform.startswith("linux"):
        command += ["-ldl", "-lm"]
    return command, entry["directory"], upstream


def build_worker(build_dir, prefix):
    if os.name == "nt":
        raise ValueError(
            "The external IQ-TREE worker build currently supports Linux and macOS."
        )
    build, prefix = Path(build_dir).resolve(), Path(prefix).resolve()
    library = build / "libiqtree.a"
    if not library.is_file():
        raise ValueError(
            "No libiqtree.a in the IQ-TREE library build directory. " + BUILD_HELP
        )
    source = Path(__file__).with_name("data_iqtree") / "worker.cpp"
    metadata = {
        "NWKIT_IQTREE_LIBRARY_SHA256": file_hash(library),
        "NWKIT_IQTREE_ADAPTER_SHA256": file_hash(source),
        "NWKIT_IQTREE_SOURCE_REVISION": "unavailable",
    }
    destination = prefix / "bin"
    destination.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(
        prefix=".nwkit-iqtree-build-", dir=destination
    ) as temporary:
        output = Path(temporary) / WORKER
        command, cwd, upstream = _compile_command(
            build, source, output, library, metadata
        )
        revision = (
            subprocess.run(
                ["git", "-C", str(upstream), "rev-parse", "HEAD"],
                capture_output=True,
                text=True,
                check=False,
            )
            if shutil.which("git")
            else None
        )
        if revision is not None and revision.returncode == 0:
            sha = revision.stdout.strip()
            if re.fullmatch(r"[0-9a-f]{40,64}", sha):
                metadata["NWKIT_IQTREE_SOURCE_REVISION"] = sha
                command, cwd, _ = _compile_command(
                    build, source, output, library, metadata
                )
        subprocess.run(command, cwd=cwd, check=True)
        capabilities = worker_capabilities(output)
        os.replace(output, destination / WORKER)
    return dict(executable=str(destination / WORKER), **capabilities)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    build = commands.add_parser(
        "build",
        help="Compile only NWKIT's adapter against an existing IQ-TREE library.",
    )
    build.add_argument("--build-dir", required=True, type=Path)
    build.add_argument("--prefix", required=True, type=Path)
    check = commands.add_parser(
        "check", help="Inspect an installed external worker without compiling."
    )
    check.add_argument("--worker")
    check.add_argument(
        "--interface", choices=["auto", "cli", "library"], default="library"
    )
    args = parser.parse_args(argv)
    try:
        result = (
            build_worker(args.build_dir, args.prefix)
            if args.command == "build"
            else select_worker(args.interface, args.worker)
        )
    except (ValueError, OSError, KeyError, subprocess.SubprocessError) as exc:
        parser.exit(1, f"IQ-TREE library setup failed: {exc}\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()

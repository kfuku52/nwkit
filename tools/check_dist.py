"""Validate reproducibility and contents of built NWKIT distributions."""

import subprocess
import sys
import tarfile
import zipfile
from pathlib import Path

from nwkit import __version__

PROJECT_ROOT = Path(__file__).resolve().parents[1]
WHEEL_NAME = f"nwkit-{__version__}-py3-none-any.whl"


def _wheel_members(path: Path) -> set[str]:
    with zipfile.ZipFile(path) as wheel:
        return set(wheel.namelist())


def main() -> int:
    direct_wheel = PROJECT_ROOT / "direct-dist" / WHEEL_NAME
    wheel = PROJECT_ROOT / "dist" / WHEEL_NAME
    sdist = PROJECT_ROOT / "dist" / f"nwkit-{__version__}.tar.gz"
    for artifact in (direct_wheel, wheel, sdist):
        if not artifact.is_file():
            raise FileNotFoundError(f"Expected distribution was not built: {artifact}")
    if direct_wheel.read_bytes() != wheel.read_bytes():
        raise RuntimeError("Direct and sdist-built wheels are not byte-for-byte equal.")

    required_wheel = {
        "nwkit/__init__.py",
        "nwkit/image_metadata.py",
        "nwkit/root.py",
        "nwkit/data_tree/apgiv.nwk",
    }
    forbidden_names = {"nwkit/_mad.py", "THIRD_PARTY_NOTICES"}
    for candidate in (direct_wheel, wheel):
        members = _wheel_members(candidate)
        missing = required_wheel - members
        if missing:
            raise RuntimeError(f"Wheel is missing required members: {sorted(missing)}")
        if any(
            member in forbidden_names or member.endswith("/THIRD_PARTY_NOTICES")
            for member in members
        ):
            raise RuntimeError(f"Wheel contains a forbidden member: {candidate}")

    with tarfile.open(sdist, "r:gz") as archive:
        sdist_members = set(archive.getnames())
    for suffix in (
        "/nwkit/__init__.py",
        "/nwkit/image_metadata.py",
        "/nwkit/root.py",
        "/nwkit/data_tree/apgiv.nwk",
        "/constraints-dev.txt",
        "/setup.py",
        "/tests/test_properties.py",
        "/tools/check_dist.py",
    ):
        if not any(member.endswith(suffix) for member in sdist_members):
            raise RuntimeError(f"Source distribution is missing {suffix}.")
    if any(
        member.endswith(("/nwkit/_mad.py", "/THIRD_PARTY_NOTICES"))
        for member in sdist_members
    ):
        raise RuntimeError("Source distribution contains a forbidden member.")

    subprocess.run(
        [sys.executable, "-m", "twine", "check", str(wheel), str(sdist)],
        check=True,
    )
    subprocess.run(["check-wheel-contents", str(wheel)], check=True)
    print("Distribution reproducibility and contents are valid.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

"""Build the checksum-verified ETE4 wheel with its Windows module-path fix.

Remove this workaround when the upstream package builds Windows extension
module names correctly (or publishes a compatible Windows wheel). This is a
CI build artifact, not an upper bound on NWKIT's runtime ETE dependency.
"""

import argparse
import hashlib
import subprocess
import sys
import tarfile
import tempfile
import urllib.request
from pathlib import Path

SOURCE_URL = "https://files.pythonhosted.org/packages/fa/53/0ceedd63e4c872f668188a5de786a00b6b8854cb285a3dc69cd107c62038/ete4-4.4.0.tar.gz"
SOURCE_SHA256 = "0306ab66c2a1685946f94a4912718a5812dfc3059070300fc25d0e2b881bb8d1"


def patch_setup(source):
    old = "name = path.replace('/', '.')[:-len('.pyx')]"
    new = r"name = path.replace('\\', '/').replace('/', '.')[:-len('.pyx')]"
    if source.count(old) != 1:
        raise ValueError(
            "ETE4 setup.py changed; revalidate or remove the Windows workaround."
        )
    return source.replace(old, new)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("wheel_dir", type=Path)
    args = parser.parse_args(argv)
    destination = args.wheel_dir.resolve()
    destination.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="nwkit-ete-windows-") as directory:
        root = Path(directory)
        archive = root / "ete4.tar.gz"
        urllib.request.urlretrieve(SOURCE_URL, archive)
        if hashlib.sha256(archive.read_bytes()).hexdigest() != SOURCE_SHA256:
            raise ValueError("ETE4 source checksum mismatch.")
        with tarfile.open(archive) as source_archive:
            source_archive.extractall(root, filter="data")
        source = root / "ete4-4.4.0"
        setup = source / "setup.py"
        setup.write_text(
            patch_setup(setup.read_text(encoding="utf-8")), encoding="utf-8"
        )
        subprocess.run(
            [
                sys.executable,
                "-m",
                "pip",
                "wheel",
                "--no-deps",
                str(source),
                "--wheel-dir",
                str(destination),
            ],
            check=True,
        )


if __name__ == "__main__":
    main()

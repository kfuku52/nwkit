import os
import shutil
import stat
from pathlib import Path

from setuptools import setup
from setuptools.command.bdist_wheel import bdist_wheel as _bdist_wheel
from setuptools.command.build_py import build_py as _build_py

PROJECT_ROOT = Path(__file__).resolve().parent
SOURCE_PACKAGE_DIR = (PROJECT_ROOT / "nwkit").resolve()
PROJECT_BUILD_DIR = (PROJECT_ROOT / "build").resolve()


class CleanBuildPy(_build_py):
    """Rebuild the package directory without retaining deleted modules."""

    def run(self):
        build_root = Path(self.build_lib).resolve()
        package_dir = build_root / "nwkit"
        try:
            resolved_package_dir = package_dir.resolve()
        except RuntimeError as exc:
            raise RuntimeError(
                "Refusing to clean an unresolved build package path: '{}'.".format(
                    package_dir
                )
            ) from exc
        if (
            resolved_package_dir.parent != build_root
            or resolved_package_dir.name != "nwkit"
        ):
            raise RuntimeError(
                "Refusing to clean a package path outside build_lib: '{}'.".format(
                    resolved_package_dir
                )
            )
        is_inside_project = (
            resolved_package_dir == PROJECT_ROOT
            or PROJECT_ROOT in resolved_package_dir.parents
        )
        is_inside_project_build = (
            resolved_package_dir == PROJECT_BUILD_DIR
            or PROJECT_BUILD_DIR in resolved_package_dir.parents
        )
        if resolved_package_dir == SOURCE_PACKAGE_DIR or (
            is_inside_project and not is_inside_project_build
        ):
            raise RuntimeError(
                "Refusing to clean the source tree as build output: '{}'.".format(
                    resolved_package_dir
                )
            )
        try:
            package_stat = os.lstat(package_dir)
        except FileNotFoundError:
            package_stat = None
        if package_stat is not None:
            if stat.S_ISLNK(package_stat.st_mode):
                raise RuntimeError(
                    "Refusing to clean a symbolic-link build package: '{}'.".format(
                        package_dir
                    )
                )
            if not stat.S_ISDIR(package_stat.st_mode):
                raise RuntimeError(
                    "Build package path is not a directory: '{}'.".format(package_dir)
                )
            shutil.rmtree(package_dir)
        super().run()


class ReproducibleBdistWheel(_bdist_wheel):
    """Match the sdist file modes in the staged package and wheel metadata."""

    def write_wheelfile(self, wheelfile_base, *args, **kwargs):
        super().write_wheelfile(wheelfile_base, *args, **kwargs)
        dist_info = Path(wheelfile_base)
        # Setuptools copies the modes of package data and license files.
        # Normalize the staged files before their permissions enter the ZIP.
        for directory in (dist_info.parent / "nwkit", dist_info):
            for path in directory.rglob("*"):
                if not path.is_symlink() and path.is_file():
                    path.chmod(0o644)


setup(cmdclass={"build_py": CleanBuildPy, "bdist_wheel": ReproducibleBdistWheel})

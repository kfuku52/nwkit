import ast
import re
from pathlib import Path

from setuptools import find_packages, setup


BASE_DIR = Path(__file__).resolve().parent
INIT_PATH = BASE_DIR / "nwkit" / "__init__.py"
README_PATH = BASE_DIR / "README.md"

match = re.search(r"__version__\s+=\s+(.*)", INIT_PATH.read_text(encoding="utf-8"))
if not match:
    raise RuntimeError("Failed to parse __version__ from nwkit/__init__.py")
version = str(ast.literal_eval(match.group(1)))


setup(
    name="nwkit",
    version=version,
    description="Tools for processing Newick trees",
    long_description=README_PATH.read_text(encoding="utf-8"),
    long_description_content_type="text/markdown",
    license="BSD-3-Clause",
    author="Kenji Fukushima",
    author_email="kfuku52@gmail.com",
    url="https://github.com/kfuku52/nwkit",
    keywords=["phylogenetics", "newick", "tree"],
    packages=find_packages(),
    python_requires=">=3.9",
    install_requires=["ete4>=4.4.0", "biopython", "pandas", "requests", "numpy"],
    extras_require={"test": ["pytest>=7"]},
    scripts=["nwkit/nwkit"],
    include_package_data=True,
    package_data={"": ["data_tree/*.nwk"]},
    classifiers=[
        "Environment :: Console",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: BSD License",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3 :: Only",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
)

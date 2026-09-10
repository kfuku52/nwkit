"""Build deterministic inputs and isolated revision snapshots for benchmark.py."""

import io
import json
import math
import shutil
import subprocess
import sys
import tarfile
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
PATCH_FILES = [
    "nwkit/asr.py",
    "nwkit/asr_figure.py",
    "nwkit/asr_output.py",
    "nwkit/pca.py",
    "tests/test_asr_output.py",
    "tests/test_pca.py",
]


def archive(root, label, revision):
    target = root / label
    target.mkdir()
    data = subprocess.check_output(["git", "archive", revision], cwd=ROOT)
    with tarfile.open(fileobj=io.BytesIO(data)) as source:
        source.extractall(target)
    for path in target.rglob("*"):
        if path.is_file():
            path.chmod(0o755 if path.stat().st_mode & 0o111 else 0o644)
    return target


def tree(n, comb=False):
    nodes = [f"t{i}:1" for i in range(n)]
    if comb:
        text = f"({nodes[0]},{nodes[1]}):1"
        for node in nodes[2:]:
            text = f"({text},{node}):1"
        return text + ";"
    while len(nodes) > 1:
        nodes = [
            f"({nodes[i]},{nodes[i + 1]}):1" if i + 1 < len(nodes) else nodes[i]
            for i in range(0, len(nodes), 2)
        ]
    return nodes[0] + ";"


def main():
    root = Path(tempfile.mkdtemp(prefix="nwkit-measured-"))
    before = archive(root, "before", "8ad9544")
    reviewed = archive(root, "reviewed", "e2465b0")
    radte = archive(root, "radte-first", "0423427")
    after = root / "after"
    shutil.copytree(reviewed, after)
    for name in PATCH_FILES:
        shutil.copyfile(ROOT / name, after / name)
    inputs = root / "inputs"
    inputs.mkdir()
    cases = []

    def add(name, args, tables=("result.tsv",), baseline=before):
        cases.append(
            dict(
                name=name,
                before=str(baseline),
                after=str(after),
                args=args,
                tables=list(tables),
            )
        )

    add("startup", ["--version"], tables=())
    for n in [30, 60, 64, 128, 256, 512, 1024, 2048]:
        (inputs / f"tree{n}.nwk").write_text(tree(n))
        (inputs / f"traits{n}.tsv").write_text(
            "leaf_name\tx\ty\tstate\tbinary\n"
            + "".join(
                f"t{i}\t{math.sin(i * 1.731):.12g}\t{math.sin(i * 0.19) + math.cos(i * 0.713):.12g}"
                f"\t{'a' if i % 3 else 'b'}\t{int(i % 7 in (0, 2, 3))}\n"
                for i in range(n)
            )
        )
    for shape, n in [("balanced", 2048), ("comb", 1600)]:
        source = inputs / f"{shape}{n}.nwk"
        source.write_text(tree(n, comb=shape == "comb"))
        add(
            f"nwk2table-{shape}-{n}",
            ["nwk2table", "-i", str(source), "-o", "result.tsv"],
        )
    for model, sizes, column in [
        ("BM", [128, 1024], "y"),
        ("OU", [128, 512], "y"),
        ("ER", [128, 1024], "state"),
        ("MV-BM", [64, 256], "x,y"),
    ]:
        for n in sizes:
            args = [
                "asr",
                "-i",
                str(inputs / f"tree{n}.nwk"),
                "--input-rooted",
                "yes",
                "--trait",
                str(inputs / f"traits{n}.tsv"),
                "--state-column",
                column,
                "--trait-type",
                "discrete" if model == "ER" else "continuous",
                "--model",
                model,
                "-o",
                "result.tsv",
            ]
            add(f"asr-{model}-{n}", args)
    for n in [128, 512]:
        add(
            f"regress-lambda-{n}",
            [
                "regress",
                "--tree",
                str(inputs / f"tree{n}.nwk"),
                "--data",
                str(inputs / f"traits{n}.tsv"),
                "--responses",
                "y",
                "--predictors",
                "x",
                "--evolution-model",
                "lambda",
                "-o",
                "result.tsv",
            ],
        )
    add(
        "regress-binomial-30",
        [
            "regress",
            "--tree",
            str(inputs / "tree30.nwk"),
            "--data",
            str(inputs / "traits30.tsv"),
            "--responses",
            "binary",
            "--predictors",
            "x",
            "--response-family",
            "binary=binomial",
            "--categorical-responses",
            "binary",
            "--coefficient-penalty",
            "none",
            "-o",
            "result.tsv",
        ],
    )
    ensemble = inputs / "ensemble.nwk"
    ensemble.write_text((tree(128) + "\n") * 4)
    add(
        "asr-ensemble-128",
        [
            "asr",
            "-i",
            str(inputs / "tree128.nwk"),
            "--trait",
            str(inputs / "traits128.tsv"),
            "--state-column",
            "y",
            "--model",
            "BM",
            "--tree-ensemble",
            str(ensemble),
            "--tree-ensemble-out",
            "ensemble.tsv",
            "-o",
            "result.tsv",
        ],
        tables=("result.tsv", "ensemble.tsv"),
        baseline=reviewed,
    )
    add(
        "pca-256",
        [
            "pca",
            "-i",
            str(inputs / "tree256.nwk"),
            "--trait",
            str(inputs / "traits256.tsv"),
            "--columns",
            "x,y",
            "-o",
            "result.tsv",
            "--eigenvalues-out",
            "eigenvalues.tsv",
        ],
        tables=("result.tsv", "eigenvalues.tsv"),
        baseline=reviewed,
    )
    sys.path.insert(0, str(reviewed / "tools"))
    sys.path.insert(0, str(reviewed))
    from benchmark_radte import simulate

    for n in [8, 16]:
        folder = inputs / f"radte{n}"
        simulate(folder, n, 128, 381, 0.1)
        add(
            f"radte-branch-{n}",
            [
                "radte",
                "--gene-tree",
                str(folder / "gene.nwk"),
                "--species-tree",
                str(folder / "species.nwk"),
                "--species-map-tsv",
                str(folder / "mapping.tsv"),
                "--reconcile",
                "lca",
                "--max-age",
                "100",
                "--seed",
                "381",
                "--out-prefix",
                "result",
            ],
            tables=("result.nodes.tsv",),
            baseline=radte,
        )
    manifest = dict(
        revisions=dict(before="8ad9544", reviewed="e2465b0", radte_first="0423427"),
        after_patch_files=PATCH_FILES,
        runs=str(root / "runs"),
        cases=cases,
    )
    (root / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    Path("/tmp/nwkit-measured-path").write_text(str(root))
    print(root)


if __name__ == "__main__":
    main()

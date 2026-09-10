"""Reproduce review findings using disposable files, without editing source inputs."""

import hashlib
import json
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))


def run(directory):
    import numpy as np

    from nwkit.cli import main
    from nwkit.evolution import build_evolutionary_process
    from nwkit.phylogenetic_pca import fit_pca
    from nwkit.util import read_tree

    traits = directory / "traits.tsv"
    traits.write_text("leaf_name\tstate\nA\t1\nB\t2\nC\t4\nD\t5\n")
    base = [
        "asr",
        "-i",
        "((A:1,B:1):1,(C:1,D:1):1);",
        "--input-rooted",
        "yes",
        "--trait",
        str(traits),
        "--state-column",
        "state",
        "--trait-type",
        "continuous",
    ]
    regimes = directory / "regimes.tsv"
    regimes.write_text("branch_id\tregime\n" + "".join(f"{i}\tr\n" for i in range(7)))
    parameters = directory / "parameters.tsv"
    original = "regime\tsigma2\nr\t1\n"
    parameters.write_text(original)
    main(
        base
        + [
            "--model",
            "BMS",
            "--regime-map",
            str(regimes),
            "--regime-parameters",
            str(parameters),
            "-o",
            str(parameters),
        ]
    )
    overwrite = {
        "input_changed": parameters.read_text() != original,
        "replacement_header": parameters.read_text().splitlines()[0],
    }
    assert overwrite["input_changed"]

    reference = directory / "reference.tsv"
    ensemble = directory / "ensemble.tsv"
    trees = directory / "ensemble.nwk"
    reference.write_text("PREVIOUS REFERENCE\n")
    ensemble.write_text("PREVIOUS ENSEMBLE\n")
    trees.write_text("((A:1,B:1):1,(C:1,E:1):1);\n")
    try:
        main(
            base
            + [
                "--model",
                "BM",
                "--sigma2",
                "1",
                "-o",
                str(reference),
                "--tree-ensemble",
                str(trees),
                "--tree-ensemble-out",
                str(ensemble),
            ]
        )
    except ValueError as exc:
        failure = str(exc)
    else:
        raise AssertionError("Expected an ensemble tip-set mismatch")
    partial = {
        "error": failure,
        "reference_changed": reference.read_text() != "PREVIOUS REFERENCE\n",
        "ensemble_changed": ensemble.read_text() != "PREVIOUS ENSEMBLE\n",
    }
    assert partial["reference_changed"] and not partial["ensemble_changed"]

    count = 200
    source = "(t0:1,t1:1):1"
    for i in range(2, count):
        source = f"({source},t{i}:1):1"
    source += ";"
    tree_path = directory / "comb.nwk"
    table_path = directory / "pca-traits.tsv"
    tree_path.write_text(source + "\n")
    table_path.write_text(
        "leaf_name\tx\ty\n"
        + "".join(f"t{i}\t{i}\t{(i * i) % 97}\n" for i in range(count))
    )
    try:
        main(
            [
                "pca",
                "-i",
                str(tree_path),
                "--input-rooted",
                "yes",
                "--trait",
                str(table_path),
                "--columns",
                "x,y",
                "-o",
                str(directory / "scores.tsv"),
            ]
        )
    except RecursionError as exc:
        pca_error = str(exc)
    else:
        raise AssertionError("Expected the recursive PCA copy to fail")
    tree = read_tree(source, 1, True, quiet=True, rooted="yes")
    covariance = build_evolutionary_process(tree, allow_zero=True).tip_covariance(
        [f"t{i}" for i in range(count)]
    )
    fit = fit_pca(covariance, np.array([[i, (i * i) % 97] for i in range(count)]))
    assert fit.status == "ok" and np.isfinite(fit.eigenvalues).all()
    return {
        "regime_parameter_overwrite": overwrite,
        "ensemble_partial_publication": partial,
        "pca_recursive_copy": {
            "tips": count,
            "error": pca_error,
            "same_input_direct_fit_status": fit.status,
            "same_input_direct_fit_eigenvalues": fit.eigenvalues.tolist(),
        },
    }


if __name__ == "__main__":
    with tempfile.TemporaryDirectory(prefix="nwkit-review-") as temporary:
        results = run(Path(temporary))
    results["head"] = subprocess.check_output(
        ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True
    ).strip()
    results["python"] = sys.version
    results["source_sha256"] = {
        name: hashlib.sha256((ROOT / name).read_bytes()).hexdigest()
        for name in ("nwkit/asr.py", "nwkit/asr_tree_ensemble.py", "nwkit/pca.py")
    }
    Path(__file__).with_name("results.json").write_text(
        json.dumps(results, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(results, ensure_ascii=False, indent=2))

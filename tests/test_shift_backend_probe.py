"""An actual original/corrected backend must fail/pass the pBIC contract."""

import csv
import hashlib
import json
import os
import subprocess
from pathlib import Path
from types import SimpleNamespace

import pytest

from nwkit.shift_backend import run_backend
from nwkit.shift_backend_probe import (
    PBIC_CONTRACT,
    R_PBIC_PROBE,
    collect_pbic_validation,
)


@pytest.fixture
def attestation(tmp_path):
    """Small local files exercise ingestion, not correctness of the R method."""
    library = tmp_path / "library"
    library.mkdir()
    (library / "DESCRIPTION").write_bytes(b"Version: 3.0.9\n")

    def write(name, rows):
        with (tmp_path / name).open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=rows[0], delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)

    write(
        "pbic-identity.tsv",
        [{"contract": PBIC_CONTRACT, "version": "3.0.9", "library": str(library)}],
    )
    write(
        "pbic-files.tsv",
        [
            {
                "path": "DESCRIPTION",
                "md5": hashlib.md5(
                    (library / "DESCRIPTION").read_bytes(), usedforsecurity=False
                ).hexdigest(),
            }
        ],
    )
    labels = []
    for root in ("OUfixedRoot", "OUrandomRoot"):
        for phase in ("fixed", "estimated"):
            for representation in ("free", "singleton"):
                for quantity in ("likelihood", "pBIC"):
                    labels.append(f"{root}/{phase}/{representation}/{quantity}")
        labels.extend(
            f"{root}/{suffix}"
            for suffix in (
                "fixed/coordinate_equivalence",
                "shared/likelihood",
                "shared/pBIC",
            )
        )
    checks = [
        {"check": label, "actual": 1, "expected": 1, "absolute_error": 0}
        for label in labels
    ]
    write("pbic-checks.tsv", checks)
    return tmp_path, checks, write


def test_validation_records_the_installed_content(attestation):
    directory, _, _ = attestation
    result = collect_pbic_validation(directory)
    assert result["status"] == "passed"
    assert result["installed_files_sha256"] == {
        "DESCRIPTION": hashlib.sha256(b"Version: 3.0.9\n").hexdigest()
    }
    json.dumps(result, allow_nan=False)


@pytest.mark.parametrize("mutation", ["missing", "duplicate", "unknown", "nan", "gap"])
def test_invalid_capability_results_are_rejected(attestation, mutation):
    directory, checks, write = attestation
    if mutation == "missing":
        checks.pop()
    elif mutation == "duplicate":
        checks[-1] = checks[0]
    elif mutation == "unknown":
        checks[0]["check"] = "different-probe"
    elif mutation == "nan":
        checks[0]["actual"] = float("nan")
    else:
        checks[0]["actual"] = 2
        checks[0]["absolute_error"] = 1
    write("pbic-checks.tsv", checks)
    with pytest.raises(ValueError, match="pBIC capability"):
        collect_pbic_validation(directory)


def test_package_replacement_invalidates_attestation(attestation):
    directory, _, _ = attestation
    (directory / "library" / "DESCRIPTION").write_bytes(b"Version: 3.0.9\nchanged")
    with pytest.raises(ValueError, match="changed during inference"):
        collect_pbic_validation(directory)


@pytest.mark.skipif(
    not os.environ.get("NWKIT_TEST_RSCRIPT"), reason="Optional R backend"
)
def test_corrected_backend_matches_dense_reference(tmp_path):
    script = tmp_path / "probe.R"
    script.write_text(
        R_PBIC_PROBE
        + "\nset.seed(240); before <- .Random.seed; nwkit_pbic_preflight(); "
        "stopifnot(identical(before, .Random.seed))\n"
    )
    subprocess.run(
        [os.environ["NWKIT_TEST_RSCRIPT"], "--vanilla", str(script)],
        cwd=tmp_path,
        check=True,
        capture_output=True,
        text=True,
    )
    result = collect_pbic_validation(tmp_path)
    assert result["contract"] == PBIC_CONTRACT
    assert len(result["checks"]) == 22
    assert max(row["absolute_error"] for row in result["checks"]) <= 2e-6


@pytest.mark.skipif(
    not os.environ.get("NWKIT_TEST_RSCRIPT")
    or not os.environ.get("NWKIT_TEST_ORIGINAL_R_LIBS"),
    reason="Requires an isolated original 3.0.9 installation",
)
def test_original_backend_is_rejected_before_reading_user_data(tmp_path, monkeypatch):
    monkeypatch.setenv("R_LIBS", os.environ["NWKIT_TEST_ORIGINAL_R_LIBS"])
    args = SimpleNamespace(
        rscript=os.environ["NWKIT_TEST_RSCRIPT"],
        criterion="pBIC",
        max_shifts=2,
        root_model="OUfixedRoot",
        search_strategy="exhaustive",
        exhaustive_max_configurations=5000,
        seed=1,
        bootstrap=0,
        bootstrap_seed=1,
        convergence=False,
    )
    # There is deliberately no tree.nwk: the scientific guard must run first.
    with pytest.raises(RuntimeError, match="OUfixedRoot/fixed/free/pBIC"):
        run_backend(tmp_path, args)
    assert not (tmp_path / "model.tsv").exists()
    assert not (tmp_path / "fit.rds").exists()


@pytest.mark.skipif(
    not os.environ.get("NWKIT_TEST_RSCRIPT"), reason="Optional R backend"
)
def test_cli_exports_capability_for_actual_fit(tmp_path):
    from nwkit.cli import main

    examples = Path(__file__).resolve().parents[1] / "examples" / "shift"
    model = tmp_path / "model.json"
    main(
        [
            "shift",
            "--selection",
            "ic",
            "--criterion",
            "pBIC",
            "-i",
            str(examples / "tree.nwk"),
            "--trait",
            str(examples / "traits.tsv"),
            "--state-column",
            "value",
            "--rscript",
            os.environ["NWKIT_TEST_RSCRIPT"],
            "--max-shifts",
            "1",
            "--model-out",
            str(model),
            "-o",
            str(tmp_path / "regimes.tsv"),
        ]
    )
    result = json.loads(model.read_text())
    validation = result["execution"]["pbic_validation"]
    assert validation["status"] == "passed"
    assert validation["version"] == result["backend_version"]
    assert validation["installed_files_sha256"]


@pytest.mark.skipif(
    not os.environ.get("NWKIT_TEST_RSCRIPT")
    or not os.environ.get("NWKIT_TEST_ORIGINAL_R_LIBS"),
    reason="Requires an isolated original 3.0.9 installation",
)
@pytest.mark.parametrize("criterion", ["pBIC", "BIC"])
def test_original_cli_rejects_pbic_without_substitution_or_output_damage(
    tmp_path, monkeypatch, criterion
):
    from nwkit.cli import main

    monkeypatch.setenv("R_LIBS", os.environ["NWKIT_TEST_ORIGINAL_R_LIBS"])
    examples = Path(__file__).resolve().parents[1] / "examples" / "shift"
    model, table = tmp_path / "model.json", tmp_path / "regimes.tsv"
    model.write_text("existing model\n")
    table.write_text("existing table\n")
    arguments = [
        "shift",
        "--selection",
        "ic",
        "--criterion",
        criterion,
        "-i",
        str(examples / "tree.nwk"),
        "--trait",
        str(examples / "traits.tsv"),
        "--state-column",
        "value",
        "--rscript",
        os.environ["NWKIT_TEST_RSCRIPT"],
        "--max-shifts",
        "0",
        "--model-out",
        str(model),
        "-o",
        str(table),
    ]
    if criterion == "pBIC":
        with pytest.raises(RuntimeError, match="pBIC capability check failed"):
            main(arguments)
        assert model.read_text() == "existing model\n"
        assert table.read_text() == "existing table\n"
    else:
        main(arguments)
        result = json.loads(model.read_text())
        assert result["criterion"] == "BIC"
        assert result["execution"]["pbic_validation"] == {
            "status": "not_applicable",
            "criterion": "BIC",
        }

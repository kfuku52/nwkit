"""Exercise adaptive decisions and real upstream PAML continuation."""

import json
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pytest

from nwkit.cli import main
from nwkit.mcmctree_runner import Sampling, _next_round, run_mcmctree
from nwkit.posterior_diagnostics import diagnose_chains


@pytest.mark.parametrize(
    "rhat,ess,draws,state,next_draws",
    [
        (1.001, 401, 20000, "passed", None),
        (1.001, 200, 20000, "unconverged", 48000),
        (1.4, 8, 20000, "unconverged", 40000),
        (1.001, 200, 60000, "unconverged", None),
    ],
)
def test_next_round_uses_rhat_ess_and_bound(
    tmp_path, rhat, ess, draws, state, next_draws
):
    path = tmp_path / "diagnostics.tsv"
    path.write_text(f"parameter\trhat\tess_bulk\tess_tail\nt_n4\t{rhat}\t410\t{ess}\n")
    result = _next_round(path, draws, Sampling(20000, 100, 20000, 60000))
    assert result["state"] == state and result.get("next_draws") == next_draws


@pytest.fixture
def iq2mc_control(tmp_path):
    if not shutil.which("iqtree3") or not shutil.which("mcmctree"):
        pytest.skip("IQ-TREE 3 and upstream PAML MCMCtree are required")
    rng = np.random.default_rng(15)
    ancestor = rng.choice(list("ACGT"), size=500)
    with (tmp_path / "alignment.fa").open("w") as output:
        for taxon in "abc":
            sequence = ancestor.copy()
            mask = rng.random(500) < 0.15
            sequence[mask] = rng.choice(list("ACGT"), size=mask.sum())
            output.write(f">{taxon}\n{''.join(sequence)}\n")
    (tmp_path / "clades.nwk").write_text(
        "((a:0.1,b:0.1)'B(0.1,0.2,0.025,0.025)':0.1,c:0.2)'B(0.3,0.4,0.025,0.025)';\n"
    )
    subprocess.run(
        [
            "iqtree3",
            "-s",
            "alignment.fa",
            "-te",
            "clades.nwk",
            "-m",
            "JC",
            "-T",
            "1",
            "--dating",
            "mcmctree",
            "--mcmc-bds",
            "1,1,0.5",
            "--mcmc-iter",
            "20,10,20",
            "--prefix",
            "iq2mc",
            "--seed",
            "37",
        ],
        cwd=tmp_path,
        check=True,
        capture_output=True,
        timeout=60,
    )
    return tmp_path / "iq2mc.mcmctree.ctl"


def test_real_checkpoints_extend_each_chain_on_the_same_grid(
    iq2mc_control, tmp_path, monkeypatch
):
    rounds = []

    def diagnose(paths, output, **_):
        generations = [
            [int(row.split()[0]) for row in path.read_text().splitlines()[1:]]
            for path in paths
        ]
        assert generations == [list(range(10, len(generations[0]) * 10 + 1, 10))] * 4
        rounds.append(len(generations[0]))
        rhat = "1.1" if rounds[-1] == 20 else "1.001"
        reports = ["parameter\trhat\tess_bulk\tess_tail"]
        for name in paths[0].read_text().splitlines()[0].split()[1:]:
            reports.append(f"{name}\t{rhat}\t450\t450")
        output.write_text("\n".join(reports) + "\n")

    monkeypatch.setattr("nwkit.mcmctree_runner.diagnose_chains", diagnose)
    out = tmp_path / "results"
    original = iq2mc_control.read_bytes()
    assert (
        main(
            [
                "mcmctree-run",
                "-i",
                str(iq2mc_control),
                "--out-dir",
                str(out),
                "--max-generations",
                "1020",
            ]
        )
        == 0
    )
    assert rounds == [20, 40]
    report = json.loads((out / "convergence.json").read_text())
    assert report["state"] == "passed"
    assert report["generations_per_chain"] == 420
    assert report["pooled_samples"] == 160
    assert "95%HPD" in (out / "FigTree.tre").read_text()
    assert (out / "diagnostics.tsv").is_file()
    assert iq2mc_control.read_bytes() == original
    root = Path(report["run"])
    resumed = sorted(root.glob("round-002/chain-*/attempt-*/stdout.log"))
    assert len(resumed) == 4 and all(
        "Reading MCMC state" in log.read_text() for log in resumed
    )
    assert all("burnin = 0" in log.with_name("run.ctl").read_text() for log in resumed)
    receipts = sorted(root.glob("round-*/chain-*/completed.json"))
    assert (
        len({json.loads(path.read_text())["identity"]["seed"] for path in receipts})
        == 8
    )
    run_mcmctree(iq2mc_control, out, max_generations=1020)
    assert sorted(root.glob("round-002/chain-*/attempt-*/stdout.log")) == resumed

    def invalid_diagnostics(paths, output, **_):
        output.write_text("parameter\trhat\tess_bulk\tess_tail\nt_n4\tNA\tNA\tNA\n")

    committed = {path: path.read_bytes() for path in out.iterdir() if path.is_file()}
    monkeypatch.setattr("nwkit.mcmctree_runner.diagnose_chains", invalid_diagnostics)
    with pytest.raises(ValueError, match="MCMCtree diagnostics failed"):
        run_mcmctree(iq2mc_control, out, max_generations=1020)
    assert all(path.read_bytes() == content for path, content in committed.items())


@pytest.fixture
def r_posterior():
    if not shutil.which("Rscript"):
        pytest.skip("R posterior is required")
    if subprocess.run(
        ["Rscript", "-e", 'stopifnot(requireNamespace("posterior", quietly=TRUE))'],
        capture_output=True,
        timeout=20,
    ).returncode:
        pytest.skip("R posterior is required")


def test_real_diagnostics_publish_budget_limited_tree(
    iq2mc_control, tmp_path, r_posterior
):
    out = tmp_path / "budgeted"
    report = run_mcmctree(iq2mc_control, out, max_generations=220)
    assert report["state"] == "unconverged"
    assert report["reason"] == "generation_limit"
    assert report["pooled_samples"] == 80
    assert (out / "mcmctree.out").is_file()


def test_parallel_diagnostics_match_serial_and_reject_invalid_draws(
    tmp_path, r_posterior
):
    rng = np.random.default_rng(310)
    paths = []
    for chain in range(4):
        path = tmp_path / f"chain-{chain}.txt"
        data = np.column_stack([np.arange(1, 201), rng.normal(size=(200, 8))])
        np.savetxt(path, data, header="Gen a b c d e f g h", comments="")
        paths.append(path)
    serial, parallel = tmp_path / "serial.tsv", tmp_path / "parallel.tsv"
    diagnose_chains(paths, serial)
    diagnose_chains(paths, parallel, jobs=4)
    assert serial.read_bytes() == parallel.read_bytes()
    diagnose_chains(paths, parallel, jobs=4, stopping_thresholds=(1.01, 1e8))
    full_rows = [line.split("\t") for line in serial.read_text().splitlines()]
    pending_rows = [line.split("\t") for line in parallel.read_text().splitlines()]
    assert pending_rows[0] == full_rows[0]
    assert [row[:4] for row in pending_rows] == [row[:4] for row in full_rows]
    assert all(row[4:] == ["NA"] * 3 for row in pending_rows[1:])
    diagnose_chains(paths, parallel, jobs=4, stopping_thresholds=(2, 1))
    assert serial.read_bytes() == parallel.read_bytes()
    lines = paths[-1].read_text().splitlines()
    fields = lines[1].split()
    fields[1] = "NA"
    lines[1] = " ".join(fields)
    paths[-1].write_text("\n".join(lines) + "\n")
    with pytest.raises(subprocess.CalledProcessError):
        diagnose_chains(paths, parallel, jobs=4)
    assert serial.read_bytes() == parallel.read_bytes()

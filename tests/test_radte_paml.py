import json
import os
import shutil
import subprocess
import time
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest
from ete4 import Tree

from nwkit.cli import main, parser
from nwkit.radte_paml import (
    _read_chain,
    chain_diagnostics,
    execute_paml,
    paml_alignment,
    paml_control,
    paml_sequence_settings,
    paml_tree,
)
from nwkit.radte_sequence import read_alignment_sequences
from tests.test_radte import cli_inputs, small_chronology
from tests.test_radte_sequence import simulated_alignment, write_alignment


def test_paml_mirrors_emit_only_one_speciation_prior():
    c = small_chronology()
    text, aliases, groups = paml_tree(c)
    assert groups == 1
    assert text.count("#1") == 2
    assert text.count("B{") == 1
    assert ">1<3" in text
    assert set(aliases) == {n.name for n in c.gene.leaves()}
    assert "D" not in text  # Original names are represented by the alias table.


def test_paml_control_and_native_option_exclusivity(tmp_path):
    inputs = cli_inputs(tmp_path)
    args = parser.parse_args(
        ["radte", *inputs, "--reconcile", "lca", "--out-prefix", str(tmp_path / "out")]
    )
    control = paml_control(args, 1, 3, model="hky", usedata=1)
    assert "duplication = 1" in control
    assert "seed = 3" in control
    assert "usedata = 1" in control
    assert "seqtype = 0" in control
    assert "model = 4" in control
    with pytest.raises(ValueError, match="require --backend"):
        main(
            [
                "radte",
                *inputs,
                "--reconcile",
                "lca",
                "--out-prefix",
                str(tmp_path / "out"),
                "--mcmctree-samples",
                "100",
            ]
        )


def test_paml_detects_alphabet_and_accepts_explicit_lg(tmp_path):
    args = parser.parse_args(
        ["radte", *cli_inputs(tmp_path), "--out-prefix", str(tmp_path / "paml")]
    )
    aliases = {"A": "g1", "B": "g2"}
    path = write_alignment(tmp_path, {"B": "AC\tGU", "A": "ACX?"})
    sequences = read_alignment_sequences(path, aliases)
    assert paml_sequence_settings(args, sequences) == ("hky", "exact")
    assert paml_alignment(sequences, aliases, "dna") == "2 4\ng2  ACGT\ng1  ACN?\n"

    path = write_alignment(tmp_path, {"B": "Lbjz", "A": "AX-?"})
    sequences = read_alignment_sequences(path, aliases)
    assert paml_sequence_settings(args, sequences) == ("lg", "approximate")
    assert paml_alignment(sequences, aliases, "protein") == "2 4\ng2  LXXX\ng1  AX-?\n"
    args.mcmctree_likelihood = "exact"
    with pytest.raises(ValueError, match="approximate for protein"):
        paml_sequence_settings(args, sequences)

    path = write_alignment(tmp_path, {"B": "ACGT", "A": "ARNT"})
    args.substitution_model = "lg"
    args.mcmctree_likelihood = None
    assert paml_sequence_settings(args, read_alignment_sequences(path, aliases)) == (
        "lg",
        "approximate",
    )


@pytest.mark.parametrize(
    ("header", "rows", "error"),
    [
        ("Gen t_n3 mu lnL", "1 2 .1 -3\n3 2 .2 -3\n", "sample generations"),
        ("Gen t_n3 mu lnL", "1 2 .1 -3\n2 2 nan -3\n", "finite values"),
        ("Gen t_n3 mu mu lnL", "1 2 .1 .1 -3\n2 2 .2 .2 -3\n", "duplicate columns"),
    ],
)
def test_paml_rejects_invalid_chain(tmp_path, header, rows, error):
    args = SimpleNamespace(mcmctree_samples=2, mcmctree_sampfreq=1, mcmctree_clock=1)
    (tmp_path / "mcmc.txt").write_text(f"{header}\n{rows}", encoding="utf-8")
    with pytest.raises(ValueError, match=error):
        _read_chain(tmp_path, args, Tree("(A,B);", parser=1))


@pytest.mark.skipif(os.name != "posix", reason="Tests POSIX process group cleanup")
def test_paml_timeout_stops_approximation_children(tmp_path, monkeypatch):
    executable = tmp_path / "mcmctree"
    executable.write_text(
        "#!/bin/sh\n(touch started; sleep 1; touch escaped) &\nwait\n", encoding="utf-8"
    )
    executable.chmod(0o700)
    wait = subprocess.Popen.wait

    def wait_after_child_starts(process, timeout=None):
        if timeout is not None:
            deadline = time.monotonic() + 5
            while not (tmp_path / "started").is_file():
                assert time.monotonic() < deadline
                time.sleep(0.01)
        return wait(process, timeout=timeout)

    monkeypatch.setattr(subprocess.Popen, "wait", wait_after_child_starts)
    args = SimpleNamespace(mcmctree_timeout=0.25)
    with pytest.raises(subprocess.TimeoutExpired):
        execute_paml(str(executable), tmp_path, args)
    time.sleep(1)
    assert not (tmp_path / "escaped").exists()


def test_chain_diagnostics_identify_separated_chains():
    rng = np.random.default_rng(11)
    a = pd.DataFrame(dict(Gen=np.arange(1000), t_n5=rng.normal(0, 1, 1000)))
    b = pd.DataFrame(dict(Gen=np.arange(1000), t_n5=rng.normal(0, 1, 1000)))
    assert chain_diagnostics([a, b])["status"] == "passed-basic-diagnostics"
    b.t_n5 += 4
    assert chain_diagnostics([a, b])["status"] == "convergence-not-established"
    assert chain_diagnostics([a])["status"] == "insufficient-chains-or-samples"


@pytest.mark.integration
@pytest.mark.skipif(
    shutil.which("mcmctree") is None,
    reason="Optional external PAML executable is not installed",
)
@pytest.mark.parametrize(
    ("model", "likelihood"), [("jc69", "exact"), ("jc69", "approximate"), ("lg", None)]
)
def test_real_mcmctree_mirror_samples_and_soft_bound_metadata(
    tmp_path, model, likelihood
):
    c = small_chronology()
    inputs = cli_inputs(tmp_path)
    sequences = simulated_alignment(c, sites=500)
    if model == "lg":
        if shutil.which("codeml") is None:
            pytest.skip("Optional external CODEML executable is not installed")
        # Map the variable DNA sites onto amino acid states for the PAML smoke test.
        sequences = {
            name: seq.translate(str.maketrans("ACGT", "ELPQ"))
            for name, seq in sequences.items()
        }
    alignment = write_alignment(tmp_path, sequences)
    prefix = str(tmp_path / "paml")
    main(
        [
            "radte",
            *inputs,
            "--reconcile",
            "lca",
            "--alignment",
            str(alignment),
            "--backend",
            "mcmctree",
            *(
                ["--substitution-model", model, "--mcmctree-likelihood", likelihood]
                if likelihood is not None
                else []
            ),
            "--out-prefix",
            prefix,
            "--figure-out",
            str(tmp_path / "paml.pdf"),
            "--mcmctree-burnin",
            "100",
            "--mcmctree-samples",
            "100",
            "--mcmctree-sampfreq",
            "2",
        ]
    )
    manifest = json.loads((tmp_path / "paml.manifest.json").read_text())
    assert manifest["method"] == "mcmctree-posterior-mean"
    metadata = manifest["sequence_model"]
    assert metadata["model"] == model
    assert metadata["likelihood"] == (likelihood or "approximate")
    if model == "lg":
        work = Path(metadata["work_directory"])
        prefit = work / "approximation"
        control = (prefit / "tmp0001.ctl").read_text()
        assert "seqtype = 2" in control
        assert "model = 2" in control
        assert "aaRatefile = lg.dat" in control
        summary = (prefit / "out.BV").read_bytes()
        for chain in ("chain-1", "chain-2"):
            assert (work / chain / "in.BV").read_bytes() == summary
            assert "usedata = 2 in.BV" in (work / chain / "mcmctree.ctl").read_text()
    assert manifest["calibration_policy"] == "PAML-soft-root-and-speciation-priors"
    assert "convergence-not-established" in manifest["diagnostics"]
    trace = pd.read_csv(tmp_path / "paml.mcmctree-trace.tsv", sep="\t")
    assert (tmp_path / "paml.pdf").read_bytes().startswith(b"%PDF")
    assert len(trace) == 200
    np.testing.assert_array_equal(trace.t_n6, trace.t_n7)
    nodes = pd.read_csv(tmp_path / "paml.nodes.tsv", sep="\t")
    assert set(nodes.bound_policy) == {"PAML-soft-prior"}
    assert (
        nodes.set_index("gene_name").estimated_age["S1"]
        == nodes.set_index("gene_name").estimated_age["S2"]
    )

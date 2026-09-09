import json
import shutil

import numpy as np
import pandas as pd
import pytest

from nwkit.cli import main, parser
from nwkit.radte_paml import chain_diagnostics, paml_control, paml_tree
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
    control = paml_control(args, 1, 3)
    assert "duplication = 1" in control
    assert "seed = 3" in control
    assert "usedata = 1" in control
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
@pytest.mark.parametrize("likelihood", ["exact", "approximate"])
def test_real_mcmctree_mirror_samples_and_soft_bound_metadata(tmp_path, likelihood):
    c = small_chronology()
    inputs = cli_inputs(tmp_path)
    alignment = write_alignment(tmp_path, simulated_alignment(c, sites=500))
    prefix = str(tmp_path / "paml")
    main(
        [
            "radte",
            *inputs,
            "--reconcile",
            "lca",
            "--alignment",
            str(alignment),
            "--substitution-model",
            "jc69",
            "--backend",
            "mcmctree",
            "--mcmctree-likelihood",
            likelihood,
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

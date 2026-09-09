import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from nwkit.cli import main
from nwkit.radte_ensemble import read_ensembles
from tests.test_radte import cli_inputs


def test_species_chronogram_samples_preserve_joint_age_dependence(tmp_path):
    gene = tmp_path / "gene.nwk"
    species = tmp_path / "species.nwk"
    mapping = tmp_path / "map.tsv"
    ensemble = tmp_path / "species-samples.nwk"
    gene.write_text(
        "(((A_1:0.1,B_1:0.1)S1:0.05,(A_2:0.1,B_2:0.1)S2:0.05)D:0.05,C_1:0.2)R;"
    )
    species.write_text("((A:10,B:10)AB:10,C:20)R;")
    mapping.write_text(
        "leaf_name\tspecies_label\nA_1\tA\nB_1\tB\nA_2\tA\nB_2\tB\nC_1\tC\n"
    )
    factors = np.linspace(0.8, 1.2, 20)
    ensemble.write_text(
        "\n".join(f"((A:{10 * f},B:{10 * f})AB:{10 * f},C:{20 * f})R;" for f in factors)
    )
    prefix = tmp_path / "result"
    main(
        [
            "radte",
            "--gene-tree",
            str(gene),
            "--species-tree",
            str(species),
            "--species-map-tsv",
            str(mapping),
            "--reconcile",
            "lca",
            "--out-prefix",
            str(prefix),
            "--species-tree-ensemble",
            str(ensemble),
            "--uncertainty",
            "input-ensemble",
        ]
    )
    nodes = pd.read_csv(str(prefix) + ".nodes.tsv", sep="\t").set_index("gene_name")
    assert nodes.loc["D", "estimated_age"] == pytest.approx(15, abs=1e-4)
    assert nodes.loc["S1", "interval_lower"] == nodes.loc["S2", "interval_lower"]
    assert nodes.sample_clade_presence.eq(1).all()
    samples = pd.read_csv(str(prefix) + ".age-samples.tsv", sep="\t").pivot(
        index="replicate", columns="shared_age_id", values="estimated_age"
    )
    root = samples[nodes.loc["R", "shared_age_id"]]
    ab = samples[nodes.loc["S1", "shared_age_id"]]
    np.testing.assert_allclose(root, 2 * ab, rtol=0, atol=1e-12)
    np.testing.assert_allclose(ab, 10 * factors)
    manifest = json.loads(Path(str(prefix) + ".manifest.json").read_text())
    assert manifest["input_ensemble"]["successes"] == 20


def test_paired_ensembles_reject_different_counts(tmp_path):
    gene, species = tmp_path / "genes.nwk", tmp_path / "species.nwk"
    gene.write_text("(A:1,B:1);\n" * 2)
    species.write_text("(A:1,B:1);\n" * 3)
    with pytest.raises(ValueError, match="equal sample counts"):
        read_ensembles(
            SimpleNamespace(
                gene_tree_ensemble=str(gene), species_tree_ensemble=str(species)
            )
        )


@pytest.mark.parametrize("option", ["--gene-tree-ensemble", "--species-tree-ensemble"])
def test_ensemble_input_cannot_be_replaced_by_an_output(tmp_path, option):
    inputs = cli_inputs(tmp_path)
    prefix = str(tmp_path / "result")
    ensemble = Path(prefix + ".species.tsv")
    ensemble.write_text("(A:10,B:10);\n" * 2)
    before = ensemble.read_bytes()
    with pytest.raises(ValueError, match="overwrite|replace"):
        main(
            [
                "radte",
                *inputs,
                "--reconcile",
                "lca",
                option,
                str(ensemble),
                "--uncertainty",
                "input-ensemble",
                "--out-prefix",
                prefix,
            ]
        )
    assert ensemble.read_bytes() == before


def test_gene_topology_samples_report_absent_reference_duplications(tmp_path):
    gene, species = tmp_path / "gene.nwk", tmp_path / "species.nwk"
    ensemble, mapping = tmp_path / "genes.nwk", tmp_path / "map.tsv"
    original = "((A_1:0.05,A_2:0.05)D_A:0.05,B_1:0.1)Root;"
    alternative = "((A_1:0.05,B_1:0.05)S:0.05,A_2:0.1)Root;"
    gene.write_text(original)
    species.write_text("(A:10,B:10)AB;")
    ensemble.write_text("\n".join([original, alternative] * 20))
    mapping.write_text("leaf_name\tspecies_label\nA_1\tA\nA_2\tA\nB_1\tB\n")
    prefix = tmp_path / "result"
    main(
        [
            "radte",
            "--gene-tree",
            str(gene),
            "--species-tree",
            str(species),
            "--species-map-tsv",
            str(mapping),
            "--reconcile",
            "lca",
            "--max-age",
            "30",
            "--gene-tree-ensemble",
            str(ensemble),
            "--uncertainty",
            "input-ensemble",
            "--out-prefix",
            str(prefix),
        ]
    )
    nodes = pd.read_csv(str(prefix) + ".nodes.tsv", sep="\t").set_index("gene_name")
    assert nodes.loc["D_A", "sample_clade_presence"] == 0.5
    assert nodes.loc["D_A", "sample_event_presence"] == 0.5
    assert pd.isna(nodes.loc["D_A", "interval_lower"])
    assert nodes.loc["D_A", "estimated_age"] == pytest.approx(5, abs=1e-4)

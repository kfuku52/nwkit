from pathlib import Path

import numpy as np
import pytest
from scipy.linalg import expm
from scipy.optimize import approx_fprime

from nwkit.radte_codon import CODONS, DISTANCES, empirical_matrix
from nwkit.radte_sequence import SequenceLikelihood, read_alignment
from tests.test_radte import small_chronology
from tests.test_radte_sequence import write_alignment


@pytest.mark.parametrize("model", ["gy94", "ecmk07", "ecmrest"])
@pytest.mark.parametrize("frequency", ["f", "f1x4", "f3x4", "fq", None])
def test_codon_likelihood_gradient_and_normalization(tmp_path, model, frequency):
    c = small_chronology()
    alignment = write_alignment(
        tmp_path,
        dict(
            A_1="TTTATGGGGAAN---",
            A_2="TTCATAGGGAAGNNN",
            B_1="CTTATTGGTAAAACA",
            B_2="CTCATCGGCAAGACC",
        ),
    )
    likelihood = SequenceLikelihood(
        c, alignment, model=model, codon_frequencies=frequency, gamma_categories=1
    )
    lengths = np.array([n.dist for n in c.edges])
    value, gradient = likelihood.value_gradient(lengths)
    assert np.isfinite(value)
    np.testing.assert_allclose(
        gradient,
        approx_fprime(lengths, lambda x: likelihood.value_gradient(x)[0], 1e-7),
        rtol=2e-4,
        atol=2e-4,
    )
    np.testing.assert_allclose(
        likelihood.transition(0.13, 1)[0], expm(likelihood.q * 0.13), atol=1e-12
    )
    assert np.sum(likelihood.pi[:, None] * likelihood.q * DISTANCES) == pytest.approx(1)
    np.testing.assert_allclose(likelihood.pi @ likelihood.q, 0, atol=1e-14)
    np.testing.assert_allclose(
        likelihood.pi[:, None] * likelihood.q,
        (likelihood.pi[:, None] * likelihood.q).T,
        atol=1e-14,
    )
    if model in {"gy94", "ecmrest"}:
        assert np.all(likelihood.q[DISTANCES > 1] == 0)
    else:
        assert np.any(likelihood.q[DISTANCES > 1] > 0)
    boot = likelihood.bootstrap(np.random.default_rng(3))
    assert boot.raw_matrix.shape == (4, 5)
    assert boot.codon_frequencies == likelihood.codon_frequencies
    assert boot.omega == likelihood.omega


@pytest.mark.parametrize(
    "sequence,error",
    [
        ("AT", "divisible"),
        ("TAA", "Stop codon"),
        ("TAR", "Stop codon"),
        ("A--", "Partial codon"),
        ("ATZ", "Unsupported"),
    ],
)
def test_bad_codons(tmp_path, sequence, error):
    alignment = write_alignment(tmp_path, {"A": sequence})
    with pytest.raises(ValueError, match=error):
        read_alignment(alignment, ["A"], "codon")


def test_high_bit_and_ambiguity(tmp_path):
    alignment = write_alignment(tmp_path, {"A": "GGGNNN---AAR"})
    matrix, states = read_alignment(alignment, ["A"], "codon")
    assert states == CODONS
    assert matrix[0, 0] == 1 << 60
    assert matrix[0, 1] == matrix[0, 2] == (1 << 61) - 1
    assert matrix[0, 3] == (1 << CODONS.index("AAA")) + (1 << CODONS.index("AAG"))


def test_original_ecm_reference_entries():
    for model, first in [("ecmk07", 16.011531), ("ecmrest", 11.192024)]:
        exchange, pi = empirical_matrix(model)
        assert exchange[0, 1] == first
        assert pi.sum() == pytest.approx(1)


def test_generax_species_labels_use_clades_and_keep_ages(tmp_path):
    from ete4 import Tree

    from nwkit.radte_inputs import remap_generax_species

    gene = Tree("(A_1:0.1[&&NHX:S=A],B_1:0.1[&&NHX:S=B])g:0[&&NHX:S=0];", parser=1)
    species = Tree("(A:10,B:10)Root;", parser=1)
    source = tmp_path / "generax.nwk"
    source.write_text("(B:0.1,A:0.2)0;")
    remap_generax_species(gene, species, source)
    assert gene.props["S"] == "Root"
    assert species["A"].dist == 10
    source.write_text("(A:1,C:1)0;")
    with pytest.raises(ValueError, match="identical rooted clades"):
        remap_generax_species(gene, species, source)


@pytest.mark.integration
@pytest.mark.parametrize("model", ["ecmk07", "ecmrest"])
def test_empirical_likelihood_matches_iqtree(tmp_path, model):
    import re
    import shutil
    import subprocess

    from ete4 import Tree

    from nwkit.radte_inputs import build_chronology
    from nwkit.reconcile import build_reconciliation_table

    iqtree = shutil.which("iqtree") or shutil.which("iqtree2")
    if not iqtree:
        pytest.skip("IQ-TREE is required for independent codon likelihood validation")
    gene = Tree("((A:0.12,B:0.15)X:0.17,(C:0.13,D:0.11)Y:0.18)Root;", parser=1)
    species = Tree("((A:5,B:5)X:5,(C:5,D:5)Y:5)Root;", parser=1)
    c = build_chronology(
        gene, species, build_reconciliation_table(gene, species, {n: n for n in "ABCD"})
    )
    rng = np.random.default_rng(123)
    alignment = write_alignment(
        tmp_path, {n: "".join(rng.choice(CODONS, 100)) for n in "ABCD"}
    )
    likelihood = SequenceLikelihood(c, alignment, model=model, gamma_categories=1)
    expected = -likelihood.value_gradient(np.array([n.dist for n in c.edges]))[0]
    # IQ-TREE's ECM implementation uses unit mean codon-event rate. Convert
    # our expected nucleotide changes per codon into its branch-length units.
    factor = -(likelihood.pi @ np.diag(likelihood.q))
    tree = tmp_path / "fixed.nwk"
    tree.write_text(
        f"((A:{0.12 * factor:.17g},B:{0.15 * factor:.17g}):{0.17 * factor:.17g},(C:{0.13 * factor:.17g},D:{0.11 * factor:.17g}):{0.18 * factor:.17g});"
    )
    prefix = tmp_path / "reference"
    result = subprocess.run(
        [
            iqtree,
            "-s",
            str(alignment),
            "--seqtype",
            "CODON",
            "-m",
            model.upper(),
            "-te",
            str(tree),
            "-blfix",
            "-T",
            "1",
            "--prefix",
            str(prefix),
            "-redo",
        ],
        capture_output=True,
        text=True,
        timeout=60,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    report = Path(str(prefix) + ".iqtree").read_text()
    observed = float(
        re.search(r"Log-likelihood of the tree:\s*([-0-9.]+)", report).group(1)
    )
    assert observed == pytest.approx(expected, abs=6e-5)


def test_gy94_estimates_shared_kappa_and_omega(tmp_path):
    from nwkit.radte_sequence_fit import fit_sequence_model

    c = small_chronology()
    rng = np.random.default_rng(35)
    ancestor = rng.choice(CODONS, 120)
    sequences = {}
    for name in ["A_1", "A_2", "B_1", "B_2"]:
        seq = ancestor.copy()
        sites = rng.choice(120, 15, replace=False)
        seq[sites] = rng.choice(CODONS, len(sites))
        sequences[name] = "".join(seq)
    alignment = write_alignment(tmp_path, sequences)
    likelihood = SequenceLikelihood(
        c, alignment, model="gy94", kappa=20, omega=10, gamma_categories=1
    )
    initial = np.array([n.dist for n in c.edges])
    before = likelihood.value_gradient(initial)[0]
    result = fit_sequence_model(likelihood, initial, fit_kappa=True, fit_omega=True)
    assert result["objective"] < before - 1
    assert likelihood.kappa != 20
    assert likelihood.omega != 10
    assert likelihood.fit_settings["fit_omega"]
    assert likelihood.bootstrap(rng).fit_settings["fit_omega"]


def test_gy94_synonymous_and_nonsynonymous_rate_ratios():
    from nwkit.radte_codon import codon_matrix

    q, _ = codon_matrix(
        "gy94", np.ones((1, 1), dtype=np.uint64), np.ones(1), 3, 0.2, "fq"
    )
    aaa, aag, aac, gaa = (CODONS.index(c) for c in ("AAA", "AAG", "AAC", "GAA"))
    # AAA/AAG both encode Lys; AAC encodes Asn and GAA encodes Glu.
    assert q[aaa, aag] / q[aaa, aac] == pytest.approx(3 / 0.2)
    assert q[aaa, gaa] / q[aaa, aac] == pytest.approx(3)


@pytest.mark.parametrize("model", ["gy94", "ecmk07", "ecmrest"])
def test_nonstandard_genetic_code_is_not_silently_misinterpreted(tmp_path, model):
    with pytest.raises(ValueError, match="standard genetic code 1"):
        SequenceLikelihood(small_chronology(), None, model=model, genetic_code=2)

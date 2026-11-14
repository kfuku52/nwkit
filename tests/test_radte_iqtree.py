"""Exercise real IQ-TREE derivatives and the dating-engine boundary."""

import shutil

import numpy as np
import pytest

from nwkit.radte_iqtree import IQTreeLikelihood, model_tokens
from nwkit.radte_sequence import SequenceLikelihood
from tests.test_radte import small_chronology
from tests.test_radte_sequence import write_alignment


@pytest.mark.parametrize(
    "model", ["MFP", "NONREV", "GY+ASC", "GY+G4+R4", "GY+G4+G4", "GY+Fblah"]
)
def test_unsupported_iqtree_models_fail_explicitly(model):
    with pytest.raises(ValueError):
        model_tokens(model)


def alignment(tmp_path):
    from nwkit.radte_codon import CODONS

    rng = np.random.default_rng(83)
    root = rng.choice(CODONS, 120)
    sequences = {}
    for name in ["A_1", "A_2", "B_1", "B_2"]:
        seq = root.copy()
        columns = rng.choice(120, 20, replace=False)
        seq[columns] = rng.choice(CODONS, 20)
        sequences[name] = "".join(seq)
    return write_alignment(tmp_path, sequences)


@pytest.fixture
def iqtree():
    executable = shutil.which("iqtree3")
    if not executable:
        pytest.skip("IQ-TREE runtime required")
    return executable


@pytest.mark.integration
@pytest.mark.parametrize(
    "model,native",
    [
        ("GY{0.5,2}+FQ+G4{1}", "gy94"),
        ("ECMrest+G4{1}", "ecmrest"),
        ("ECMK07+G4{1}", "ecmk07"),
    ],
)
def test_iqtree_likelihood_and_derivatives_match_native(
    tmp_path, iqtree, model, native
):
    c = small_chronology()
    path = alignment(tmp_path)
    exact = IQTreeLikelihood(c, path, model, executable=iqtree)
    reference = SequenceLikelihood(
        c, path, model=native, codon_frequencies="fq" if native == "gy94" else "model"
    )
    lengths = np.array([n.dist for n in c.edges])
    # ECMK07 normalizes multi-nucleotide codon events differently in IQ-TREE.
    scale = -reference.pi @ np.diag(reference.q)
    value, gradient = exact.value_gradient(lengths * scale)
    expected, score = reference.value_gradient(lengths)
    assert value == pytest.approx(expected, abs=2e-5)
    np.testing.assert_allclose(gradient * scale, score, atol=3e-4, rtol=1e-5)
    _, _, g, h, mapping, _ = exact.evaluate(lengths)
    # IQ2MC exports observed diagonals but OPG off-diagonals. Only check
    # diagonals here; NWKIT builds full curvature from IQ-TREE scores.
    for j in range(len(h)):
        edge = int(np.flatnonzero(mapping == j)[0])
        plus, minus = lengths.copy(), lengths.copy()
        plus[edge] += 1e-4
        minus[edge] -= 1e-4
        gp = exact.evaluate(plus)[2]
        gm = exact.evaluate(minus)[2]
        assert ((gp[j] - gm[j]) / 2e-4) == pytest.approx(
            h[j] if h.ndim == 1 else h[j, j], abs=2, rel=0.01
        )
    assert exact.frozen_model == model


@pytest.mark.integration
@pytest.mark.parametrize(
    "model",
    [
        "GY+F3X4+R4",
        "MG+F3X4+G4",
        "GY+F3X4+I+G4",
        "GY+F3X4+I+R4",
        "GTR+G4",
        "HKY+G4",
        "LG+R4",
    ],
)
def test_iqtree_prefit_freezes_model(tmp_path, iqtree, model):
    c = small_chronology()
    path = alignment(tmp_path)
    if model.startswith("LG"):
        path = write_alignment(
            tmp_path,
            {
                name: "ACDEFGHIKLMNPQRSTVWY" * 4 + tail
                for name, tail in zip(
                    ["A_1", "A_2", "B_1", "B_2"], ["AA", "AC", "CD", "DE"], strict=True
                )
            },
        )
    exact = IQTreeLikelihood(c, path, model, executable=iqtree)
    assert "{" in exact.frozen_model
    assert np.isfinite(exact.value_gradient(np.ones(len(c.edges)) * 0.15)[0])
    assert exact.prefit_nll == pytest.approx(
        exact.value_gradient(exact.initial_lengths)[0], abs=1e-5
    )


@pytest.mark.integration
@pytest.mark.parametrize(
    "gene_text",
    [
        "(A_1:0.1,(B_1:0.12,(A_2:0.2,B_2:0.21)S:0.13)X:0.17)D;",
        "((B_1:0.12,(A_2:0.2,B_2:0.21)S:0.13)X:0.17,A_1:0.1)D;",
    ],
)
def test_tiny_lengths_and_single_tip_root_preserve_scores(tmp_path, iqtree, gene_text):
    c = small_chronology(gene_text=gene_text)
    path = alignment(tmp_path)
    model = "GY{0.5,2}+FQ+G4{1}"
    exact = IQTreeLikelihood(c, path, model, executable=iqtree)
    native = SequenceLikelihood(c, path, model="gy94", codon_frequencies="fq")
    lengths = np.array([node.dist for node in c.edges])
    lengths[exact.edge_id[c.gene["A_2"]]] = 1e-12
    expected, score = native.value_gradient(lengths)
    value, gradient = exact.value_gradient(lengths)
    assert value == pytest.approx(expected, abs=3e-5)
    np.testing.assert_allclose(gradient, score, rtol=2e-5, atol=1e-3)
    assert c.gene.name == "D"


@pytest.mark.integration
def test_partial_gy_parameters_and_bootstrap(tmp_path, iqtree):
    c = small_chronology()
    exact = IQTreeLikelihood(
        c, alignment(tmp_path), "GY{0.5}+F3X4+G4{1}", executable=iqtree
    )
    assert exact.frozen_model.startswith("GY{0.5,")
    boot = exact.bootstrap(np.random.default_rng(3))
    assert boot.model == exact.model
    assert all(len(boot.sequences[n]) == len(exact.sequences[n]) for n in exact.names)
    assert any(boot.sequences[n] != exact.sequences[n] for n in exact.names)
    assert np.isfinite(boot.value_gradient(boot.initial_lengths)[0])


@pytest.mark.parametrize(
    "options,error",
    [
        (["--iqtree-model", "GY+R4"], "require --sequence-engine"),
        (
            [
                "--sequence-engine",
                "iqtree",
                "--iqtree-model",
                "GY+R4",
                "--gamma-categories",
                "4",
            ],
            "complete model",
        ),
        (
            [
                "--sequence-engine",
                "iqtree",
                "--substitution-model",
                "gy94",
                "--omega",
                "0.5",
            ],
            "both",
        ),
        (
            [
                "--sequence-engine",
                "iqtree",
                "--iqtree-model",
                "GY+G4",
                "--iqtree-threads",
                "0",
            ],
            "threads",
        ),
        (
            [
                "--sequence-engine",
                "iqtree",
                "--iqtree-model",
                "GY+G4",
                "--genetic-code",
                "0",
            ],
            "genetic code",
        ),
        (
            [
                "--sequence-engine",
                "iqtree",
                "--iqtree-model",
                "GY+G4",
                "--iqtree-executable",
                "/no-such-iqtree",
            ],
            "not found",
        ),
    ],
)
def test_bad_iqtree_controls_preserve_outputs(tmp_path, options, error):
    from nwkit.cli import main
    from tests.test_radte import cli_inputs

    old = tmp_path / "result.dated.nwk"
    old.write_text("previous validated result")
    with pytest.raises(ValueError, match=error):
        main(
            [
                "radte",
                *cli_inputs(tmp_path),
                "--reconcile",
                "lca",
                "--alignment",
                str(alignment(tmp_path)),
                "--out-prefix",
                str(tmp_path / "result"),
                *options,
            ]
        )
    assert old.read_text() == "previous validated result"
    assert not (tmp_path / "result.manifest.json").exists()


@pytest.mark.parametrize(
    "options,expected",
    [
        (["--substitution-model", "gy94"], "GY+F3X4+G4"),
        (["--substitution-model", "ecmk07", "--gamma-categories", "1"], "ECMK07"),
        (
            ["--substitution-model", "hky", "--kappa", "3", "--gamma-shape", "0.5"],
            "HKY{3.0}+G4{0.5}",
        ),
        (
            ["--substitution-model", "gtr", "--gtr-exchangeabilities", "1,2,3,4,5,6"],
            "GTR{1,2,3,4,5,6}+G4",
        ),
        (["--substitution-model", "lg-f", "--gamma-categories", "1"], "LG+F"),
    ],
)
def test_existing_controls_translate_to_iqtree(tmp_path, options, expected):
    from nwkit.cli import parser
    from nwkit.radte_iqtree import requested_model
    from tests.test_radte import cli_inputs

    args = parser.parse_args(
        [
            "radte",
            *cli_inputs(tmp_path),
            "--alignment",
            str(alignment(tmp_path)),
            "--sequence-engine",
            "iqtree",
            "--out-prefix",
            str(tmp_path / "out"),
            *options,
        ]
    )
    assert requested_model(args) == expected


@pytest.mark.integration
def test_three_sequence_fixed_topology(tmp_path, iqtree):
    c = small_chronology(gene_text="((A_1:0.1,B_1:0.2)S:0.1,A_2:0.3)D;")
    path = write_alignment(
        tmp_path, {"A_1": "ACGTTGCA", "B_1": "TCGATGCA", "A_2": "AGGTTGGA"}
    )
    exact = IQTreeLikelihood(c, path, "JC", executable=iqtree)
    native = SequenceLikelihood(c, path, model="jc69", gamma_categories=1)
    lengths = np.array([n.dist for n in c.edges])
    value, score = exact.value_gradient(lengths)
    expected, gradient = native.value_gradient(lengths)
    assert value == pytest.approx(expected, abs=1e-6)
    np.testing.assert_allclose(score, gradient, atol=1e-5, rtol=1e-5)


def test_nonfinite_dependency_derivatives_are_not_accepted(tmp_path):
    from nwkit.radte_iqtree import read_export

    prefix = tmp_path / "broken"
    (tmp_path / "broken.mcmctree.hessian").write_text(
        "4\n(A:1,B:2,(C:3,D:4):5);\n1 2 5 3 4\n1 2 inf 3 4\nHessian\n"
        + "1 0 0 0 1\n" * 5
    )
    with pytest.raises(ValueError, match="Nonfinite IQ-TREE gradient"):
        read_export(prefix)


@pytest.mark.integration
def test_identical_sequences_keep_distinct_gene_edges(tmp_path, iqtree):
    c = small_chronology()
    path = write_alignment(
        tmp_path,
        {"A_1": "ACGTTGCA", "B_1": "TCGATGCA", "A_2": "ACGTTGCA", "B_2": "AGGTTGGA"},
    )
    exact = IQTreeLikelihood(c, path, "JC", executable=iqtree)
    value, gradient = exact.value_gradient(np.full(len(c.edges), 0.15))
    assert len(exact.names) == 4
    assert len(gradient) == 6
    assert np.isfinite(value) and np.isfinite(gradient).all()


def test_missing_iq2mc_export_is_explicit(tmp_path):
    from nwkit.radte_iqtree import read_export

    with pytest.raises(ValueError, match="required IQ2MC"):
        read_export(tmp_path / "missing")
    (tmp_path / "short.mcmctree.hessian").write_text("4\n(A,B,C,D);\n")
    with pytest.raises(ValueError, match="Incomplete"):
        read_export(tmp_path / "short")


@pytest.mark.integration
@pytest.mark.parametrize(
    "model", ["GY{0.5,2}+FQ+G4{1}", "ECMK07+G4{1}", "GY+F3X4+I+R4"]
)
def test_standard_cli_cache_reuses_complete_export(tmp_path, iqtree, model):
    c = small_chronology()
    exact = IQTreeLikelihood(
        c, alignment(tmp_path), model, executable=iqtree, interface="cli"
    )
    try:
        assert "version 3." in exact.version
        for scale in [1.0, 1.03, 0.95]:
            lengths = np.array([n.dist for n in c.edges]) * scale
            before = exact.evaluations
            value, gradient = exact.value_gradient(lengths)
            full = exact.evaluate(lengths)
            assert exact.evaluations == before + 1
            assert full[1] == value
            np.testing.assert_array_equal(full[2][full[4]], gradient)
            assert np.isfinite(full[3]).all()
    finally:
        exact.close()


@pytest.mark.parametrize(
    "version", ["IQ-TREE version 2.4.0", "unrecognized executable"]
)
def test_requires_iqtree3_before_likelihood_work(tmp_path, monkeypatch, version):
    monkeypatch.setattr(shutil, "which", lambda executable: "/test/iqtree3")
    calls = []

    def command(self, args):
        calls.append(args)
        return version

    monkeypatch.setattr(IQTreeLikelihood, "_command", command)
    with pytest.raises(ValueError, match="IQ-TREE 3 or later"):
        IQTreeLikelihood(
            small_chronology(), alignment(tmp_path), "GY+FQ", interface="cli"
        )
    assert calls == [["--version"]]


@pytest.mark.integration
def test_official_cli_tiny_codon_branches_match_matrix_exponential(tmp_path, iqtree):
    from scipy.linalg import expm

    c = small_chronology()
    path = alignment(tmp_path)
    reference = SequenceLikelihood(
        c, path, model="gy94", codon_frequencies="fq", gamma_categories=1
    )

    def transition(length, rate):
        p = expm(reference.q * length * rate)
        return p, rate * reference.q @ p

    reference.transition = transition
    exact = IQTreeLikelihood(c, path, "GY{0.5,2}+FQ", executable=iqtree)
    ordinary = np.array([n.dist for n in c.edges])
    tiny = ordinary.copy()
    for i, node in enumerate(c.edges):
        if node.name == "B_1":
            tiny[i] = 8.050103916839662e-11
        elif node.name == "B_2":
            tiny[i] = 6.291831685936013e-11
    try:
        for lengths in (ordinary, tiny, ordinary * 1.02):
            value, gradient = exact.value_gradient(lengths)
            expected, score = reference.value_gradient(lengths)
            # IQ2MC text scores have the same limited export precision as the
            # ordinary-branch comparison above, unlike a binary session protocol.
            assert value == pytest.approx(expected, abs=2e-5)
            np.testing.assert_allclose(gradient, score, rtol=1e-5, atol=3e-4)
    finally:
        exact.close()


@pytest.fixture
def library_worker():
    from nwkit.iqtree_library import find_worker

    found = find_worker()
    if found is None:
        pytest.skip("Externally built IQ-TREE library worker required")
    return found["executable"]


@pytest.mark.integration
@pytest.mark.parametrize(
    "model",
    [
        "GY{0.5,2}+FQ+G4{1}",
        "GY+F3X4+R4",
        "ECMK07+G4{1}",
        "HKY+I+G4",
    ],
)
def test_library_reuses_loaded_model_and_matches_official_cli(
    tmp_path, iqtree, library_worker, model
):
    c = small_chronology()
    path = alignment(tmp_path)
    exact = IQTreeLikelihood(
        c, path, model, executable=iqtree, interface="library", worker=library_worker
    )
    reference = IQTreeLikelihood(c, path, model, executable=iqtree, interface="cli")
    try:
        assert exact.interface == "library-worker-v1"
        pid = exact.worker.process.pid
        exact.alignment.unlink()  # A new IQ-TREE invocation could no longer load this input.
        base = np.array([n.dist for n in c.edges])
        for lengths in [base, base * [1.2, 0.7, 2, 0.8, 1.1, 0.9], base]:
            exact.cache.clear()
            value, score = exact.value_gradient(lengths)
            expected, gradient = reference.value_gradient(lengths)
            assert value == pytest.approx(expected, abs=2e-5)
            np.testing.assert_allclose(score, gradient, rtol=1e-5, atol=3e-4)
            assert exact.worker.process.pid == pid
        assert exact.worker.process.poll() is None
    finally:
        exact.close()
        reference.close()
    assert exact.worker.process.poll() == 0


@pytest.mark.integration
def test_library_two_threads_matches_single_thread_cli(
    tmp_path, iqtree, library_worker
):
    c = small_chronology()
    path = alignment(tmp_path)
    exact = IQTreeLikelihood(
        c,
        path,
        "GY{0.5,2}+FQ+G4{1}",
        executable=iqtree,
        threads=2,
        interface="library",
        worker=library_worker,
    )
    reference = IQTreeLikelihood(
        c, path, "GY{0.5,2}+FQ+G4{1}", executable=iqtree, interface="cli"
    )
    try:
        lengths = np.array([n.dist for n in c.edges]) * 1.1
        value, gradient = exact.value_gradient(lengths)
        expected, score = reference.value_gradient(lengths)
        assert value == pytest.approx(expected, abs=2e-5)
        np.testing.assert_allclose(gradient, score, rtol=1e-5, atol=3e-4)
    finally:
        exact.close()
        reference.close()


@pytest.mark.integration
def test_library_bootstrap_keeps_the_selected_interface(
    tmp_path, iqtree, library_worker
):
    exact = IQTreeLikelihood(
        small_chronology(),
        alignment(tmp_path),
        "GY{0.5,2}+FQ",
        executable=iqtree,
        interface="library",
        worker=library_worker,
    )
    try:
        replicate = exact.bootstrap(np.random.default_rng(41))
        try:
            assert replicate.interface == exact.interface == "library-worker-v1"
            assert replicate.worker.process.pid != exact.worker.process.pid
            assert (
                replicate.worker_info["library_sha256"]
                == exact.worker_info["library_sha256"]
            )
            assert np.isfinite(replicate.value_gradient(replicate.initial_lengths)[0])
        finally:
            replicate.close()
        assert exact.worker.process.poll() is None
    finally:
        exact.close()

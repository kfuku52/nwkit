import json

import numpy as np
import pytest
from ete4 import Tree
from scipy.linalg import expm
from scipy.optimize import approx_fprime

from nwkit.cli import main
from nwkit.radte import _likelihood_data, read_likelihood_summary
from nwkit.radte_inputs import build_chronology
from nwkit.radte_model import DatingProblem, fit_dates
from nwkit.radte_sequence import (
    QuadraticLikelihood,
    SequenceLikelihood,
    build_quadratic,
    gamma_rates,
    read_alignment,
)
from nwkit.reconcile import build_reconciliation_table
from tests.test_radte import cli_inputs, small_chronology


def write_alignment(tmp_path, sequences):
    path = tmp_path / "alignment.fasta"
    path.write_text("".join(f">{name}\n{seq}\n" for name, seq in sequences.items()))
    return path


@pytest.mark.parametrize("categories", [1, 4])
def test_two_tip_jc69_likelihood_matches_closed_form(tmp_path, categories):
    g = Tree("(A:0.1,B:0.2)Root;", parser=1)
    s = Tree("(A:10,B:10)Root;", parser=1)
    events = build_reconciliation_table(g, s, {"A": "A", "B": "B"})
    c = build_chronology(g, s, events)
    alignment = write_alignment(tmp_path, {"A": "AAAA", "B": "AGGN"})
    likelihood = SequenceLikelihood(
        c, alignment, model="jc69", gamma_categories=categories
    )
    value, gradient = likelihood.value_gradient(np.array([0.1, 0.2]))
    p_same = np.mean([0.25 + 0.75 * np.exp(-4 * 0.3 * r / 3) for r in likelihood.rates])
    p_diff = np.mean([0.25 - 0.25 * np.exp(-4 * 0.3 * r / 3) for r in likelihood.rates])
    expected = -np.log(0.25 * p_same) - 2 * np.log(0.25 * p_diff) - np.log(0.25)
    assert value == pytest.approx(expected, abs=1e-12)
    assert gradient[0] == pytest.approx(gradient[1], abs=1e-12)
    shifted, _ = likelihood.value_gradient(np.array([0.29, 0.01]))
    assert shifted == pytest.approx(value, abs=1e-12)


@pytest.mark.parametrize(
    "model", ["jc69", "hky", "gtr", "f81", "poisson", "lg", "lg-f"]
)
def test_pruning_gradient_matches_finite_difference(tmp_path, model):
    c = small_chronology()
    sequences = {
        "A_1": "ACGTAAGT",
        "A_2": "ACTTACTT",
        "B_1": "GCGTCACT",
        "B_2": "GTGTAACC",
    }
    alignment = write_alignment(tmp_path, sequences)
    likelihood = SequenceLikelihood(c, alignment, model=model)
    lengths = np.array([n.dist for n in c.edges])
    _, gradient = likelihood.value_gradient(lengths)
    numeric = approx_fprime(lengths, lambda x: likelihood.value_gradient(x)[0], 1e-7)
    np.testing.assert_allclose(gradient, numeric, atol=5e-5, rtol=5e-5)
    # Direct matrix exponential is an independent check of the eigen routine.
    p, dp = likelihood.transition(0.2, 1.3)
    np.testing.assert_allclose(p, expm(likelihood.q * 0.26), atol=1e-13)
    np.testing.assert_allclose(dp, likelihood.q @ p * 1.3, atol=1e-13)


def test_sequence_date_gradient_and_root_split_invariance(tmp_path):
    c = small_chronology()
    alignment = write_alignment(
        tmp_path,
        {"A_1": "ACGTACGT", "A_2": "ACGTTCGT", "B_1": "GCGTACCT", "B_2": "GCGTTCCT"},
    )
    likelihood = SequenceLikelihood(c, alignment, model="jc69", gamma_categories=1)
    problem = DatingProblem(c, likelihood=likelihood, rate_sd=0.4)
    x = problem.initial_parameters()
    value, gradient = problem.value_gradient(x)
    numeric = approx_fprime(x, lambda x: problem.value_gradient(x)[0], 1e-7)
    np.testing.assert_allclose(gradient, numeric, atol=2e-5, rtol=2e-5)
    assert np.isfinite(value)
    first, _ = fit_dates(c, likelihood=likelihood, starts=3)
    root_edges = [n for n in c.edges if n.up is c.gene]
    root_edges[0].dist, root_edges[1].dist = 0.25, 0.05
    second, _ = fit_dates(c, likelihood=likelihood, starts=3)
    assert second.log_rate_sd == pytest.approx(first.log_rate_sd, abs=1e-10)
    np.testing.assert_allclose(first.ages, second.ages, atol=2e-4)


def test_gamma_categories_have_unit_mean_and_validate_shape():
    for shape in [0.1, 0.5, 1, 10]:
        rates = gamma_rates(shape, 4)
        assert rates.mean() == pytest.approx(1)
        assert np.all(np.diff(rates) > 0)
    with pytest.raises(ValueError):
        gamma_rates(0, 4)


@pytest.mark.parametrize(
    "sequences",
    [
        {"A": "AC", "B": "ACC"},
        {"A": "AC", "X": "AC"},
        {"A": "AZ", "B": "AC"},
    ],
)
def test_alignment_validation(tmp_path, sequences):
    path = write_alignment(tmp_path, sequences)
    with pytest.raises(ValueError):
        read_alignment(path, ["A", "B"])


def test_quadratic_root_pair_and_gradient():
    likelihood = QuadraticLikelihood(
        np.log([0.3, 0.1]),
        np.array([0.1, -0.2]),
        np.array([[4.0, 0.5], [0.5, 2.0]]),
        12.0,
        np.array([0, 0, 1]),
    )
    lengths = np.array([0.1, 0.2, 0.12])
    value, gradient = likelihood.value_gradient(lengths)
    np.testing.assert_allclose(
        gradient,
        approx_fprime(lengths, lambda x: likelihood.value_gradient(x)[0], 1e-8),
        atol=1e-6,
    )
    assert likelihood.value_gradient(np.array([0.2, 0.1, 0.12]))[0] == value


def simulated_alignment(c, sites=2000, seed=12):
    """Independent JC69 simulator, used to test the pruning implementation."""
    rng = np.random.default_rng(seed)
    ancestors = {c.gene: rng.integers(0, 4, sites)}
    for node in c.gene.traverse():
        if node is c.gene:
            continue
        parent = ancestors[node.up]
        # JC69: retain state with exp(-4b/3), otherwise redraw uniformly.
        redraw = rng.random(sites) > np.exp(-4 * node.dist / 3)
        values = parent.copy()
        values[redraw] = rng.integers(0, 4, redraw.sum())
        ancestors[node] = values
    return {
        str(n.name): "".join(np.array(list("ACGT"))[ancestors[n]])
        for n in c.gene.leaves()
    }


@pytest.mark.parametrize("model", ["hky", "gtr"])
def test_substitution_prefit_improves_likelihood_and_bootstrap_is_isolated(
    tmp_path, model
):
    from nwkit.radte_sequence_fit import default_sequence_model, fit_sequence_model

    c = small_chronology()
    alignment = write_alignment(tmp_path, simulated_alignment(c, sites=2000))
    assert default_sequence_model(alignment) == "gtr"
    likelihood = SequenceLikelihood(c, alignment, model=model, gamma_categories=1)
    lengths = np.array([n.dist for n in c.edges])
    before = likelihood.value_gradient(lengths)[0]
    fit = fit_sequence_model(
        likelihood, lengths, fit_kappa=model == "hky", fit_gtr=model == "gtr"
    )
    assert fit["status"] == "estimated-unclocked-conditional-model"
    assert likelihood.value_gradient(likelihood.initial_lengths)[0] <= before + 1e-7
    sample = likelihood.bootstrap(np.random.default_rng(1))
    assert sample.fit_settings == likelihood.fit_settings
    if model == "gtr":
        original = likelihood.exchangeabilities.copy()
        sample.exchangeabilities[0] *= 2
        np.testing.assert_array_equal(likelihood.exchangeabilities, original)


def test_lg_equilibrium_and_full_protein_alphabet(tmp_path):
    from nwkit.radte_sequence_fit import default_sequence_model

    alignment = write_alignment(
        tmp_path,
        {name: "ARNDCQEGHILKMFPSTWYV" for name in ["A_1", "A_2", "B_1", "B_2"]},
    )
    assert default_sequence_model(alignment) == "lg"
    likelihood = SequenceLikelihood(small_chronology(), alignment, model="lg")
    np.testing.assert_allclose(likelihood.pi @ likelihood.q, 0, atol=1e-12)
    assert -np.dot(likelihood.pi, np.diag(likelihood.q)) == pytest.approx(1)


def test_fixed_model_prefit_removes_input_length_dependence_of_rate_variance(tmp_path):
    from nwkit.radte_sequence_fit import fit_sequence_model

    c = small_chronology()
    alignment = write_alignment(tmp_path, simulated_alignment(c, sites=3000))
    fits = []
    for factor in [1.0, 1.5]:
        likelihood = SequenceLikelihood(c, alignment, model="jc69", gamma_categories=1)
        fit_sequence_model(likelihood, np.array([n.dist for n in c.edges]) * factor)
        fit, _ = fit_dates(c, likelihood=likelihood)
        fits.append(fit)
    np.testing.assert_allclose(fits[0].ages, fits[1].ages, rtol=1e-4, atol=1e-5)
    assert fits[0].log_rate_sd == pytest.approx(fits[1].log_rate_sd, rel=1e-4)


def test_quadratic_prefit_has_identifiable_dimension_and_local_agreement(tmp_path):
    c = small_chronology()
    alignment = write_alignment(tmp_path, simulated_alignment(c))
    exact = SequenceLikelihood(c, alignment, model="jc69", gamma_categories=1)
    initial = np.array([n.dist for n in c.edges])
    quadratic = build_quadratic(exact, initial)
    assert len(quadratic.center) == len(c.edges) - 1
    fractions = (
        initial / np.bincount(quadratic.mapping, weights=initial)[quadratic.mapping]
    )
    optimum = np.exp(quadratic.center[quadratic.mapping]) * fractions
    valid, value_error, gradient_error = quadratic.check(optimum)
    assert valid and value_error < 1e-7 and gradient_error < 1e-7


def test_joint_map_summary_retains_variance_and_ages_with_changed_input_lengths(
    tmp_path,
):
    from nwkit.radte_sequence_fit import fit_sequence_model

    c = small_chronology()
    alignment = write_alignment(tmp_path, simulated_alignment(c, sites=5000))
    exact = SequenceLikelihood(c, alignment, model="jc69", gamma_categories=1)
    fit_sequence_model(exact, np.array([n.dist for n in c.edges]))
    quadratic = build_quadratic(exact, exact.initial_lengths)
    original, _ = fit_dates(c, likelihood=quadratic, inference="joint-map")
    path = tmp_path / "likelihood.json"
    path.write_text(json.dumps(_likelihood_data(quadratic, c, {})))
    # These lengths imply zero rate variance if accidentally used as data.
    for node in c.edges:
        node.dist = 0.1
    summary, _ = read_likelihood_summary(path, c)
    restored, _ = fit_dates(c, likelihood=summary, inference="joint-map")
    assert original.log_rate_sd > 0.3
    assert restored.log_rate_sd == pytest.approx(original.log_rate_sd, abs=1e-10)
    np.testing.assert_allclose(restored.ages, original.ages, atol=1e-6)
    assert restored.objective == pytest.approx(original.objective, abs=1e-7)


def test_alignment_cli_and_summary_roundtrip(tmp_path):
    c = small_chronology()
    inputs = cli_inputs(tmp_path)
    alignment = write_alignment(tmp_path, simulated_alignment(c, sites=5000))
    prefix = str(tmp_path / "sequence")
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
            "--gamma-categories",
            "1",
            "--rate-sd",
            "0.4",
            "--out-prefix",
            prefix,
        ]
    )
    manifest = json.loads(open(prefix + ".manifest.json").read())
    assert manifest["method"] == "sequence-marginal-quadratic"
    summary = json.loads(open(prefix + ".likelihood.json").read())
    if summary["schema"] == "nwkit-radte-quadratic-v1":
        likelihood, _ = read_likelihood_summary(prefix + ".likelihood.json", c)
        assert likelihood.exact is None
        main(
            [
                "radte",
                *inputs,
                "--reconcile",
                "lca",
                "--likelihood-summary",
                prefix + ".likelihood.json",
                "--rate-sd",
                "0.4",
                "--out-prefix",
                str(tmp_path / "summary"),
            ]
        )
        first = json.loads(open(prefix + ".manifest.json").read())
        second = json.loads((tmp_path / "summary.manifest.json").read_text())
        assert second["objective"] == pytest.approx(first["objective"], abs=1e-7)


def test_extremely_short_branches_have_finite_likelihood(tmp_path):
    c = small_chronology()
    alignment = write_alignment(
        tmp_path, {"A_1": "A", "B_1": "G", "A_2": "C", "B_2": "T"}
    )
    likelihood = SequenceLikelihood(c, alignment, model="jc69", gamma_categories=1)
    value, gradient = likelihood.value_gradient(np.full(len(c.edges), 1e-20))
    assert np.isfinite(value) and np.isfinite(gradient).all()

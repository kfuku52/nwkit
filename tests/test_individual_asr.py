"""Joint covariance: independent dense GLS, Rphylopars, missingness and CLI."""

from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from scipy.optimize._numdiff import approx_derivative

from nwkit.cli import main, parser
from nwkit.evolution import build_evolutionary_process
from nwkit.individual_asr import read_individual_data
from nwkit.individual_covariance import (
    _Objective,
    conditional_vector,
    covariance_components,
    fit_individual_covariance,
    profile_gaussian,
)
from nwkit.util import read_tree

EXAMPLE = Path(__file__).resolve().parents[1] / "examples" / "individual_asr"


def arguments(*extra):
    return [
        "asr",
        "--infile",
        str(EXAMPLE / "tree.nwk"),
        "--trait",
        str(EXAMPLE / "individuals.tsv"),
        "--state-column",
        "x,y",
        "--model",
        "MV-BM",
        "--within-species-covariance",
        "full",
        *extra,
    ]


@pytest.fixture(scope="module")
def problem():
    args = parser.parse_args(arguments())
    tree = read_tree(args.infile, 0, False)
    process = build_evolutionary_process(tree)
    data = read_individual_data(args, sorted(tree.leaf_names()), ("x", "y"))
    return tree, process, data, process.tip_covariance(data.names)


def fit_problem(problem, **kwargs):
    _, _, data, c = problem
    return fit_individual_covariance(
        c,
        data.values,
        data.species,
        data.individuals,
        data.coordinates,
        dimension=2,
        **kwargs,
    )


@pytest.fixture(scope="module")
def fitted(problem):
    return fit_problem(problem)


def direct_profile(c, data, sigma, within, method):
    # Deliberately assemble scalar covariances with independent indexing loops.
    n = len(data.values)
    v = np.empty((n, n))
    for i in range(n):
        for j in range(n):
            a, b = data.coordinates[i], data.coordinates[j]
            v[i, j] = c[data.species[i], data.species[j]] * sigma[a, b]
            if data.individuals[i] == data.individuals[j]:
                v[i, j] += within[a, b]
    x = np.eye(2)[data.coordinates]
    inverse = np.linalg.inv(v)
    beta_cov = np.linalg.inv(x.T @ inverse @ x)
    beta = beta_cov @ x.T @ inverse @ data.values
    residual = data.values - x @ beta
    value = (
        len(residual) * np.log(2 * np.pi)
        + np.linalg.slogdet(v)[1]
        + residual @ inverse @ residual
    )
    if method == "REML":
        value += np.linalg.slogdet(x.T @ inverse @ x)[1] - 2 * np.log(2 * np.pi)
    return -0.5 * value, beta, beta_cov, v, x


@pytest.mark.parametrize("method", ["ML", "REML"])
@pytest.mark.parametrize("within", ["full", "diagonal"])
def test_original_unit_likelihood_and_root(problem, method, within):
    fit = fit_problem(problem, method=method, within=within)
    _, _, data, c = problem
    ll, beta, beta_cov, _, _ = direct_profile(c, data, fit.sigma, fit.within, method)
    assert fit.log_likelihood == pytest.approx(ll, abs=1e-8)
    assert fit.root_mean == pytest.approx(beta, abs=1e-9)
    assert fit.root_covariance == pytest.approx(beta_cov, abs=1e-9)
    assert fit.optimizer_success and fit.fit_status == "ok"
    if within == "diagonal":
        assert fit.within[0, 1] == 0


@pytest.mark.parametrize("method", ["ML", "REML"])
@pytest.mark.parametrize("full", [True, False])
def test_analytic_gradient(problem, method, full):
    _, _, data, c = problem
    x = np.eye(2)[data.coordinates]
    obj = _Objective(
        data.values,
        x,
        c[np.ix_(data.species, data.species)],
        data.individuals[:, None] == data.individuals[None, :],
        data.coordinates,
        method,
        full,
    )
    params = np.array(
        [0.1, -0.3, 0.2, -0.7, 0.2, -0.8] if full else [0.1, -0.3, 0.2, -0.7, -0.8]
    )
    gradient = obj(params)[1]
    numeric = approx_derivative(lambda b: obj(b)[0], params).ravel()
    assert gradient == pytest.approx(numeric, rel=2e-6, abs=2e-6)


def test_species_and_individual_conditioning_against_dense_oracle(problem, fitted):
    tree, process, data, c = problem
    _, beta, bc, v, x = direct_profile(c, data, fitted.sigma, fitted.within, "REML")
    inverse = np.linalg.inv(v)
    for individual in (None, 1):
        s = 0
        cross = c[s, data.species][None, :] * fitted.sigma[:, data.coordinates]
        prior = c[s, s] * fitted.sigma
        same = None if individual is None else data.individuals == individual
        if same is not None:
            cross += same[None, :] * fitted.within[:, data.coordinates]
            prior = prior + fitted.within
        mean = beta + cross @ inverse @ (data.values - x @ beta)
        remainder = np.eye(2) - cross @ inverse @ x
        covariance = prior - cross @ inverse @ cross.T + remainder @ bc @ remainder.T
        actual = conditional_vector(
            fitted, data.coordinates, c[s, data.species], c[s, s], same_individual=same
        )
        assert actual[0] == pytest.approx(mean, abs=1e-8)
        assert actual[1] == pytest.approx(covariance, abs=1e-8)
    root = conditional_vector(fitted, data.coordinates, np.zeros(len(data.values)), 0)
    assert root[0] == pytest.approx(fitted.root_mean)
    assert root[1] == pytest.approx(fitted.root_covariance)


def test_scale_translation_and_branch_units(problem, fitted):
    _, _, data, c = problem
    scale = np.array([1e8, 1e-7])
    offset = np.array([1e9, 2e-6])
    second = fit_individual_covariance(
        c * 17,
        data.values * scale[data.coordinates] + offset[data.coordinates],
        data.species,
        data.individuals,
        data.coordinates,
        dimension=2,
    )
    assert second.sigma / scale[:, None] / scale[None, :] * 17 == pytest.approx(
        fitted.sigma, abs=2e-5
    )
    assert second.within / scale[:, None] / scale[None, :] == pytest.approx(
        fitted.within, abs=2e-5
    )
    assert (second.root_mean - offset) / scale == pytest.approx(
        fitted.root_mean, abs=2e-5
    )
    shift = np.log(scale[data.coordinates]).sum() - np.log(scale).sum()
    assert second.log_likelihood == pytest.approx(
        fitted.log_likelihood - shift, abs=1e-7
    )


def test_observation_order_does_not_change_fit(problem, fitted):
    _, _, data, c = problem
    order = np.random.default_rng(5).permutation(len(data.values))
    other = fit_individual_covariance(
        c,
        data.values[order],
        data.species[order],
        data.individuals[order],
        data.coordinates[order],
        dimension=2,
    )
    assert other.log_likelihood == pytest.approx(fitted.log_likelihood, abs=1e-7)
    assert other.sigma == pytest.approx(fitted.sigma, abs=2e-5)


def test_no_replication_rejected(problem):
    _, _, data, c = problem
    keep = np.array([data.keys[i][1] == "1" for i in data.individuals])
    with pytest.raises(ValueError, match="replicate degrees"):
        fit_individual_covariance(
            c,
            data.values[keep],
            data.species[keep],
            data.individuals[keep],
            data.coordinates[keep],
            dimension=2,
        )


def test_zero_within_variation_rejected(problem):
    _, _, data, c = problem
    values = data.species.astype(float) + data.coordinates
    with pytest.raises(ValueError, match="Zero within-species"):
        fit_individual_covariance(
            c, values, data.species, data.individuals, data.coordinates, dimension=2
        )


def test_trait_pair_without_paired_individuals_rejected(problem):
    _, _, data, c = problem
    with pytest.raises(ValueError, match="paired within-species"):
        fit_individual_covariance(
            c,
            data.values,
            data.species,
            data.individuals * 2 + data.coordinates,
            data.coordinates,
            dimension=2,
        )
    fit = fit_individual_covariance(
        c,
        data.values,
        data.species,
        data.individuals * 2 + data.coordinates,
        data.coordinates,
        dimension=2,
        within="diagonal",
    )
    assert fit.optimizer_success


def test_singular_paired_within_contrasts_rejected(problem):
    _, _, data, c = problem
    species = np.repeat(np.arange(8), 6)
    individuals = np.repeat(np.arange(24), 2)
    coordinates = np.tile([0, 1], 24)
    values = individuals.astype(float) * (1 + coordinates) + species
    with pytest.raises(ValueError, match="Singular paired"):
        fit_individual_covariance(
            c, values, species, individuals, coordinates, dimension=2
        )


def test_covariance_design_not_identifiable(problem):
    _, _, data, c = problem
    with pytest.raises(ValueError, match="identifiable"):
        fit_individual_covariance(
            np.ones_like(c),
            data.values,
            data.species,
            data.individuals,
            data.coordinates,
            dimension=2,
        )


def test_cli_outputs_are_distinct_estimands(tmp_path):
    paths = {
        key: tmp_path / (key + ".tsv")
        for key in ("outfile", "model-out", "individual-out", "covariance-out")
    }
    tree_out = tmp_path / "nodes.nwk"
    extra = [item for key, path in paths.items() for item in ("--" + key, str(path))]
    assert main(arguments(*extra, "--tree-out", str(tree_out), "--target", "all")) in (
        0,
        None,
    )
    nodes = pd.read_csv(paths["outfile"], sep="\t")
    individual = pd.read_csv(paths["individual-out"], sep="\t")
    model = pd.read_csv(paths["model-out"], sep="\t").iloc[0]
    assert len(nodes) == 30
    assert set(nodes.loc[nodes.node_class == "leaf", "estimand"]) == {
        "latent_species_mean"
    }
    assert (nodes.loc[nodes.node_class == "leaf", "variance"] > 0).all()
    assert len(individual) == 56 and individual.is_imputed.sum() == 1
    assert (individual.loc[~individual.is_imputed, "variance"] == 0).all()
    assert (individual.loc[individual.is_imputed, "variance"] > 0).all()
    assert (
        not model.parameter_uncertainty_included
        and model.root_mean_uncertainty_included
    )
    assert tree_out.is_file()


@pytest.mark.parametrize(
    "extra,match",
    [
        (["--model", "MV-OU"], "requires --model MV-BM"),
        (["--standard-error-column", "sx,sy"], "does not support"),
        (["--bootstrap-out", "a"], "does not support"),
        (["--figure-out", "a.png"], "does not support"),
    ],
)
def test_unsupported_cli_options(extra, match, capsys):
    with pytest.raises(ValueError, match=match):
        main(arguments(*extra))


def test_output_alias_preserves_input(tmp_path, capsys):
    source = tmp_path / "data.tsv"
    source.write_bytes((EXAMPLE / "individuals.tsv").read_bytes())
    before = source.read_bytes()
    with pytest.raises(ValueError, match="overwrite"):
        main(arguments("--trait", str(source), "--individual-out", str(source)))
    assert source.read_bytes() == before


@pytest.mark.parametrize(
    "mutation,match",
    [
        ("duplicate", "Duplicate"),
        ("blank_id", "nonempty"),
        ("unknown_trait", "absent"),
        ("infinite", "non-finite"),
    ],
)
def test_invalid_long_input(tmp_path, mutation, match, capsys):
    table = pd.read_csv(EXAMPLE / "individuals.tsv", sep="\t", dtype=str)
    if mutation == "duplicate":
        table = pd.concat([table, table.iloc[:1]])
    elif mutation == "blank_id":
        table.loc[0, "individual_id"] = ""
    elif mutation == "unknown_trait":
        table.loc[0, "trait"] = "typo"
    else:
        table.loc[0, "value"] = "inf"
    path = tmp_path / "bad.tsv"
    table.to_csv(path, sep="\t", index=False)
    with pytest.raises(ValueError, match=match):
        main(arguments("--trait", str(path)))


def test_fixed_profile_and_components(problem):
    _, _, data, c = problem
    sigma = np.array([[1, 0.3], [0.3, 0.7]])
    w = np.diag([0.2, 0.4])
    v = covariance_components(
        c, data.species, data.individuals, data.coordinates, sigma, w
    )
    expected, beta, *_ = direct_profile(c, data, sigma, w, "REML")
    fit = profile_gaussian(data.values, np.eye(2)[data.coordinates], v, "REML")
    assert -fit.value == pytest.approx(expected, abs=1e-9)
    assert fit.beta == pytest.approx(beta)


@pytest.mark.parametrize(
    "method,within,ll,sigma,w",
    [
        (
            "REML",
            "full",
            -68.33088365421,
            [[1.5104815013644, 0.8834425281471], [0.8834425281471, 0.6130486745009]],
            [[0.5832552992086, 0.2269636412353], [0.2269636412353, 0.3928968475086]],
        ),
        (
            "ML",
            "full",
            -68.55056414398,
            [[1.2702069495857, 0.7580971186198], [0.7580971186198, 0.5175814723749]],
            [[0.5899079053504, 0.2258513291621], [0.2258513291621, 0.3977991120154]],
        ),
        (
            "REML",
            "diagonal",
            -70.79712203657,
            [[1.4826651857146, 0.9456275376191], [0.9456275376191, 0.6476783438966]],
            [[0.5855340746582, 0], [0, 0.3903777834634]],
        ),
        (
            "ML",
            "diagonal",
            -70.8057825924197,
            [
                [1.239495611380254, 0.822559843503292],
                [0.822559843503292, 0.552298621470945],
            ],
            [[0.602323605260564, 0], [0, 0.400818806459041]],
        ),
    ],
)
def test_rphylopars_0310_reference(problem, method, within, ll, sigma, w):
    # examples/individual_asr/reference.R: same tree, individual grouping,
    # missing coordinate, unequal replication, common W, and likelihood convention.
    fit = fit_problem(problem, method=method, within=within)
    assert fit.log_likelihood == pytest.approx(ll, abs=2e-7)
    assert fit.sigma == pytest.approx(np.array(sigma), abs=2e-5)
    assert fit.within == pytest.approx(np.array(w), abs=2e-5)


def test_rphylopars_species_and_root_reference(problem, fitted):
    _, _, data, c = problem
    mean, _ = conditional_vector(fitted, data.coordinates, c[0, data.species], c[0, 0])
    assert mean == pytest.approx([2.953263237475, 6.034505470594], abs=2e-6)
    assert fitted.root_mean == pytest.approx([3.756005608840, 6.684008910934], abs=2e-6)


def test_missing_species_trait_and_entire_individual(tmp_path):
    table = pd.read_csv(EXAMPLE / "individuals.tsv", sep="\t", dtype=str)
    table.loc[(table.leaf_name == "H") & (table.trait == "y"), "value"] = "NA"
    table = pd.concat(
        [
            table,
            pd.DataFrame(
                [{"leaf_name": "A", "individual_id": "NA", "trait": "x", "value": "NA"}]
            ),
        ]
    )
    source = tmp_path / "missing.tsv"
    table.to_csv(source, sep="\t", index=False)
    nodes_path, individuals_path = tmp_path / "nodes.tsv", tmp_path / "individuals.tsv"
    main(
        arguments(
            "--trait",
            str(source),
            "--target",
            "missing-leaf",
            "--outfile",
            str(nodes_path),
            "--individual-out",
            str(individuals_path),
        )
    )
    nodes = pd.read_csv(nodes_path, sep="\t")
    predictions = pd.read_csv(individuals_path, sep="\t", keep_default_na=False)
    assert set(nodes.name) == {"H"}
    assert list(nodes.is_imputed) == [False, True]
    assert (predictions.loc[predictions.individual_id == "NA", "variance"] > 0).all()
    assert len(predictions.loc[predictions.individual_id == "NA"]) == 2


def test_output_transaction_recovers_on_tree_write_failure(tmp_path, monkeypatch):
    import nwkit.individual_asr as module

    paths = [tmp_path / name for name in ("nodes.tsv", "individuals.tsv", "tree.nwk")]
    for path in paths:
        path.write_text("original")

    def fail(*args, **kwargs):
        raise OSError("tree export failed")

    monkeypatch.setattr(module, "write_multivariate_tree", fail)
    with pytest.raises(OSError, match="tree export failed"):
        main(
            arguments(
                "--outfile",
                str(paths[0]),
                "--individual-out",
                str(paths[1]),
                "--tree-out",
                str(paths[2]),
            )
        )
    assert all(path.read_text() == "original" for path in paths)


@pytest.mark.slow
def test_simulated_covariances_recover_separate_components():
    # Independent matrix-normal simulation with 24 species, four individuals,
    # and a fixed, ignorable missing-coordinate pattern. Verify average recovery
    # with broad Monte Carlo tolerances, not a fortuitous single fitted dataset.
    sigma = np.array([[1.1, 0.3], [0.3, 0.8]])
    within = np.array([[0.35, 0.12], [0.12, 0.5]])
    n, reps, p = 24, 4, 2
    clade = np.repeat(np.arange(6), 4)
    c = np.eye(n) * 0.8 + 0.7 * (clade[:, None] == clade[None, :])
    species = np.repeat(np.arange(n), reps * p)
    individuals = np.repeat(np.arange(n * reps), p)
    coordinates = np.tile(np.arange(p), n * reps)
    observed = np.arange(len(species)) % 17 != 0
    estimates = []
    for seed in range(12):
        rng = np.random.default_rng(seed + 301)
        latent = (
            np.linalg.cholesky(c)
            @ rng.normal(size=(n, p))
            @ np.linalg.cholesky(sigma).T
        )
        errors = rng.normal(size=(n, reps, p)) @ np.linalg.cholesky(within).T
        values = (latent[:, None, :] + errors).reshape(-1)
        fit = fit_individual_covariance(
            c,
            values[observed],
            species[observed],
            individuals[observed],
            coordinates[observed],
            dimension=p,
        )
        assert fit.fit_status == "ok"
        estimates.append([fit.sigma, fit.within])
    means = np.mean(estimates, axis=0)
    assert means[0] == pytest.approx(sigma, abs=0.22)
    assert means[1] == pytest.approx(within, abs=0.065)


def test_three_trait_partial_observations():
    rng = np.random.default_rng(74)
    n, reps, p = 10, 5, 3
    c = np.eye(n) + 0.4
    species = np.repeat(np.arange(n), reps * p)
    individuals = np.repeat(np.arange(n * reps), p)
    coordinates = np.tile(np.arange(p), n * reps)
    latent = rng.normal(size=(n, p))
    values = (latent[:, None, :] + rng.normal(size=(n, reps, p)) * 0.4).reshape(-1)
    keep = np.arange(len(values)) % 13 != 0
    fit = fit_individual_covariance(
        c,
        values[keep],
        species[keep],
        individuals[keep],
        coordinates[keep],
        dimension=p,
    )
    assert fit.sigma.shape == (3, 3) and fit.fit_status == "ok"
    assert np.linalg.eigvalsh(fit.within).min() > 0


def test_whole_evolutionary_component_near_zero_is_flagged():
    species = np.repeat(np.arange(8), 8)
    individuals = np.repeat(np.arange(32), 2)
    coordinates = np.tile([0, 1], 32)
    values = np.tile(
        np.array([[-1, -0.5], [1, 0.5], [-0.3, 1], [0.3, -1]]).reshape(-1), 8
    )
    fit = fit_individual_covariance(
        np.eye(8), values, species, individuals, coordinates, dimension=2
    )
    assert fit.fit_status == "boundary_covariance"
    assert np.max(np.abs(fit.sigma)) < 1e-7
    assert np.linalg.eigvalsh(fit.within).min() > 0.1


def test_stdout_written_once_after_success(capsys):
    from io import StringIO

    main(arguments("--target", "all"))
    text = capsys.readouterr().out
    assert text.count("branch_id\t") == 1
    assert len(pd.read_csv(StringIO(text), sep="\t")) == 30


def test_stdout_remains_empty_on_export_failure(tmp_path, monkeypatch, capsys):
    import nwkit.individual_asr as module

    def fail(*args, **kwargs):
        raise OSError("tree export failed")

    monkeypatch.setattr(module, "write_multivariate_tree", fail)
    with pytest.raises(OSError, match="tree export failed"):
        main(arguments("--tree-out", str(tmp_path / "tree.nwk")))
    assert capsys.readouterr().out == ""


def test_unobserved_long_branch_does_not_change_fit(problem, fitted):
    _, _, data, c = problem
    augmented = np.zeros((9, 9))
    augmented[:8, :8] = c
    augmented[8, 8] = 1e20
    other = fit_individual_covariance(
        augmented,
        data.values,
        data.species,
        data.individuals,
        data.coordinates,
        dimension=2,
    )
    assert other.log_likelihood == pytest.approx(fitted.log_likelihood, abs=1e-8)
    assert other.sigma == pytest.approx(fitted.sigma, abs=1e-8)

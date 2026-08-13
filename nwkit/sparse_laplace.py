"""Sparse latent-Gaussian helpers for large phylogenetic models."""

from dataclasses import dataclass
from typing import Mapping

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import SuperLU, eigsh, splu


@dataclass(frozen=True)
class SparseCovarianceModel:
    """A normalized tip covariance represented by sparse latent variables."""

    precision: sparse.csc_matrix
    tip_loading: sparse.csr_matrix
    logdet_covariance: float
    sampling_parent: np.ndarray
    sampling_transition: np.ndarray
    sampling_variance: np.ndarray
    covariance_scale: float = 1.0
    sampling_precision_factor: sparse.csr_matrix | None = None

    @property
    def n_tips(self) -> int:
        return int(self.tip_loading.shape[0])

    @property
    def n_states(self) -> int:
        return int(self.precision.shape[0])

    def sample(self, rng: np.random.Generator, variance: float = 1.0) -> np.ndarray:
        """Draw one tip vector without materializing its covariance."""
        scale = float(np.sqrt(variance))
        if len(self.sampling_parent) == self.n_states:
            states = np.empty(self.n_states, dtype=float)
            innovations = rng.normal(size=self.n_states)
            for index in range(self.n_states):
                parent = int(self.sampling_parent[index])
                parent_value = 0.0 if parent < 0 else states[parent]
                states[index] = (
                    self.sampling_transition[index] * parent_value
                    + scale
                    * np.sqrt(self.sampling_variance[index])
                    * innovations[index]
                )
        else:
            precision_factor = self.precision_factor()
            factor = factor_sparse_positive_definite(self.precision)
            perturbation = precision_factor.T @ rng.normal(
                size=precision_factor.shape[0]
            )
            states = scale * factor.solve(np.asarray(perturbation).reshape(-1))
        return np.asarray(self.tip_loading @ states, dtype=float)

    def precision_factor(self) -> sparse.csr_matrix:
        """Return a sparse B satisfying precision = B.T @ B."""
        if self.sampling_precision_factor is not None:
            return sparse.csr_matrix(self.sampling_precision_factor, dtype=float)
        if (
            len(self.sampling_parent) != self.n_states
            or len(self.sampling_transition) != self.n_states
            or len(self.sampling_variance) != self.n_states
        ):
            raise ValueError("This sparse covariance has no exact sampling factor.")
        if (
            not np.isfinite(self.sampling_transition).all()
            or not np.isfinite(self.sampling_variance).all()
            or np.any(self.sampling_variance <= 0.0)
        ):
            raise ValueError("Sparse covariance sampling parameters are invalid.")
        rows = np.repeat(np.arange(self.n_states, dtype=int), 2)
        columns = np.empty(self.n_states * 2, dtype=int)
        values = np.empty(self.n_states * 2, dtype=float)
        inverse_sd = 1.0 / np.sqrt(self.sampling_variance)
        columns[0::2] = np.arange(self.n_states, dtype=int)
        values[0::2] = inverse_sd
        columns[1::2] = np.maximum(self.sampling_parent, 0)
        values[1::2] = np.where(
            self.sampling_parent < 0,
            0.0,
            -self.sampling_transition * inverse_sd,
        )
        factor = sparse.csr_matrix(
            (values, (rows, columns)), shape=(self.n_states, self.n_states)
        )
        factor.eliminate_zeros()
        return factor

    def materialize(self) -> np.ndarray:
        """Materialize the tip covariance for small reference calculations."""
        factor = factor_sparse_positive_definite(self.precision)
        solved = factor.solve(self.tip_loading.T.toarray())
        covariance = np.asarray(self.tip_loading @ solved, dtype=float)
        return (covariance + covariance.T) / 2.0


@dataclass(frozen=True)
class SparsePositiveDefiniteFactor:
    """Sparse LU factorization used as an SPD solve and log-determinant."""

    factor: SuperLU
    logdet: float

    def solve(self, values: np.ndarray) -> np.ndarray:
        return np.asarray(self.factor.solve(np.asarray(values, dtype=float)))


@dataclass(frozen=True)
class SparseLatentModel:
    """Combined independent covariance components for a Laplace calculation."""

    precision: sparse.csc_matrix
    loading: sparse.csr_matrix
    sampling_precision_factor: sparse.csr_matrix | None
    prior_logdet: float
    component_slices: Mapping[str, slice]
    component_loadings: Mapping[str, sparse.csr_matrix]

    @property
    def n_states(self) -> int:
        return int(self.precision.shape[0])


@dataclass(frozen=True)
class SparseLatentSampler:
    """Reusable exact sampler for one sparse latent-Gaussian model."""

    loading: sparse.csr_matrix
    precision_factor: sparse.csr_matrix
    solver: SparsePositiveDefiniteFactor

    def sample(self, rng: np.random.Generator) -> np.ndarray:
        perturbation = self.precision_factor.T @ rng.normal(
            size=self.precision_factor.shape[0]
        )
        states = self.solver.solve(np.asarray(perturbation).reshape(-1))
        return np.asarray(self.loading @ states).reshape(-1)


@dataclass(frozen=True)
class ContinuousPredictorUncertainty:
    """Low-rank species uncertainty lifted to gene tips by integer mapping."""

    factor: np.ndarray | sparse.spmatrix
    observation_index: np.ndarray
    row_scale: np.ndarray | None = None


@dataclass(frozen=True)
class GroupedPredictorUncertainty:
    """Independent per-species factor uncertainty shared by its gene tips."""

    factors: tuple[np.ndarray, ...]
    observation_index: np.ndarray


@dataclass(frozen=True)
class JointPredictorUncertainty:
    """Correlated multicolumn uncertainty with common latent factor columns."""

    factors: tuple[np.ndarray | sparse.spmatrix, ...]
    row_scale: np.ndarray | None = None


@dataclass(frozen=True)
class GmrfPredictorUncertainty:
    """Species-level predictor uncertainty retained as a sparse GMRF."""

    model: SparseCovarianceModel
    observation_index: np.ndarray
    row_scale: np.ndarray | None = None


def factor_sparse_positive_definite(
    matrix: sparse.spmatrix,
) -> SparsePositiveDefiniteFactor:
    """Factor a symmetric positive-definite sparse matrix with SuperLU.

    SuperLU itself accepts nonsingular indefinite matrices.  Check symmetry and
    the smallest algebraic eigenvalue explicitly before using its absolute
    determinant as an SPD log-determinant.
    """
    values = sparse.csc_matrix(matrix, dtype=float)
    if values.shape[0] != values.shape[1] or values.shape[0] == 0:
        raise np.linalg.LinAlgError("Sparse SPD matrix must be non-empty and square.")
    difference = values - values.T
    scale = max(1.0, float(np.max(np.abs(values.data), initial=0.0)))
    tolerance = np.finfo(float).eps * scale * max(1, values.shape[0]) * 100.0
    if difference.nnz and float(np.max(np.abs(difference.data))) > tolerance:
        raise np.linalg.LinAlgError("Sparse SPD matrix must be symmetric.")
    symmetric = ((values + values.T) * 0.5).tocsc()
    if values.shape[0] == 1:
        smallest = float(symmetric[0, 0])
    else:
        try:
            smallest = float(
                eigsh(
                    symmetric,
                    k=1,
                    which="SA",
                    return_eigenvectors=False,
                    tol=1e-7,
                )[0]
            )
        except Exception as exc:  # pragma: no cover - backend-specific failures
            raise np.linalg.LinAlgError("Sparse SPD eigenvalue check failed.") from exc
    if not np.isfinite(smallest) or smallest <= tolerance:
        raise np.linalg.LinAlgError("Sparse matrix is not positive definite.")
    return factor_sparse_nonsingular(symmetric)


def factor_sparse_nonsingular(
    matrix: sparse.spmatrix,
) -> SparsePositiveDefiniteFactor:
    """Factor a general nonsingular sparse matrix and return log(abs(det))."""
    values = sparse.csc_matrix(matrix, dtype=float)
    if values.shape[0] != values.shape[1] or values.shape[0] == 0:
        raise np.linalg.LinAlgError("Sparse matrix must be non-empty and square.")
    factor = splu(values, permc_spec="COLAMD")
    diagonal = np.asarray(factor.U.diagonal(), dtype=float)
    if not len(diagonal) or np.any(diagonal == 0.0) or not np.isfinite(diagonal).all():
        raise np.linalg.LinAlgError("Sparse matrix factorization is singular.")
    logdet = float(np.sum(np.log(np.abs(diagonal))))
    if not np.isfinite(logdet):
        raise np.linalg.LinAlgError("Sparse matrix log-determinant is non-finite.")
    return SparsePositiveDefiniteFactor(factor, logdet)


def combine_sparse_covariance_models(
    components: Mapping[str, tuple[float, SparseCovarianceModel]],
    *,
    random_dimension: int = 1,
) -> SparseLatentModel:
    """Combine covariance components in independent latent coordinates."""
    identity = sparse.eye(random_dimension, format="csc")
    precisions = []
    sampling_factors = []
    loadings = []
    component_slices = {}
    component_loadings = {}
    prior_logdet = 0.0
    position = 0
    for name, (variance, model) in components.items():
        if not np.isfinite(variance) or variance <= 0.0:
            raise ValueError("Sparse covariance-component variances must be positive.")
        precision = sparse.kron(model.precision / variance, identity, format="csc")
        sampling_factor = sparse.kron(
            model.precision_factor() / np.sqrt(variance), identity, format="csr"
        )
        loading = sparse.kron(model.tip_loading, identity, format="csr")
        width = precision.shape[0]
        component_slices[name] = slice(position, position + width)
        component_loadings[name] = loading
        position += width
        precisions.append(precision)
        sampling_factors.append(sampling_factor)
        loadings.append(loading)
        prior_logdet += random_dimension * (
            model.logdet_covariance + model.n_states * np.log(variance)
        )
    if not precisions:
        raise ValueError("At least one sparse covariance component is required.")
    return SparseLatentModel(
        precision=sparse.block_diag(precisions, format="csc"),
        loading=sparse.hstack(loadings, format="csr"),
        sampling_precision_factor=sparse.block_diag(sampling_factors, format="csr"),
        prior_logdet=float(prior_logdet),
        component_slices=component_slices,
        component_loadings=component_loadings,
    )


def append_identity_latent_components(
    model: SparseLatentModel,
    loadings: Mapping[str, sparse.spmatrix],
) -> SparseLatentModel:
    """Append standard-normal latent components with supplied loadings."""
    if not loadings:
        return model
    precision_blocks = [model.precision]
    sampling_blocks = (
        None
        if model.sampling_precision_factor is None
        else [model.sampling_precision_factor]
    )
    loading_blocks = [model.loading]
    slices = dict(model.component_slices)
    component_loadings = dict(model.component_loadings)
    position = model.n_states
    for name, loading_value in loadings.items():
        loading = sparse.csr_matrix(loading_value, dtype=float)
        if loading.shape[0] != model.loading.shape[0]:
            raise ValueError("Sparse latent-component loading has the wrong height.")
        width = loading.shape[1]
        if width == 0:
            continue
        precision_blocks.append(sparse.eye(width, format="csc"))
        if sampling_blocks is not None:
            sampling_blocks.append(sparse.eye(width, format="csr"))
        loading_blocks.append(loading)
        slices[name] = slice(position, position + width)
        component_loadings[name] = loading
        position += width
    return SparseLatentModel(
        precision=sparse.block_diag(precision_blocks, format="csc"),
        loading=sparse.hstack(loading_blocks, format="csr"),
        sampling_precision_factor=(
            None
            if sampling_blocks is None
            else sparse.block_diag(sampling_blocks, format="csr")
        ),
        prior_logdet=model.prior_logdet,
        component_slices=slices,
        component_loadings=component_loadings,
    )


def append_latent_component(
    model: SparseLatentModel,
    name: str,
    precision: sparse.spmatrix,
    loading: sparse.spmatrix,
    logdet_covariance: float,
    sampling_precision_factor: sparse.spmatrix | None = None,
) -> SparseLatentModel:
    """Append one general latent GMRF component."""
    loading = sparse.csr_matrix(loading, dtype=float)
    precision = sparse.csc_matrix(precision, dtype=float)
    if loading.shape != (model.loading.shape[0], precision.shape[0]) or (
        precision.shape[0] != precision.shape[1]
    ):
        raise ValueError("Appended sparse latent component has invalid dimensions.")
    position = model.n_states
    slices = dict(model.component_slices)
    slices[name] = slice(position, position + precision.shape[0])
    loadings = dict(model.component_loadings)
    loadings[name] = loading
    combined_sampling_factor = None
    if (
        model.sampling_precision_factor is not None
        and sampling_precision_factor is not None
    ):
        sampling_factor = sparse.csr_matrix(sampling_precision_factor, dtype=float)
        if sampling_factor.shape[1] != precision.shape[0]:
            raise ValueError("Appended sparse sampling factor has invalid dimensions.")
        combined_sampling_factor = sparse.block_diag(
            [model.sampling_precision_factor, sampling_factor], format="csr"
        )
    return SparseLatentModel(
        precision=sparse.block_diag([model.precision, precision], format="csc"),
        loading=sparse.hstack([model.loading, loading], format="csr"),
        sampling_precision_factor=combined_sampling_factor,
        prior_logdet=model.prior_logdet + float(logdet_covariance),
        component_slices=slices,
        component_loadings=loadings,
    )


def prepare_sparse_latent_sampler(model: SparseLatentModel) -> SparseLatentSampler:
    """Factor a sparse latent model once for repeated exact draws."""
    if model.sampling_precision_factor is None:
        raise ValueError("Sparse latent model has no exact sampling factor.")
    return SparseLatentSampler(
        loading=model.loading,
        precision_factor=model.sampling_precision_factor,
        solver=factor_sparse_positive_definite(model.precision),
    )


def continuous_predictor_loading(
    uncertainty: ContinuousPredictorUncertainty,
    slopes: np.ndarray,
) -> sparse.csr_matrix:
    """Project one continuous predictor's low-rank factor through its slopes."""
    factor_value = uncertainty.factor
    indices = np.asarray(uncertainty.observation_index, dtype=int)
    slopes = np.asarray(slopes, dtype=float).reshape(-1)
    if sparse.issparse(factor_value):
        factor = sparse.csr_matrix(factor_value, dtype=float)
        factor_rows = factor.shape[0]
        factor_valid = factor.ndim == 2 and np.isfinite(factor.data).all()
    else:
        factor = np.asarray(factor_value, dtype=float)
        factor_rows = len(factor)
        factor_valid = factor.ndim in {1, 2} and np.isfinite(factor).all()
    if (
        not factor_valid
        or indices.ndim != 1
        or np.any(indices < 0)
        or np.any(indices >= factor_rows)
        or not np.isfinite(slopes).all()
    ):
        raise ValueError("Continuous predictor factor is malformed.")
    if isinstance(factor, np.ndarray) and factor.ndim == 1:
        selected = sparse.csr_matrix(
            (
                factor[indices],
                (np.arange(len(indices), dtype=int), indices),
            ),
            shape=(len(indices), len(factor)),
        )
    else:
        selected = sparse.csr_matrix(factor[indices])
    if uncertainty.row_scale is not None:
        row_scale = np.asarray(uncertainty.row_scale, dtype=float)
        if row_scale.shape != (len(indices),) or not np.isfinite(row_scale).all():
            raise ValueError("Continuous predictor row scale is malformed.")
        selected = selected.multiply(row_scale[:, None]).tocsr()
    return sparse.kron(selected, sparse.csr_matrix(slopes[:, None]), format="csr")


def grouped_predictor_loading(
    uncertainty: GroupedPredictorUncertainty,
    slopes: np.ndarray,
) -> sparse.csr_matrix:
    """Project per-species categorical factors without gene-tip covariance."""
    indices = np.asarray(uncertainty.observation_index, dtype=int)
    slopes = np.asarray(slopes, dtype=float)
    if (
        indices.ndim != 1
        or np.any(indices < 0)
        or np.any(indices >= len(uncertainty.factors))
    ):
        raise ValueError("Grouped predictor factor indices are malformed.")
    ranks = np.asarray([factor.shape[1] for factor in uncertainty.factors], dtype=int)
    offsets = np.concatenate([[0], np.cumsum(ranks)])
    response_dimension = slopes.shape[1]
    rows = []
    columns = []
    values = []
    for observation, species_index in enumerate(indices):
        factor = np.asarray(uncertainty.factors[int(species_index)], dtype=float)
        if factor.ndim != 2 or factor.shape[0] != slopes.shape[0]:
            raise ValueError("Grouped predictor factor is malformed.")
        block = slopes.T @ factor
        for response in range(response_dimension):
            nonzero = np.flatnonzero(block[response])
            rows.extend([observation * response_dimension + response] * len(nonzero))
            columns.extend((offsets[species_index] + nonzero).tolist())
            values.extend(block[response, nonzero].tolist())
    return sparse.csr_matrix(
        (values, (rows, columns)),
        shape=(len(indices) * response_dimension, int(offsets[-1])),
    )


def joint_predictor_loading(
    uncertainty: JointPredictorUncertainty,
    slopes: np.ndarray,
) -> sparse.csr_matrix:
    """Project common-latent multicolumn uncertainty through its slopes."""
    slopes = np.asarray(slopes, dtype=float).reshape(-1)
    if len(uncertainty.factors) != len(slopes) or not len(slopes):
        raise ValueError("Joint predictor uncertainty has the wrong dimension.")
    factors = [sparse.csr_matrix(factor, dtype=float) for factor in uncertainty.factors]
    shape = factors[0].shape
    if any(factor.shape != shape for factor in factors):
        raise ValueError("Joint predictor factors must have matching dimensions.")
    loading = sum(
        (float(slope) * factor for slope, factor in zip(slopes, factors, strict=True)),
        start=sparse.csr_matrix(shape),
    ).tocsr()
    if uncertainty.row_scale is not None:
        row_scale = np.asarray(uncertainty.row_scale, dtype=float)
        if row_scale.shape != (shape[0],) or not np.isfinite(row_scale).all():
            raise ValueError("Joint predictor row scale is malformed.")
        loading = loading.multiply(row_scale[:, None]).tocsr()
    return loading


def gmrf_predictor_loading(
    uncertainty: GmrfPredictorUncertainty,
    slopes: np.ndarray,
) -> sparse.csr_matrix:
    """Lift a sparse species posterior GMRF to gene-tip linear predictors."""
    indices = np.asarray(uncertainty.observation_index, dtype=int)
    slopes = np.asarray(slopes, dtype=float).reshape(-1)
    if (
        indices.ndim != 1
        or np.any(indices < 0)
        or np.any(indices >= uncertainty.model.n_tips)
    ):
        raise ValueError("GMRF predictor indices are malformed.")
    selected = uncertainty.model.tip_loading[indices]
    return sparse.kron(selected, sparse.csr_matrix(slopes[:, None]), format="csr")


def condition_sparse_tip_model(
    prior: SparseCovarianceModel,
    prior_variance: float,
    sampling_variance: np.ndarray,
) -> SparseCovarianceModel:
    """Condition a tree GMRF on independent noisy tip observations."""
    sampling = np.asarray(sampling_variance, dtype=float)
    if (
        sampling.shape != (prior.n_tips,)
        or not np.isfinite(sampling).all()
        or np.any(sampling <= 0.0)
        or not np.isfinite(prior_variance)
        or prior_variance <= 0.0
    ):
        raise ValueError(
            "Sparse predictor conditioning requires positive diagonal sampling "
            "variance and evolutionary variance."
        )
    observation_precision = prior.tip_loading.T @ sparse.diags(1.0 / sampling)
    posterior_precision = (
        prior.precision / prior_variance + observation_precision @ prior.tip_loading
    ).tocsc()
    posterior_sampling_factor = sparse.vstack(
        [
            prior.precision_factor() / np.sqrt(prior_variance),
            sparse.diags(1.0 / np.sqrt(sampling), format="csr") @ prior.tip_loading,
        ],
        format="csr",
    )
    logdet_covariance = -factor_sparse_positive_definite(posterior_precision).logdet
    return SparseCovarianceModel(
        precision=posterior_precision,
        tip_loading=prior.tip_loading,
        logdet_covariance=logdet_covariance,
        sampling_parent=np.empty(0, dtype=int),
        sampling_transition=np.empty(0, dtype=float),
        sampling_variance=np.empty(0, dtype=float),
        sampling_precision_factor=posterior_sampling_factor,
    )


def sparse_group_covariance(labels) -> SparseCovarianceModel:
    """Represent a same-group random intercept without an n-by-n matrix."""
    group_index: dict[str, int] = {}
    indices = []
    for label in labels:
        key = str(label)
        if key not in group_index:
            group_index[key] = len(group_index)
        indices.append(group_index[key])
    rows = np.arange(len(indices), dtype=int)
    loading = sparse.csr_matrix(
        (np.ones(len(indices)), (rows, np.asarray(indices, dtype=int))),
        shape=(len(indices), len(group_index)),
    )
    n_states = len(group_index)
    return SparseCovarianceModel(
        precision=sparse.eye(n_states, format="csc"),
        tip_loading=loading,
        logdet_covariance=0.0,
        sampling_parent=np.full(n_states, -1, dtype=int),
        sampling_transition=np.zeros(n_states, dtype=float),
        sampling_variance=np.ones(n_states, dtype=float),
    )

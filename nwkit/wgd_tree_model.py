"""Fixed species-colored gene-topology DL/WGD likelihood with exact branch flow.

Rates and the candidate event are supplied, not fitted here. Gene branch
lengths are not used. The ancestral-copy prior is a unit pure-birth Yule stem,
with rate log(root_mean), hence a positive geometric count prior at the root.
"""

import math
from dataclasses import dataclass
from functools import lru_cache

import numpy as np
from scipy.special import logsumexp

from nwkit.wgd_count_model import (
    CountTree,
    MultiplicationEvent,
    _birth_death_log_probabilities,
    _validate_fixed_events,
)


@dataclass(frozen=True)
class GeneTopology:
    children: tuple[tuple[int, ...], ...]
    species: tuple[str, ...]
    node_ids: tuple[str, ...]

    def __post_init__(self):
        n = len(self.children)
        if n < 1 or len(self.species) != n or len(self.node_ids) != n:
            raise ValueError(
                "Gene topology arrays must have matching nonempty lengths."
            )
        seen: list[int] = []
        for node, children in enumerate(self.children):
            if len(children) not in {0, 2} or any(
                not 0 <= child < node for child in children
            ):
                raise ValueError("Gene topology must be binary and in postorder.")
            if (
                children
                and self.species[node]
                or not children
                and not self.species[node]
            ):
                raise ValueError("Only gene tips have species labels.")
            seen.extend(children)
        if sorted(seen) != list(range(n - 1)) or len(set(self.node_ids)) != n:
            raise ValueError(
                "Gene topology must be a tree with unique node identifiers."
            )

    @property
    def signatures(self):
        result: list[tuple] = []
        for node, children in enumerate(self.children):
            result.append(
                ("tip", self.species[node])
                if not children
                else ("node", *sorted((result[children[0]], result[children[1]])))
            )
        return result


@dataclass(frozen=True)
class TopologyResult:
    log_likelihood: float
    origin_probabilities: np.ndarray
    origin_node_ids: tuple[str, ...]
    precision_bits: int
    log_likelihood_error: float
    origin_probability_error: float


def lineage_survival(survival, duplication, loss, length):
    """Backward survival under linear DL; avoid subtracting extinction from one."""
    if length == 0 or survival == 0:
        return survival
    difference = duplication - loss
    if difference == 0:
        return survival / (1 + duplication * survival * length)
    if difference > 0:
        exponent = math.exp(-difference * length)
        denominator = (
            duplication * survival * -math.expm1(-difference * length)
            + difference * exponent
        )
        return survival * difference / denominator
    exponent = math.exp(difference * length)
    denominator = -difference + duplication * survival * -math.expm1(
        difference * length
    )
    return survival * -difference * exponent / denominator


@lru_cache(maxsize=4096)
def lineage_log_parameters(duplication, loss, length):
    """Log extinction, geometric parameter, survival and geometric success."""
    return _birth_death_log_probabilities(duplication, loss, length)


def lineage_log_complements(log_survival, log_extinction, duplication, loss, length):
    """Propagate both small complements without forming one minus the other."""
    log_a, log_b, log_survives, log_success = lineage_log_parameters(
        duplication, loss, length
    )
    denominator = float(np.logaddexp(log_success, log_b + log_survival))
    return (
        min(0.0, log_survives + log_survival - denominator),
        min(
            0.0,
            float(
                np.logaddexp(
                    log_a, log_survives + log_success + log_extinction - denominator
                )
            ),
        ),
    )


def _root_clade_masks(tree, first_child):
    root_side = [0] * len(tree.parents)
    for node, parent in enumerate(tree.parents[1:], 1):
        root_side[node] = (
            (1 if node == first_child else 2) if parent == 0 else root_side[parent]
        )
    return {
        name: root_side[node]
        for name, node in zip(tree.tip_names, tree.tip_nodes, strict=True)
    }


def _requires_multiple_root_copies(gene_tree, species_side):
    masks: list[int] = []
    for name, children in zip(gene_tree.species, gene_tree.children, strict=True):
        masks.append(
            (masks[children[0]] | masks[children[1]])
            if children
            else species_side[name]
        )
    return sorted(masks[child] for child in gene_tree.children[-1]) != [1, 2]


def _fixed_doublings(species_tree, events):
    by_node = _validate_fixed_events(species_tree, tuple(events))
    if any(event.multiplicity != 2 for event in by_node.values()):
        raise ValueError(
            "Fixed-topology likelihood requires all fixed background events to be doublings."
        )
    return by_node


def _integrate_product(coefficients, left, right, factors, dtype):
    i_values = np.flatnonzero(np.isfinite(left[:, 0]))
    j_values = np.flatnonzero(np.isfinite(right[:, 0]))
    # With a single coefficient on one side, all products target distinct
    # powers. Otherwise retain the original reduction order within a power.
    if len(j_values) == 1:
        _add_polynomial_products(
            coefficients, left, right, i_values, j_values[0], factors, dtype
        )
    else:
        for i in i_values:
            _add_polynomial_products(
                coefficients, left, right, i, j_values, factors, dtype
            )


def _add_polynomial_products(coefficients, left, right, i, j, factors, dtype):
    powers = i + j + 1
    factor = factors[powers].astype(dtype)
    coefficients[powers, 0] = np.logaddexp(
        coefficients[powers, 0], factor + left[i, 0] + right[j, 0]
    )
    if coefficients.shape[1] > 1:
        coefficients[powers, 1:] = np.logaddexp(
            coefficients[powers, 1:],
            factor[:, None]
            + np.logaddexp(
                left[i, 1:] + np.asarray(right[j, 0])[..., None],
                np.asarray(left[i, 0])[..., None] + right[j, 1:],
            ),
        )


def _integration_factors(tip_count):
    log_power = np.array([0.0, *(math.log(k) for k in range(1, tip_count))])
    return np.stack((-log_power, math.log(2) - log_power))


def _gene_splits(gene, targets):
    nodes = np.asarray(targets, dtype=np.intp)
    children = np.asarray([gene.children[node] for node in nodes], dtype=np.intp)
    return nodes, children.reshape(-1, 2).T


def _log_root_survival_normalizer(excess_mean, survival, log_survival):
    if survival >= np.finfo(float).tiny or excess_mean == 0:
        return math.log1p(excess_mean * survival)
    return math.log1p(math.exp(math.log(excess_mean) + log_survival))


def _add_species_split(result, left, right, nodes, first, second):
    result[nodes, 0] = np.logaddexp(result[nodes, 0], left[first, 0] + right[second, 0])
    if result.shape[1] > 1:
        result[nodes, 1:] = np.logaddexp(
            result[nodes, 1:],
            np.logaddexp(
                left[first, 1:] + right[second, 0, None],
                left[first, 0, None] + right[second, 1:],
            ),
        )


class TopologyLikelihood:
    def __init__(
        self,
        species_tree: CountTree,
        gene_tree: GeneTopology,
        *,
        detection=None,
        branch_groups=None,
        rate_scales=None,
        ascertainment="observed",
        max_tips=128,
        fixed_events=(),
    ):
        self.tree, self.gene = species_tree, gene_tree
        self._fixed_by_node = _fixed_doublings(species_tree, fixed_events)
        self.children = species_tree.children
        if any(len(children) not in {0, 2} for children in self.children):
            raise ValueError("Fixed-topology DL/WGD requires a binary species tree.")
        tips = [species for species in gene_tree.species if species]
        if len(tips) > max_tips or max_tips < 1:
            raise ValueError("Gene topology exceeds the supplied numerical tip limit.")
        if set(tips) - set(species_tree.tip_names):
            raise ValueError("Every gene-tree tip must map to a species-tree tip.")
        self.detection = (
            np.ones(len(species_tree.tip_names))
            if detection is None
            else np.asarray(detection, dtype=float)
        )
        if (
            self.detection.shape != (len(species_tree.tip_names),)
            or np.any(~np.isfinite(self.detection))
            or np.any(self.detection < 0)
            or np.any(self.detection > 1)
        ):
            raise ValueError("Known detection probabilities must be in [0,1].")
        if any(
            name in tips and self.detection[i] == 0
            for i, name in enumerate(species_tree.tip_names)
        ):
            raise ValueError("Observed genes cannot have zero detection probability.")
        self.groups = (
            tuple(0 for _ in species_tree.parents)
            if branch_groups is None
            else tuple(branch_groups)
        )
        if len(self.groups) != len(species_tree.parents) or any(
            not isinstance(group, int) or group < 0 for group in self.groups
        ):
            raise ValueError(
                "Branch-rate group indices must be nonnegative and cover nodes."
            )
        self.num_groups = max(self.groups[1:]) + 1
        self.scales = (
            np.ones(1) if rate_scales is None else np.asarray(rate_scales, dtype=float)
        )
        if (
            self.scales.ndim != 1
            or not len(self.scales)
            or np.any(~np.isfinite(self.scales))
            or np.any(self.scales <= 0)
        ):
            raise ValueError("Family-rate scales must be finite and positive.")
        species_side = _root_clade_masks(self.tree, self.children[0][0])
        both_root_clades = len({species_side[tip] for tip in tips}) == 2
        if ascertainment not in {"observed", "root-clades"}:
            raise ValueError("Unknown topology ascertainment rule.")
        if ascertainment == "root-clades" and not both_root_clades:
            raise ValueError(
                "Topology does not satisfy root-clade observation selection."
            )
        self.ascertainment = ascertainment
        self.requires_multiple_root_copies = (
            both_root_clades and _requires_multiple_root_copies(gene_tree, species_side)
        )
        signatures = gene_tree.signatures
        self.coefficients = np.array(
            [
                0
                if not children
                else 1
                if signatures[children[0]] == signatures[children[1]]
                else 2
                for children in gene_tree.children
            ]
        )
        self.targets = tuple(
            node for node, children in enumerate(gene_tree.children) if children
        )
        self.split_nodes, self.split_children = _gene_splits(gene_tree, self.targets)
        self.tip_counts: list[int] = []
        for children in gene_tree.children:
            self.tip_counts.append(
                1 if not children else sum(self.tip_counts[child] for child in children)
            )
        self.integration_factors = _integration_factors(len(tips))

    def _edge(self, survival, logs, duplication, loss, length, dtype, log_survival):
        if not length:
            return survival, logs.copy()
        _, log_b, log_survives, log_success = lineage_log_parameters(
            duplication, loss, length
        )
        denominator = float(np.logaddexp(log_success, log_b + log_survival))
        new_survival = math.exp(min(0.0, log_survives + log_survival - denominator))
        log_h = log_survives + log_success - 2 * denominator
        if duplication == 0:
            return new_survival, logs.astype(dtype) + log_h
        if np.all(np.isneginf(logs[:, 0])):
            # Duplication multiplies existing subtree probabilities; an edge
            # with no possible observed descendants cannot create any.
            return new_survival, logs.copy()
        log_z = log_b - denominator
        columns = logs.shape[1]
        polynomials: dict[int, np.ndarray] = {}
        result = np.full(logs.shape, -np.inf, dtype=dtype)
        # L_g = H F_g(z), where F'_g = c_g F_left F_right. Every F is
        # a positive polynomial of degree at most the observed tip count - 1.
        for node, children in enumerate(self.gene.children):
            if not children:
                coefficients = logs[node : node + 1].astype(dtype, copy=False)
                polynomials[node] = coefficients
                result[node] = dtype(log_h) + coefficients[0]
                continue
            coefficients = np.full(
                (self.tip_counts[node], columns), -np.inf, dtype=dtype
            )
            coefficients[0] = logs[node]
            left, right = (polynomials.pop(child) for child in children)
            factors = self.integration_factors[self.coefficients[node] - 1]
            _integrate_product(coefficients, left, right, factors, dtype)
            polynomials[node] = coefficients
            powers = np.arange(len(coefficients), dtype=dtype) * dtype(log_z)
            powers[0] = 0
            result[node] = dtype(log_h) + logsumexp(
                coefficients + powers[:, None], axis=0
            )
        return new_survival, result

    def _jump(self, survival, logs, event, log_extinction, *, mark=True):
        q = event.retention
        if q == 0:
            return survival, logs.copy()
        result = np.full_like(logs, -np.inf)
        log_carry = float(
            np.logaddexp(
                math.log1p(-q) if q < 1 else -np.inf, math.log(2 * q) + log_extinction
            )
        )
        if np.isfinite(log_carry):
            result = logs + log_carry
        nodes = self.split_nodes
        left, right = self.split_children
        factor_by_multiplicity = np.array([math.log(q), math.log(2 * q)])
        factor = factor_by_multiplicity[self.coefficients[nodes] - 1]
        split = factor + logs[left, 0] + logs[right, 0]
        result[nodes, 0] = np.logaddexp(result[nodes, 0], split)
        if logs.shape[1] > 1:
            terms = np.stack(
                (
                    result[nodes, 1:],
                    factor[:, None] + logs[left, 1:] + logs[right, 0, None],
                    factor[:, None] + logs[left, 0, None] + logs[right, 1:],
                ),
                axis=-1,
            )
            markers = logsumexp(terms, axis=-1)
            # Only the marker for this node gains a new origin at the pulse.
            if mark:
                diagonal = np.arange(len(nodes))
                markers[diagonal, diagonal] = logsumexp(
                    np.column_stack((terms[diagonal, diagonal], split)), axis=1
                )
            result[nodes, 1:] = markers
        return survival * (1 + q - q * survival), result

    def _category(self, rates, root_mean, event, scale, dtype, *, origins=True):
        nodes = len(self.tree.parents)
        survivors = np.zeros(nodes)
        log_survivors = np.full(nodes, -np.inf)
        log_extinctions = np.zeros(nodes)
        columns = 1 + len(self.targets) if event is not None and origins else 1
        probabilities: dict[int, np.ndarray] = {}
        tips = dict(
            zip(self.tree.tip_nodes, range(len(self.tree.tip_nodes)), strict=True)
        )
        for node in reversed(range(nodes)):
            if node in tips:
                i = tips[node]
                probabilities[node] = np.full(
                    (len(self.gene.children), columns), -np.inf
                )
                survivors[node] = self.detection[i]
                log_survivors[node] = (
                    math.log(self.detection[i]) if self.detection[i] else -np.inf
                )
                log_extinctions[node] = (
                    math.log1p(-self.detection[i]) if self.detection[i] < 1 else -np.inf
                )
                for gene, name in enumerate(self.gene.species):
                    if name == self.tree.tip_names[i]:
                        probabilities[node][gene, 0] = math.log(self.detection[i])
            else:
                left, right = self.children[node]
                log_survivors[node] = np.logaddexp(
                    log_survivors[left], log_extinctions[left] + log_survivors[right]
                )
                log_extinctions[node] = log_extinctions[left] + log_extinctions[right]
                survivors[node] = math.exp(log_survivors[node])
                probabilities[node] = np.logaddexp(
                    probabilities[left] + log_extinctions[right],
                    probabilities[right] + log_extinctions[left],
                )
                first, second = self.split_children
                for a_species, b_species in ((left, right), (right, left)):
                    _add_species_split(
                        probabilities[node],
                        probabilities[a_species],
                        probabilities[b_species],
                        self.split_nodes,
                        first,
                        second,
                    )
                del probabilities[left], probabilities[right]
            if node:
                duplication, loss = rates[self.groups[node]] * scale
                length = self.tree.lengths[node]
                pulse = (
                    event
                    if event is not None and node == event.node
                    else self._fixed_by_node.get(node)
                )
                if pulse is not None:
                    survivors[node], probabilities[node] = self._edge(
                        survivors[node],
                        probabilities[node],
                        duplication,
                        loss,
                        length * (1 - pulse.fraction),
                        dtype,
                        log_survivors[node],
                    )
                    log_survivors[node], log_extinctions[node] = (
                        lineage_log_complements(
                            log_survivors[node],
                            log_extinctions[node],
                            duplication,
                            loss,
                            length * (1 - pulse.fraction),
                        )
                    )
                    survivors[node], probabilities[node] = self._jump(
                        survivors[node],
                        probabilities[node],
                        pulse,
                        log_extinctions[node],
                        mark=pulse is event,
                    )
                    q = pulse.retention
                    old_extinction = log_extinctions[node]
                    log_survivors[node] = min(
                        0.0,
                        log_survivors[node] + math.log1p(q * math.exp(old_extinction)),
                    )
                    log_extinctions[node] += np.logaddexp(
                        math.log1p(-q) if q < 1 else -np.inf,
                        math.log(q) + old_extinction if q else -np.inf,
                    )
                    length *= pulse.fraction
                survivors[node], probabilities[node] = self._edge(
                    survivors[node],
                    probabilities[node],
                    duplication,
                    loss,
                    length,
                    dtype,
                    log_survivors[node],
                )
                log_survivors[node], log_extinctions[node] = lineage_log_complements(
                    log_survivors[node],
                    log_extinctions[node],
                    duplication,
                    loss,
                    length,
                )
                survivors[node] = math.exp(log_survivors[node])
        root_survival = survivors[0]
        r = root_mean - 1
        root_normalizer = _log_root_survival_normalizer(
            r, root_survival, log_survivors[0]
        )
        if self.ascertainment == "observed":
            selection_log = math.log(root_mean) + log_survivors[0] - root_normalizer
        else:
            left, right = self.children[0]
            with np.errstate(divide="ignore"):
                log_r = np.log(r)
            # Keep 1 + 2*r + r**2 * P(either clade) in range even when r**2
            # overflows or the single-copy survival underflows.
            log_factor = np.logaddexp(
                math.log1p(r), log_r + np.logaddexp(0, log_r + log_survivors[0])
            )
            selection_log = (
                math.log(root_mean)
                + log_survivors[left]
                + log_survivors[right]
                + log_factor
                - _log_root_survival_normalizer(r, survivors[left], log_survivors[left])
                - _log_root_survival_normalizer(
                    r, survivors[right], log_survivors[right]
                )
                - root_normalizer
            )
        if not np.isfinite(selection_log):
            raise ValueError(
                "Topology observation selection is unrepresentable numerically."
            )
        _, root = self._edge(
            survivors[0],
            probabilities[0],
            math.log(root_mean),
            0,
            1,
            dtype,
            log_survivors[0],
        )
        return root[-1], selection_log

    def _evaluate(self, rates, root_mean, event, dtype, *, origins=True):
        categories = [
            self._category(rates, root_mean, event, scale, dtype, origins=origins)
            for scale in self.scales
        ]
        likelihoods = np.array([item[0][0] for item in categories])
        selections = np.array([item[1] for item in categories])
        raw = logsumexp(likelihoods)
        if not np.isfinite(raw):
            raise ValueError(
                "Topology has zero or unrepresentable probability under this model."
            )
        likelihood = float(raw - logsumexp(selections))
        probabilities = (
            np.zeros(len(self.targets) if origins else 0)
            if event is None or not origins
            else np.exp(
                logsumexp(np.array([item[0][1:] for item in categories]), axis=0) - raw
            )
        )
        if np.any(probabilities > 1 + 1e-7) or np.any(~np.isfinite(probabilities)):
            raise ValueError(
                "Conditional origin probabilities exceed numerical bounds."
            )
        return likelihood, np.minimum(probabilities, 1)

    def evaluate(
        self,
        rates,
        root_mean,
        event: MultiplicationEvent | None = None,
        *,
        tolerance=1e-6,
        origin_tolerance=1e-5,
        origins=True,
    ):
        """Evaluate one fixed topology, optionally omitting all origin markers."""
        rates = np.asarray(rates, dtype=float)
        if (
            rates.shape != (self.num_groups, 2)
            or np.any(~np.isfinite(rates))
            or np.any(rates < 0)
        ):
            raise ValueError(
                "Topology rates must be finite nonnegative group-by-2 values."
            )
        if not np.isfinite(root_mean) or root_mean < 1:
            raise ValueError("Topology root mean must be finite and >=1.")
        if root_mean == 1 and self.requires_multiple_root_copies:
            raise ValueError(
                "Topology crosses the species-root split and is impossible with one ancestral root copy."
            )
        if event is not None and (
            event.multiplicity != 2 or not 0 < event.node < len(self.tree.parents)
        ):
            raise ValueError(
                "Fixed-topology model supports one doubling on a non-root branch."
            )
        if event is not None and event.node in self._fixed_by_node:
            raise ValueError(
                "A topology candidate cannot overlap a fixed genome event."
            )
        if any(
            not np.isfinite(value) or value <= 0
            for value in (tolerance, origin_tolerance)
        ):
            raise ValueError(
                "Topology numerical tolerances must be finite and positive."
            )
        # At zero retention the pulse is the identity in the unbounded
        # topology model. There is no marked origin, including when other
        # fixed background doublings are present.
        if event is not None and event.retention == 0:
            event = None
        first = self._evaluate(rates, root_mean, event, np.float64, origins=origins)
        second = self._evaluate(rates, root_mean, event, np.longdouble, origins=origins)
        error = abs(second[0] - first[0])
        origin_error = float(np.max(abs(second[1] - first[1]))) if len(second[1]) else 0
        if error > tolerance or origin_error > origin_tolerance:
            raise ValueError(
                "Fixed-topology likelihood or origin assignment failed arithmetic precision check."
            )
        return TopologyResult(
            second[0],
            np.asarray(second[1], dtype=float),
            tuple(self.gene.node_ids[node] for node in self.targets) if origins else (),
            np.finfo(np.longdouble).nmant + 1,
            error,
            origin_error,
        )

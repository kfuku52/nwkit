"""Fixed species-colored gene-topology DL/WGD likelihood with exact branch flow.

Rates and the candidate event are supplied, not fitted here. Gene branch
lengths are not used. The ancestral-copy prior is a unit pure-birth Yule stem,
with rate log(root_mean), hence a positive geometric count prior at the root.
"""

import math
from dataclasses import dataclass

import numpy as np
from scipy.special import logsumexp

from nwkit.wgd_count_model import CountTree, MultiplicationEvent


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


def lineage_log_parameters(duplication, loss, length):
    """Log extinction, geometric parameter, survival and geometric success."""
    if length == 0 or duplication == loss == 0:
        return -np.inf, -np.inf, 0.0, 0.0
    difference = abs(duplication - loss)
    if difference == 0:
        log_product = math.log(duplication) + math.log(length)
        denominator = float(np.logaddexp(0, log_product))
        log_a = log_b = log_product - denominator
        log_survives = log_success = -denominator
    else:
        change = -math.expm1(-difference * length)
        log_change = (
            math.log(change) if change else math.log(difference) + math.log(length)
        )
        denominator = float(
            np.logaddexp(
                math.log(difference),
                math.log(min(duplication, loss)) + log_change
                if min(duplication, loss)
                else -np.inf,
            )
        )
        complement = math.log(difference) - denominator
        log_a = math.log(loss) + log_change - denominator if loss else -np.inf
        log_b = (
            math.log(duplication) + log_change - denominator if duplication else -np.inf
        )
        log_survives = complement - (difference * length if duplication < loss else 0)
        log_success = complement - (difference * length if duplication > loss else 0)
    return log_a, log_b, log_survives, log_success


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
    ):
        self.tree, self.gene = species_tree, gene_tree
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
        self._validate_ascertainment(ascertainment, tips)
        self.ascertainment = ascertainment
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
        self.tip_counts: list[int] = []
        for children in gene_tree.children:
            self.tip_counts.append(
                1 if not children else sum(self.tip_counts[child] for child in children)
            )

    def _validate_ascertainment(self, ascertainment, tips):
        if ascertainment not in {"observed", "root-clades"}:
            raise ValueError("Unknown topology ascertainment rule.")
        if ascertainment == "root-clades":
            descendant_species: list[set[str]] = [set() for _ in self.tree.parents]
            for node, name in zip(
                self.tree.tip_nodes, self.tree.tip_names, strict=True
            ):
                descendant_species[node].add(name)
            for node in reversed(range(len(self.tree.parents))):
                for child in self.children[node]:
                    descendant_species[node].update(descendant_species[child])
            if any(
                not (descendant_species[child] & set(tips))
                for child in self.children[0]
            ):
                raise ValueError(
                    "Topology does not satisfy root-clade observation selection."
                )

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
        log_z = log_b - denominator
        columns = logs.shape[1]
        polynomials: list[np.ndarray] = []
        result = np.full(logs.shape, -np.inf, dtype=dtype)
        # L_g = H F_g(z), where F'_g = c_g F_left F_right. Every F is
        # a positive polynomial of degree at most the observed tip count - 1.
        for node, children in enumerate(self.gene.children):
            coefficients = np.full(
                (self.tip_counts[node], columns), -np.inf, dtype=dtype
            )
            coefficients[0] = logs[node]
            if children:
                left, right = (polynomials[child] for child in children)
                for i in np.flatnonzero(np.isfinite(left[:, 0])):
                    for j in np.flatnonzero(np.isfinite(right[:, 0])):
                        power = i + j + 1
                        factor = dtype(
                            math.log(self.coefficients[node]) - math.log(power)
                        )
                        coefficients[power, 0] = np.logaddexp(
                            coefficients[power, 0], factor + left[i, 0] + right[j, 0]
                        )
                        coefficients[power, 1:] = np.logaddexp(
                            coefficients[power, 1:],
                            factor
                            + np.logaddexp(
                                left[i, 1:] + right[j, 0], left[i, 0] + right[j, 1:]
                            ),
                        )
            polynomials.append(coefficients)
            powers = np.arange(len(coefficients), dtype=dtype) * dtype(log_z)
            powers[0] = 0
            result[node] = dtype(log_h) + logsumexp(
                coefficients + powers[:, None], axis=0
            )
        return new_survival, result

    def _jump(self, survival, logs, event, log_extinction):
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
        for node, children in enumerate(self.gene.children):
            if not children:
                continue
            left, right = children
            factor = math.log(q * self.coefficients[node])
            split = factor + logs[left, 0] + logs[right, 0]
            result[node, 0] = np.logaddexp(result[node, 0], split)
            for column, target in enumerate(self.targets, 1):
                terms = [
                    result[node, column],
                    factor + logs[left, column] + logs[right, 0],
                    factor + logs[left, 0] + logs[right, column],
                ]
                if node == target:
                    terms.append(split)
                result[node, column] = logsumexp(terms)
        return survival * (1 + q - q * survival), result

    def _category(self, rates, root_mean, event, scale, dtype):
        nodes = len(self.tree.parents)
        survivors = np.zeros(nodes)
        log_survivors = np.full(nodes, -np.inf)
        log_extinctions = np.zeros(nodes)
        columns = 1 + len(self.targets) if event is not None else 1
        probabilities = [
            np.full((len(self.gene.children), columns), -np.inf) for _ in range(nodes)
        ]
        tips = dict(
            zip(self.tree.tip_nodes, range(len(self.tree.tip_nodes)), strict=True)
        )
        for node in reversed(range(nodes)):
            if node in tips:
                i = tips[node]
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
                for gene, children in enumerate(self.gene.children):
                    if not children:
                        continue
                    first, second = children
                    for a_species, b_species in ((left, right), (right, left)):
                        a_logs, b_logs = (
                            probabilities[a_species],
                            probabilities[b_species],
                        )
                        probabilities[node][gene, 0] = np.logaddexp(
                            probabilities[node][gene, 0],
                            a_logs[first, 0] + b_logs[second, 0],
                        )
                        probabilities[node][gene, 1:] = np.logaddexp(
                            probabilities[node][gene, 1:],
                            np.logaddexp(
                                a_logs[first, 1:] + b_logs[second, 0],
                                a_logs[first, 0] + b_logs[second, 1:],
                            ),
                        )
            if node:
                duplication, loss = rates[self.groups[node]] * scale
                length = self.tree.lengths[node]
                if event is not None and node == event.node:
                    survivors[node], probabilities[node] = self._edge(
                        survivors[node],
                        probabilities[node],
                        duplication,
                        loss,
                        length * (1 - event.fraction),
                        dtype,
                        log_survivors[node],
                    )
                    log_survivors[node], log_extinctions[node] = (
                        lineage_log_complements(
                            log_survivors[node],
                            log_extinctions[node],
                            duplication,
                            loss,
                            length * (1 - event.fraction),
                        )
                    )
                    survivors[node], probabilities[node] = self._jump(
                        survivors[node],
                        probabilities[node],
                        event,
                        log_extinctions[node],
                    )
                    q = event.retention
                    old_extinction = log_extinctions[node]
                    log_survivors[node] = min(
                        0.0,
                        log_survivors[node] + math.log1p(q * math.exp(old_extinction)),
                    )
                    log_extinctions[node] += np.logaddexp(
                        math.log1p(-q) if q < 1 else -np.inf,
                        math.log(q) + old_extinction if q else -np.inf,
                    )
                    length *= event.fraction
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
        if self.ascertainment == "observed":
            selection_log = (
                math.log(root_mean) + log_survivors[0] - math.log1p(r * root_survival)
            )
        else:
            left, right = self.children[0]
            selection_log = (
                math.log(root_mean)
                + log_survivors[left]
                + log_survivors[right]
                + math.log1p(2 * r + r * r * root_survival)
                - math.log1p(r * survivors[left])
                - math.log1p(r * survivors[right])
                - math.log1p(r * root_survival)
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

    def _evaluate(self, rates, root_mean, event, dtype):
        categories = [
            self._category(rates, root_mean, event, scale, dtype)
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
            np.zeros(len(self.targets))
            if event is None
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
    ):
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
        if event is not None and (
            event.multiplicity != 2 or not 0 < event.node < len(self.tree.parents)
        ):
            raise ValueError(
                "Fixed-topology model supports one doubling on a non-root branch."
            )
        if any(
            not np.isfinite(value) or value <= 0
            for value in (tolerance, origin_tolerance)
        ):
            raise ValueError(
                "Topology numerical tolerances must be finite and positive."
            )
        first = self._evaluate(rates, root_mean, event, np.float64)
        second = self._evaluate(rates, root_mean, event, np.longdouble)
        error = abs(second[0] - first[0])
        origin_error = (
            float(np.max(abs(second[1] - first[1]))) if len(self.targets) else 0
        )
        if error > tolerance or origin_error > origin_tolerance:
            raise ValueError(
                "Fixed-topology likelihood or origin assignment failed arithmetic precision check."
            )
        return TopologyResult(
            second[0],
            np.asarray(second[1], dtype=float),
            tuple(self.gene.node_ids[node] for node in self.targets),
            np.finfo(np.longdouble).nmant + 1,
            error,
            origin_error,
        )

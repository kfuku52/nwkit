"""Independent topology-only conditional reference, using dense CTMC exponentials.

No production coalescent transitions, count messages, samplers or detection
functions are imported. This does not simulate branch lengths or sequences.
"""

import math
from collections import defaultdict

import numpy as np
from scipy.linalg import expm


class ConditionalReference:
    def __init__(self, root, ne, max_work=100000):
        if (
            isinstance(ne, bool)
            or not math.isfinite(ne)
            or ne <= 0
            or not math.isfinite(2 * ne)
        ):
            raise ValueError("Independent reference requires finite positive Ne.")
        if type(max_work) is not int or max_work < 1:
            raise ValueError(
                "Independent reference work cap must be a positive integer."
            )
        self.root, self.inside, self.outside, self.matrices = root, {}, {}, {}
        self.order, self.work, self.limit = [], 0, max_work
        pending = [(root, None, False)]
        while pending:
            node, parent, visited = pending.pop()
            if not visited:
                pending.append((node, parent, True))
                pending.extend((child, node, False) for child in node.children)
                continue
            self.order.append(node)
            counts = {int(node.kind == "tip"): 1.0}
            if node.children:
                counts = {0: 1.0}
                for child in node.children:
                    combined = defaultdict(float)
                    for k, p in counts.items():
                        for j, q in self.outside[child].items():
                            self.charge(1)
                            combined[k + j] += p * q
                    counts = dict(combined)
            self.inside[node] = counts
            n = max(counts)
            duration = 0 if parent is None else (parent.age - node.age) / (2 * ne)
            if not math.isfinite(duration) or duration < 0:
                raise ValueError("Invalid reference branch duration.")
            self.charge((n + 1) ** 3)
            generator = np.zeros((n + 1, n + 1))
            for k in range(2, n + 1):
                generator[k, k - 1] = k * (k - 1) / 2
                generator[k, k] = -generator[k, k - 1]
            # Pure-death structural zeros are exact; expm can leak roundoff above
            # the diagonal. Never clip a negative feasible transition.
            matrix = np.tril(expm(generator * duration))
            matrix[1:, 0] = 0
            if np.any(~np.isfinite(matrix)) or np.any(matrix < 0):
                raise ArithmeticError("Invalid reference CTMC transition matrix.")
            self.matrices[node] = matrix
            outputs = defaultdict(float)
            for k, p in counts.items():
                for j in range(k + 1):
                    if matrix[k, j] > 0:
                        outputs[j] += p * matrix[k, j]
            if node.daughter and n:
                if outputs[1] <= 0:
                    raise ArithmeticError("Reference daughter probability underflowed.")
                outputs = {1: 1.0}
            self.outside[node] = dict(outputs)

    def charge(self, work):
        self.work += work
        if self.work > self.limit:
            raise ValueError("Independent conditional reference exceeded work cap.")

    @staticmethod
    def choose(weights, rng):
        total = math.fsum(weights.values())
        if not math.isfinite(total) or total <= 0:
            raise ArithmeticError("Invalid independent conditional weights.")
        target = rng.random() * total
        partial = 0.0
        last = None
        for key, weight in weights.items():
            if weight > 0:
                last = key
                partial += weight
                if target < partial:
                    return key
        return last

    def sample(self, rng):
        requested = {self.root: self.choose(self.inside[self.root], rng)}
        for node in reversed(self.order):
            end = requested[node]
            start = self.choose(
                {
                    k: p * self.matrices[node][k, end]
                    for k, p in self.inside[node].items()
                },
                rng,
            )
            if len(node.children) == 1:
                requested[node.children[0]] = start
            elif len(node.children) == 2:
                a, b = node.children
                left = self.choose(
                    {
                        k: p * self.outside[b].get(start - k, 0)
                        for k, p in self.outside[a].items()
                    },
                    rng,
                )
                requested[a], requested[b] = left, start - left
            elif node.children:
                raise ValueError("Reference requires binary locus histories.")
        forests = {}
        for i, node in enumerate(self.order):
            forest = [gene for child in node.children for gene in forests.pop(child)]
            if node.kind == "tip":
                forest = [("tip", node.species, i)]
            end = 1 if node is self.root and forest else requested[node]
            while len(forest) > end:
                pairs = [(a, b) for a in range(len(forest)) for b in range(a)]
                a, b = pairs[int(rng.integers(len(pairs)))]
                pair = sorted((forest.pop(a), forest.pop(b)))
                forest.append(("node", *pair))
            forests[node] = forest
        return forests[self.root][0] if forests[self.root] else None


def detect(gene, probabilities, rng):
    if gene is None:
        return None
    if gene[0] == "tip":
        return ("tip", gene[1]) if rng.random() < probabilities[gene[1]] else None
    left, right = (detect(child, probabilities, rng) for child in gene[1:])
    if left is None or right is None:
        return right if left is None else left
    return ("node", *sorted((left, right)))


def tips(gene):
    return (
        0 if gene is None else 1 if gene[0] == "tip" else sum(tips(c) for c in gene[1:])
    )


def reference_selected(population, point, config, rng, *, count):
    from pilot import reference_locus

    result, attempts = [], 0
    while len(result) < count:
        attempts += 1
        if attempts > config["max_attempts"]:
            raise ValueError("Independent conditional selection exceeded attempt cap.")
        locus = reference_locus(population, point, config, rng)
        sampler = ConditionalReference(locus, point.ne, config["max_coalescent_states"])
        observation = detect(sampler.sample(rng), config["detection"], rng)
        if 2 <= tips(observation) <= config["max_observed_tips"]:
            result.append(observation)
    return result, attempts

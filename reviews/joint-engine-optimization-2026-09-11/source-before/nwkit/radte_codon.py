"""Standard-code GY94 and Kosiol et al. (2007) empirical codon models."""

from functools import lru_cache
from importlib.resources import files
from itertools import product

import numpy as np

# NCBI translation table 1, enumerated in T,C,A,G order.
_CODE = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
TRANSLATION = dict(zip(map("".join, product("TCAG", repeat=3)), _CODE, strict=True))
CODONS = tuple(c for c, aa in TRANSLATION.items() if aa != "*")
CODON_MODELS = {"gy94", "ecmk07", "ecmrest"}
DISTANCES = np.array(
    [[sum(x != y for x, y in zip(a, b, strict=True)) for b in CODONS] for a in CODONS]
)


def encode_codons(seqs, names, codes):
    length = len(seqs[names[0]])
    if length % 3:
        raise ValueError("Codon alignment length must be divisible by three.")
    matrix = np.empty((len(names), length // 3), dtype=np.uint64)
    ids = {c: i for i, c in enumerate(CODONS)}
    for i, name in enumerate(names):
        for j in range(length // 3):
            token = seqs[name][3 * j : 3 * j + 3]
            if "-" in token and token != "---":
                raise ValueError(
                    f"Partial codon gap in {name} at codon {j + 1}: {token}"
                )
            if any(c not in codes for c in token):
                raise ValueError(
                    f"Unsupported codon in {name} at codon {j + 1}: {token}"
                )
            possible = map("".join, product(*(codes[c] for c in token)))
            mask = sum(1 << ids[c] for c in possible if c in ids)
            if not mask:
                raise ValueError(f"Stop codon in {name} at codon {j + 1}: {token}")
            matrix[i, j] = mask
    return matrix, CODONS


@lru_cache(maxsize=2)
def empirical_matrix(model):
    name = "ECMunrest" if model == "ecmk07" else "ECMrest"
    text = files("nwkit").joinpath(f"data_model/{name}.txt").read_text(encoding="utf-8")
    tokens = " ".join(
        line for line in text.splitlines() if not line.startswith("#")
    ).split()
    n = len(CODONS)
    end = n * (n - 1) // 2
    exchange = np.zeros((n, n))
    exchange[np.tril_indices(n, -1)] = np.asarray(tokens[:end], dtype=float)
    exchange += exchange.T
    frequencies = np.asarray(tokens[end : end + n], dtype=float)
    order = tokens[end + n : end + 2 * n]
    if set(order) != set(CODONS) or len(set(order)) != n:
        raise ValueError("Invalid packaged ECM codon order.")
    indices = [order.index(c) for c in CODONS]
    return exchange[np.ix_(indices, indices)], frequencies[indices] / frequencies.sum()


def codon_matrix(model, patterns, weights, kappa, omega, frequency):
    counts = np.array(
        [
            np.sum(weights * np.sum(patterns == np.uint64(1 << i), axis=0))
            for i in range(len(CODONS))
        ],
        dtype=float,
    )
    if frequency == "f":
        pi = counts + 0.5
    elif frequency in {"f1x4", "f3x4"}:
        # Count known bases even when another position in the codon is ambiguous.
        nt = np.full((3, 4), 0.5)
        for pos in range(3):
            for base_index, base in enumerate("ACGT"):
                mask = np.uint64(
                    sum(1 << i for i, c in enumerate(CODONS) if c[pos] == base)
                )
                known = (patterns & mask) == patterns
                nt[pos, base_index] += np.sum(weights * np.sum(known, axis=0))
        if frequency == "f1x4":
            nt[:] = nt.sum(axis=0)
        nt /= nt.sum(axis=1, keepdims=True)
        pi = np.array(
            [np.prod([nt[p, "ACGT".index(b)] for p, b in enumerate(c)]) for c in CODONS]
        )
    elif frequency == "fq":
        pi = np.ones(len(CODONS))
    elif frequency == "model" and model != "gy94":
        pi = empirical_matrix(model)[1].copy()
    else:
        raise ValueError(
            "Invalid codon frequencies; GY94 requires f, f1x4, f3x4 or fq."
        )
    pi /= pi.sum()
    if model == "gy94":
        if not np.isfinite([kappa, omega]).all() or min(kappa, omega) <= 0:
            raise ValueError("GY94 kappa and omega must be finite and positive.")
        exchange = np.zeros((len(CODONS), len(CODONS)))
        for i, a in enumerate(CODONS):
            for j, b in enumerate(CODONS):
                if DISTANCES[i, j] != 1:
                    continue
                x, y = next((x, y) for x, y in zip(a, b, strict=True) if x != y)
                transition = {x, y} in ({"A", "G"}, {"C", "T"})
                exchange[i, j] = (kappa if transition else 1) * (
                    omega if TRANSLATION[a] != TRANSLATION[b] else 1
                )
    else:
        exchange = empirical_matrix(model)[0]
    q = exchange * pi[None, :]
    np.fill_diagonal(q, 0)
    # Branch length = expected nucleotide changes per codon site. A multi-base
    # ECM transition contributes its Hamming distance, not one event.
    q /= np.sum(pi[:, None] * q * DISTANCES)
    np.fill_diagonal(q, -q.sum(axis=1))
    return q, pi

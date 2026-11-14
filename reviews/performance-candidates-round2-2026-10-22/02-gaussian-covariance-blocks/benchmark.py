"""LCA-descendant-block covariance prototype; production sources are untouched."""

import argparse
import hashlib
import inspect
import json
import sys
import textwrap
from contextlib import contextmanager
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from support import environment, measurement_lock, paired, peak_bytes

import nwkit
import nwkit.gaussian_tree as gt
from nwkit.evolution import build_evolutionary_process
from nwkit.ordinary_regression import _fit_ordinary_gaussian
from nwkit.util import read_tree

ORIGINAL = gt.GaussianTreeProcess.covariance
WORKSPACE_BYTES = 4 * 1024**2


def _fill_blocks(covariance, lca, requested, variances, loading, factors, stable):
    """Every requested pair is written at its LCA, with bounded product tiles."""
    local = [[] for _ in lca.nodes]
    for output, index in enumerate(requested):
        local[index].append(output)
    descendants = {}
    output_loading = loading[np.asarray(requested, dtype=int)]
    for ancestor in range(len(lca.nodes) - 1, -1, -1):
        own = np.asarray(local[ancestor], dtype=int)
        groups = [own] if len(own) else []
        for child in lca.nodes[ancestor].children:
            group = descendants.pop(lca.index_by_node[child])
            if len(group):
                groups.append(group)
        merged = np.empty(sum(map(len, groups)), dtype=int)
        if len(own):
            own_value = variances[ancestor]
            if factors is not None:
                # Repeated requests for one node are off-diagonal pairs in
                # the original API: retain its loading-product arithmetic.
                with np.errstate(over="ignore", under="ignore", invalid="ignore"):
                    own_value = (
                        loading[ancestor] * loading[ancestor] * factors[ancestor]
                    )
                if not np.isfinite(own_value) or (
                    own_value == 0 and variances[ancestor] > 0
                ):
                    own_value = stable(ancestor, ancestor, ancestor)
            covariance[np.ix_(own, own)] = own_value
            covariance[own, own] = variances[ancestor]
        start = 0
        for group in groups:
            if start:
                # Two product temporaries may coexist; tile both axes.
                # Include the finite/zero masks as well as the two products.
                max_entries = max(1, WORKSPACE_BYTES // 24)
                columns = min(len(group), max_entries)
                rows = max(1, min(start, max_entries // columns))
                for second in range(0, len(group), columns):
                    right = group[second : second + columns]
                    for first in range(0, start, rows):
                        left = merged[first : min(start, first + rows)]
                        if factors is None:
                            covariance[np.ix_(left, right)] = variances[ancestor]
                            covariance[np.ix_(right, left)] = variances[ancestor]
                            continue
                        with np.errstate(
                            over="ignore", under="ignore", invalid="ignore"
                        ):
                            block = (
                                output_loading[left, None]
                                * output_loading[right][None, :]
                                * factors[ancestor]
                            )
                        bad = ~np.isfinite(block) | (
                            (block == 0) & (variances[ancestor] > 0)
                        )
                        for row, column in zip(*np.nonzero(bad), strict=True):
                            one, two = sorted((int(left[row]), int(right[column])))
                            block[row, column] = stable(
                                requested[one], requested[two], ancestor
                            )
                        covariance[np.ix_(left, right)] = block
                        covariance[np.ix_(right, left)] = block.T
                        del block, bad
            merged[start : start + len(group)] = group
            start += len(group)
        descendants[ancestor] = merged


def _candidate_function(minimum_nodes=64):
    source = textwrap.dedent(inspect.getsource(ORIGINAL))
    source = source.replace("def covariance(", "def candidate_covariance(", 1)
    source = source.replace(
        "    nodes = tuple(nodes)\n",
        f"    nodes = tuple(nodes)\n    if len(nodes) < {minimum_nodes}:\n        return ORIGINAL(self, nodes)\n",
        1,
    )
    begin = source.index("    covariance = np.empty")
    end = source.index("    if not np.isfinite(covariance).all():", begin)
    replacement = """    covariance = np.empty((len(nodes), len(nodes)), dtype=float)
    factors = None
    if not unit_loadings:
        if not direct_loadings:
            return ORIGINAL(self, nodes)
        with np.errstate(over="ignore", under="ignore", divide="ignore"):
            factors = variances / (cumulative_loading * cumulative_loading)
        if not np.isfinite(factors).all():
            return ORIGINAL(self, nodes)
    _fill_blocks(covariance, lca, requested, variances, cumulative_loading,
                 factors, stable_covariance)
"""
    source = source[:begin] + replacement + source[end:]
    scope = dict(gt.__dict__)
    scope.update(ORIGINAL=ORIGINAL, _fill_blocks=_fill_blocks)
    exec(compile(source, "<candidate02_covariance>", "exec"), scope)
    return scope["candidate_covariance"]


CANDIDATE = _candidate_function()


@contextmanager
def patched():
    gt.GaussianTreeProcess.covariance = CANDIDATE
    try:
        yield
    finally:
        gt.GaussianTreeProcess.covariance = ORIGINAL


def fixture(size, topology="balanced", model="brownian", parameter=None):
    leaves = [f"T{i}:{0.2 + (i % 11) / 20:.12g}" for i in range(size)]
    if topology == "star":
        newick = "(" + ",".join(leaves) + ")R;"
    elif topology == "pectinate":
        newick = f"({leaves[0]},{leaves[1]})I0:0.2"
        for index, leaf in enumerate(leaves[2:], 1):
            newick = f"({newick},{leaf})I{index}:0.2"
        newick = newick.rsplit(":", 1)[0] + ";"
    else:
        current = leaves
        counter = 0
        while len(current) > 1:
            next_nodes = []
            for index in range(0, len(current), 2):
                if index + 1 == len(current):
                    next_nodes.append(current[index])
                else:
                    next_nodes.append(
                        f"({current[index]},{current[index + 1]})I{counter}:0.3"
                    )
                    counter += 1
            current = next_nodes
        newick = current[0].rsplit(":", 1)[0] + ";"
    tree = read_tree(newick, "1", True, quiet=True, rooted="yes")
    process = build_evolutionary_process(tree, model=model, parameter=parameter)
    return tree, process, tuple(tree.leaves())


def compare_arrays(before, after):
    np.testing.assert_array_equal(before, after)
    return {
        "exact": True,
        "max_absolute_difference": float(np.max(np.abs(before - after), initial=0)),
        "shape": list(before.shape),
    }


def compare_fit(before, after):
    differences = {}
    for key in before:
        left, right = before[key], after[key]
        if isinstance(left, dict):
            assert left == right, key
        elif isinstance(left, (np.ndarray, np.number, float, int)):
            np.testing.assert_array_equal(left, right, err_msg=key)
            differences[key] = "exact"
        else:
            assert left == right, key
    assert set(before) == set(after)
    return {"all_fields_exact": True, "fields": sorted(before)}


def edge_checks():
    candidate = _candidate_function(0)
    results = []
    tree, process, leaves = fixture(48)
    nodes = tuple(tree.traverse())
    rng = np.random.default_rng(871)
    # Requested internal nodes, repeats, order, root and absent observations.
    for label, selected in (
        (
            "arbitrary_repeated",
            tuple(nodes[i] for i in rng.integers(0, len(nodes), 97)),
        ),
        ("reversed_all", nodes[::-1]),
        ("only_root_repeated", (tree,) * 40),
        ("empty", ()),
    ):
        results.append(
            {
                "case": label,
                **compare_arrays(
                    ORIGINAL(process, selected), candidate(process, selected)
                ),
            }
        )
    for label, slopes, variance, root_variance in (
        ("signed", lambda i: -0.8 if i % 2 else 0.7, 0.2, 0.6),
        ("zero", lambda i: 0 if i % 5 == 0 else 0.8, 0.2, 0.6),
        ("underflow", lambda i: 1e-180, 0.2, 0.6),
        ("extreme_variance", lambda i: 0.8, 1e280, 1e280),
        ("small_variance", lambda i: 0.8, 1e-280, 1e-280),
    ):
        modified = gt.GaussianTreeProcess(
            tree,
            {
                node: gt.GaussianTransition(slopes(i), 0.0, variance)
                for i, node in enumerate(nodes[1:])
            },
            gt.GaussianRootPrior("gaussian", 0.0, root_variance),
            "custom",
        )
        with np.errstate(over="ignore", under="ignore", invalid="ignore"):
            result = compare_arrays(
                ORIGINAL(modified, leaves), candidate(modified, leaves)
            )
        results.append({"case": label, **result})
    for trial in range(8):
        transitions = {
            node: gt.GaussianTransition(
                float(rng.choice([-1, 1]) * rng.uniform(0.65, 1.25)),
                float(rng.normal()),
                float(rng.uniform(0.05, 3)),
            )
            for node in nodes[1:]
        }
        modified = gt.GaussianTreeProcess(
            tree, transitions, gt.GaussianRootPrior("gaussian", 0, 0.7), "affine"
        )
        selected = tuple(nodes[index] for index in rng.integers(0, len(nodes), 150))
        results.append(
            {
                "case": f"random_signed_affine_duplicates_{trial}",
                **compare_arrays(
                    ORIGINAL(modified, selected), candidate(modified, selected)
                ),
            }
        )
    _, ou, _ = fixture(48, model="ou", parameter=0.43)
    selected = list(ou.tree.traverse()) * 2
    rng.shuffle(selected)
    results.append(
        {
            "case": "ou_shuffled_duplicates",
            **compare_arrays(ORIGINAL(ou, selected), candidate(ou, selected)),
        }
    )
    outsider = fixture(2)[0]
    for label, modified, selected in (
        ("outside_tree", process, (outsider,) * 33),
        (
            "flat_root",
            gt.GaussianTreeProcess(
                tree,
                process.transitions,
                gt.GaussianRootPrior("flat", variance=None),
                "custom",
            ),
            leaves,
        ),
        (
            "overflow",
            gt.GaussianTreeProcess(
                tree,
                {node: gt.GaussianTransition(1e150, 0, 1e300) for node in nodes[1:]},
                gt.GaussianRootPrior("gaussian", 0, 1e300),
                "custom",
            ),
            leaves,
        ),
    ):
        failures = []
        for function in (ORIGINAL, candidate):
            try:
                with np.errstate(over="ignore", under="ignore", invalid="ignore"):
                    function(modified, selected)
            except Exception as error:
                failures.append((type(error).__name__, str(error)))
            else:
                failures.append(None)
        assert failures[0] is not None and failures[0] == failures[1], (label, failures)
        results.append({"case": label, "exception": list(failures[0]), "exact": True})
    return results


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output", type=Path, default=Path(__file__).with_name("results.json")
    )
    args = parser.parse_args()
    with measurement_lock():
        results = {
            "environment": {**environment(), "nwkit_source": nwkit.__file__},
            "baseline_source_sha256": hashlib.sha256(
                inspect.getsource(ORIGINAL).encode()
            ).hexdigest(),
            "scope": "whole covariance calls including index/setup; and whole fixed Brownian PGLS fit",
            "workspace_bytes": WORKSPACE_BYTES,
            "edge_checks": edge_checks(),
            "workloads": [],
        }
        for size, topology, model, parameter in (
            (16, "balanced", "brownian", None),
            (64, "balanced", "brownian", None),
            (128, "balanced", "brownian", None),
            (512, "balanced", "brownian", None),
            (512, "balanced", "ou", 0.4),
            (256, "pectinate", "brownian", None),
            (256, "pectinate", "ou", 0.01),
            (512, "star", "brownian", None),
        ):
            _, process, leaves = fixture(size, topology, model, parameter)

            def baseline(process=process, leaves=leaves):
                return ORIGINAL(process, leaves)

            def candidate(process=process, leaves=leaves):
                return CANDIDATE(process, leaves)

            results["workloads"].append(
                {
                    "size": size,
                    "topology": topology,
                    "model": model,
                    "equivalence": compare_arrays(baseline(), candidate()),
                    "timing": paired(baseline, candidate),
                    "baseline_peak_bytes": peak_bytes(baseline),
                    "candidate_peak_bytes": peak_bytes(candidate),
                }
            )
        for size in (64, 256):
            tree, _, leaves = fixture(size)
            names = [str(leaf.name) for leaf in leaves]
            rng = np.random.default_rng(123 + size)
            x = np.column_stack([np.ones(size), rng.normal(size=size)])
            y = x @ np.asarray([0.5, -0.8]) + rng.normal(size=size)

            def fit(y=y, x=x, size=size, tree=tree, names=names):
                return _fit_ordinary_gaussian(
                    y,
                    x,
                    np.zeros(size),
                    tree,
                    names,
                    evolution_model="brownian",
                    evolution_parameter=None,
                    branch_length="original",
                    custom_covariance=None,
                    reml=True,
                )

            def candidate_fit(fit=fit):
                with patched():
                    return fit()

            before, after = fit(), candidate_fit()
            results["workloads"].append(
                {
                    "case": "whole_ordinary_brownian_fit",
                    "size": size,
                    "equivalence": compare_fit(before, after),
                    "timing": paired(fit, candidate_fit),
                    "baseline_peak_bytes": peak_bytes(fit),
                    "candidate_peak_bytes": peak_bytes(candidate_fit),
                }
            )
        args.output.write_text(json.dumps(results, indent=2, allow_nan=False) + "\n")
        print(
            json.dumps(
                [
                    {
                        "case": row["case"]
                        if "case" in row
                        else f"{row['model']}_{row['topology']}",
                        "size": row["size"],
                        "speedup": row["timing"]["speedup"],
                    }
                    for row in results["workloads"]
                ],
                indent=2,
            )
        )


if __name__ == "__main__":
    main()

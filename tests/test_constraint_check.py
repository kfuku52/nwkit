import json
import random

import pytest
from ete4 import Tree

from nwkit.cli import main
from nwkit.constraint_check import check_constraint


def tree(text):
    return Tree(text, parser=1)


def test_partial_constraint_and_free_tip_retention():
    constraint = tree("((A,B),C,D);")
    for inferred in ("(A,(B,(C,D)),X);", "(C,(D,(A,B)),X);", "(((X,B),A),(C,D));"):
        result = check_constraint(tree(inferred), constraint, list("ABCDX"))
        assert (
            result["ok"]
            and result["required_splits"] == result["displayed_splits"] == 1
        )
    result = check_constraint(tree("((A,C),B,D);"), constraint, list("ABCDX"))
    assert (
        not result["ok"]
        and result["missing_tips"] == ["X"]
        and result["free_tips"] == ["X"]
    )
    lo, hi = result["violated_intervals"][0]
    assert set(result["constraint_tip_order"][lo : hi + 1]) == {"C", "D"}
    result = check_constraint(tree("(A,B,X);"), constraint, list("ABCD"))
    assert result["extra_tips"] == ["X"] and result["missing_constraint_tips"] == [
        "C",
        "D",
    ]
    assert result["displayed_splits"] is None
    assert (
        check_constraint(tree("(A,B,C);"), tree("A;"), list("ABC"))["required_splits"]
        == 0
    )
    with pytest.raises(ValueError):
        check_constraint(tree("(A,A,C);"), constraint, list("ABCD"))


def test_against_independent_split_enumeration():
    rng = random.Random(17)

    def random_tree(names):
        parts = list(names)
        while len(parts) > 1:
            rng.shuffle(parts)
            n = rng.randint(2, min(4, len(parts)))
            parts[:n] = ["(" + ",".join(parts[:n]) + ")"]
        return tree(parts[0] + ";")

    def splits(t, tips):
        groups = set()
        for node in t.traverse():
            group = tips.intersection(node.leaf_names())
            if 1 < len(group) < len(tips) - 1:
                groups.add(frozenset(group if "A" not in group else tips - group))
        return groups

    for _ in range(80):
        fixed = set("ABCDEFG")
        c, t = random_tree(sorted(fixed)), random_tree("ABCDEFGHI")
        required = splits(c, fixed)
        result = check_constraint(t, c, list("ABCDEFGHI"))
        assert result["required_splits"] == len(required)
        assert result["displayed_splits"] == len(required & splits(t, fixed))


def test_cli_match_provenance_and_unmatched(tmp_path):
    backbone, labels, mapping, report, out = [
        tmp_path / s for s in ("backbone", "labels", "map", "report", "out")
    ]
    backbone.write_text("(Poales,Fabales);")
    labels.write_text("A_a\nB_b\nC_c\n")
    mapping.write_text(
        "leaf_name\ttaxonomy_query\tsource\tanchor_rank\nA_a\tPoales\tfamily_to_order\torder\nB_b\tFabales\trecord\torder\n"
    )
    assert main(
        [
            "constrain",
            "--backbone",
            "user",
            "--infile",
            str(backbone),
            "--species-list",
            str(labels),
            "--species-map-tsv",
            str(mapping),
            "--report",
            str(report),
            "--outfile",
            str(out),
        ]
    ) in (0, None)
    import csv

    with report.open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert [(r["status"], r["anchor"], r["source"]) for r in rows] == [
        ("constrained", "Poales", "family_to_order"),
        ("constrained", "Fabales", "record"),
        ("free", "", ""),
    ]
    result = tmp_path / "result"
    with pytest.raises(ValueError, match="failed"):
        main(
            [
                "check-constraint",
                "--infile",
                str(out),
                "--constraint",
                str(out),
                "--species-list",
                str(labels),
                "--outfile",
                str(result),
            ]
        )
    assert json.loads(result.read_text())["missing_tips"] == ["C_c"]


def test_ncbi_resolution_and_invalid_explicit_taxid(tmp_path, monkeypatch):
    import csv

    class Taxonomy:
        def get_name_translator(self, names):
            known = {"Rosa": [10], "Homo sapiens": [20], "Unclear name": [10, 20]}
            return {q: known[q] for q in names if q in known}

        def get_lineage(self, tid):
            if tid not in (10, 20):
                raise ValueError(f"Could not find taxid: {tid}")
            return [1, tid]

        def get_taxid_translator(self, ids):
            return {i: {1: "root", 10: "Rosa", 20: "Homo sapiens"}[i] for i in ids}

        def get_rank(self, ids):
            return {i: {1: "no rank", 10: "genus", 20: "species"}[i] for i in ids}

    monkeypatch.setattr(
        "nwkit.constrain.get_ete_ncbitaxa", lambda args=None: Taxonomy()
    )
    labels, report, out = [tmp_path / s for s in ("labels", "report", "out")]
    labels.write_text("Rosa_unknown\nHomo_sapiens\nUnclear_name\nAbsent_label\n")
    opts = [
        "constrain",
        "--backbone",
        "ncbi",
        "--report",
        str(report),
        "--outfile",
        str(out),
    ]
    assert main([*opts, "--species-list", str(labels)]) in (None, 0)
    with report.open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert [(r["taxid"], r["match_method"], r["anchor_rank"]) for r in rows] == [
        ("10", "genus_fallback", "genus"),
        ("20", "exact_name", "species"),
        ("", "ambiguous_name", ""),
        ("", "unmatched", ""),
    ]
    assert set(tree(out.read_text()).leaf_names()) == {"Rosa_unknown", "Homo_sapiens"}
    before = out.read_bytes()
    labels.write_text("leaf_name\ttaxid\nRosa\t10\nUnknown\t333\n")
    with pytest.raises(ValueError, match="taxid"):
        main([*opts, "--taxid-tsv", str(labels)])
    assert out.read_bytes() == before
    with report.open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert [(r["query_source"], r["reason"]) for r in rows] == [
        ("taxid_tsv", ""),
        ("taxid_tsv", "unknown_taxid"),
    ]


def test_check_diagnostics_cannot_overwrite_constraint(tmp_path):
    c = tmp_path / "c"
    c.write_text("(A,B,C);")
    inferred = tmp_path / "inferred"
    inferred.write_text("(A,B,C);")
    labels = tmp_path / "l"
    labels.write_text("A\nB\nC\n")
    with pytest.raises(ValueError, match=r"distinct from input paths \(constraint\)"):
        main(
            [
                "check-constraint",
                "--constraint",
                str(c),
                "--species-list",
                str(labels),
                "--infile",
                str(inferred),
                "--audit",
                str(c),
            ]
        )
    assert c.read_text() == "(A,B,C);"

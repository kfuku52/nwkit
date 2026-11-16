import csv
import hashlib
import io
import json
import os
import sys
from argparse import Namespace
from concurrent.futures import ThreadPoolExecutor
from dataclasses import replace
from itertools import combinations
from threading import Barrier, Event

import pytest

from nwkit import output_transaction, vertlife
from nwkit.cli import main
from nwkit.rooting_state import is_rooted
from nwkit.util import read_tree

DATASET_ID = "amphibians-2018-consensus"
REFERENCE = (
    "(Homo_sapiens:12,(((Xenopus_laevis:2.000000123456789,Xenopus_tropicalis:2):3,"
    "Ambystoma_mexicanum:5):4,(Bufo_bufo:6,Rana_temporaria:6):3):3);"
)
SELECTED = ("Xenopus_laevis", "Xenopus_tropicalis", "Bufo_bufo", "Rana_temporaria")


@pytest.fixture
def reference(tmp_path, monkeypatch):
    return _write_reference(tmp_path / "reference.nwk", monkeypatch, REFERENCE)


def _write_reference(path, monkeypatch, newick):
    data = newick.encode()
    path.write_bytes(data)
    dataset = replace(
        vertlife.DATASETS[DATASET_ID], sha256=hashlib.sha256(data).hexdigest()
    )
    monkeypatch.setitem(vertlife.DATASETS, DATASET_ID, dataset)
    return path


def _source_options(reference):
    return ["--vertlife-dataset", DATASET_ID, "--vertlife-file", str(reference)]


def _constrain_options(species):
    return ["constrain", "--backbone", "vertlife", "--species-list", str(species)]


def _read_report(path):
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_constrain_dated_and_topology(reference, tmp_path):
    species = tmp_path / "species.txt"
    species.write_text(
        "\n".join(("Xenopus laevis", *SELECTED[1:])), encoding="utf-8-sig"
    )
    dated = tmp_path / "dated.nwk"
    topology = tmp_path / "topology.nwk"
    report = tmp_path / "matches.tsv"
    options = _constrain_options(species) + _source_options(reference)
    main(
        options
        + ["--branch-length", "original", "--report", str(report), "-o", str(dated)]
    )
    main(options + ["-o", str(topology)])
    source = read_tree(str(reference), 0, True, quiet=True)
    tree = read_tree(str(dated), 0, True, quiet=True)
    assert is_rooted(tree)
    assert tree.dist is None
    assert set(tree.leaf_names()) == set(SELECTED)
    assert ":" not in topology.read_text(encoding="utf-8")
    for left, right in combinations(SELECTED, 2):
        assert tree.get_distance(left, right) == pytest.approx(
            source.get_distance(left, right), rel=0, abs=2e-15
        )
    rows = _read_report(report)
    assert rows[0]["vertlife_leaf_name"] == "Xenopus_laevis"
    assert rows[0]["source_sha256"] == vertlife.DATASETS[DATASET_ID].sha256


def test_constrain_dated_rejects_missing_and_duplicate_species(reference, tmp_path):
    species = tmp_path / "species.txt"
    species.write_text(
        "Xenopus_laevis_1\nXenopus_laevis_2\nBufo_bufo\nMissing_name\n",
        encoding="utf-8",
    )
    tree = tmp_path / "dated.nwk"
    tree.write_text("previous output", encoding="utf-8")
    report = tmp_path / "matches.tsv"
    args = _constrain_options(species) + _source_options(reference)
    args += ["--branch-length", "original", "--report", str(report), "-o", str(tree)]
    with pytest.raises(ValueError, match="Unmatched VertLife labels: Missing_name"):
        main(args)
    assert tree.read_text(encoding="utf-8") == "previous output"
    assert [row["status"] for row in _read_report(report)] == [
        "duplicate_species",
        "duplicate_species",
        "exact",
        "unmatched",
    ]
    species.write_text(
        "Xenopus_laevis_1\nXenopus_laevis_2\nBufo_bufo\n", encoding="utf-8"
    )
    with pytest.raises(ValueError, match="one label per species"):
        main(args)
    main(args + ["--collapse", "yes"])
    assert set(read_tree(str(tree), 0, True, quiet=True).leaf_names()) == {
        "Xenopus_laevis",
        "Bufo_bufo",
    }


def test_constrain_species_stdin_and_explicit_mapping(reference, tmp_path, monkeypatch):
    monkeypatch.setattr(sys, "stdin", io.StringIO("\ufeffalias\nBufo_bufo\n"))
    mapping = tmp_path / "mapping.tsv"
    mapping.write_text(
        "leaf_name\tspecies_label\nalias\tXenopus_laevis\n", encoding="utf-8"
    )
    tree = tmp_path / "tree.nwk"
    main(
        _constrain_options("-")
        + _source_options(reference)
        + ["--species-map-tsv", str(mapping), "-o", str(tree)]
    )
    assert set(read_tree(str(tree), 9, True, quiet=True).leaf_names()) == {
        "alias",
        "Bufo_bufo",
    }


def test_constrain_cache_reused_and_audited(reference, tmp_path, monkeypatch):
    downloads = []

    class Response:
        def __enter__(self):
            return self

        def __exit__(self, *args):
            pass

        def raise_for_status(self):
            pass

        def iter_content(self, *, chunk_size):
            yield reference.read_bytes()

    def download(url, *, stream, timeout):
        downloads.append(url)
        return Response()

    monkeypatch.setattr(vertlife.requests, "get", download)
    species = tmp_path / "species.txt"
    species.write_text("\n".join(SELECTED), encoding="utf-8")
    audit = tmp_path / "audit.json"
    output = tmp_path / "output.nwk"
    args = _constrain_options(species) + [
        "--vertlife-dataset",
        DATASET_ID,
        "--download-dir",
        str(tmp_path / "cache"),
        "--audit",
        str(audit),
        "-o",
        str(output),
    ]
    main(args)
    first = output.read_bytes()
    log = json.loads(audit.read_text(encoding="utf-8"))
    asset = next(item for item in log["inputs"] if item["argument"] == "vertlife_cache")
    assert asset["sha256"] == vertlife.DATASETS[DATASET_ID].sha256
    main(args)
    assert output.read_bytes() == first
    assert downloads == [vertlife.DATASETS[DATASET_ID].url]


def test_constrain_rejects_corrupt_source_and_protects_it(reference, tmp_path):
    species = tmp_path / "species.txt"
    species.write_text("\n".join(SELECTED), encoding="utf-8")
    args = _constrain_options(species) + _source_options(reference)
    with pytest.raises(ValueError, match="must not overwrite input"):
        main(args + ["--report", str(reference)])
    reference.write_text("(broken);", encoding="utf-8")
    output = tmp_path / "tree.nwk"
    with pytest.raises(ValueError, match="checksum verification failed"):
        main(args + ["-o", str(output)])
    assert not output.exists()


def test_concurrent_cold_cache_shares_verified_download(
    reference, tmp_path, monkeypatch
):
    count = []
    start = Barrier(4)

    def download(dataset):
        count.append(dataset.id)
        # Keep the first download open while the other readers reach the cache.
        Event().wait(0.1)
        return reference.read_bytes()

    def load():
        args = Namespace(
            vertlife_dataset=DATASET_ID, download_dir=str(tmp_path / "cache")
        )
        start.wait(timeout=5)
        return vertlife.load_vertlife(args)[1].write(parser=0)

    monkeypatch.setattr(vertlife, "_download", download)
    with ThreadPoolExecutor(max_workers=4) as pool:
        results = [pool.submit(load) for _ in range(4)]
        assert len({future.result(timeout=10) for future in results}) == 1
    assert count == [DATASET_ID]
    assert not list((tmp_path / "cache").rglob("*.lock"))


def test_cache_rejects_audit_path_on_download_lock(tmp_path):
    args = Namespace(vertlife_dataset=DATASET_ID, download_dir=str(tmp_path / "cache"))
    lock = vertlife.vertlife_cache_lock_paths(args)[0][1]
    species = tmp_path / "species.txt"
    species.write_text("\n".join(SELECTED), encoding="utf-8")
    with pytest.raises(
        ValueError, match="must not overwrite input VertLife download lock"
    ):
        main(
            _constrain_options(species)
            + [
                "--vertlife-dataset",
                DATASET_ID,
                "--download-dir",
                args.download_dir,
                "--audit",
                str(lock),
            ]
        )
    assert not lock.exists()


def test_failed_cache_download_can_be_retried(reference, tmp_path, monkeypatch):
    args = Namespace(vertlife_dataset=DATASET_ID, download_dir=str(tmp_path / "cache"))
    cache = vertlife.vertlife_input_path(args)
    results = iter([b"wrong content", reference.read_bytes()])
    monkeypatch.setattr(vertlife, "_download", lambda dataset: next(results))
    with pytest.raises(ValueError, match="checksum verification failed"):
        vertlife.load_vertlife(args)
    assert not cache.exists()
    assert not list((tmp_path / "cache").rglob("*.lock"))
    vertlife.load_vertlife(args)
    assert cache.read_bytes() == reference.read_bytes()


@pytest.mark.skipif(not hasattr(os, "mkfifo"), reason="POSIX FIFO required")
def test_vertlife_rejects_named_pipe(tmp_path):
    path = tmp_path / "pipe.tre"
    os.mkfifo(path)
    with pytest.raises(ValueError, match="VertLife source must be a regular file"):
        vertlife.load_vertlife(
            Namespace(vertlife_dataset=DATASET_ID, vertlife_file=str(path))
        )


def test_root_uses_reference_split(reference, tmp_path):
    source = (
        "((Xenopus_laevis:1.12345678901234,Xenopus_tropicalis:2)Xenopus_laevis:3"
        "[&&NHX:support=0.9876543210987],"
        "Bufo_bufo:4.2345678901234,Rana_temporaria:5);"
    )
    output = tmp_path / "rooted.nwk"
    main(
        [
            "root",
            "--method",
            "taxonomy",
            "--taxonomy-source",
            "vertlife",
            "-i",
            source,
            "-o",
            str(output),
        ]
        + _source_options(reference)
    )
    tree = read_tree(str(output), 1, True, quiet=True)
    assert is_rooted(tree)
    assert {frozenset(child.leaf_names()) for child in tree.children} == {
        frozenset(SELECTED[:2]),
        frozenset(SELECTED[2:]),
    }
    assert any(
        node.name == "Xenopus_laevis" and node.support == 0.9876543210987
        for node in tree.traverse()
        if not node.is_leaf
    )
    original = read_tree(source, 1, True, quiet=True)
    before = {leaf.name: leaf for leaf in original.leaves()}
    after = {leaf.name: leaf for leaf in tree.leaves()}
    for left, right in combinations(SELECTED, 2):
        assert tree.get_distance(after[left], after[right]) == pytest.approx(
            original.get_distance(before[left], before[right]), rel=0, abs=4e-15
        )


@pytest.mark.parametrize(
    "newick",
    [
        "(Xenopus_laevis_r1,(Xenopus_laevis_r2,(Bufo_bufo,Rana_temporaria)));",
        "(Xenopus_laevis_r1,Xenopus_laevis_r2,(Bufo_bufo,Rana_temporaria));",
        "((Bufo_bufo,Rana_temporaria),Xenopus_laevis_r1,Xenopus_laevis_r2);",
    ],
)
def test_root_accepts_monophyly_across_arbitrary_input_root(
    reference, tmp_path, newick
):
    output = tmp_path / "rooted.nwk"
    main(
        [
            "root",
            "--method",
            "taxonomy",
            "--taxonomy-source",
            "vertlife",
            "-i",
            "[&U]" + newick,
            "-o",
            str(output),
        ]
        + _source_options(reference)
    )
    rooted = read_tree(str(output), 0, True, quiet=True)
    assert {frozenset(child.leaf_names()) for child in rooted.children} == {
        frozenset(("Xenopus_laevis_r1", "Xenopus_laevis_r2")),
        frozenset(SELECTED[2:]),
    }


def test_root_rejects_species_split_across_distinct_clades(reference):
    source = "((Xenopus_laevis_r1,Bufo_bufo),(Xenopus_laevis_r2,Rana_temporaria));"
    with pytest.raises(ValueError, match="monophyletic on the unrooted input topology"):
        main(
            [
                "root",
                "--method",
                "taxonomy",
                "--taxonomy-source",
                "vertlife",
                "-i",
                source,
            ]
            + _source_options(reference)
        )


def test_root_tries_configured_fallback_when_local_reference_is_missing(
    tmp_path, monkeypatch
):
    labels = set(SELECTED)
    monkeypatch.setattr(
        "nwkit.root._build_opentree_reference_tree",
        lambda tree, args: (set(SELECTED[:2]), labels, set()),
    )
    tree = "((Xenopus_laevis,Xenopus_tropicalis),Bufo_bufo,Rana_temporaria);"
    output = tmp_path / "rooted.nwk"
    main(
        [
            "root",
            "--method",
            "taxonomy",
            "--taxonomy-source",
            "vertlife,opentree",
            "-i",
            tree,
            "-o",
            str(output),
        ]
        + _source_options(tmp_path / "missing.nwk")
    )
    rooted = read_tree(str(output), 0, True, quiet=True)
    assert {frozenset(child.leaf_names()) for child in rooted.children} == {
        frozenset(SELECTED[:2]),
        frozenset(SELECTED[2:]),
    }


def test_root_reports_ambiguous_partial_split(reference):
    source = "((Xenopus_laevis,Xenopus_tropicalis),(Bufo_bufo,Rana_temporaria),Unknown_species);"
    with pytest.raises(
        ValueError, match="interfere with the vertlife-derived root position"
    ):
        main(
            [
                "root",
                "--method",
                "taxonomy",
                "--taxonomy-source",
                "vertlife",
                "-i",
                source,
            ]
            + _source_options(reference)
        )


def test_mcmctree_reports_ages_without_altering_input_tree(reference, tmp_path):
    source = (
        "[&R](((Xenopus_laevis:2.12345678901,Xenopus_tropicalis:2)Xenopus:1"
        "[&&NHX:support=0.9876543210987],"
        "Bufo_bufo:3)Other:3,Rana_temporaria:6)Root;"
    )
    tree_path = tmp_path / "unchanged.nwk"
    report = tmp_path / "ages.tsv"
    main(
        ["mcmctree", "-i", source, "--report", str(report), "-o", str(tree_path)]
        + _source_options(reference)
    )
    before = read_tree(source, 1, True, quiet=True)
    after = read_tree(str(tree_path), 1, True, quiet=True)
    assert [node.props for node in before.traverse() if not node.is_root] == [
        node.props for node in after.traverse() if not node.is_root
    ]
    assert after.name == "Root"
    rows = {row["node_name"]: row for row in _read_report(report)}
    assert rows["Root"]["status"] == "child_clades_differ"
    assert rows["Root"]["age_ma"] == ""
    assert rows["Other"]["status"] == "clade_absent"
    assert rows["Xenopus"]["status"] == "matched"
    assert float(rows["Xenopus"]["age_ma"]) == pytest.approx(2.0000000617283944)
    assert float(rows["Xenopus"]["max_tip_path_myr"]) > 2


def test_mcmctree_keeps_names_carried_alongside_numeric_support(reference, tmp_path):
    source = (
        "((Xenopus_laevis,Xenopus_tropicalis)0.9876543210987[&&NHX:name=Xenopus],"
        "(Bufo_bufo,Rana_temporaria));"
    )
    tree_path, report = (tmp_path / name for name in ("tree.nwk", "ages.tsv"))
    main(
        [
            "mcmctree",
            "-i",
            source,
            "-f",
            "0",
            "-o",
            str(tree_path),
            "--report",
            str(report),
        ]
        + _source_options(reference)
    )
    clade = read_tree(str(tree_path), 1, True, quiet=True).children[0]
    assert clade.name == "Xenopus"
    assert clade.support == 0.9876543210987


def test_mcmctree_serializes_numeric_nhx_fields_as_numbers(reference, tmp_path):
    source = (
        "((Xenopus_laevis[&&NHX:dist=0.1234567890123],Xenopus_tropicalis)"
        "[&&NHX:support=0.9876543210987],(Bufo_bufo,Rana_temporaria));"
    )
    tree_path, report = (tmp_path / name for name in ("tree.nwk", "ages.tsv"))
    main(
        [
            "mcmctree",
            "-i",
            source,
            "-f",
            "0",
            "-o",
            str(tree_path),
            "--report",
            str(report),
        ]
        + _source_options(reference)
    )
    clade = read_tree(str(tree_path), 0, True, quiet=True).children[0]
    assert clade.support == 0.9876543210987
    assert clade.children[0].dist == 0.1234567890123


def test_mcmctree_does_not_match_noncontiguous_reference_clades(reference, tmp_path):
    report = tmp_path / "ages.tsv"
    source = "((Xenopus_laevis,Bufo_bufo),(Xenopus_tropicalis,Rana_temporaria));"
    main(
        ["mcmctree", "-i", source, "--report", str(report)] + _source_options(reference)
    )
    rows = _read_report(report)
    assert [row["status"] for row in rows] == [
        "child_clades_differ",
        "clade_absent",
        "clade_absent",
    ]
    assert all(row["age_ma"] == "" for row in rows)


def test_mcmctree_reports_rooted_polytomy_from_constrain(tmp_path, monkeypatch):
    reference = _write_reference(
        tmp_path / "reference.nwk",
        monkeypatch,
        "(Homo_sapiens:12,(Bufo_bufo:5,Rana_temporaria:5,Ambystoma_mexicanum:5):7);",
    )
    species = tmp_path / "species.txt"
    species.write_text(
        "Bufo_bufo\nRana_temporaria\nAmbystoma_mexicanum\n", encoding="utf-8"
    )
    tree, copy, report = (
        tmp_path / name for name in ("dated.nwk", "copy.nwk", "ages.tsv")
    )
    source = _source_options(reference)
    main(
        _constrain_options(species)
        + source
        + ["--branch-length", "original", "-o", str(tree)]
    )
    main(
        ["mcmctree", "-i", str(tree), "--report", str(report), "-o", str(copy)] + source
    )
    [row] = _read_report(report)
    assert (row["branch_id"], row["status"], row["age_ma"]) == ("0", "matched", "5.0")
    assert tree.read_bytes() == copy.read_bytes()


def test_mcmctree_rolls_back_tree_and_report_on_publish_failure(
    reference, tmp_path, monkeypatch
):
    tree, report = (tmp_path / name for name in ("tree.nwk", "ages.tsv"))
    tree.write_text("old tree", encoding="utf-8")
    report.write_text("old report", encoding="utf-8")
    real_replace = output_transaction.replace_output

    def fail_on_tree(source, target):
        real_replace(source, target)
        if ".stage." in str(source) and os.path.realpath(target) == str(tree.resolve()):
            raise OSError("injected tree installation failure")

    monkeypatch.setattr(output_transaction, "replace_output", fail_on_tree)
    with pytest.raises(OSError, match="injected tree"):
        main(
            [
                "mcmctree",
                "-i",
                "((Xenopus_laevis,Xenopus_tropicalis),(Bufo_bufo,Rana_temporaria));",
                "--report",
                str(report),
                "-o",
                str(tree),
            ]
            + _source_options(reference)
        )
    assert tree.read_text(encoding="utf-8") == "old tree"
    assert report.read_text(encoding="utf-8") == "old report"


@pytest.mark.parametrize(
    "conflicting",
    [["--angiocal", "v1.0"], ["--timetree", "point"], ["--posterior", "mcmc.txt"]],
)
def test_mcmctree_rejects_combining_calibration_modes(reference, tmp_path, conflicting):
    with pytest.raises(ValueError, match="cannot be combined"):
        main(
            ["mcmctree", "--report", str(tmp_path / "ages.tsv")]
            + _source_options(reference)
            + conflicting
        )

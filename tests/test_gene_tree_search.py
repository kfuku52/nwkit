"""Small exhaustive references, coupled anomalies, and the supplied BMI1 input."""

import hashlib
import io
import json
import os
import subprocess
import sys
import time
from pathlib import Path

import pandas as pd
import pytest
from ete4 import Tree

from nwkit.cli import main
from nwkit.gene_tree_search import _best_evaluation
from nwkit.gene_tree_search_generax import (
    Evaluation,
    _parse_result,
    _run_generax,
    evaluate_candidates,
    read_alignment,
    tree_text,
)
from nwkit.gene_tree_search_model import (
    DLContext,
    Proposal,
    _remove_and_suppress,
    discover_proposals,
    generate_candidates,
    topology_key,
)
from nwkit.reconcile import build_reconciliation_table
from nwkit.species_parser import get_species_parser

DATA = Path(__file__).parent / "data" / "gene_tree_search" / "bmi1"


def context(species, gene):
    parser = get_species_parser()
    return DLContext(
        species, {n.name: parser.parse(n.name).species_label for n in gene.leaves()}
    )


def test_deep_tree_serialization_and_candidate_copy_preserve_all_tips():
    species = Tree("(A_a,B_b);", parser=1)
    gene = Tree()
    cursor = gene
    for index in range(1200):
        cursor.add_child(name=f"A_a_{index}")
        cursor = cursor.add_child()
    cursor.name = "B_b_0"
    ctx = context(species, gene)
    candidates, _ = generate_candidates(gene, ctx, [], max_candidates=1)
    assert len(candidates) == 1
    assert candidates[0].tree is not gene
    assert topology_key(candidates[0].tree) == topology_key(gene)
    text = tree_text(candidates[0].tree)
    reconstructed = Tree(text, parser=1)
    assert set(reconstructed.leaf_names()) == set(gene.leaf_names())
    assert all(
        len(node.children) == 2 for node in reconstructed.traverse() if not node.is_leaf
    )


@pytest.fixture
def coupled():
    species = Tree("(A_a,(B_b,C_c));", parser=1)
    gene = Tree("(((A_a_1,B_b_1),C_c_1),(A_a_2,B_b_2));", parser=1)
    return species, gene, context(species, gene)


def test_overlap_cover_detects_two_tips_before_single_move_search(coupled):
    _, gene, ctx = coupled
    proposals, coverage = discover_proposals(gene, ctx, max_proposals=100)
    covers = {p.tips for p in proposals if "species_overlap_cover" in p.reasons}
    assert {
        ("A_a_1", "B_b_1"),
        ("A_a_1", "B_b_2"),
        ("A_a_2", "B_b_1"),
        ("A_a_2", "B_b_2"),
    } <= covers
    # At this duplication, neither singleton can remove the species overlap.
    for tip in ("A_a_1", "B_b_2"):
        backbone, _ = _remove_and_suppress(gene, (tip,))
        left, right = backbone.children
        assert {ctx.mapping[n.name] for n in left.leaves()} & {
            ctx.mapping[n.name] for n in right.leaves()
        }
    assert not coverage["set_states_truncated"]


def test_joint_endpoints_keep_all_tips_and_can_bypass_non_improving_singletons(coupled):
    _, gene, ctx = coupled
    proposed = [Proposal(("A_a_1", "B_b_2"), {"species_overlap_cover"})]
    candidates, coverage = generate_candidates(
        gene, ctx, proposed, beam_width=10000, max_candidates=10000
    )
    assert coverage["beam_states_discarded"] == 0
    assert all(set(c.tree.leaf_names()) == set(gene.leaf_names()) for c in candidates)
    assert all(
        all(n.is_leaf or len(n.children) == 2 for n in c.tree.traverse())
        for c in candidates
    )
    assert any(
        c.components == 2 and c.cost.total < candidates[0].cost.total
        for c in candidates
    )


def _all_rooted_trees(names):
    if len(names) == 1:
        return [names[0]]
    from itertools import combinations

    output = []
    # Include the first name on one side to eliminate mirrored partitions.
    for size in range(1, len(names)):
        for rest in combinations(names[1:], size - 1):
            left = (names[0], *rest)
            right = tuple(n for n in names if n not in left)
            for a in _all_rooted_trees(left):
                for b in _all_rooted_trees(right):
                    output.append(f"({a},{b})")
    return output


def test_exhaustive_small_tree_search_and_shared_reconciliation_reference(coupled):
    species, gene, ctx = coupled
    proposals, _ = discover_proposals(gene, ctx, max_proposals=100)
    candidates, coverage = generate_candidates(
        gene, ctx, proposals, beam_width=10000, max_candidates=10000
    )
    all_trees = [
        Tree(text + ";", parser=1)
        for text in _all_rooted_trees(tuple(sorted(gene.leaf_names())))
    ]
    assert len(all_trees) == 105
    assert {topology_key(c.tree) for c in candidates} == {
        topology_key(t) for t in all_trees
    }
    assert coverage["candidate_topologies_truncated"] == 0
    unrooted, report = generate_candidates(
        gene, ctx, proposals, beam_width=10000, max_candidates=10000, rooted=False
    )
    assert len(unrooted) == 15
    assert {topology_key(c.tree, rooted=False) for c in unrooted} == {
        topology_key(t, rooted=False) for t in all_trees
    }
    assert report["topology_equivalence"] == "unrooted"
    for candidate in candidates:
        table = build_reconciliation_table(
            candidate.tree.copy(), species.copy(), ctx.mapping, event_source="lca"
        )
        assert candidate.cost.duplications == (table.event_type == "duplication").sum()
        assert candidate.cost.losses == table.implied_losses.sum()


def test_real_ancient_duplication_is_not_eliminated_by_parsimony():
    species = Tree("(A_a,(B_b,C_c));", parser=1)
    gene = Tree("((A_a_1,(B_b_1,C_c_1)),(A_a_2,(B_b_2,C_c_2)));", parser=1)
    ctx = context(species, gene)
    proposals, _ = discover_proposals(gene, ctx, max_proposals=100)
    candidates, _ = generate_candidates(
        gene, ctx, proposals, beam_width=100, max_candidates=1000
    )
    assert ctx.cost(gene).duplications == 1 and ctx.cost(gene).losses == 0
    assert min(c.cost.total for c in candidates) == 1


def test_clade_move_counts_tips_separately_from_components(coupled):
    _, gene, ctx = coupled
    candidates, _ = generate_candidates(
        gene, ctx, [Proposal(("A_a_2", "B_b_2"))], beam_width=100
    )
    assert any(len(c.moved_tips) == 2 and c.components == 1 for c in candidates)


def test_limits_and_child_order_are_auditable(coupled):
    _, gene, ctx = coupled
    proposals, coverage = discover_proposals(
        gene, ctx, max_set_states=1, max_proposals=2
    )
    assert coverage["set_states_truncated"]
    assert coverage["proposal_sets_truncated"] > 0
    candidates, coverage = generate_candidates(
        gene, ctx, proposals, beam_width=1, max_candidates=2
    )
    assert len(candidates) <= 2 and coverage["beam_states_discarded"] > 0
    mirrored = gene.copy()
    for node in mirrored.traverse():
        node.children.reverse()
    original, _ = discover_proposals(gene, ctx, max_proposals=100)
    reversed_order, _ = discover_proposals(mirrored, ctx, max_proposals=100)
    assert [(p.tips, p.reasons) for p in original] == [
        (p.tips, p.reasons) for p in reversed_order
    ]


def test_bmi1_is_detected_generally_and_fixture_has_provenance():
    manifest = json.loads((DATA / "PROVENANCE.json").read_text())
    for name, record in manifest["files"].items():
        assert (
            hashlib.sha256((DATA / name).read_bytes()).hexdigest() == record["sha256"]
        )
    gene = Tree((DATA / "gene.nwk").read_text(), parser=1)
    species = Tree((DATA / "species.nwk").read_text(), parser=1)
    ctx = context(species, gene)
    proposals, _ = discover_proposals(gene, ctx, max_proposals=32)
    assert proposals[0].tips == ("Nymphaea_colorata_GeneID116261861",)
    assert proposals[0].diagnostic_gain == 6
    # Removing the dataset's names must give the same detection; no gene-ID rule.
    renamed = gene.copy()
    mapping = {}
    for number, leaf in enumerate(renamed.leaves()):
        mapping[f"tip_{number}"] = ctx.mapping[leaf.name]
        leaf.name = f"tip_{number}"
    anonymous, _ = discover_proposals(
        renamed, DLContext(species, mapping), max_proposals=32
    )
    assert len(anonymous[0].tips) == 1 and anonymous[0].diagnostic_gain == 6
    records = read_alignment(DATA / "alignment.fa.gz", gene.leaf_names())
    assert len(records) == 79 and {len(s) for s in records.values()} == {3273}


@pytest.mark.parametrize(
    "text,match",
    [
        (">A\nAC\n>B\nA\n", "common"),
        (">A\nAC\n>A\nAC\n", "unique"),
        (">A\nAC\n>C\nAC\n", "mismatch"),
        (">A\n--\n>B\nAC\n", "observed"),
    ],
)
def test_alignment_failures(text, match, tmp_path):
    path = tmp_path / "alignment.fa"
    path.write_text(text)
    with pytest.raises(ValueError, match=match):
        read_alignment(path, ("A", "B"))


def test_gene_rax_eval_rejects_changed_topology(coupled, tmp_path):
    _, gene, ctx = coupled
    candidates, _ = generate_candidates(gene, ctx, [Proposal(("A_a_1",))], beam_width=3)
    directory = tmp_path / "results" / "baseline"
    directory.mkdir(parents=True)
    (directory / "stats.txt").write_text("-100 -10\n0.1 0.1\n")
    (directory / "geneTree.newick").write_text(tree_text(candidates[1].tree))
    with pytest.raises(ValueError, match="changed topology"):
        _parse_result(tmp_path, candidates[0], root_policy="keep")


def _cli_args(tmp_path, coupled):
    species, gene, _ = coupled
    (tmp_path / "gene.nwk").write_text(tree_text(gene))
    (tmp_path / "species.nwk").write_text(tree_text(species))
    return [
        "gene-tree-search",
        "-i",
        str(tmp_path / "gene.nwk"),
        "--species-tree",
        str(tmp_path / "species.nwk"),
        "--tree-out",
        str(tmp_path / "best.nwk"),
        "--report-out",
        str(tmp_path / "report.json"),
        "-o",
        str(tmp_path / "scores.tsv"),
    ]


@pytest.mark.parametrize("token", [False, True])
def test_proposal_only_never_repairs_input_on_event_count(coupled, tmp_path, token):
    args = _cli_args(tmp_path, coupled)
    args += [
        "--rooting-token",
        "yes" if token else "no",
        "--candidates-out",
        str(tmp_path / "candidates.tsv"),
    ]
    main(args)
    best = Tree((tmp_path / "best.nwk").read_text().removeprefix("[&R]"), parser=1)
    assert topology_key(best) == topology_key(coupled[1])
    assert (tmp_path / "best.nwk").read_text().startswith("[&R]") is token
    import csv

    with (tmp_path / "candidates.tsv").open() as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    assert rows and all(row["newick"].startswith("[&R]") is token for row in rows)
    assert (
        json.loads((tmp_path / "report.json").read_text())["evaluated_topologies"] == 0
    )


def test_likelihood_gate_and_evaluation_failure_preserve_existing_outputs(
    coupled, tmp_path, monkeypatch
):
    args = _cli_args(tmp_path, coupled)
    alignment = tmp_path / "alignment.fa"
    alignment.write_text("".join(f">{n}\nACGT\n" for n in coupled[1].leaf_names()))
    args += [
        "--evaluation",
        "generax",
        "--alignment",
        str(alignment),
        "--subst-model",
        "GTR+G4",
        "--workdir",
        str(tmp_path / "work"),
    ]

    def score(candidates, *_args, **_kwargs):
        # Every event-reducing candidate is worse under the sequence+DL model.
        return {
            c.id: Evaluation(
                -100 if c.id == "baseline" else -200, -10, tree_text(c.tree)
            )
            for c in candidates
        }, {}

    monkeypatch.setattr("nwkit.gene_tree_search.evaluate_candidates", score)
    main(args)
    frame = pd.read_csv(tmp_path / "scores.tsv", sep="\t")
    assert frame.loc[frame.selected, "candidate_id"].tolist() == ["baseline"]
    before = {p: p.read_bytes() for p in tmp_path.glob("*.tsv")}

    def fail(*_args, **_kwargs):
        raise RuntimeError("backend failed")

    monkeypatch.setattr("nwkit.gene_tree_search.evaluate_candidates", fail)
    with pytest.raises(RuntimeError, match="backend failed"):
        main(args)
    assert all(p.read_bytes() == value for p, value in before.items())


def test_output_cannot_overwrite_an_input(coupled, tmp_path):
    args = _cli_args(tmp_path, coupled)
    args[args.index("-o") + 1] = str(tmp_path / "species.nwk")
    with pytest.raises(ValueError, match="replace|input"):
        main(args)


def test_round_protocol_keeps_alignment_and_fits_rates_independently(
    coupled, tmp_path, monkeypatch
):
    species, gene, ctx = coupled
    candidates, _ = generate_candidates(gene, ctx, [Proposal(("A_a_1",))], beam_width=2)
    alignment = tmp_path / "alignment.fa"
    alignment.write_text("".join(f">{n}\nACGT\n" for n in gene.leaf_names()))
    calls = []

    def run(argv, **kwargs):
        assert kwargs["timeout"] == 3600
        calls.append(argv)
        prefix = Path(kwargs["cwd"]) / argv[argv.index("--prefix") + 1]
        families = (
            Path(kwargs["cwd"]) / argv[argv.index("--families") + 1]
        ).read_text()
        assert "subst_model = GTR+G4" in families
        assert f"alignment = {Path('..') / 'alignment.fasta'}" in families
        assert f"mapping = {Path('..') / 'mapping.tsv'}" in families
        for candidate in candidates:
            directory = prefix / "results" / candidate.id
            directory.mkdir(parents=True)
            source = Path(kwargs["cwd"]) / (candidate.id + ".nwk")
            (directory / "geneTree.newick").write_text(source.read_text())
            (directory / "stats.txt").write_text(
                "-100.1 -10.12\nReconciliation rates = 0.2 0.2"
            )

    monkeypatch.setattr("nwkit.gene_tree_search_generax._run_generax", run)
    scores, metadata = evaluate_candidates(
        candidates,
        species,
        ctx.mapping,
        alignment,
        tmp_path / "work with spaces # and 'quotes'",
        subst_model="GTR+G4",
        rounds=2,
        root_policy="keep",
    )
    assert len(calls) == 2 and len(metadata["rounds"]) == 2
    assert all(
        "--per-family-rates" in call and "--enforce-gene-tree-root" in call
        for call in calls
    )
    assert all(call[call.index("--strategy") + 1] == "EVAL" for call in calls)
    assert all(call[call.index("--rec-weight") + 1] == "1.0" for call in calls)
    assert metadata["alignment_tips"] == 5 and metadata["alignment_sites"] == 4
    assert set(scores) == {c.id for c in candidates}
    assert scores["baseline"].rounding_bound == pytest.approx(0.055)


def test_rounds_retain_best_baseline_and_warm_start_best_fit(
    coupled, tmp_path, monkeypatch
):
    species, gene, ctx = coupled
    candidates, _ = generate_candidates(gene, ctx, [Proposal(("A_a_1",))], beam_width=2)
    alignment = tmp_path / "alignment.fa"
    alignment.write_text("".join(f">{n}\nACGT\n" for n in gene.leaf_names()))
    starts = []
    first_fits = {}

    def evaluate(candidates, *_args, **kwargs):
        starts.append(kwargs["starting_text"])
        number = len(starts)
        fits = {}
        for c in candidates:
            tree = c.tree.copy()
            next(tree.leaves()).dist = number / 10
            value = (-100 if number == 1 else -110) if c.id == "baseline" else -105
            fits[c.id] = Evaluation(value, 0, tree_text(tree))
        if number == 1:
            first_fits.update(fits)
        return fits, {"likelihoods": {}}

    monkeypatch.setattr("nwkit.gene_tree_search_generax._evaluate_round", evaluate)
    scores, metadata = evaluate_candidates(
        candidates,
        species,
        ctx.mapping,
        alignment,
        tmp_path / "work",
        subst_model="GTR+G4",
        rounds=3,
    )
    assert scores["baseline"].joint == -100
    assert _best_evaluation(scores) == "baseline"
    assert metadata["best_round_per_topology"]["baseline"] == 1
    assert starts[2] == {name: fit.optimized_tree for name, fit in first_fits.items()}


def test_resolved_improvement_wins_even_if_highest_point_estimate_is_ambiguous():
    scores = {
        "baseline": Evaluation(-100, 0, "", 0.1),
        "uncertain": Evaluation(-98, 0, "", 3),
        "resolved": Evaluation(-99, 0, "", 0.1),
    }
    assert _best_evaluation(scores) == "resolved"


def test_deduplication_uses_shortest_move_metadata_even_for_reserved_endpoints(coupled):
    _, gene, ctx = coupled
    proposals = [
        Proposal(("A_a_1", "B_b_1")),
        Proposal(("A_a_2",)),
    ]
    minimum = {}
    for proposal in proposals:
        generated, _ = generate_candidates(
            gene,
            ctx,
            [proposal],
            beam_width=1000,
            max_candidates=1000,
            rooted=False,
        )
        for c in generated:
            identity = topology_key(c.tree, rooted=False)
            minimum[identity] = min(minimum.get(identity, 100), len(c.moved_tips))
    combined, _ = generate_candidates(
        gene,
        ctx,
        proposals,
        beam_width=1000,
        max_candidates=1000,
        rooted=False,
    )
    # The same unrooted topology is reached through a two-tip clade move or a
    # singleton from the opposite side. A reserved endpoint must use the final
    # canonical singleton, even if its supplied root has a higher heuristic cost.
    assert all(
        len(c.moved_tips) == minimum[topology_key(c.tree, rooted=False)]
        for c in combined
    )


@pytest.mark.parametrize("source", ("inline", "stdin"))
def test_tree_text_inputs_have_hashes_and_work_with_reports(
    coupled, tmp_path, monkeypatch, source
):
    args = _cli_args(tmp_path, coupled)
    text = tree_text(coupled[1])
    args[args.index("-i") + 1] = text if source == "inline" else "-"
    args[args.index("--species-tree") + 1] = tree_text(coupled[0])
    monkeypatch.setattr(sys, "stdin", io.StringIO(text))
    main(args)
    record = json.loads((tmp_path / "report.json").read_text())["inputs"]["infile"]
    assert record == {
        "source": source,
        "sha256": hashlib.sha256(text.encode()).hexdigest(),
    }


@pytest.mark.parametrize("option", ("--sets-out", "--report-out", "--candidates-out"))
def test_audit_cannot_overlap_search_outputs(coupled, tmp_path, option):
    args = _cli_args(tmp_path, coupled)
    path = tmp_path / "protected"
    path.write_text("old content")
    if option in args:
        args[args.index(option) + 1] = str(path)
    else:
        args += [option, str(path)]
    args += ["--audit", str(path)]
    with pytest.raises(ValueError, match="audit|Audit"):
        main(args)
    assert path.read_text() == "old content"


def test_missing_output_directory_fails_before_search(coupled, tmp_path, monkeypatch):
    args = _cli_args(tmp_path, coupled)
    args[args.index("-o") + 1] = str(tmp_path / "absent" / "scores.tsv")

    def fail(*_args, **_kwargs):
        pytest.fail("search must not start for an invalid output path")

    monkeypatch.setattr("nwkit.gene_tree_search.discover_proposals", fail)
    with pytest.raises(ValueError, match="directory must already exist"):
        main(args)


@pytest.mark.parametrize("model", ("GTR+G4", "HKY{1/2}+G4", "DNA010010+G"))
def test_alignment_unknown_dna_n_is_rejected(coupled, tmp_path, model):
    alignment = tmp_path / "dna.fa"
    alignment.write_text("\ufeff>A\nNNNN\n>B\nACGT\n", encoding="utf-8")
    with pytest.raises(ValueError, match="no observed"):
        read_alignment(alignment, ("A", "B"), subst_model=model)


def test_protein_asparagine_is_observed(tmp_path):
    alignment = tmp_path / "protein.fa"
    alignment.write_text("\ufeff>A\nNNNN\n>B\nACGT\n", encoding="utf-8")
    assert read_alignment(alignment, ("A", "B"), subst_model="LG+G4")["A"] == "NNNN"


@pytest.mark.parametrize(
    "failure", ("negative_branch", "nonfinite_score", "multifurcation")
)
def test_invalid_fitted_results_are_rejected(coupled, tmp_path, failure):
    _, gene, ctx = coupled
    candidate = generate_candidates(gene, ctx, [])[0][0]
    tree = gene.copy()
    stats = "-100 -10\n"
    if failure == "negative_branch":
        next(tree.leaves()).dist = -1
    elif failure == "nonfinite_score":
        stats = "nan -10\n"
    else:
        child = tree.children[0]
        grandchildren = list(child.children)
        child.detach()
        for node in grandchildren:
            node.detach()
            tree.add_child(node)
    directory = tmp_path / "results" / "baseline"
    directory.mkdir(parents=True)
    (directory / "stats.txt").write_text(stats)
    (directory / "geneTree.newick").write_text(tree_text(tree))
    with pytest.raises(ValueError, match="Negative|Nonfinite|bifurcating|rooted"):
        _parse_result(tmp_path, candidate, root_policy="optimize")


@pytest.mark.skipif(os.name != "posix", reason="Local MPI process groups require POSIX")
@pytest.mark.parametrize("failure", ("timeout", "nonzero_exit"))
def test_failed_launcher_stops_children(tmp_path, failure):
    ready = tmp_path / "ready"
    marker = tmp_path / "orphan"
    child = (
        "import pathlib,time; "
        f"pathlib.Path({str(ready)!r}).touch(); "
        f"time.sleep(1); pathlib.Path({str(marker)!r}).touch()"
    )
    parent = (
        "import subprocess,sys,pathlib,time\n"
        f"subprocess.Popen([sys.executable,'-c',{child!r}])\n"
        f"while not pathlib.Path({str(ready)!r}).exists(): time.sleep(0.01)\n"
        + ("time.sleep(20)" if failure == "timeout" else "raise SystemExit(7)")
    )
    with (tmp_path / "log").open("w") as log:
        error = (
            subprocess.TimeoutExpired
            if failure == "timeout"
            else subprocess.CalledProcessError
        )
        with pytest.raises(error):
            _run_generax(
                [sys.executable, "-c", parent], cwd=tmp_path, log=log, timeout=0.5
            )
    assert ready.exists(), "Launcher must have started its child before timeout"
    time.sleep(1)
    assert not marker.exists()


def test_failed_launcher_is_fatal(tmp_path):
    with (tmp_path / "log").open("w") as log:
        with pytest.raises(subprocess.CalledProcessError) as failure:
            _run_generax(
                [sys.executable, "-c", "raise SystemExit(7)"],
                cwd=tmp_path,
                log=log,
                timeout=10,
            )
    assert failure.value.returncode == 7


def test_two_tip_candidate_can_win_when_every_singleton_is_worse(
    coupled, tmp_path, monkeypatch
):
    args = _cli_args(tmp_path, coupled)
    alignment = tmp_path / "alignment.fa"
    alignment.write_text("".join(f">{n}\nACGT\n" for n in coupled[1].leaf_names()))
    args += [
        "--evaluation",
        "generax",
        "--alignment",
        str(alignment),
        "--subst-model",
        "GTR+G4",
        "--workdir",
        str(tmp_path / "work"),
        "--max-evaluations",
        "128",
    ]

    def score(candidates, *_args, **_kwargs):
        assert any(len(c.moved_tips) == 2 for c in candidates)
        return {
            c.id: Evaluation(
                -100 if not c.moved_tips else -200 if len(c.moved_tips) == 1 else -50,
                -10,
                tree_text(c.tree),
            )
            for c in candidates
        }, {}

    monkeypatch.setattr("nwkit.gene_tree_search.evaluate_candidates", score)
    main(args)
    frame = pd.read_csv(tmp_path / "scores.tsv", sep="\t")
    assert frame.loc[frame.selected, "num_moved_tips"].item() >= 2


def test_rounding_ambiguity_does_not_replace_baseline(coupled, tmp_path, monkeypatch):
    args = _cli_args(tmp_path, coupled)
    alignment = tmp_path / "alignment.fa"
    alignment.write_text("".join(f">{n}\nACGT\n" for n in coupled[1].leaf_names()))
    args += [
        "--evaluation",
        "generax",
        "--alignment",
        str(alignment),
        "--subst-model",
        "GTR+G4",
        "--workdir",
        str(tmp_path / "work"),
    ]

    def score(candidates, *_args, **_kwargs):
        return {
            c.id: Evaluation(
                -100 if c.id == "baseline" else -99.9,
                -10,
                tree_text(c.tree),
                rounding_bound=0.5,
            )
            for c in candidates
        }, {}

    monkeypatch.setattr("nwkit.gene_tree_search.evaluate_candidates", score)
    main(args)
    assert (
        json.loads((tmp_path / "report.json").read_text())["selected_candidate"]
        == "baseline"
    )


@pytest.mark.slow
@pytest.mark.parametrize("root_policy", ("keep", "optimize"))
@pytest.mark.parametrize("rec_model", ("UndatedDL", "UndatedDTL"))
def test_real_generax_on_bmi1_when_available(tmp_path, root_policy, rec_model):
    import shutil

    if shutil.which("generax") is None:
        pytest.skip("GeneRax executable unavailable")
    gene = Tree((DATA / "gene.nwk").read_text(), parser=1)
    species = Tree((DATA / "species.nwk").read_text(), parser=1)
    ctx = context(species, gene)
    proposals, _ = discover_proposals(gene, ctx, max_proposals=1)
    candidates, _ = generate_candidates(
        gene, ctx, proposals, beam_width=2, max_candidates=3
    )
    scores, metadata = evaluate_candidates(
        candidates,
        species,
        ctx.mapping,
        DATA / "alignment.fa.gz",
        tmp_path / "GeneRax with spaces # and 'quotes'",
        subst_model="GTR+G4",
        rounds=2,
        rec_model=rec_model,
        root_policy=root_policy,
    )
    assert set(scores) == {c.id for c in candidates}
    assert len(metadata["rounds"]) == 2
    assert metadata["alignment_tips"] == 79 and metadata["alignment_sites"] == 3273

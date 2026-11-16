import json
from io import StringIO

import pandas as pd
import pytest

from nwkit.cli import main
from nwkit.util import read_tree
from nwkit.wgd_count import _count_tree


def model(tree):
    design, _, _ = _count_tree(read_tree(str(tree), "auto", True))
    fit = {"rates": [[0.1, 0.2]], "root_mean": 1.2}
    return {
        "species_tree": {
            "parents": list(design.parents),
            "lengths": list(design.lengths),
        },
        "species_event_ids": list(design.clade_ids),
        "tip_names": list(design.tip_names),
        "branch_rate_groups": [0] * len(design.parents),
        "detection_probabilities": [1] * len(design.tip_names),
        "family_rate_scales": [1],
        "background": fit,
        "candidates": [
            {
                "event": {
                    **fit,
                    "event": {
                        "node": 1,
                        "retention": 0.8,
                        "fraction": 0.5,
                        "multiplicity": 2,
                    },
                }
            }
        ],
    }


def inputs(tmp_path):
    species = tmp_path / "species.nwk"
    species.write_text("(A_a:1,B_b:1);")
    gene = tmp_path / "gene.nwk"
    gene.write_text("(A_a_1:1,A_a_2:1);")
    fit = tmp_path / "fit.json"
    fit.write_text(json.dumps(model(species)))
    return gene, species, fit


def command(gene, species, fit):
    return [
        "wgd-tree",
        "-i",
        str(gene),
        "--species-tree",
        str(species),
        "--count-model",
        str(fit),
        "--species-parser",
        "taxonomic",
    ]


def test_cli_fixed_tree_likelihood_and_conditional_assignments(tmp_path):
    gene, species, fit = inputs(tmp_path)
    output, report = tmp_path / "origins.tsv", tmp_path / "likelihood.json"
    main(command(gene, species, fit) + ["-o", str(output), "--model-out", str(report)])
    rows = pd.read_csv(output, sep="\t")
    assert len(rows) == 1
    assert 0 < rows.conditional_wgd_probability.iloc[0] < 1
    assert (
        rows.probability_meaning.iloc[0]
        == "latent_node_origin_given_topology_parameters_and_one_event"
    )
    assert rows.parameter_source.iloc[0] == "supplied_count_fit"
    assert (
        json.loads(report.read_text())["method"]
        == "native-fixed-colored-topology-DL-WGD-v1"
    )


def test_species_child_reordering_is_valid_but_length_changes_are_not(tmp_path):
    gene, species, fit = inputs(tmp_path)
    first, second = tmp_path / "first.tsv", tmp_path / "second.tsv"
    main(command(gene, species, fit) + ["-o", str(first)])
    species.write_text("(B_b:1,A_a:1);")
    main(command(gene, species, fit) + ["-o", str(second)])
    assert first.read_text() == second.read_text()
    species.write_text("(B_b:1,A_a:2);")
    with pytest.raises(ValueError, match="branch lengths"):
        main(command(gene, species, fit) + ["-o", str(second)])
    assert first.read_text() == second.read_text()


def test_cli_reuses_zero_retention_only_for_exact_background_parameters(
    tmp_path, monkeypatch
):
    from nwkit.wgd_tree_model import TopologyLikelihood

    gene, species, fit = inputs(tmp_path)
    record = json.loads(fit.read_text())
    (first,) = record["candidates"]
    first["event"]["event"]["retention"] = 0
    second = {
        "event": {
            **first["event"],
            "rates": [[0.15, 0.2]],
            "event": {**first["event"]["event"], "node": 2},
        }
    }
    record["candidates"].append(second)
    fit.write_text(json.dumps(record))
    evaluated = []
    original = TopologyLikelihood.evaluate

    def evaluate(self, rates, root, event=None, **options):
        evaluated.append(event.node if event is not None else None)
        return original(self, rates, root, event, **options)

    monkeypatch.setattr(TopologyLikelihood, "evaluate", evaluate)
    output = tmp_path / "origins.tsv"
    main(command(gene, species, fit) + ["-o", str(output)])
    rows = pd.read_csv(output, sep="\t")
    assert evaluated == [None, 2]
    assert list(rows.conditional_wgd_probability) == [0.0, 0.0]
    assert rows.conditional_log_likelihood_difference.iloc[0] == 0
    assert rows.conditional_log_likelihood_difference.iloc[1] != 0


def test_cli_fixed_background_event_follows_clade_not_preorder_number(tmp_path):
    gene, species, fit = inputs(tmp_path)
    record = json.loads(fit.read_text())
    record["fixed_background_events"] = [
        dict(
            species_event_id=record["species_event_ids"][2],
            branch_id=2,
            retention=0.6,
            event_fraction=0.3,
            multiplicity=2,
        )
    ]
    fit.write_text(json.dumps(record))
    first, second, report = (
        tmp_path / name for name in ("one.tsv", "two.tsv", "report.json")
    )
    main(command(gene, species, fit) + ["-o", str(first), "--model-out", str(report)])
    assert (
        json.loads(report.read_text())["fixed_background_events"]
        == record["fixed_background_events"]
    )
    assert (
        "latent_node_origin_at_candidate_given_topology_and_fixed_background"
        in first.read_text()
    )
    species.write_text("(B_b:1,A_a:1);")
    main(command(gene, species, fit) + ["-o", str(second)])
    assert first.read_bytes() == second.read_bytes()


def test_outputs_do_not_overwrite_inputs_and_transfers_are_not_absorbed(tmp_path):
    gene, species, fit = inputs(tmp_path)
    original = fit.read_text()
    with pytest.raises(ValueError):
        main(command(gene, species, fit) + ["-o", str(fit)])
    assert fit.read_text() == original
    gene.write_text("(A_a_1:1,A_a_2:1)[&&NHX:H=Y];")
    with pytest.raises(ValueError, match="transfers"):
        main(command(gene, species, fit))


def test_unknown_gene_species_is_not_lost_or_converted_to_zero(tmp_path):
    gene, species, fit = inputs(tmp_path)
    gene.write_text("(A_a_1:1,Unknown_x_1:1);")
    with pytest.raises(ValueError, match="every|Every"):
        main(command(gene, species, fit))


def test_count_model_can_be_the_single_stdin_input(tmp_path, monkeypatch, capsys):
    gene, species, fit = inputs(tmp_path)
    arguments = command(gene, species, fit)
    arguments[arguments.index("--count-model") + 1] = "-"
    monkeypatch.setattr("sys.stdin", StringIO(fit.read_text()))
    main(arguments + ["-o", "-"])
    table = pd.read_csv(StringIO(capsys.readouterr().out), sep="\t")
    assert len(table) == 1
    assert 0 < table.conditional_wgd_probability.iloc[0] < 1


def test_count_model_and_gene_tree_cannot_both_own_stdin(tmp_path):
    gene, species, fit = inputs(tmp_path)
    arguments = command(gene, species, fit)
    arguments[arguments.index("--count-model") + 1] = "-"
    arguments[arguments.index("-i") + 1] = "-"
    with pytest.raises(ValueError, match="only one input"):
        main(arguments)


@pytest.mark.parametrize("declaration, option", [("[&U]", "auto"), ("", "no")])
def test_species_rooting_is_respected_and_failed_inputs_preserve_outputs(
    tmp_path, declaration, option
):
    gene, species, fit = inputs(tmp_path)
    output, report = tmp_path / "origins.tsv", tmp_path / "likelihood.json"
    output.write_text("existing origins\n")
    report.write_text("existing report\n")
    species.write_text(declaration + "(A_a:1,B_b:1);")
    with pytest.raises(ValueError, match="rooted"):
        main(
            command(gene, species, fit)
            + [
                "--species-tree-rooted",
                option,
                "-o",
                str(output),
                "--model-out",
                str(report),
            ]
        )
    assert output.read_text() == "existing origins\n"
    assert report.read_text() == "existing report\n"

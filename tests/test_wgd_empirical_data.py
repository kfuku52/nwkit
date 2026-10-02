import importlib.util
from pathlib import Path

import pandas as pd
import pytest

spec = importlib.util.spec_from_file_location(
    "prepare_wgd_empirical_counts",
    Path(__file__).resolve().parents[1] / "tools" / "prepare_wgd_empirical_counts.py",
)
study = importlib.util.module_from_spec(spec)
spec.loader.exec_module(study)
prepare_counts = study.prepare_counts


def inputs(tmp_path):
    counts, tree = tmp_path / "raw.tsv", tmp_path / "tree.nwk"
    counts.write_text("A\tB\tD\n1\t0\t1\n0\t0\t2\n5\t1\t1\n1\t1\t0\n1\t1\t1\n")
    tree.write_text("(D:1,(A:1,B:1):1);")
    return counts, tree


def test_study_adapter_audits_root_selection_without_copy_number_filtering(tmp_path):
    counts, tree = inputs(tmp_path)
    output = tmp_path / "study"
    metadata = prepare_counts(counts, tree, output, 0, 620)
    assert metadata["raw_families"] == 5
    assert metadata["root_spanning_families"] == 3
    assert metadata["selected_families"] == 3
    assert metadata["selected_max_copy_number"] == 5
    table = pd.read_csv(output / "counts.tsv", sep="\t")
    assert list(table.family_id) == [
        "source_row_00001",
        "source_row_00003",
        "source_row_00005",
    ]
    assert metadata["known_event_locations_used"] is False
    audit = pd.read_csv(output / "family_selection.tsv", sep="\t")
    assert list(audit.selected) == [1, 0, 1, 0, 1]
    assert all(
        Path(path).exists()
        for path in [metadata["source_counts"], metadata["source_tree"]]
    )


def test_uniform_study_sampling_is_reproducible_and_existing_evidence_is_preserved(
    tmp_path,
):
    counts, tree = inputs(tmp_path)
    first, second = tmp_path / "first", tmp_path / "second"
    prepare_counts(counts, tree, first, 2, 620)
    prepare_counts(counts, tree, second, 2, 620)
    assert (first / "counts.tsv").read_bytes() == (second / "counts.tsv").read_bytes()
    assert (first / "family_selection.tsv").read_bytes() == (
        second / "family_selection.tsv"
    ).read_bytes()
    with pytest.raises(ValueError, match="new"):
        prepare_counts(counts, tree, first, 2, 620)
    assert counts.read_text().startswith("A\tB\tD\n")


@pytest.mark.parametrize("bad", ["NA", "-1", "0.5", "inf", str(2**63)])
def test_complete_count_study_rejects_missing_and_invalid_values(tmp_path, bad):
    counts, tree = inputs(tmp_path)
    counts.write_text(f"A\tB\tD\n{bad}\t1\t1\n")
    with pytest.raises(ValueError):
        prepare_counts(counts, tree, tmp_path / "output", 0, 620)
    assert not (tmp_path / "output").exists()


@pytest.mark.parametrize("count", [2**53 + 1, 2**63 - 1])
def test_study_adapter_preserves_integer_counts_without_float_rounding(tmp_path, count):
    counts, tree = inputs(tmp_path)
    counts.write_text(f"A\tB\tD\n{count}\t1\t1\n")
    output = tmp_path / "output"
    metadata = prepare_counts(counts, tree, output, 0, 620)
    result = pd.read_csv(output / "counts.tsv", sep="\t", dtype=str)
    assert result.A.iloc[0] == str(count)
    assert metadata["selected_max_copy_number"] == count

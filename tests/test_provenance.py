import hashlib
import io
import json
import os
import stat
import sys
from argparse import Namespace

import pytest

from nwkit import provenance as provenance_mod
from nwkit.cli import main
from nwkit.provenance import _argument_dict, _TeeTextWriter


def _sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_argument_dict_omits_private_runtime_caches():
    args = Namespace(
        command="constrain",
        species_parser="legacy",
        _nwkit_species_parser_cache=object(),
        _nwkit_species_parser_cache_key=("legacy", "pattern", None),
        handler=object(),
    )
    assert _argument_dict(args) == {
        "command": "constrain",
        "species_parser": "legacy",
    }


def test_global_audit_records_arguments_hashes_and_messages(tmp_path):
    infile = tmp_path / "input.nwk"
    outfile = tmp_path / "output.nwk"
    audit = tmp_path / "audit.jsonl"
    infile.write_text("((A:1,B:1):1,C:1);")
    main(
        [
            "label",
            "-i",
            str(infile),
            "-o",
            str(outfile),
            "--target",
            "intnode",
            "--audit",
            str(audit),
        ]
    )
    record = json.loads(audit.read_text().strip())
    assert record["schema"] == "nwkit-audit-v1"
    assert record["status"] == "ok"
    assert record["command"] == "label"
    assert record["inputs"][0]["sha256"] == _sha256(infile)
    assert record["outputs"][0]["sha256"] == _sha256(outfile)
    assert record["primary_input"]["kind"] == "newick"
    assert record["primary_input"]["first_tree_tip_count"] == 3
    assert any(
        "Number of labeled target nodes" in message for message in record["messages"]
    )


def test_audit_does_not_treat_non_path_option_values_as_input_files(
    monkeypatch, tmp_path
):
    infile = tmp_path / "input.nwk"
    audit = tmp_path / "audit.jsonl"
    infile.write_text("(A:1,B:1);")
    (tmp_path / "auto").write_text("not an input\n")
    monkeypatch.chdir(tmp_path)

    main(
        [
            "info",
            "-i",
            str(infile),
            "--format",
            "auto",
            "--audit",
            str(audit),
        ]
    )

    record = json.loads(audit.read_text())
    assert {item["path"] for item in record["inputs"]} == {str(infile.resolve())}


def test_audit_hashes_auxiliary_standard_input(monkeypatch, tmp_path):
    infile = tmp_path / "input.nwk"
    outfile = tmp_path / "output.nwk"
    report = tmp_path / "sample.tsv"
    audit = tmp_path / "audit.jsonl"
    infile.write_text("((A:1,B:1):1,C:1);")
    trait_text = "leaf_name\tscore\nA\t1\nB\t2\nC\t3\n"
    monkeypatch.setattr(sys, "stdin", io.StringIO(trait_text))

    main(
        [
            "sample",
            "-i",
            str(infile),
            "--trait",
            "-",
            "-n",
            "1",
            "-o",
            str(outfile),
            "--report",
            str(report),
            "--audit",
            str(audit),
        ]
    )

    record = json.loads(audit.read_text().strip())
    assert record["stdin"]["argument"] == "trait"
    assert record["stdin"]["sha256"] == hashlib.sha256(trait_text.encode()).hexdigest()
    assert record["primary_input"]["first_tree_tip_count"] == 3


def test_audit_records_skim_group_table_outputs(tmp_path):
    infile = tmp_path / "input.nwk"
    trait = tmp_path / "trait.tsv"
    outfile = tmp_path / "output.nwk"
    prefix = tmp_path / "groups"
    audit = tmp_path / "audit.jsonl"
    infile.write_text("((A:1,B:1):1,(C:1,D:1):1);")
    trait.write_text("leaf_name\tgroup\nA\tx\nB\tx\nC\ty\nD\ty\n")

    main(
        [
            "skim",
            "-i",
            str(infile),
            "--trait",
            str(trait),
            "--group-by",
            "group",
            "-o",
            str(outfile),
            "--group-table-prefix",
            str(prefix),
            "--audit",
            str(audit),
        ]
    )

    record = json.loads(audit.read_text().strip())
    output_paths = {item["path"] for item in record["outputs"]}
    assert str(outfile.resolve()) in output_paths
    assert str((tmp_path / "groups.all.tsv").resolve()) in output_paths
    assert str((tmp_path / "groups.sampled.tsv").resolve()) in output_paths


def test_audit_records_replicate_contrast_inputs_and_auxiliary_outputs(tmp_path):
    infile = tmp_path / "gene.nwk"
    trait = tmp_path / "expression.tsv"
    outfile = tmp_path / "contrasts.tsv"
    covariance = tmp_path / "sampling-covariance.tsv"
    summary = tmp_path / "tip-summary.tsv"
    audit = tmp_path / "audit.jsonl"
    infile.write_text("((A:1,B:1):1,C:2);")
    trait.write_text(
        "leaf_name\tsample\texpression\n"
        "A\ta1\t1\nA\ta2\t3\n"
        "B\tb1\t4\nB\tb2\t6\n"
        "C\tc1\t7\nC\tc2\t9\n"
    )

    main(
        [
            "contrast",
            "--infile",
            str(infile),
            "--trait",
            str(trait),
            "--columns",
            "expression",
            "--biological-id",
            "sample",
            "--sampling-covariance-out",
            str(covariance),
            "--tip-summary-out",
            str(summary),
            "--outfile",
            str(outfile),
            "--audit",
            str(audit),
        ]
    )

    record = json.loads(audit.read_text())
    input_paths = {item["path"] for item in record["inputs"]}
    output_paths = {item["path"] for item in record["outputs"]}
    assert input_paths == {str(infile.resolve()), str(trait.resolve())}
    assert output_paths == {
        str(outfile.resolve()),
        str(covariance.resolve()),
        str(summary.resolve()),
    }


def test_audit_records_pgls_sampling_covariance_and_random_effect_output(tmp_path):
    response = tmp_path / "response.tsv"
    predictor = tmp_path / "predictor.tsv"
    covariance = tmp_path / "sampling-covariance.tsv"
    outfile = tmp_path / "pgls.tsv"
    random_effects = tmp_path / "random-effects.tsv"
    audit = tmp_path / "audit.jsonl"
    response.write_text(
        "tree_id\tgene_clade_id\tlineage_clade_id\tevent_type\teligible\t"
        "coverage_status\tspecies_event_id\tspecies_event_taxa\t"
        "species_numerator_event_id\tspecies_denominator_event_id\ttrait\t"
        "evolution_model\tevolution_parameter_name\tevolution_parameter\t"
        "branch_length_mode\t"
        "raw_contrast\tcontrast_variance\tsampling_variance\n"
        "OG1\tg1\tl1\tspeciation\tyes\tcomplete\te1\tt1\tn1\td1\texpression\tbrownian\t\t\toriginal\t2\t1\t0.1\n"
        "OG1\tg2\tl1\tspeciation\tyes\tcomplete\te2\tt2\tn2\td2\texpression\tbrownian\t\t\toriginal\t4\t1\t0.1\n"
        "OG1\tg3\tl1\tspeciation\tyes\tcomplete\te3\tt3\tn3\td3\texpression\tbrownian\t\t\toriginal\t6\t1\t0.1\n"
    )
    predictor.write_text(
        "tree_id\tbranch_clade_id\tdescendant_taxa\tnumerator_clade_id\t"
        "denominator_clade_id\ttrait\tevolution_model\t"
        "evolution_parameter_name\tevolution_parameter\tbranch_length_mode\t"
        "raw_contrast\n"
        "species\te1\tt1\tn1\td1\tbody_size\tbrownian\t\t\toriginal\t1\n"
        "species\te2\tt2\tn2\td2\tbody_size\tbrownian\t\t\toriginal\t2\n"
        "species\te3\tt3\tn3\td3\tbody_size\tbrownian\t\t\toriginal\t3\n"
    )
    covariance.write_text(
        "tree_id\ttrait\tcontrast_id_1\tcontrast_id_2\tsampling_covariance\n"
        "OG1\texpression\tg1\tg1\t0.1\n"
        "OG1\texpression\tg1\tg2\t0\n"
        "OG1\texpression\tg1\tg3\t0\n"
        "OG1\texpression\tg2\tg2\t0.1\n"
        "OG1\texpression\tg2\tg3\t0\n"
        "OG1\texpression\tg3\tg3\t0.1\n"
    )

    main(
        [
            "regress",
            "--response-contrasts",
            str(response),
            "--predictor-contrasts",
            str(predictor),
            "--response-sampling-covariance",
            str(covariance),
            "--responses",
            "expression",
            "--predictors",
            "body_size",
            "--random-effects-out",
            str(random_effects),
            "--outfile",
            str(outfile),
            "--audit",
            str(audit),
        ]
    )

    record = json.loads(audit.read_text())
    input_paths = {item["path"] for item in record["inputs"]}
    output_paths = {item["path"] for item in record["outputs"]}
    assert input_paths == {
        str(response.resolve()),
        str(predictor.resolve()),
        str(covariance.resolve()),
    }
    assert output_paths == {str(outfile.resolve()), str(random_effects.resolve())}


def test_audit_records_compose_manifest_dependencies_and_all_roles(tmp_path):
    target = tmp_path / "target.nwk"
    source = tmp_path / "source.nwk"
    manifest = tmp_path / "compose.json"
    outfile = tmp_path / "output.nwk"
    audit = tmp_path / "audit.jsonl"
    target.write_text("((A:1,B:1):1,C:1);")
    source.write_text("((A:1,B:1)AB:1,C:1);")
    manifest.write_text(json.dumps({"name": "source.nwk", "support": "source.nwk"}))
    main(
        [
            "compose",
            "-i",
            str(target),
            "--manifest",
            str(manifest),
            "-o",
            str(outfile),
            "--audit",
            str(audit),
        ]
    )
    record = json.loads(audit.read_text().strip())
    source_record = next(
        item for item in record["inputs"] if item["path"] == str(source.resolve())
    )
    assert source_record["arguments"] == ["manifest:name", "manifest:support"]


def test_audit_records_draw_tip_image_manifest_and_assets(tmp_path):
    Image = pytest.importorskip("PIL.Image")
    infile = tmp_path / "input.nwk"
    image = tmp_path / "tip.png"
    manifest = tmp_path / "manifest.tsv"
    outfile = tmp_path / "tree.svg"
    audit = tmp_path / "audit.jsonl"
    infile.write_text("(A:1,B:1);")
    Image.new("RGBA", (8, 8), (40, 120, 200, 255)).save(image)
    manifest.write_text("leaf_name\tlocal_path\nA\ttip.png\nB\ttip.png\n")

    main(
        [
            "draw",
            "-i",
            str(infile),
            "-o",
            str(outfile),
            "--species-overlap-node-plot",
            "no",
            "--tip-image-manifest",
            str(manifest),
            "--audit",
            str(audit),
        ]
    )

    record = json.loads(audit.read_text().strip())
    records_by_path = {item["path"]: item for item in record["inputs"]}
    assert records_by_path[str(manifest.resolve())]["argument"] == "tip_image_manifest"
    assert records_by_path[str(image.resolve())]["arguments"] == [
        "tip_image_manifest:asset"
    ]


def test_audit_records_densitree_tree_collection(tmp_path):
    infile = tmp_path / "reference.nwk"
    samples = tmp_path / "posterior.trees"
    outfile = tmp_path / "densitree.svg"
    audit = tmp_path / "audit.jsonl"
    infile.write_text("((A:1,B:1):1,(C:1,D:1):1);")
    samples.write_text("((A:1,B:1):1,(C:1,D:1):1);\n((A:1,C:1):1,(B:1,D:1):1);\n")

    main(
        [
            "draw",
            "-i",
            str(infile),
            "-o",
            str(outfile),
            "--species-overlap-node-plot",
            "no",
            "--densitree-trees",
            str(samples),
            "--densitree",
            "all",
            "--audit",
            str(audit),
        ]
    )

    record = json.loads(audit.read_text().strip())
    records_by_path = {item["path"]: item for item in record["inputs"]}
    assert records_by_path[str(samples.resolve())]["argument"] == "densitree_trees"


def test_stderr_capture_is_bounded():
    writer = _TeeTextWriter(io.StringIO(), capture=True, max_lines=3, max_line_chars=10)
    writer.write("warning one\nline two\nline three\nline four\n")
    assert writer.captured_lines == ["warning on", "line two", "line three"]
    assert writer.warning_lines == ["warning on"]
    writer.write("final warning")
    assert writer.captured_warning_lines == ["warning on", "final warn"]


def test_audit_path_cannot_overwrite_primary_output(tmp_path):
    infile = tmp_path / "input.nwk"
    output = tmp_path / "same.out"
    infile.write_text("(A:1,B:1);")
    with pytest.raises(ValueError, match="Output paths must be distinct"):
        main(["label", "-i", str(infile), "-o", str(output), "--audit", str(output)])


def test_audit_path_cannot_modify_an_input_file(tmp_path):
    infile = tmp_path / "input.nwk"
    outfile = tmp_path / "output.nwk"
    original = "(A:1,B:1);"
    infile.write_text(original)
    infile.chmod(0o640)

    with pytest.raises(ValueError, match="distinct from input"):
        main(
            [
                "label",
                "-i",
                str(infile),
                "-o",
                str(outfile),
                "--audit",
                str(infile),
            ]
        )

    assert infile.read_text() == original
    if os.name != "nt":
        assert stat.S_IMODE(infile.stat().st_mode) == 0o640
    assert not outfile.exists()


def test_audit_path_cannot_create_a_declared_missing_input(tmp_path):
    missing_infile = tmp_path / "future-input.nwk"

    with pytest.raises(ValueError, match="distinct from input"):
        main(
            [
                "info",
                "-i",
                str(missing_infile),
                "--audit",
                str(missing_infile),
            ]
        )

    assert not missing_infile.exists()


def test_missing_input_audit_collision_respects_case_insensitive_filesystems(
    monkeypatch,
    tmp_path,
):
    monkeypatch.setattr(
        "nwkit.util._filesystem_is_case_insensitive",
        lambda path: True,
    )
    missing_infile = tmp_path / "Future-Input.nwk"
    audit = tmp_path / "future-input.nwk"

    with pytest.raises(ValueError, match="distinct from input"):
        main(
            [
                "info",
                "-i",
                str(missing_infile),
                "--audit",
                str(audit),
            ]
        )

    assert not audit.exists()


def test_audit_hard_link_cannot_modify_an_input_file(tmp_path):
    infile = tmp_path / "input.nwk"
    audit_link = tmp_path / "audit.jsonl"
    outfile = tmp_path / "output.nwk"
    original = "(A:1,B:1);"
    infile.write_text(original)
    infile.chmod(0o640)
    os.link(infile, audit_link)

    with pytest.raises(ValueError, match="distinct from input"):
        main(
            [
                "label",
                "-i",
                str(infile),
                "-o",
                str(outfile),
                "--audit",
                str(audit_link),
            ]
        )

    assert infile.read_text() == original
    if os.name != "nt":
        assert stat.S_IMODE(infile.stat().st_mode) == 0o640
    assert not outfile.exists()


def test_audit_path_cannot_modify_draw_manifest_asset(tmp_path):
    infile = tmp_path / "input.nwk"
    manifest = tmp_path / "manifest.tsv"
    asset = tmp_path / "tip.svg"
    outfile = tmp_path / "output.svg"
    infile.write_text("(A:1,B:1);")
    asset.write_text('<svg xmlns="http://www.w3.org/2000/svg"/>')
    asset.chmod(0o640)
    original_asset = asset.read_bytes()
    manifest.write_text("leaf_name\tlocal_path\nA\ttip.svg\n")

    with pytest.raises(ValueError, match="distinct from input"):
        main(
            [
                "draw",
                "-i",
                str(infile),
                "-o",
                str(outfile),
                "--species-overlap-node-plot",
                "no",
                "--tip-image-manifest",
                str(manifest),
                "--audit",
                str(asset),
            ]
        )

    assert asset.read_bytes() == original_asset
    if os.name != "nt":
        assert stat.S_IMODE(asset.stat().st_mode) == 0o640
    assert not outfile.exists()


def test_audit_snapshots_input_before_in_place_output(tmp_path):
    tree_path = tmp_path / "tree.nwk"
    audit = tmp_path / "audit.jsonl"
    original = "((A:1,B:1):1,C:1);"
    tree_path.write_text(original)
    original_digest = hashlib.sha256(original.encode()).hexdigest()

    main(
        [
            "label",
            "-i",
            str(tree_path),
            "-o",
            str(tree_path),
            "--target",
            "intnode",
            "--audit",
            str(audit),
        ]
    )

    record = json.loads(audit.read_text())
    input_record = next(
        item for item in record["inputs"] if item["path"] == str(tree_path.resolve())
    )
    output_record = next(
        item for item in record["outputs"] if item["path"] == str(tree_path.resolve())
    )
    assert input_record["sha256"] == original_digest
    assert output_record["sha256"] == _sha256(tree_path)
    assert output_record["sha256"] != original_digest


@pytest.mark.skipif(
    os.name == "nt",
    reason="POSIX file modes are unavailable on Windows",
)
def test_new_audit_file_is_private(tmp_path):
    infile = tmp_path / "input.nwk"
    audit = tmp_path / "audit.jsonl"
    infile.write_text("(A:1,B:1);")

    main(["info", "-i", str(infile), "--audit", str(audit)])

    assert stat.S_IMODE(audit.stat().st_mode) == 0o600


@pytest.mark.skipif(
    os.name == "nt",
    reason="POSIX file modes are unavailable on Windows",
)
def test_existing_audit_permissions_are_tightened(tmp_path):
    infile = tmp_path / "input.nwk"
    audit = tmp_path / "audit.jsonl"
    infile.write_text("(A:1,B:1);")
    audit.write_text('{"existing": true}\n')
    audit.chmod(0o666)

    main(["info", "-i", str(infile), "--audit", str(audit)])

    assert stat.S_IMODE(audit.stat().st_mode) == 0o600
    assert len(audit.read_text().splitlines()) == 2


def test_audit_symlink_is_rejected_before_handler_runs(tmp_path):
    infile = tmp_path / "input.nwk"
    outfile = tmp_path / "output.nwk"
    audit_target = tmp_path / "audit-target.jsonl"
    audit_link = tmp_path / "audit.jsonl"
    infile.write_text("(A:1,B:1);")
    audit_target.write_text("unchanged\n")
    audit_link.symlink_to(audit_target)

    with pytest.raises(ValueError, match="symbolic link"):
        main(
            [
                "label",
                "-i",
                str(infile),
                "-o",
                str(outfile),
                "--audit",
                str(audit_link),
            ]
        )

    assert not outfile.exists()
    assert audit_target.read_text() == "unchanged\n"


def test_audit_directory_is_rejected_before_handler_runs(tmp_path):
    infile = tmp_path / "input.nwk"
    outfile = tmp_path / "output.nwk"
    audit_dir = tmp_path / "audit"
    infile.write_text("(A:1,B:1);")
    audit_dir.mkdir()

    with pytest.raises(ValueError, match="regular file"):
        main(
            [
                "label",
                "-i",
                str(infile),
                "-o",
                str(outfile),
                "--audit",
                str(audit_dir),
            ]
        )

    assert not outfile.exists()


@pytest.mark.skipif(not hasattr(os, "mkfifo"), reason="FIFOs are not available")
def test_audit_directory_hash_rejects_fifo_without_opening_it(tmp_path):
    output_dir = tmp_path / "output"
    output_dir.mkdir()
    fifo_path = output_dir / "blocking.fifo"
    os.mkfifo(fifo_path)

    with pytest.raises(ValueError, match="only supports regular files"):
        provenance_mod._sha256_directory(output_dir)


@pytest.mark.skipif(not hasattr(os, "mkfifo"), reason="FIFOs are not available")
def test_audit_rejects_declared_fifo_input_before_handler_runs(tmp_path):
    fifo_path = tmp_path / "input.fifo"
    audit = tmp_path / "audit.jsonl"
    os.mkfifo(fifo_path)
    handler_called = False

    def handler(args):
        nonlocal handler_called
        handler_called = True

    with pytest.raises(ValueError, match="regular files"):
        provenance_mod.run_with_audit(
            args=Namespace(
                command="test",
                audit=str(audit),
                infile=str(fifo_path),
            ),
            argv=["test"],
            handler=handler,
        )

    assert handler_called is False
    assert not audit.exists()


def test_requested_audit_write_failure_is_not_silenced(monkeypatch, tmp_path):
    infile = tmp_path / "input.nwk"
    audit = tmp_path / "audit.jsonl"
    infile.write_text("(A:1,B:1);")

    def fail_write(path, record):
        raise OSError("audit disk failed")

    monkeypatch.setattr(provenance_mod, "_write_audit", fail_write)
    with pytest.raises(OSError, match="audit disk failed"):
        main(["info", "-i", str(infile), "--audit", str(audit)])


def test_audit_failure_does_not_mask_original_handler_error(
    monkeypatch, tmp_path, capsys
):
    audit = tmp_path / "audit.jsonl"
    args = Namespace(command="test", audit=str(audit), infile=None)

    def fail_handler(parsed_args):
        raise ValueError("original handler failure")

    def fail_write(path, record):
        raise OSError("audit disk failed")

    monkeypatch.setattr(provenance_mod, "_write_audit", fail_write)

    with pytest.raises(ValueError, match="original handler failure"):
        provenance_mod.run_with_audit(
            args=args,
            argv=["test"],
            handler=fail_handler,
        )

    assert "audit disk failed" in capsys.readouterr().err


def test_audit_spools_standard_input_with_bounded_reads(monkeypatch, tmp_path):
    class GuardedInput(io.StringIO):
        def read(self, size=-1):
            assert size > 0
            return super().read(size)

    infile = tmp_path / "input.nwk"
    outfile = tmp_path / "output.nwk"
    report = tmp_path / "report.tsv"
    audit = tmp_path / "audit.jsonl"
    infile.write_text("((A:1,B:1):1,C:1);")
    trait_text = "leaf_name\tscore\nA\t1\nB\t2\nC\t3\n"
    monkeypatch.setattr(sys, "stdin", GuardedInput(trait_text))

    main(
        [
            "sample",
            "-i",
            str(infile),
            "--trait",
            "-",
            "-n",
            "1",
            "-o",
            str(outfile),
            "--report",
            str(report),
            "--audit",
            str(audit),
        ]
    )

    record = json.loads(audit.read_text())
    assert record["stdin"]["bytes"] == len(trait_text.encode())
    assert record["stdin"]["sha256"] == hashlib.sha256(trait_text.encode()).hexdigest()


@pytest.mark.parametrize(
    "derived_name", ["manifest.tsv", "ATTRIBUTION.md", "unmatched.tsv"]
)
def test_audit_rejects_image_derived_output_collisions(tmp_path, derived_name):
    infile = tmp_path / "input.nwk"
    out_dir = tmp_path / "images"
    infile.write_text("(A:1,B:1);")

    with pytest.raises(ValueError, match="Output paths must be distinct"):
        main(
            [
                "image",
                "-i",
                str(infile),
                "--out-dir",
                str(out_dir),
                "--audit",
                str(out_dir / derived_name),
            ]
        )

    assert not out_dir.exists()


def test_audit_rejects_path_nested_under_image_output_directory(tmp_path):
    infile = tmp_path / "input.nwk"
    out_dir = tmp_path / "images"
    infile.write_text("(A:1,B:1);")

    with pytest.raises(ValueError, match="must not be inside output directory"):
        main(
            [
                "image",
                "-i",
                str(infile),
                "--out-dir",
                str(out_dir),
                "--audit",
                str(out_dir / "audit.jsonl"),
            ]
        )

    assert not out_dir.exists()

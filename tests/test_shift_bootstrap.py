"""Bootstrap mapping, failure accounting and real parametric refit checks."""

import csv
import json
import os

import pytest

from nwkit.cli import main
from tests import test_shift as shift_support

shift_inputs = shift_support.shift_inputs


def bootstrap_backend(directory, args):
    result = shift_support.fake_backend(directory, args)
    (directory / "bootstrap.tsv").write_text("attempted\tsuccessful\tfailed\n4\t3\t1\n")
    # Complementary root-child shifts yield the same unlabeled tip partition.
    (directory / "bootstrap-configurations.tsv").write_text(
        "success_index\tclades\n1\tt0/t1\n2\tt2/t3\n3\t\n"
    )
    with (directory / "bootstrap-failures.tsv").open("w", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t")
        writer.writerow(["message", "count"])
        writer.writerow(['fit failed: "example"\twith newline\nhere', 1])
    return result


def test_bootstrap_support_and_partition_equivalence(
    shift_inputs, tmp_path, monkeypatch
):
    monkeypatch.setattr("nwkit.shift.run_backend", bootstrap_backend)
    main([*shift_inputs, "--bootstrap", "4"])
    boot = json.loads((tmp_path / "model.json").read_text())["bootstrap"]
    assert (boot["attempted"], boot["successful"], boot["failed"]) == (4, 3, 1)
    assert boot["successful_configurations"] == [[1], [2], []]
    assert [r["frequency"] for r in boot["configuration_frequencies"]] == [1 / 3] * 3
    assert boot["tip_partition_frequencies"][0] == {
        "groups": [["A", "B"], ["C", "D"]],
        "count": 2,
        "frequency": 2 / 3,
    }
    assert [r["count"] for r in boot["edge_inclusion"]] == [1, 1, 0, 0, 0, 0]
    assert boot["failure_messages"][0]["message"].endswith("newline\nhere")


@pytest.mark.parametrize(
    "filename,text",
    [
        ("bootstrap.tsv", "attempted\tsuccessful\tfailed\n4\t0\t4\n"),
        ("bootstrap.tsv", "attempted\tsuccessful\tfailed\n5\t3\t2\n"),
        ("bootstrap-configurations.tsv", "success_index\tclades\n1\tt0/t1;t0/t1\n"),
        ("bootstrap-configurations.tsv", "success_index\tclades\n1\tt0/t1/t2/t3\n"),
        ("bootstrap-configurations.tsv", "success_index\tclades\n1\tunknown\n"),
        ("bootstrap-configurations.tsv", "success_index\tclades\n1\t\n"),
        ("bootstrap-failures.tsv", "message\tcount\nerror\t2\n"),
    ],
)
def test_invalid_bootstrap_preserves_outputs(
    shift_inputs, tmp_path, monkeypatch, filename, text
):
    def backend(directory, args):
        result = bootstrap_backend(directory, args)
        (directory / filename).write_text(text)
        return result

    monkeypatch.setattr("nwkit.shift.run_backend", backend)
    model = tmp_path / "model.json"
    model.write_text("original")
    with pytest.raises(ValueError, match="[Bb]ootstrap"):
        main([*shift_inputs, "--bootstrap", "4"])
    assert model.read_text() == "original"


@pytest.mark.parametrize(
    "option,value",
    [("--bootstrap", "-1"), ("--bootstrap", "2147483648"), ("--bootstrap-seed", "-1")],
)
def test_bootstrap_options_preflight(shift_inputs, monkeypatch, option, value):
    def backend(*args):
        pytest.fail("Invalid options must fail before R runs")

    monkeypatch.setattr("nwkit.shift.run_backend", backend)
    with pytest.raises(ValueError, match="bootstrap"):
        main([*shift_inputs, option, value])


@pytest.mark.integration
@pytest.mark.skipif(
    not os.environ.get("NWKIT_TEST_RSCRIPT"), reason="Optional R backend"
)
@pytest.mark.parametrize("with_error", [False, True])
def test_real_bootstrap_reproducible(shift_inputs, tmp_path, with_error):
    extra = []
    if with_error:
        (tmp_path / "trait.tsv").write_text(
            "leaf_name\tx\tse\nD\t4\t0.1\nB\t2\t0.2\nA\t1\t0\nC\t3\t0.3\n"
        )
        extra = ["--standard-error-column", "se"]
    command = [
        *shift_inputs,
        "--rscript",
        os.environ["NWKIT_TEST_RSCRIPT"],
        "--max-shifts",
        "1",
        "--search-strategy",
        "exhaustive",
        "--criterion",
        "BIC",
        "--bootstrap",
        "3",
        "--bootstrap-seed",
        "43",
        *extra,
    ]
    main(command)
    first = json.loads((tmp_path / "model.json").read_text())["bootstrap"]
    main(command)
    second = json.loads((tmp_path / "model.json").read_text())["bootstrap"]
    assert first == second
    assert first["successful"] == 3
    assert first["failed"] == 0
    assert sum(row["count"] for row in first["configuration_frequencies"]) == 3


def test_missing_failure_message_remains_accounted_for(
    shift_inputs, tmp_path, monkeypatch, capsys
):
    def backend(directory, args):
        result = bootstrap_backend(directory, args)
        (directory / "bootstrap-failures.tsv").write_text("message\tcount\n")
        return result

    monkeypatch.setattr("nwkit.shift.run_backend", backend)
    main([*shift_inputs, "--bootstrap", "4"])
    boot = json.loads((tmp_path / "model.json").read_text())["bootstrap"]
    assert boot["failures_without_message"] == 1
    assert boot["successful"] == 3
    assert "1/4 bootstrap refits failed" in capsys.readouterr().err


@pytest.mark.parametrize(
    "text", ['message\tcount\n"unterminated\t1\n', 'message\tcount\n"error"junk\t1\n']
)
def test_malformed_diagnostic_quotes_rejected(shift_inputs, monkeypatch, text):
    def backend(directory, args):
        result = bootstrap_backend(directory, args)
        (directory / "bootstrap-failures.tsv").write_text(text)
        return result

    monkeypatch.setattr("nwkit.shift.run_backend", backend)
    with pytest.raises(ValueError, match="quoting"):
        main([*shift_inputs, "--bootstrap", "4"])


@pytest.mark.integration
@pytest.mark.skipif(
    not os.environ.get("NWKIT_TEST_RSCRIPT"), reason="Optional R backend"
)
def test_all_small_tree_partitions_against_public_r_api(tmp_path):
    import subprocess

    from nwkit.shift_bootstrap import _partition
    from nwkit.util import assign_branch_ids, read_tree

    tree = read_tree("((A:1,B:1):1,(C:1,D:1):1);", "auto", True, quiet=True)
    ids = assign_branch_ids(tree)
    mapping = {"/".join(sorted(n.leaf_names())): ids[n] for n in ids}
    script = tmp_path / "partitions.R"
    script.write_text(r"""
tr <- ape::read.tree(text="((A:1,B:1):1,(C:1,D:1):1);")
keys <- vapply(tr$edge[,2], function(child) {
    tips <- if(child <= length(tr$tip.label)) tr$tip.label[child] else
        ape::extract.clade(tr, child)$tip.label
    paste(sort(tips), collapse="/")
}, "")
configs <- c(list(integer()), as.list(seq_len(nrow(tr$edge))),
             combn(seq_len(nrow(tr$edge)), 2, simplify=FALSE))
rows <- lapply(configs, function(config) {
    partition <- kfl1ou::shift_tip_partition(tr, config)
    data.frame(clades=paste(keys[config], collapse=";"),
        groups=paste(sort(vapply(split(names(partition), partition),
            function(tips) paste(sort(tips), collapse="/"), "")), collapse=";"))
})
write.table(do.call(rbind, rows), "partitions.tsv", sep="\t", quote=FALSE, row.names=FALSE)
""")
    subprocess.run(
        [os.environ["NWKIT_TEST_RSCRIPT"], "--vanilla", str(script)],
        cwd=tmp_path,
        check=True,
        capture_output=True,
    )
    with (tmp_path / "partitions.tsv").open() as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    assert len(rows) == 22
    for row in rows:
        selected = {mapping[key] for key in row["clades"].split(";") if key}
        expected = tuple(tuple(group.split("/")) for group in row["groups"].split(";"))
        assert _partition(tree, ids, selected) == expected


@pytest.mark.integration
@pytest.mark.skipif(
    not os.environ.get("NWKIT_TEST_RSCRIPT"), reason="Optional R backend"
)
@pytest.mark.parametrize("root_model", ["OUfixedRoot", "OUrandomRoot"])
def test_support_matches_direct_public_bootstrap(shift_inputs, tmp_path, root_model):
    import subprocess

    (tmp_path / "trait.tsv").write_text(
        "leaf_name\tx\tse\nD\t4\t0.1\nB\t2\t0.2\nA\t1\t0\nC\t3\t0.3\n"
    )
    main(
        [
            *shift_inputs,
            "--rscript",
            os.environ["NWKIT_TEST_RSCRIPT"],
            "--max-shifts",
            "1",
            "--criterion",
            "BIC",
            "--search-strategy",
            "exhaustive",
            "--root-model",
            root_model,
            "--standard-error-column",
            "se",
            "--bootstrap",
            "5",
            "--bootstrap-seed",
            "87",
            "--fit-out",
            str(tmp_path / "fit.rds"),
        ]
    )
    model = json.loads((tmp_path / "model.json").read_text())
    script = tmp_path / "reference.R"
    script.write_text(r"""
fit <- readRDS("fit.rds")
stopifnot(fit$l1ou.options$max.nShifts == 1,
          fit$l1ou.options$criterion == "BIC",
          !is.null(fit$l1ou.options$input_error))
boot <- kfl1ou::l1ou_bootstrap_support(fit, nItrs=5, seed=87, type="parametric")
tr <- fit$tree
keys <- vapply(tr$edge[,2], function(child) {
    tips <- if(child <= length(tr$tip.label)) tr$tip.label[child] else
        ape::extract.clade(tr, child)$tip.label
    paste(sort(tips), collapse="/")
}, "")
write.table(data.frame(clades=keys, rate=boot$detection.rate), "rates.tsv",
            sep="\t", row.names=FALSE, quote=FALSE)
""")
    subprocess.run(
        [os.environ["NWKIT_TEST_RSCRIPT"], "--vanilla", str(script)],
        cwd=tmp_path,
        check=True,
        capture_output=True,
    )
    from nwkit.util import assign_branch_ids, read_tree

    tree = read_tree(model["analysis_tree_tokens"], "auto", True, quiet=True)
    ids = assign_branch_ids(tree)
    mapping = {"/".join(sorted(n.leaf_names())): ids[n] for n in ids}
    rates = {
        row["branch_id"]: row["frequency"]
        for row in model["bootstrap"]["edge_inclusion"]
    }
    with (tmp_path / "rates.tsv").open() as stream:
        for row in csv.DictReader(stream, delimiter="\t"):
            assert rates[mapping[row["clades"]]] == pytest.approx(float(row["rate"]))


@pytest.mark.integration
@pytest.mark.skipif(
    not os.environ.get("NWKIT_TEST_RSCRIPT"), reason="Optional R backend"
)
def test_actual_r_bootstrap_failure_preserves_all_outputs(
    shift_inputs, tmp_path, monkeypatch, capsys
):
    from nwkit import shift_backend

    script = shift_backend.R_SCRIPT.replace(
        "boot <- if (convergence.enabled)",
        'stop("all bootstrap replicates failed.")\n    boot <- if (convergence.enabled)',
    )
    monkeypatch.setattr(shift_backend, "R_SCRIPT", script)
    paths = [tmp_path / "model.json", tmp_path / "fit.rds", tmp_path / "regimes.tsv"]
    for path in paths:
        path.write_bytes(b"original")
    with pytest.raises(RuntimeError, match="all bootstrap replicates failed"):
        main(
            [
                *shift_inputs,
                "--rscript",
                os.environ["NWKIT_TEST_RSCRIPT"],
                "--max-shifts",
                "0",
                "--bootstrap",
                "1",
                "--fit-out",
                str(paths[1]),
                "-o",
                str(paths[2]),
            ]
        )
    assert all(path.read_bytes() == b"original" for path in paths)
    assert capsys.readouterr().out == ""

import io

import pandas as pd
import pytest

from nwkit.cli import main
from nwkit.tree_outputs import write_tree_with_tables
from nwkit.util import read_tree
from tests.helpers import make_args


@pytest.mark.parametrize(
    "command", ["sample", "skim", "annotate", "transfer", "compose"]
)
def test_failed_tree_output_preserves_companion_files(tmp_path, command):
    table = tmp_path / "input.tsv"
    table.write_text("leaf_name\ttrait\nA\tx\nB\tx\n")
    source = tmp_path / "source.nwk"
    source.write_text("(A:1,B:1)R;")
    report = tmp_path / "report.tsv"
    outputs = [report]
    options = {
        "sample": ["--n", "1", "--report", str(report)],
        "skim": ["--group-table-prefix", str(tmp_path / "group")],
        "annotate": ["--table", str(table), "--report", str(report)],
        "transfer": [
            "--infile2",
            str(source),
            "--format2",
            "1",
            "--name",
            "yes",
            "--report",
            str(report),
        ],
        "compose": [
            "--name-source",
            str(source),
            "--source-format",
            "1",
            "--report",
            str(report),
        ],
    }
    if command == "skim":
        outputs = [tmp_path / "group.all.tsv", tmp_path / "group.sampled.tsv"]
    for output in outputs:
        output.write_text("previous output\n")
    with pytest.raises(OSError):
        main(
            [
                command,
                "--infile",
                str(source),
                "--format",
                "1",
                "--outfile",
                str(tmp_path / "missing" / "out.nwk"),
                *options[command],
            ]
        )
    assert all(output.read_text() == "previous output\n" for output in outputs)


def test_tree_and_table_rollback_after_install_failure(tmp_path, monkeypatch):
    from nwkit import output_transaction as transaction

    tree_path, table_path = tmp_path / "tree.nwk", tmp_path / "table.tsv"
    for path in (tree_path, table_path):
        path.write_text("previous\n")
    original = transaction.replace_output
    calls = 0

    def fail_second(source, target):
        nonlocal calls
        calls += 1
        if calls == 2:
            raise OSError("injected installation failure")
        return original(source, target)

    monkeypatch.setattr(transaction, "replace_output", fail_second)
    tree = read_tree("(A:1,B:1)R;", "1", True, quiet=True)
    with pytest.raises(OSError, match="injected"):
        write_tree_with_tables(
            tree,
            make_args(outfile=str(tree_path)),
            format=1,
            tables=[(table_path, pd.DataFrame({"x": [1]}))],
        )
    assert tree_path.read_text() == table_path.read_text() == "previous\n"
    assert sorted(path.name for path in tmp_path.iterdir()) == ["table.tsv", "tree.nwk"]


def test_stream_failure_restores_table(tmp_path):
    class BrokenStream(io.StringIO):
        def write(self, text):
            raise OSError("broken stream")

    path = tmp_path / "table.tsv"
    path.write_text("previous\n")
    tree = read_tree("(A:1,B:1)R;", "1", True, quiet=True)
    with pytest.raises(OSError, match="broken stream"):
        write_tree_with_tables(
            tree,
            make_args(outfile=BrokenStream()),
            format=1,
            tables=[(path, pd.DataFrame({"x": [1]}))],
        )
    assert path.read_text() == "previous\n"


def test_successful_stream_and_table(tmp_path, capsys):
    tree = read_tree("(A:1,B:1)R;", "1", True, quiet=True)
    path = tmp_path / "table.tsv"
    write_tree_with_tables(
        tree,
        make_args(outfile="-"),
        format=1,
        tables=[(path, pd.DataFrame({"x": [1]}))],
    )
    assert capsys.readouterr().out == "(A:1,B:1)R;\n"
    assert path.read_text() == "x\n1\n"


def test_standalone_tree_preserves_existing_output_on_partial_write(
    tmp_path, monkeypatch
):
    import os

    from nwkit import output_transaction as transaction
    from nwkit.util import write_tree

    path = tmp_path / "tree.nwk"
    path.write_text("previous output\n")
    original_fdopen = os.fdopen

    class FailingWriter:
        def __init__(self, handle):
            self.handle = handle

        def __enter__(self):
            return self

        def __exit__(self, *args):
            self.handle.close()

        def write(self, text):
            self.handle.write(text[:5])
            self.handle.flush()
            raise OSError("injected partial write")

    def fdopen(descriptor, mode, **kwargs):
        handle = original_fdopen(descriptor, mode, **kwargs)
        return FailingWriter(handle) if mode == "w" else handle

    monkeypatch.setattr(transaction.os, "fdopen", fdopen)
    tree = read_tree("(A:1,B:1);", "1", True, quiet=True)
    with pytest.raises(OSError, match="injected partial write"):
        write_tree(tree, make_args(outfile=str(path)), format=1)
    assert path.read_text() == "previous output\n"
    assert list(tmp_path.iterdir()) == [path]


def test_tree_serializer_reuses_enclosing_transaction_stage(tmp_path):
    from nwkit.output_transaction import output_transaction
    from nwkit.util import write_tree

    path = tmp_path / "tree.nwk"
    tree = read_tree("(A:1,B:1);", "1", True, quiet=True)
    with output_transaction([path]) as staged:
        write_tree(tree, make_args(outfile=staged[path]), format=1)
        assert not path.exists()
    assert path.read_text() == "(A:1,B:1);"
    # Leaving the context must also restore normal standalone publication.
    write_tree(tree, make_args(outfile=str(path)), format=1)
    assert path.read_text() == "(A:1,B:1);"

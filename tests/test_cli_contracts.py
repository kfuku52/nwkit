"""Exercise every command's real parser and handler with small offline inputs."""

import argparse
import json

import pandas as pd
import pytest

from nwkit.cli import main, parser
from nwkit.util import read_tree

pytestmark = pytest.mark.integration

CASES = {
    "annotate": ["--table", "{data}"],
    "asr": ["--trait", "{data}", "--state-column", "state", "--rate", "0.2"],
    "asrcompare": [
        "--trait",
        "{data}",
        "--state-column",
        "state",
        "--models",
        "ER,ARD",
        "--figure-out",
        "{pdf}",
    ],
    "constrain": ["--backbone", "user", "--species-list", "{species}"],
    "collapse": ["--max-dist", "0.5"],
    "compose": ["--length-source", "{other}"],
    "cladefreq": [],
    "consensus": [],
    "contrast": ["--trait", "{data}", "--columns", "x"],
    "diff": ["-i2", "{other}"],
    "dist": ["-i2", "{other}"],
    "draw": ["--species-overlap-node-plot", "no"],
    "drop": [],
    "info": [],
    "image": [
        "--source",
        "wikimedia",
        "--out-dir",
        "{images}",
        "--download-dir",
        "{downloads}",
    ],
    "intersection": ["-i2", "{other}"],
    "label": [],
    "rename": ["--pattern", "Genus_a", "--replacement", "Renamed_a"],
    "reconcile": ["--species-tree", "{tree}"],
    "regress": [
        "--tree",
        "{tree}",
        "--data",
        "{data}",
        "--responses",
        "y",
        "--predictors",
        "x",
        "--evolution-model",
        "independent",
    ],
    "mark": ["--insert-text", "marked"],
    "mcmctree": [
        "--left-species",
        "Genus_a",
        "--right-species",
        "Genus_b",
        "--lower-bound",
        "1",
    ],
    "monophyly": [],
    "nhx2nwk": [],
    "nwk2table": [],
    "printlabel": [],
    "prune": ["--pattern", "Genus_a"],
    "rescale": ["--factor", "2"],
    "root": ["--method", "midpoint"],
    "rootcompare": ["--methods", "midpoint", "--figure-out", "{pdf}"],
    "sanitize": [],
    "shuffle": ["--seed", "1"],
    "skim": ["--seed", "1"],
    "sample": [],
    "subtree": ["--leaves", "Genus_a,Genus_b"],
    "transfer": ["-i2", "{other}"],
    "validate": [],
    "table2nwk": [],
    "help": ["asr"],
}

TABLE_COLUMNS = {
    "asr": "map_state",
    "asrcompare": "model",
    "cladefreq": "frequency",
    "contrast": "standardized_contrast",
    "diff": "status",
    "dist": "metric",
    "reconcile": "event_type",
    "regress": "coefficient",
    "monophyly": "status",
    "nwk2table": "branch_id",
    "rootcompare": "method",
    "validate": "tree_id",
}


def test_every_registered_command_has_a_handler_contract():
    subcommands = next(
        action
        for action in parser._actions
        if isinstance(action, argparse._SubParsersAction)
    )
    assert set(subcommands.choices) == set(CASES)


@pytest.mark.parametrize("command", CASES)
def test_command_parser_reaches_its_real_handler(
    command, tmp_path, monkeypatch, capsys
):
    names = [f"Genus_{letter}" for letter in "abcde"]
    tree = tmp_path / "tree.nwk"
    tree.write_text("(((Genus_a:1,Genus_b:1):1,Genus_c:2):1,(Genus_d:1,Genus_e:1):2);")
    other = tmp_path / "other.nwk"
    other.write_text(tree.read_text())
    data = tmp_path / "data.tsv"
    pd.DataFrame(
        {
            "leaf_name": names,
            "x": [1, 2, 3, 4, 5],
            "y": [2, 5, 5, 9, 8],
            "state": ["x", "x", "y", "y", "y"],
        }
    ).to_csv(data, sep="\t", index=False)
    species = tmp_path / "species.txt"
    species.write_text("\n".join(names) + "\n")
    paths = {
        "tree": tree,
        "other": other,
        "data": data,
        "species": species,
        "images": tmp_path / "images",
        "downloads": tmp_path / "downloads",
        "pdf": tmp_path / "roots.pdf",
    }
    arguments = [command, *(argument.format(**paths) for argument in CASES[command])]
    if command == "help":
        with pytest.raises(SystemExit) as exc:
            main(arguments)
        assert exc.value.code == 0
        assert "--root-prior" in capsys.readouterr().out
        return
    if command == "image":
        # The only external boundary in these cases: a provider with no hits.
        monkeypatch.setattr(
            "nwkit.image.provider_json_request", lambda *args, **kwargs: {}
        )
        arguments += ["-i", str(tree)]
        main(arguments)
        manifests = list(paths["images"].glob("*.tsv"))
        assert manifests
        assert any("leaf_name" in path.read_text() for path in manifests)
        return
    output = tmp_path / ("drawing.svg" if command == "draw" else "result.txt")
    if command == "table2nwk":
        table = tmp_path / "tree.tsv"
        main(["nwk2table", "-i", str(tree), "-o", str(table)])
        arguments += ["-i", str(table)]
    elif command != "regress":
        arguments += ["-i", str(tree)]
    arguments += ["-o", str(output)]
    main(arguments)
    assert output.is_file() and output.stat().st_size > 0
    if command in TABLE_COLUMNS:
        frame = pd.read_csv(output, sep="\t")
        assert TABLE_COLUMNS[command] in frame.columns
        if command in {"asrcompare", "rootcompare"}:
            assert paths["pdf"].read_bytes().startswith(b"%PDF")
    elif command == "draw":
        assert "<svg" in output.read_text()
    elif command not in {"info", "printlabel", "mcmctree"}:
        result = read_tree(str(output), "auto", True, quiet=True)
        assert list(result.leaves())


@pytest.mark.parametrize(
    "arguments", [["--version"], ["--help"], ["regress", "--help"]]
)
def test_help_and_version_do_not_import_numerical_libraries(arguments):
    import subprocess
    import sys

    code = """import json, sys
from nwkit.cli import main
try:
    main(sys.argv[1:])
except SystemExit as exc:
    assert exc.code in (None, 0)
print(json.dumps(sorted(set(sys.modules) & {'numpy', 'pandas', 'scipy', 'ete4', 'matplotlib'})))
"""
    result = subprocess.run(
        [sys.executable, "-c", code, *arguments],
        check=True,
        text=True,
        capture_output=True,
    )
    assert json.loads(result.stdout.splitlines()[-1]) == []

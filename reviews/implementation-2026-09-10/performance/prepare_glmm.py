"""Generate two known difficult GLMM cases after prepare.py; argument is its root."""

import hashlib
import json
import sys
from pathlib import Path

import pandas as pd


def main():
    root = Path(sys.argv[1]).resolve()
    sys.path.insert(0, str(root / "reviewed" / "tools"))
    sys.path.insert(0, str(root / "reviewed"))
    from regression_calibration_design import cases, generate, seed_for

    manifest = json.loads((root / "manifest.json").read_text())
    extra = []
    for name in ["binomial-n60-p0.05", "negative-binomial-n30-boundary"]:
        case = next(c for c in cases() if c.name == name)
        data = generate(case, seed_for(20260910, name, 1, "data"))
        folder = root / "inputs" / name
        folder.mkdir(exist_ok=True)
        (folder / "tree.nwk").write_text(data["tree"])
        names = [f"S{i}" for i in range(len(data["y"]))]
        pd.DataFrame(dict(leaf_name=names, x=data["x"][:, 0], y=data["y"])).to_csv(
            folder / "traits.tsv", sep="\t", index=False, float_format="%.17g"
        )
        frame = pd.DataFrame(data["covariance"], index=names, columns=names)
        frame.index.name = "leaf_name"
        frame.to_csv(folder / "covariance.tsv", sep="\t", float_format="%.17g")
        args = [
            "regress",
            "--tree",
            str(folder / "tree.nwk"),
            "--data",
            str(folder / "traits.tsv"),
            "--responses",
            "y",
            "--predictors",
            "x",
            "--response-family",
            "y=" + case.family,
            "--coefficient-penalty",
            "none",
            "--evolution-model",
            "custom",
            "--evolution-covariance",
            str(folder / "covariance.tsv"),
            "-o",
            "result.tsv",
        ]
        if case.family == "binomial":
            args += ["--categorical-responses", "y"]
        extra.append(
            dict(
                name=name,
                before=str(root / "before"),
                after=str(root / "after"),
                args=args,
                tables=["result.tsv"],
            )
        )
    manifest["cases"] = extra
    manifest["input_sha256"] = {
        str(p.relative_to(root / "inputs")): hashlib.sha256(p.read_bytes()).hexdigest()
        for p in (root / "inputs").rglob("*")
        if p.is_file()
    }
    (root / "glmm-manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")


if __name__ == "__main__":
    main()

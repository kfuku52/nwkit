"""Run the frozen paired OU simulation protocol against a corrected backend."""

import argparse
import csv
import hashlib
import json
import os
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
import shutil
import subprocess

from shift_alpha_design import cases, generate_case, protocol
from shift_joint_candidates import enumerate_candidates

from nwkit.util import assign_branch_ids, read_tree


def write_table(path, rows, fields=None):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields or list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def prepare_case(root, case, specification):
    folder = root / case["case_id"]
    folder.mkdir()
    newick, truth = generate_case(case)
    (folder / "tree.nwk").write_text(newick + "\n")
    (folder / "truth.json").write_text(json.dumps(truth, indent=2) + "\n")
    (folder / "case.json").write_text(json.dumps(case, indent=2) + "\n")
    write_table(folder / "traits.tsv", [
        {"leaf_name": name, "value": value, "se": case["standard_error"]}
        for name, value in zip(truth["tip_names"], truth["observations"], strict=True)
    ])
    tree = read_tree(newick, "auto", True, quiet=True)
    ids = assign_branch_ids(tree)
    write_table(folder / "branches.tsv", [
        {"branch_id": branch, "clade": "/".join(sorted(node.leaf_names()))}
        for node, branch in ids.items() if branch
    ])
    candidates, _ = enumerate_candidates(tree)
    write_table(folder / "candidates.tsv", [
        {"candidate_id": row["candidate_id"],
         "shifts": ";".join(map(str, row["shift_branch_ids"])),
         "groups": "|".join(";".join(map(str, group)) for group in row["groups"])}
        for row in candidates
    ])
    height = max(tree.get_distance(tree, node) for node in tree.leaves())
    write_table(folder / "settings.tsv", [
        {"floor_id": floor["floor_id"], "root_model": case["root_model"],
         "lower": floor["alpha_height"] / height,
         "upper": specification["alpha_upper_height"] / height,
         "starting": specification["alpha_start_height"] / height}
        for floor in specification["floors"]
    ])
    return folder


def run_case(folder, rscript, backend, environment):
    with (folder / "backend.log").open("w") as log:
        result = subprocess.run([rscript, str(backend), "run"], cwd=folder,
                                env=environment, stdout=log, stderr=subprocess.STDOUT,
                                check=False)
    return {"case_id": folder.name, "returncode": result.returncode}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--rscript", required=True)
    parser.add_argument("--replicates", type=int, default=50)
    parser.add_argument("--extension-replicates", type=int, default=10)
    parser.add_argument("--seed", type=int, default=20260911)
    args = parser.parse_args()
    root = args.output.resolve()
    root.mkdir(parents=True, exist_ok=False)
    specification = protocol(args.replicates, args.extension_replicates, args.seed)
    (root / "protocol.json").write_text(json.dumps(specification, indent=2) + "\n")
    source = Path(__file__).resolve().parent
    snapshot = root / "source"
    snapshot.mkdir()
    for name in ("shift_alpha_design.py", "shift_alpha_backend.R", "validate_shift_alpha.py",
                 "shift_simulation_cases.py", "shift_joint_candidates.py"):
        shutil.copy2(source / name, snapshot / name)
    backend = snapshot / "shift_alpha_backend.R"
    environment = dict(os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1",
                       MKL_NUM_THREADS="1", VECLIB_MAXIMUM_THREADS="1")
    subprocess.run([args.rscript, str(backend), "probe", str(root / "backend-probe.tsv")],
                   env=environment, check=True)
    library = Path(next(csv.DictReader((root / "backend-probe.tsv").open(), delimiter="\t"))["library"])
    hashes = {str(path.relative_to(library)): hashlib.sha256(path.read_bytes()).hexdigest()
              for path in sorted(library.rglob("*")) if path.is_file()}
    (root / "installed-backend-sha256.json").write_text(json.dumps(hashes, indent=2) + "\n")
    folders = [prepare_case(root, case, specification) for case in cases(specification)]
    print(f"Frozen {len(folders)} datasets before fitting", flush=True)
    with ThreadPoolExecutor(max_workers=specification["workers"]) as pool:
        futures = [pool.submit(run_case, folder, args.rscript, backend, environment) for folder in folders]
        with (root / "jobs.jsonl").open("w") as stream:
            for i, future in enumerate(as_completed(futures), 1):
                result = future.result()
                stream.write(json.dumps(result) + "\n")
                stream.flush()
                print(f"{i}/{len(folders)} {result}", flush=True)


if __name__ == "__main__":
    main()

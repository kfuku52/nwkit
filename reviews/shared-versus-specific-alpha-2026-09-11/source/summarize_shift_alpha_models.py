"""Summarize the paired native estimated-alpha comparison without hiding failures."""

import argparse
import hashlib
import json
import math
from pathlib import Path


def key(row):
    job = row["job"]
    return job["traits"], job["truth"], job["scenario"], job["replicate"]


def collect(path):
    return json.loads(path.read_text()) if path.exists() else []


def paired_rows(timing, pilot):
    pairs = {}
    for row in timing + pilot:
        identity = (row["job"]["kind"], *key(row))
        modes = pairs.setdefault(identity, {})
        if row["job"]["mode"] in modes:
            raise ValueError("Duplicate model result for a paired input")
        modes[row["job"]["mode"]] = row
    for modes in pairs.values():
        if len(modes) == 2 and all("generating" in r for r in modes.values()):
            if len({r["generating"]["data_sha256"] for r in modes.values()}) != 1:
                raise ValueError("Paired model inputs differ")
    return pairs


def check_nesting(pairs):
    nesting_checks = []
    for identity, modes in pairs.items():
        if (
            identity[0] != "pilot"
            or len(modes) != 2
            or any(r["status"] != "complete" for r in modes.values())
        ):
            continue
        by_mode = {
            mode: {
                tuple(r["shift_branch_ids"]): r["log_likelihood"]
                for r in row["records"]
                if r.get("log_likelihood") is not None
            }
            for mode, row in modes.items()
        }
        for layout in by_mode["shared"].keys() & by_mode["trait-specific"].keys():
            gain = by_mode["trait-specific"][layout] - by_mode["shared"][layout]
            nesting_checks.append(
                {
                    "case": list(identity),
                    "layout": list(layout),
                    "specific_minus_shared_ll": gain,
                    "violation": gain < -1e-5,
                }
            )
    return nesting_checks


def summarize_timing(pairs, summary, lines):
    for p in [2, 5, 10]:
        for truth in ["shared", "different"]:
            modes = pairs.get(("timing", p, truth, "opposed", 0), {})
            record = {"traits": p, "truth": truth}
            labels = []
            for mode in ["shared", "trait-specific"]:
                row = modes.get(mode, {"status": "missing"})
                record[mode] = {"status": row["status"]}
                if row["status"] == "complete":
                    record[mode].update(
                        seconds=row["seconds"],
                        median_seconds=row["median_seconds"],
                        alpha=row["fits"][-1]["alpha"],
                        log_likelihood=row["fits"][-1]["log_likelihood"],
                        BIC=row["fits"][-1]["criteria"]["BIC"]["score"],
                        engine=row["fits"][-1]["engine"],
                        optimizer=row["fits"][-1]["optimizer"],
                    )
                    ll = [r["log_likelihood"] for r in row["fits"]]
                    if max(ll) - min(ll) > 1e-7:
                        raise ValueError("Repeated fixed-layout fits disagree")
                    labels.append(f"{row['median_seconds']:.4f}")
                else:
                    record[mode].update(stage=row.get("stage"), error=row.get("error"))
                    labels.append(
                        ">120 (timeout)"
                        if row["status"] == "timeout"
                        else row["status"]
                    )
            ratio = gain = None
            if all(
                record[m]["status"] == "complete" for m in ["shared", "trait-specific"]
            ):
                ratio = (
                    record["trait-specific"]["median_seconds"]
                    / record["shared"]["median_seconds"]
                )
                gain = (
                    record["trait-specific"]["log_likelihood"]
                    - record["shared"]["log_likelihood"]
                )
                record["shared_submodel_nesting_violation"] = gain < -1e-5
            record.update(ratio=ratio, log_likelihood_gain=gain)
            if gain is not None:
                record["BIC_specific_minus_shared"] = (
                    record["trait-specific"]["BIC"] - record["shared"]["BIC"]
                )
            summary["timing"].append(record)
            ratio_text = "—" if ratio is None else f"{ratio:.1f}×"
            gain_text = "—" if gain is None else f"{gain:.6f}"
            lines.append(
                f"| {p} | {truth} | {labels[0]} | {labels[1]} | {ratio_text} | {gain_text} |"
            )


def summarize_pilot(pilot, summary, lines):
    for truth in ["shared", "different"]:
        for scenario in ["null", "aligned", "opposed"]:
            for criterion in ["AIC", "BIC"]:
                for mode in ["shared", "trait-specific"]:
                    rows = [
                        r
                        for r in pilot
                        if r["job"]["truth"] == truth
                        and r["job"]["scenario"] == scenario
                        and r["job"]["mode"] == mode
                    ]
                    completed = [r for r in rows if r["status"] == "complete"]
                    counts = {
                        name: sum(r["selected"][criterion][name] for r in completed)
                        for name in ["any_shift", "exact_branch", "false_branch"]
                    }
                    failures = sum(r["status"] == "failed" for r in rows)
                    timeouts = sum(r["status"] == "timeout" for r in rows)
                    rec = {
                        "truth": truth,
                        "scenario": scenario,
                        "criterion": criterion,
                        "mode": mode,
                        "complete": len(completed),
                        "planned": 10,
                        "failed": failures,
                        "timeout": timeouts,
                        **counts,
                    }
                    summary["pilot"].append(rec)
                    n = len(completed)
                    labels = [
                        f"{counts[name]}/{n}"
                        for name in ["any_shift", "exact_branch", "false_branch"]
                    ]
                    lines.append(
                        f"| {truth} | {scenario} | {criterion} | {mode} | {n}/10 | {' | '.join(labels)} | {failures}/{timeouts} |"
                    )


def paired_pilot_counts(pairs):
    records = []
    for truth in ["shared", "different"]:
        for scenario in ["null", "aligned", "opposed"]:
            completed = [
                modes
                for identity, modes in pairs.items()
                if identity[0] == "pilot"
                and identity[2:4] == (truth, scenario)
                and len(modes) == 2
                and all(r["status"] == "complete" for r in modes.values())
            ]
            for criterion in ["AIC", "BIC"]:
                result = {
                    "truth": truth,
                    "scenario": scenario,
                    "criterion": criterion,
                    "complete_pairs": len(completed),
                }
                for metric in ["any_shift", "exact_branch", "false_branch"]:
                    values = [
                        (
                            m["shared"]["selected"][criterion][metric],
                            m["trait-specific"]["selected"][criterion][metric],
                        )
                        for m in completed
                    ]
                    result[metric] = {
                        "shared": sum(a for a, b in values),
                        "trait-specific": sum(b for a, b in values),
                        "shared_only": sum(a and not b for a, b in values),
                        "specific_only": sum(b and not a for a, b in values),
                    }
                records.append(result)
    return records


def validate_finished(root, part, rows):
    manifest = json.loads((root / part / "manifest.json").read_text())
    expected = {json.dumps(job, sort_keys=True) for job in manifest["jobs"]}
    actual = {json.dumps(row["job"], sort_keys=True) for row in rows}
    if expected != actual or len(rows) != len(expected):
        raise ValueError(f"Incomplete or duplicate {part} outcomes")
    fingerprint = hashlib.sha256(
        json.dumps(manifest["implementation"], sort_keys=True, allow_nan=False).encode()
    ).hexdigest()
    for row in rows:
        if row["implementation_sha256"] != fingerprint:
            raise ValueError("Runtime/source fingerprint differs across jobs")
        if "generating" not in row:
            raise ValueError("Formal experiment contains a fixture failure")
        parameters = row["generating"]["parameters"]
        alpha, diffusion = parameters["alpha"], parameters["diffusion_covariance"]
        for i, rate in enumerate(alpha):
            variance = diffusion[i][i] * -math.expm1(-2 * rate) / (2 * rate)
            if abs(variance - 1) > 1e-12:
                raise ValueError("Generating marginal tip variance differs from one")
            for j in range(i):
                correlation = diffusion[i][j] / math.sqrt(
                    diffusion[i][i] * diffusion[j][j]
                )
                if abs(correlation - 0.8) > 1e-12:
                    raise ValueError(
                        "Generating diffusion correlation differs from protocol"
                    )


def summarize_errors(root, summary, lines, require_complete):
    path = root / "timing_errors/results.json"
    if not path.exists():
        return
    rows = collect(path)
    if require_complete:
        validate_finished(root, "timing_errors", rows)
    pairs = paired_rows(rows, [])
    records = []
    lines.extend(
        [
            "\n## Observation-error timing supplement\n",
            "Two traits, 100 tips, true shifted layout. Known sampling SE 0.1; generating extra measurement variance 0.04. Both models estimate the extra diagonal measurement variances as well as alpha/full process covariance. Serial warmup plus three fits, 120-second limit per fit. Both models use general vector pruning here.\n",
            "| Generating alpha | Shared, s | Trait-specific, s | Ratio | All evaluated modes succeeded |\n|---|---:|---:|---:|---|",
        ]
    )
    for truth in ["shared", "different"]:
        modes = pairs.get(("timing_errors", 2, truth, "opposed", 0), {})
        record = {"truth": truth}
        labels = []
        checked = []
        for mode in ["shared", "trait-specific"]:
            row = modes.get(mode, {"status": "missing"})
            record[mode] = {"status": row["status"]}
            if row["status"] == "complete":
                checked.append(row["all_evaluated_modes_succeeded"])
                record[mode].update(
                    median_seconds=row["median_seconds"],
                    seconds=row["seconds"],
                    fits=row["fits"],
                    all_evaluated_modes_succeeded=checked[-1],
                )
                if (
                    max(f["log_likelihood"] for f in row["fits"])
                    - min(f["log_likelihood"] for f in row["fits"])
                    > 1e-7
                ):
                    raise ValueError("Repeated noisy fits disagree")
                labels.append(f"{row['median_seconds']:.4f}")
            else:
                checked.append(False)
                labels.append(
                    ">120 (timeout)" if row["status"] == "timeout" else row["status"]
                )
                record[mode].update(stage=row.get("stage"), error=row.get("error"))
        ratio = None
        if all(checked):
            ratio = (
                record["trait-specific"]["median_seconds"]
                / record["shared"]["median_seconds"]
            )
            record["specific_minus_shared_ll"] = (
                record["trait-specific"]["fits"][-1]["log_likelihood"]
                - record["shared"]["fits"][-1]["log_likelihood"]
            )
        record["ratio"] = ratio
        records.append(record)
        ratio_text = "—" if ratio is None else f"{ratio:.2f}×"
        lines.append(
            f"| {truth} | {labels[0]} | {labels[1]} | {ratio_text} | {checked[0]} / {checked[1]} |"
        )
    summary["timing_errors"] = records
    lines.append(
        "\nThis supplement does not test detection power with observation errors. Run `python tools/benchmark_shift_alpha_errors.py --output results/timing_errors` only after parallel jobs have finished.\n"
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    parser.add_argument("--require-complete", action="store_true")
    args = parser.parse_args()
    root = args.directory
    timing = collect(root / "timing/results.json")
    pilot = collect(root / "pilot/results.json")
    if args.require_complete:
        validate_finished(root, "timing", timing)
        validate_finished(root, "pilot", pilot)
    pairs = paired_rows(timing, pilot)
    summary = {"timing": [], "pilot": [], "paired_inputs_equal": True}
    summary["matched_layout_nesting_checks"] = check_nesting(pairs)
    summary["paired_pilot"] = paired_pilot_counts(pairs)
    lines = [
        "# Estimated shared versus trait-specific alpha: results\n",
        "See [PROTOCOL.md](PROTOCOL.md) for the frozen model, budgets, seeds and limitations. Both fitted models estimate alpha and full evolutionary covariance. This is a Docker experiment on 100 tips, not a SIF validation.\n",
        "## Fixed-layout timing\n",
        "Serial warmup plus three repeated fits on the true shifted layout. Times exclude imports and simulation. A timeout is a failed attempt to obtain a complete timing series; it is not a measured median.\n",
        "| Traits | Generating alpha | Shared, s | Trait-specific, s | Ratio | LL gain, specific − shared |\n|---:|---|---:|---:|---:|---:|",
    ]
    summarize_timing(pairs, summary, lines)
    lines.extend(
        [
            "\n`shared` generating alpha means all ones; `different` means geometric spacing from 0.25 to 4. Likelihood gain is not a speed-correctness check between identical models: the models differ. A materially negative gain instead warns that numerical optimization failed to recover the nested shared submodel. Alpha estimates and optimizer diagnostics are retained in JSON.\n",
            "## Approximate-search pilot\n",
            "Ten paired datasets per cell; full covariance, four candidate branches and five refits. AIC/BIC selection is not a calibrated significance test. Reported detection fractions are conditional on completed runs; failures/timeouts remain visible. These small cells cannot establish general detection performance.\n",
            "| True alpha | Scenario | Criterion | Model | Complete / planned | Any shift | Exact branch | False branch | Failed / timeout |\n|---|---|---|---|---:|---:|---:|---:|---:|",
        ]
    )
    summarize_pilot(pilot, summary, lines)
    summary["complete_pilot_pairs"] = sum(
        len(modes) == 2 and all(r["status"] == "complete" for r in modes.values())
        for identity, modes in pairs.items()
        if identity[0] == "pilot"
    )
    lines.extend(
        [
            "\nCompleted paired searches: "
            + str(summary["complete_pilot_pairs"])
            + "/60. The per-job JSON records the actual candidate sets, nuisance estimates and numerical-mode diagnostics. Conditional success fractions must not be compared without the failure counts.\n",
            "## Reproduction\n",
            "Use the pinned local image recorded in the protocol, one BLAS/OpenMP thread, and run timing before parallel search jobs:\n",
            "```sh\npython tools/benchmark_shift_alpha_models.py --part timing --traits 2,5,10 --repeats 3 --timeout 120 --output results/timing\npython tools/benchmark_shift_alpha_models.py --part pilot --traits 2 --replicates 10 --workers 4 --timeout 180 --output results/pilot\npython tools/summarize_shift_alpha_models.py results\n```\n",
        ]
    )
    summarize_errors(root, summary, lines, args.require_complete)
    (root / "summary.json").write_text(
        json.dumps(summary, indent=2, allow_nan=False) + "\n"
    )
    (root / "README.md").write_text("\n".join(lines))


if __name__ == "__main__":
    main()

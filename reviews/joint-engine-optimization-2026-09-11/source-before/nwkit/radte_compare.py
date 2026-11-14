"""Compare three verified, saved single-family dating analyses without refitting."""

import hashlib
import io
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages

from nwkit.output_transaction import output_transaction, validate_output_targets
from nwkit.radte_compare_plot import comparison_figure as _comparison_figure
from nwkit.radte_compare_plot import components_figure as _components_figure
from nwkit.result_plot_data import load_dating_plot_data, radte_plot_protected_paths
from nwkit.util import validate_outputs_do_not_replace_inputs

MODES = ("fixed", "bounded", "ensemble")


def comparison_paths(prefix):
    if not isinstance(prefix, str) or not prefix.strip() or prefix == "-":
        raise ValueError("Comparison requires a filesystem output prefix.")
    return {
        "figure": prefix + ".pdf",
        "comparison": prefix + ".comparison.tsv",
        "components": prefix + ".uncertainty-components.tsv",
        "manifest": prefix + ".manifest.json",
    }


def _validate_runs(runs):
    reference = runs["fixed"]
    ignored = {
        "--species-node-bounds-tsv",
        "--species-node-intervals-tsv",
        "--species-tree-ensemble",
        "--gene-tree-ensemble",
    }
    evidence = {
        k: v for k, v in reference.manifest["input_sha256"].items() if k not in ignored
    }
    settings = (
        "max_age",
        "rate_sd",
        "rate_correlation",
        "substitution_model",
        "inference",
        "likelihood",
        "kappa",
        "gamma_shape",
        "gamma_categories",
        "gtr_exchangeabilities",
    )
    for mode, data in runs.items():
        if data.manifest.get("calibration_policy") != "hard-all-events":
            raise ValueError(
                "Comparison currently requires native hard-bound analyses."
            )
        if {
            k: v for k, v in data.manifest["input_sha256"].items() if k not in ignored
        } != evidence:
            raise ValueError(
                "Comparison runs must use the same reference input evidence."
            )
        if set(data.events) != set(reference.events):
            raise ValueError("Comparison runs must retain the reference gene clades.")
        if any(
            data.events[k]["event_type"] != reference.events[k]["event_type"]
            or data.events[k]["species_event_id"]
            != reference.events[k]["species_event_id"]
            for k in data.events
        ):
            raise ValueError("Comparison reconciliation events must agree.")
        if data.manifest["method"] != reference.manifest["method"] or any(
            data.manifest["options"].get(k) != reference.manifest["options"].get(k)
            for k in settings
        ):
            raise ValueError(
                "Comparison runs must use the same estimator and clock settings."
            )
        if data.manifest["interval_level"] != reference.manifest["interval_level"]:
            raise ValueError("Comparison interval levels must agree.")
        evidence_columns = (
            "input_interval_lower",
            "input_interval_upper",
            "input_interval_level",
            "input_interval_kind",
            "input_interval_source",
        )
        if any(
            data.species_rows[k].get(col) != reference.species_rows[k].get(col)
            for k in data.species_rows
            for col in evidence_columns
        ):
            raise ValueError(
                "Comparison runs must carry the same external species interval evidence."
            )
        if (
            mode != "ensemble"
            and data.manifest["options"].get("uncertainty") == "input-ensemble"
        ):
            raise ValueError(
                "Fixed and bounded modes must be conditional fits, not input ensembles."
            )
        variable = any(r["age_min"] != r["age_max"] for r in data.species_rows.values())
        if mode in {"fixed", "ensemble"} and variable:
            raise ValueError(
                "Fixed and ensemble reference runs must use fixed species ages."
            )
        if mode == "bounded" and not variable:
            raise ValueError("Bounded run requires at least one variable species age.")
    metadata = runs["ensemble"].manifest.get("input_ensemble")
    if not metadata or not runs["ensemble"].manifest["options"].get(
        "species_tree_ensemble"
    ):
        raise ValueError("Ensemble run requires joint species-chronogram samples.")
    if runs["ensemble"].manifest["options"].get("gene_tree_ensemble"):
        raise ValueError(
            "This comparison isolates species-age uncertainty; omit a gene-tree ensemble."
        )
    expected_method = {
        "marginal-lognormal": "branch-marginal",
        "sequence-marginal-quadratic": "sequence-marginal",
        "sequence-empirical-bayes-map": "sequence-joint-map",
    }[reference.manifest["method"]]
    if set(metadata.get("sample_methods", [])) != {expected_method}:
        raise ValueError(
            "Comparison ensemble estimators must match the reference estimator."
        )
    if len(set(metadata.get("sample_methods", []))) != 1:
        raise ValueError("Comparison cannot combine mixed ensemble estimators.")


def comparison_table(runs):
    _validate_runs(runs)
    records = []
    for mode, data in runs.items():
        for kind, rows in (
            ("species", data.species_rows),
            ("duplication", data.events),
        ):
            for key, row in rows.items():
                if kind == "duplication" and row["event_type"] != "duplication":
                    continue
                if kind == "species" and float(row["age"]) == 0:
                    continue
                records.append(
                    dict(
                        mode=mode,
                        kind=kind,
                        clade_id=key,
                        shared_age_id="S:" + key
                        if kind == "species"
                        else row["shared_age_id"],
                        node=row.get("node", "")
                        if kind == "species"
                        else row.get("gene_name", ""),
                        estimated_age=row["estimated_age"],
                        interval_lower=row.get("interval_lower"),
                        interval_upper=row.get("interval_upper"),
                        interval_status=row.get("interval_status", "unavailable"),
                        age_min=row["age_min"],
                        age_max=row["age_max"],
                        estimation_status=row.get(
                            "estimation_status", "represented-gene-event"
                        ),
                        input_age=float(row["age"]) if kind == "species" else np.nan,
                        input_interval_lower=row.get("input_interval_lower"),
                        input_interval_upper=row.get("input_interval_upper"),
                        input_interval_level=row.get("input_interval_level"),
                        input_interval_kind=row.get("input_interval_kind", ""),
                        input_interval_source=row.get("input_interval_source", ""),
                    )
                )
    return pd.DataFrame(records)


def _input_hashes(inputs):
    return {
        name: hashlib.sha256(Path(path).read_bytes()).hexdigest()
        for name, path in inputs
        if Path(path).is_file()
    }


def _write_comparison_pdf(path, runs, components, level):
    existing = set(plt.get_fignums())
    try:
        with (
            plt.rc_context({"font.family": "DejaVu Sans", "font.size": 9}),
            PdfPages(path) as pdf,
        ):
            fig = _comparison_figure(runs, level)
            pdf.savefig(fig)
            plt.close(fig)
            fig = _components_figure(components, runs["fixed"])
            pdf.savefig(fig)
            plt.close(fig)
    finally:
        for number in set(plt.get_fignums()) - existing:
            plt.close(number)


def compare_main(args):
    prefixes = {mode: getattr(args, mode + "_prefix") for mode in MODES}
    paths = comparison_paths(args.out_prefix)
    inputs = [
        (mode + ":" + key, path)
        for mode, prefix in prefixes.items()
        for key, path in radte_plot_protected_paths(prefix).items()
    ]
    inputs.append(("species-tree", args.species_tree))
    validate_outputs_do_not_replace_inputs(
        inputs, list(paths.items()), label="Comparison output"
    )
    validate_output_targets(paths.values(), follow_symlinks=False)
    source_hashes = _input_hashes(inputs)
    runs = {
        mode: load_dating_plot_data(
            prefix, args.species_tree, species_rooted=args.species_tree_rooted
        )
        for mode, prefix in prefixes.items()
    }
    table = comparison_table(runs)
    components_path = prefixes["ensemble"] + ".uncertainty-components.tsv"
    contents = Path(components_path).read_bytes()
    if hashlib.sha256(contents).hexdigest() != runs["ensemble"].manifest[
        "output_sha256"
    ].get("uncertainty_components"):
        raise ValueError(
            "Uncertainty component table does not match the ensemble manifest."
        )
    components = pd.read_csv(io.BytesIO(contents), sep="\t")
    if _input_hashes(inputs) != source_hashes:
        raise ValueError(
            "Comparison inputs changed during reading; retry after source runs finish."
        )
    Path(args.out_prefix).absolute().parent.mkdir(parents=True, exist_ok=True)
    with output_transaction(list(paths.values()), follow_symlinks=False) as staged:
        staged.write_text(
            paths["comparison"],
            lambda f: table.to_csv(f, sep="\t", index=False, na_rep="NA"),
        )
        staged.write_text(
            paths["components"],
            lambda f: components.to_csv(f, sep="\t", index=False, na_rep="NA"),
        )
        _write_comparison_pdf(
            staged[paths["figure"]],
            runs,
            components,
            runs["fixed"].manifest["interval_level"],
        )
        manifest = dict(
            schema="nwkit-radte-comparison-v1",
            status="complete",
            input_sha256=source_hashes,
            prefixes=prefixes,
            output_sha256={
                k: hashlib.sha256(Path(staged[p]).read_bytes()).hexdigest()
                for k, p in paths.items()
                if k != "manifest"
            },
        )
        staged.write_text(
            paths["manifest"],
            lambda f: json.dump(manifest, f, indent=2, allow_nan=False),
        )

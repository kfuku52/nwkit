"""Concise scientific captions; complete metadata stays in the source document."""


def model_label(rec):
    model = rec["model"]
    c = model.get("conditions", {})
    c = c if isinstance(c, dict) else {}
    rates = ", ".join(
        f"{name}={c[name]:.5g}"
        for name in ("duplication", "transfer", "loss", "eta")
        if isinstance(c.get(name), (int, float))
    )
    rho = c.get("detection", {})
    rho = (
        rho
        if isinstance(rho, dict)
        and all(
            isinstance(value, (int, float)) and not isinstance(value, bool)
            for value in rho.values()
        )
        else {}
    )
    observation = (
        f"; detection={min(rho.values()):.4g}..{max(rho.values()):.4g} ({len(rho)} species)"
        if rho
        else ""
    )
    conditions = "; ".join(
        f"{key}={c[key]}" for key in ("origin", "ascertainment", "time") if key in c
    )
    return f"{rec['id']}: {rec['history_kind']}; selection={rec['selection'].get('method', 'unreported')}; model={model.get('id', model.get('source_format', 'unreported'))}; {rates}; {conditions}{observation}"


def source_label(record):
    source = record["provenance"]
    if source.get("description"):
        return "Source: " + str(source["description"])
    label = "Source: " + str(source.get("software", source.get("format", "unreported")))
    if "refitted" in source:
        label += f"; refitted={source['refitted']}"
    if source.get("saved_run_sha256"):
        label += (
            "; saved run SHA-256="
            + source["saved_run_sha256"][:16]
            + "... (full hash in source JSON)"
        )
    return label


def summary_label(rec):
    summary = rec.get("summary", {})
    values = ", ".join(
        f"{key}={value:.6g}"
        for key, value in summary.items()
        if isinstance(value, (float, int))
    )
    return f"{rec['id']} source summary: {values or 'retained in source JSON'} (summary values are not placed events)."


def runtime_label(panel):
    values = []
    for kind in ("inference", "history"):
        runtime = panel[kind + "_runtime"]
        seconds = runtime["seconds"]
        duration = f"{seconds:.5g}s" if seconds is not None else "unavailable"
        date = runtime["measured_at"] or "unreported date"
        values.append(f"{kind} runtime={duration}, {date}")
    return (
        panel["label"]
        + ": "
        + "; ".join(values)
        + ". Timing conditions and execution hashes are retained in the manifest."
    )

"""Run and extend independent upstream MCMCtree chains."""

from __future__ import annotations

import csv
import hashlib
import json
import math
import os
import re
import shutil
import signal
import struct
import subprocess
import sys
import tempfile
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path

from nwkit.file_paths import validate_outputs_do_not_replace_inputs
from nwkit.output_transaction import output_transaction
from nwkit.posterior_diagnostics import diagnose_chains

CHAINS = 4
RHAT_MAX = 1.01
ESS_MIN = 400
MAX_PAML_GENERATIONS = 2**31 - 1


def _digest(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024**2), b""):
            digest.update(block)
    return digest.hexdigest()


def _save_json(path: Path, data: dict) -> None:
    pending = path.with_suffix(".pending")
    pending.write_text(
        json.dumps(data, indent=2, allow_nan=False) + "\n", encoding="utf-8"
    )
    pending.replace(path)


@dataclass(frozen=True)
class Sampling:
    burnin: int
    sampfreq: int
    minimum_draws: int
    maximum_draws: int

    def generations(self, draws: int) -> int:
        return self.burnin + draws * self.sampfreq


def paml_control_values(text, hessian=None):
    """Translate an IQ2MC control template into upstream MCMCtree options."""
    values = {}
    for line in text.replace("FossilErrprint", "FossilErr\nprint").splitlines():
        line = re.split(r"[*#]", line, maxsplit=1)[0].strip()
        if not line:
            continue
        match = re.fullmatch(r"([A-Za-z][A-Za-z0-9_]*)\s*=\s*(.+)", line)
        if not match or match[1] in values:
            raise ValueError(f"Malformed or duplicate control setting: {line}")
        values[match[1]] = match[2].strip()

    iq2mc_hessian = values.pop("hessianfile", None)
    usedata = values.get("usedata", "").split()
    if usedata[:1] == ["2"]:
        input_path = usedata[1] if len(usedata) > 1 else iq2mc_hessian or "in.BV"
        values["usedata"] = " ".join(["2", hessian or input_path, *usedata[2:]])

    iq2mc_checkpoint = values.pop("ckpfile", None)
    checkpoint = values.get("checkpoint", "0").split()
    if iq2mc_checkpoint and len(checkpoint) < 3:
        probability = checkpoint[1] if len(checkpoint) == 2 else "0.1"
        values["checkpoint"] = f"{checkpoint[0]} {probability} {iq2mc_checkpoint}"

    if "BDparas" in values and values["BDparas"].split()[-1].upper() not in {"C", "M"}:
        # The IQ2MC MCMCtree fork used the conditional birth-death prior.
        values["BDparas"] += " C"
    return values


def _load_inputs(control: Path, max_generations: int | None):
    raw_control = control.read_bytes()
    values = paml_control_values(raw_control.decode("utf-8-sig"))
    try:
        burnin, frequency, minimum = (
            int(values[key]) for key in ("burnin", "sampfreq", "nsample")
        )
    except (ValueError, KeyError):
        raise ValueError(
            "Specify integer burnin, sampfreq and nsample in the control file."
        ) from None
    if burnin < 0 or frequency < 1 or minimum < 8:
        raise ValueError(
            "Use nonnegative burnin, positive sampfreq and at least eight retained samples."
        )
    if max_generations is None:
        max_generations = min(
            burnin + frequency * min(64 * minimum, MAX_PAML_GENERATIONS // CHAINS),
            MAX_PAML_GENERATIONS,
        )
    if not burnin + frequency * minimum <= max_generations <= MAX_PAML_GENERATIONS:
        raise ValueError(
            "The generation limit must include the initial run and fit in PAML's integer range."
        )
    maximum = (max_generations - burnin) // frequency
    if maximum > MAX_PAML_GENERATIONS // CHAINS:
        raise ValueError(
            "The pooled MCMCtree summary would exceed PAML's integer range."
        )
    if maximum > minimum and frequency * minimum <= 100:
        raise ValueError(
            "Upstream PAML needs more than 100 sampling generations to write a checkpoint."
        )

    usedata = values.get("usedata", "").split()
    if usedata[:1] != ["2"] or len(usedata) > 3:
        raise ValueError(
            "Adaptive MCMCtree requires approximate likelihood (usedata = 2)."
        )
    if float(values.get("fossilerror", "0").split()[0]) != 0:
        raise ValueError(
            "PAML checkpoints do not preserve the fossil-error indicators."
        )
    if float(values.get("BayesFactorBeta", "1")) != 1:
        raise ValueError(
            "Adaptive MCMCtree requires the ordinary posterior (BayesFactorBeta = 1)."
        )
    inputs = {
        "alignment": (control.parent / values.get("seqfile", "")).resolve(),
        "tree": (control.parent / values.get("treefile", "")).resolve(),
        "hessian": (control.parent / usedata[1]).resolve(),
    }
    for kind, path in inputs.items():
        if not path.is_file() or not path.stat().st_size:
            raise ValueError(f"Missing or empty MCMCtree {kind}: {path}")
    values["usedata"] = " ".join(["2", "{hessian}", *usedata[2:]])
    return values, inputs, Sampling(burnin, frequency, minimum, maximum), raw_control


def _next_round(report: Path, draws: int, sampling: Sampling) -> dict:
    with report.open(encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        raise ValueError("No MCMCtree diagnostic parameters.")
    worst_ess = math.inf
    high_rhat = []
    failed = []
    for row in rows:
        rhat, bulk, tail = (float(row[key]) for key in ("rhat", "ess_bulk", "ess_tail"))
        if not all(math.isfinite(x) and x > 0 for x in (rhat, bulk, tail)):
            raise ValueError(
                f"Undefined convergence diagnostics for {row['parameter']}."
            )
        if rhat >= RHAT_MAX:
            high_rhat.append(row["parameter"])
        if rhat >= RHAT_MAX or min(bulk, tail) < ESS_MIN:
            failed.append(row["parameter"])
        worst_ess = min(worst_ess, bulk, tail)
    if not failed:
        return {
            "state": "passed",
            "reason": "convergence_criteria_met",
            "failed_parameters": [],
        }

    decision = {"state": "unconverged", "failed_parameters": failed}
    if draws == sampling.maximum_draws:
        return {**decision, "reason": "generation_limit"}
    if high_rhat:
        next_draws = 2 * draws
        estimate = None
        reason = "rhat_above_threshold"
    else:
        # Budget the next check conservatively when ESS is the only concern.
        projected = (1.2 * ESS_MIN / worst_ess) * draws
        estimate = (
            sampling.generations(math.ceil(projected))
            if math.isfinite(projected)
            else None
        )
        next_draws = math.ceil(min(projected, 4 * draws))
        reason = "ess_below_threshold"
    delta = max(next_draws - draws, (100 // sampling.sampfreq) + 1)
    target = min(draws + delta, sampling.maximum_draws)
    return {
        **decision,
        "reason": reason,
        "next_draws": target,
        "estimated_total_generations_per_chain": estimate,
    }


@contextmanager
def _interruptible(stop: threading.Event):
    def interrupt(_signum, _frame):
        stop.set()

    signals = (
        (signal.SIGINT, signal.SIGTERM)
        if threading.current_thread() is threading.main_thread()
        else ()
    )
    handlers = {sig: signal.signal(sig, interrupt) for sig in signals}
    try:
        yield
    finally:
        for sig, handler in handlers.items():
            signal.signal(sig, handler)


def _execute(binary: Path, directory: Path, stop: threading.Event) -> None:
    if stop.is_set():
        raise InterruptedError("MCMCtree interrupted.")
    with (
        (directory / "stdout.log").open("w") as out,
        (directory / "stderr.log").open("w") as err,
    ):
        with subprocess.Popen(
            [str(binary), "run.ctl"],
            cwd=directory,
            stdout=out,
            stderr=err,
            env={**os.environ, "OMP_NUM_THREADS": "1"},
        ) as child:
            while True:
                try:
                    code = child.wait(timeout=0.25)
                    break
                except subprocess.TimeoutExpired:
                    if stop.is_set():
                        child.terminate()
                        try:
                            child.wait(timeout=3)
                        except subprocess.TimeoutExpired:
                            child.kill()
                            child.wait()
                        raise InterruptedError("MCMCtree interrupted.") from None
    if code:
        raise ValueError(
            f"MCMCtree exited {code}. See {directory}/stdout.log and stderr.log."
        )


def _write_control(
    directory: Path, inputs: Path, values: dict[str, str], **overrides
) -> None:
    prefix = Path(os.path.relpath(inputs, directory)).as_posix()
    config = {
        **values,
        "seqfile": f"{prefix}/alignment",
        "treefile": f"{prefix}/tree",
        "usedata": values["usedata"].replace("{hessian}", f"{prefix}/hessian"),
        "outfile": "mcmctree.out",
        "mcmcfile": "mcmc.txt",
        **overrides,
    }
    (directory / "run.ctl").write_text(
        "".join(f"{key} = {value}\n" for key, value in config.items()), encoding="utf-8"
    )


def _validated_sample_digest(path: Path, frequency: int, draws: int) -> str:
    """Validate each new segment and fingerprint it in the same pass."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        raw_header = handle.readline()
        digest.update(raw_header)
        header = raw_header.split()
        if len(header) < 2 or header[0] != b"Gen" or len(set(header)) != len(header):
            raise ValueError(f"Invalid MCMCtree sample header in {path}.")
        count = 0
        for index, line in enumerate(handle):
            digest.update(line)
            row = line.split()
            if len(row) != len(header) or not all(
                math.isfinite(float(x)) for x in row[1:]
            ):
                raise ValueError(f"Invalid MCMCtree sample row in {path}.")
            if index == 0 and frequency > 1 and row[0] == b"1":
                continue
            count += 1
            if int(row[0]) != count * frequency:
                raise ValueError(f"Irregular MCMCtree generation grid in {path}.")
        if count != draws:
            raise ValueError(
                f"Incomplete MCMCtree samples in {path}: expected {draws}, read {count}."
            )
    return digest.hexdigest()


@dataclass(frozen=True)
class Segment:
    directory: Path
    draws: int
    checkpoint_hash: str | None


def _cached_segment(
    round_dir: Path, identity: dict, outputs: set[str]
) -> Segment | None:
    try:
        receipt = json.loads((round_dir / "completed.json").read_text(encoding="utf-8"))
        name = receipt["attempt"]
        if receipt["identity"] != identity or not re.fullmatch(r"attempt-\w+", name):
            return None
        attempt = round_dir / name
        hashes = receipt["sha256"]
        if hashes.keys() != outputs or any(
            _digest(attempt / key) != hashes[key] for key in outputs
        ):
            return None
        return Segment(attempt, identity["draws"], hashes.get("mcmctree.ckpt"))
    except (OSError, KeyError, TypeError, ValueError):
        return None


def _run_segment(
    chain: Path,
    inputs: Path,
    binary: Path,
    values: dict[str, str],
    sampling: Sampling,
    draws: int,
    seed: int,
    previous: Segment | None,
    stop: threading.Event,
) -> Segment:
    chain.mkdir(parents=True, exist_ok=True)
    generations = draws * sampling.sampfreq
    writes_checkpoint = generations > 100
    outputs = {"mcmc.txt"}
    if writes_checkpoint:
        outputs.add("mcmctree.ckpt")
    identity = {
        "draws": draws,
        "seed": seed,
        "previous_checkpoint": previous.checkpoint_hash if previous else None,
    }
    if (cached := _cached_segment(chain, identity, outputs)) is not None:
        return cached
    attempt = Path(tempfile.mkdtemp(prefix="attempt-", dir=chain))
    if previous:
        if previous.checkpoint_hash is None:
            raise ValueError(
                "Cannot extend a chain without a verified end-of-segment checkpoint."
            )
        shutil.copyfile(previous.directory / "mcmctree.ckpt", attempt / "mcmctree.ckpt")
    _write_control(
        attempt,
        inputs,
        values,
        seed=seed,
        burnin=0 if previous else sampling.burnin,
        sampfreq=sampling.sampfreq,
        nsample=draws,
        print=1,
        checkpoint=f"{2 if previous else 1} 1 mcmctree.ckpt",
    )
    _execute(binary, attempt, stop)
    hashes = {
        "mcmc.txt": _validated_sample_digest(
            attempt / "mcmc.txt", sampling.sampfreq, draws
        )
    }
    if writes_checkpoint:
        with (attempt / "mcmctree.ckpt").open("rb") as handle:
            (iteration,) = struct.unpack("=i", handle.read(4))
        if iteration + 1 != generations:
            raise ValueError(
                f"MCMCtree checkpoint is not at the final generation: {attempt}"
            )
        hashes["mcmctree.ckpt"] = _digest(attempt / "mcmctree.ckpt")
    _save_json(
        chain / "completed.json",
        {"attempt": attempt.name, "identity": identity, "sha256": hashes},
    )
    return Segment(attempt, draws, hashes.get("mcmctree.ckpt"))


def _append_segment(
    segment: Segment,
    destination: Path,
    frequency: int,
    offset: int,
    header: bytes | None,
) -> bytes:
    """Append only new samples to the invocation's temporary diagnostic trace."""
    with (segment.directory / "mcmc.txt").open("rb") as stream:
        names = b"\t".join(stream.readline().split()) + b"\n"
        if header is None:
            with destination.open("wb") as handle:
                handle.write(names)
                first = stream.readline()
                if frequency == 1 or first.split(None, 1)[0] != b"1":
                    handle.write(first)
                shutil.copyfileobj(stream, handle)
            return names
        if names != header:
            raise ValueError("MCMCtree segment sample headers differ.")
        with destination.open("ab") as handle:
            for line in stream:
                generation, row = line.split(None, 1)
                if generation != b"1" or frequency == 1:
                    handle.write(b"%d\t%s" % (offset + int(generation), row))
    return names


def _summarize(
    binary: Path, values, inputs: Path, read_paths, draws, seed, scratch, stop
):
    with (scratch / "mcmc.txt").open("wb") as output:
        header = None
        for path in read_paths:
            with path.open("rb") as handle:
                names = handle.readline()
                if header is None:
                    header = names
                    output.write(header)
                elif names != header:
                    raise ValueError(
                        "MCMCtree independent chain sample headers differ."
                    )
                # PAML print=-1 skips Gen when summarizing. Each input trace already
                # has the validated number of draws and a complete sample grid.
                shutil.copyfileobj(handle, output)
    count = CHAINS * draws
    _write_control(
        scratch,
        inputs,
        values,
        seed=seed,
        sampfreq=1,
        nsample=count,
        print=-1,
        checkpoint=0,
    )
    _execute(binary, scratch, stop)
    from nwkit.convert import convert_tree_text

    return convert_tree_text(
        (scratch / "mcmctree.out").read_text(encoding="utf-8"), target="figtree"
    )


def _run_rounds(
    root, scratch, binary, values, sampling, seed, jobs, destinations, stop
):
    inputs = root / "inputs"
    completed: list[Segment | None] = [None] * CHAINS
    headers: list[bytes | None] = [None] * CHAINS
    paths = [scratch / f"chain-{index + 1:02}.txt" for index in range(CHAINS)]
    history: list[dict] = []
    target = sampling.minimum_draws
    sampled = 0
    number = 0
    started = time.monotonic()
    while True:
        if stop.is_set():
            raise InterruptedError(
                "MCMCtree interrupted. Completed segments were retained."
            )
        number += 1
        round_dir = root / f"round-{number:03}"
        delta = target - sampled
        with ThreadPoolExecutor(max_workers=min(jobs, CHAINS)) as executor:
            futures = {
                executor.submit(
                    _run_segment,
                    round_dir / f"chain-{index + 1:02}",
                    inputs,
                    binary,
                    values,
                    sampling,
                    delta,
                    1
                    + (seed - 1 + (number - 1) * CHAINS + index)
                    % (MAX_PAML_GENERATIONS - 1),
                    previous,
                    stop,
                ): index
                for index, previous in enumerate(completed)
            }
            try:
                for future in as_completed(futures):
                    index = futures[future]
                    segment = future.result()
                    headers[index] = _append_segment(
                        segment,
                        paths[index],
                        sampling.sampfreq,
                        sampled * sampling.sampfreq,
                        headers[index],
                    )
                    completed[index] = segment
            except Exception:
                stop.set()
                raise
        sampled = target
        print(
            f"MCMCtree completed {sampling.generations(sampled):,} generations per chain.",
            file=sys.stderr,
        )
        diagnostics = round_dir / "diagnostics.tsv"
        try:
            with (round_dir / "diagnostics.log").open("w") as log:
                diagnose_chains(
                    paths,
                    diagnostics,
                    log=log,
                    jobs=jobs,
                    stopping_thresholds=(RHAT_MAX, ESS_MIN)
                    if sampled < sampling.maximum_draws
                    else None,
                )
            decision = _next_round(diagnostics, sampled, sampling)
        except (OSError, ValueError, subprocess.CalledProcessError) as error:
            _save_json(
                root / "progress.json",
                {
                    "state": "diagnostic_failed",
                    "reason": str(error),
                    "rounds": history,
                },
            )
            raise ValueError(
                f"MCMCtree diagnostics failed. See {round_dir}/diagnostics.log: {error}"
            ) from error
        if stop.is_set():
            raise InterruptedError(
                "MCMCtree interrupted. Completed segments were retained."
            )
        history.append(
            {
                **decision,
                "draws_per_chain": sampled,
                "generations_per_chain": sampling.generations(sampled),
                "elapsed_seconds": time.monotonic() - started,
            }
        )
        _save_json(root / "progress.json", {"rounds": history})
        if "next_draws" in decision:
            target = decision["next_draws"]
            print(
                f"MCMCtree will check again at {sampling.generations(target):,} "
                f"generations per chain ({decision['reason']}).",
                file=sys.stderr,
            )
            continue

        tree = _summarize(binary, values, inputs, paths, sampled, seed, scratch, stop)
        status = {
            **decision,
            "schema_version": 1,
            "completed_chains": CHAINS,
            "pooled_samples": sampled * CHAINS,
            "generations_per_chain": sampling.generations(sampled),
            "max_generations_per_chain": sampling.generations(sampling.maximum_draws),
            "rhat_max": RHAT_MAX,
            "ess_min": ESS_MIN,
            "run": str(root),
            "diagnostics": str(round_dir),
            "summary": str(root / "mcmctree.out"),
            "rounds": history,
            "tree_sha256": hashlib.sha256(tree.encode("utf-8")).hexdigest(),
        }
        retained_summary = root / "mcmctree.out"
        with output_transaction(
            [*destinations.values(), retained_summary], follow_symlinks=False
        ) as staging:
            Path(staging[destinations["tree"]]).write_bytes(tree.encode("utf-8"))
            shutil.copyfile(scratch / "mcmctree.out", staging[destinations["summary"]])
            shutil.copyfile(scratch / "mcmctree.out", staging[retained_summary])
            shutil.copyfile(diagnostics, staging[destinations["diagnostics"]])
            Path(staging[destinations["status"]]).write_text(
                json.dumps(status, indent=2, allow_nan=False) + "\n",
                encoding="utf-8",
            )
        return status


def run_mcmctree(
    control: Path,
    out_dir: Path,
    *,
    max_generations: int | None = None,
    run_store: Path | None = None,
    seed: int = 1729,
    jobs: int = CHAINS,
) -> dict:
    """Run a calibrated IQ2MC/PAML control and publish a checked output bundle."""
    if (
        not isinstance(jobs, int)
        or isinstance(jobs, bool)
        or not 1 <= jobs <= MAX_PAML_GENERATIONS
    ):
        raise ValueError("Use a positive integer worker budget below 2147483648.")
    if not 1 <= seed < MAX_PAML_GENERATIONS:
        raise ValueError("Use a positive base seed below 2147483647.")
    executable = shutil.which("mcmctree")
    if not executable:
        raise ValueError("Upstream PAML's mcmctree executable is required on PATH.")
    binary = Path(executable).resolve()
    control = control.resolve()
    values, inputs, sampling, raw_control = _load_inputs(control, max_generations)
    out_dir = out_dir.resolve()
    destinations = {
        "tree": out_dir / "FigTree.tre",
        "summary": out_dir / "mcmctree.out",
        "status": out_dir / "convergence.json",
        "diagnostics": out_dir / "diagnostics.tsv",
    }
    validate_outputs_do_not_replace_inputs(
        [("control", control), *inputs.items()], destinations.items()
    )
    input_hashes = {label: _digest(path) for label, path in inputs.items()}
    contract = {
        "inputs": {"control": hashlib.sha256(raw_control).hexdigest(), **input_hashes},
        "binary_sha256": _digest(binary),
        "runner_sha256": _digest(Path(__file__)),
        "diagnostics_sha256": _digest(
            Path(__file__).with_name("data_r") / "mcmctree_diagnostics.R"
        ),
        "max_generations": sampling.generations(sampling.maximum_draws),
        "seed": seed,
    }
    run_id = hashlib.sha256(
        json.dumps(contract, sort_keys=True).encode("utf-8")
    ).hexdigest()
    root = (run_store.resolve() if run_store is not None else out_dir / "runs") / run_id
    out_dir.mkdir(parents=True, exist_ok=True)
    root.mkdir(parents=True, exist_ok=True)
    from nwkit.util import acquire_exclusive_lock

    stop = threading.Event()
    with (
        _interruptible(stop),
        acquire_exclusive_lock(root / ".lock", lock_label="MCMCtree run"),
    ):
        saved_inputs = root / "inputs"
        saved_inputs.mkdir(exist_ok=True)
        for label, source in inputs.items():
            destination = saved_inputs / label
            if not destination.is_file() or _digest(destination) != input_hashes[label]:
                pending = destination.with_suffix(".pending")
                shutil.copyfile(source, pending)
                if _digest(pending) != input_hashes[label]:
                    raise ValueError(f"MCMCtree input changed while copying: {source}")
                pending.replace(destination)
        (saved_inputs / "initial.ctl").write_bytes(raw_control)
        _save_json(root / "manifest.json", contract)
        with tempfile.TemporaryDirectory(prefix="posterior-", dir=root) as scratch:
            return _run_rounds(
                root,
                Path(scratch),
                binary,
                values,
                sampling,
                seed,
                jobs,
                destinations,
                stop,
            )


def mcmctree_run_main(args):
    status = run_mcmctree(
        args.control,
        args.out_dir,
        max_generations=args.max_generations,
        run_store=args.run_store,
        seed=args.seed,
        jobs=args.jobs,
    )
    if status["state"] != "passed":
        print(
            f"MCMCtree {status['state']}. The output dates are provisional.",
            file=sys.stderr,
        )
    return 0

---
name: validate-change
description: Select and run NWKIT checks for a concrete code, CLI, or documentation change. Use for affected verification, not performance measurement or research calibration runs.
---

# Validate a NWKIT change

## Inputs

Use the working diff (including untracked intended files), the behavior being
changed, and the requested phase: local iteration or delivery. Run from the
repository root. Read [DEVELOPMENT.md](../../../DEVELOPMENT.md) for the current
environment preflight, check modes, and change-to-test table; do not maintain a
second command matrix here.

## Procedure

1. Identify the producer and consumers of the changed behavior. Read the
   affected command guide through [docs/README.md](../../../docs/README.md).
   For interface changes, compare the actual `python -m nwkit COMMAND --help`
   with the guide and shared CLI/TSV conventions.
2. Run the environment preflight. If the existing environment cannot start or
   import runtime packages, use the documented isolated setup. Preserve the
   existing environment and record the actual failure; never treat missing
   imports or an unavailable external runtime as a passing check.
3. Select tests using the development table and inspect the selected cases.
   Include consumers of shared helpers, numerical reference/invariance cases
   for scientific changes, and failure/rollback cases for output changes.
   Collect an uncertain selection with `python -m pytest --collect-only -q`
   followed by the selected paths/options before running it.
4. Run the selected `quick` or `test` command. For a CLI/docs starting example,
   use the small offline CLI check in DEVELOPMENT.md: it exercises real parsers
   and handlers with temporary outputs. Add the changed command's tests; it is
   not sufficient evidence for a numerical model change.
5. Inspect failures and skips (`-rs`), not just the process status. An empty
   selection is a failed validation attempt. Fix the cause within task scope;
   do not loosen tolerances, remove coverage, or switch models to pass. For
   missing native tools, report the unverified backend and continue independent
   checks. Network/build failures leave the corresponding delivery gate open.
6. For delivery, follow [RELEASING.md](../../../RELEASING.md) and the existing
   prepare-github-push skill. Preserve existing build artifacts before running
   modes that clear them. Do not turn focused results into a full-suite claim.

## Result and validation

Report the changed surface, why the selection covers it, the interpreter,
commands, pass/skip/failure results, and remaining unverified behavior. Distinguish
executed checks from static inspection. Keep exploratory outputs in temporary
directories; no permanent report or new test harness is required. Success means
the applicable commands pass and limitations are explicit, not merely that a
test process exits zero with optional integrations skipped.

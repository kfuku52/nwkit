<!-- BEGIN KF AGENT POLICY: source=https://github.com/kfuku52/kf-agent-policy; version=10; sha256=82e3c0eb467582a414d9a6b2feaaaf6f5c8ae330d30f2e3efbf8c303155d0e2e -->
# Common agent policy

Repository-specific instructions override these defaults.

- Follow the user's task scope within higher-priority instructions and execution
  permissions. Complete implementation through affected verification and a result
  report; a plan or investigation ends with its requested deliverable. Continue
  authorized work without repeated approval; identify actual blocking boundaries.
- Inspect the worktree and preserve unrelated changes. Refresh remote information
  when needed; do not merge, rebase, or switch branches merely to inspect it.
- Prefer the default branch when starting work without an established branch.
  Preserve an existing task branch; follow explicit user branch instructions.
  Never create or switch branches solely for a commit, push, release, or PR.
- Change or recommend branch protection only when explicitly asked. Honor explicit
  repository-specific direct-push exceptions; otherwise report a rejected push
  without bypassing protection or inventing a branch or PR.
- Unpublished implementation details may be redesigned; preserve existing public
  APIs, file formats, and saved-data compatibility unless a breaking change is
  authorized. Update affected producers, consumers, tests, examples, and docs.
- Fix verified root causes; do not hide failures with fallbacks or weaker checks.
  Document unavoidable workarounds and their removal conditions.
- Read relevant docs and run the repository's check entrypoint for the change and
  phase. Verify affected behavior; report checks run and omitted. Repeat or broaden
  successful checks only for new changes, failures, or unresolved concerns.
- For library metadata, require demonstrated incompatibility for exact pins or
  upper bounds; keep reproducibility locks separate.
- When editing READMEs, keep them concise with useful visuals inline; put extended
  guides in linked documentation.
- For GitHub push/release work, use `prepare-github-push` in `.agents/skills/`.
  Local-only commits need no version bump; GitHub pushes require one.
- For software performance work, use `benchmark-performance` in `.agents/skills/`.
  Performance claims require comparable measurements and equivalent output.
- For GitHub Actions edits, use `optimize-github-actions` in `.agents/skills/`.
  Preserve required coverage; never run untrusted PR code on self-hosted runners.
<!-- END KF AGENT POLICY -->

# Working in NWKIT

## Start here

- Read [DEVELOPMENT.md](DEVELOPMENT.md) for environment setup and the
  change-to-test table. [README.md](README.md) describes the command surface;
  [docs/README.md](docs/README.md) routes to command guides and research evidence.
- For CLI, tree I/O, or table changes, read
  [CLI/TSV conventions](docs/guides/CLI_TSV_CONVENTIONS.md) first.
- `nwkit/cli.py` defines parsing and dispatch; command modules expose
  `*_main` handlers. Shared tree/table I/O is in `nwkit/util.py`, conventions in
  `nwkit/conventions.py`, and staged writes in `nwkit/output_transaction.py`.
  Follow imports into the numerical model modules for scientific changes.
- Before a push, read [RELEASING.md](RELEASING.md) and use
  `.agents/skills/prepare-github-push/SKILL.md`. For selecting affected checks,
  use `.agents/skills/validate-change/SKILL.md`.

## Run and verify

Run from the repository root with an activated Python >=3.10 environment.
Use the setup/import preflight in DEVELOPMENT.md before trusting an existing
`.venv`; do not replace someone else's environment to repair it.

```sh
python -m nwkit --help
python tools/check.py quick -- tests/test_cli.py tests/test_cli_contracts.py tests/test_interface_conventions.py
```

The second command runs lint, formatting checks, incremental mypy, and small
offline CLI/interface cases using temporary test outputs. It is a starting
check, not sufficient validation for every change. Select affected tests from
DEVELOPMENT.md; `test` runs only pytest, `quick` adds static checks, and `release`
is the delivery gate. Report skipped optional backends separately from passes.

## Preserve scientific and output contracts

- Consult the affected guide before changing root assumptions, time/trait
  units, ML/REML likelihoods, missing-data handling, uncertainty interpretation,
  or random-seed behavior. ASR and native SHIFT have different tree requirements;
  do not apply one model's assumptions to another.
- Preserve documented CLI/Python interfaces, Newick/NHX properties, TSV columns,
  saved-model schemas and original-unit parameters. Include relevant consumer,
  round-trip and failure/rollback tests when changing shared serialization.
- Treat `examples/`, `reviews/`, and `docs/validation/` as evidence with their own
  inputs, seeds and protocols. Do not regenerate results or adjust scientific
  parameters/tolerances merely to make a check pass. Read the local study guide
  before intentionally revising an experiment.
- Keep bundled reference data in `nwkit/data_*` intact unless the task requires
  it. Avoid hand-editing build outputs, caches, `.venv`, `*.egg-info`, or local
  `output/` and `tmp/` contents. Put trial outputs in a fresh temporary directory;
  `dist`/`release` clear build directories, so preserve existing artifacts first.

At completion, review the diff and report changed behavior, exact checks and
results, skipped/unrun checks with reasons, and any remaining limitations.

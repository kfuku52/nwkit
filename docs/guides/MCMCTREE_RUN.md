# Adaptive MCMCtree runs

Use `nwkit mcmctree-tree` to prepare calibrations or convert an existing posterior
into a dated tree. `nwkit mcmctree` remains available as a compatibility alias.

Run IQ-TREE with `--dating mcmctree` on your calibrated topology. Keep the generated
control, rooted tree, dummy alignment, and Hessian together. NWKIT runs upstream
PAML's `mcmctree` on PATH and uses R's `posterior` package via `Rscript`.

```sh
nwkit mcmctree-run --control iq2mc.mcmctree.ctl --out-dir dating --jobs 4
```

The `--jobs` worker budget applies to both sampling and diagnostics. Sampling
runs up to four chains concurrently. R can use additional workers for
diagnostic columns on systems that support forked processes. On Windows, R
computes diagnostics in one process. The runner appends each new segment to
temporary diagnostic traces and retains the original segments for restart.
Intermediate checks compute R-hat and ESS. Their TSV leaves the MCSE columns
as `NA` until the criteria are met or the sampling budget is reached.

The control's `burnin`, `sampfreq`, and `nsample` define the initial run for each of
four independently seeded chains. For `20000`, `100`, and `20000`, that is
2,020,000 generations per chain. Calibrations, the Hessian, and clock and rate
priors stay fixed across rounds. Use the ordinary approximate-likelihood posterior
(`usedata = 2`) with fossil-error indicators disabled.

The first round can finish the run when every parameter has rank-normalized
R-hat below 1.01 and bulk and tail ESS of at least 400. If a check falls short
with all R-hats below 1.02 and both ESS measures at least 300, the runner fixes a
shorter interval of about one eighth of the current sampling length and
keeps it while the diagnostics stay in that range. Farther from the criteria,
high R-hat doubles the accumulated sampling generations. For low ESS alone,
the runner projects from the worst ESS with a 20% margin and a check no later
than four times the current sampling length. The interval respects PAML's
checkpoint requirements and the remaining budget.

Specify `--burnin-percentages 2.5,5,10,20,40` to try excluding part of each
trace's beginning. A full trace that passes has priority. Otherwise, the
runner checks the percentages in increasing order, using the same cutoff for all
chains. Each percentage applies to the total generations including the
original PAML burn-in. The cutoff rounds down onto the retained sample grid.
When a cutoff first passes, the runner adds draws and checks the same generation
boundary again. It pools trimmed draws only after confirmation. If the
cutoff remains unconfirmed at the limit, the full trace supplies the result
and its diagnostic status. Omit the option to use only the full trace.

Use `--max-generations` for an explicit per-chain limit including burn-in. The
default is the burn-in plus 64 times the initial sampling generations within
PAML's integer limits. The limit rounds down to a complete sampling interval.
When it is reached, `convergence.json` marks unmet criteria and records the
output as provisional. A failed chain or invalid diagnostic leaves prior
published files intact.

Each extension starts at the verified final checkpoint of its own chain with
a distinct reproducible seed, since PAML saves the sampler state without its
random number stream. The runner checks the saved iteration and sample grid,
excludes extra off-grid rows, and pools the regularly spaced draws at the end.
Rerun the same command to reuse complete segments after interruption. Input or
sampling-setting changes get separate run directories. Use `--run-store` to
keep restart evidence on a separate persistent path.

The output bundle contains `FigTree.tre`, `mcmctree.out`, `diagnostics.tsv`, and
`convergence.json`. Times are in the control file's units. The JSON records
each round, the final status, and the path to the independent samples and
checkpoint evidence. Its `sample_window` records the effective burn-in and the
first and last retained generation used for both the final diagnostics and
the dated tree. Python callers can use
`nwkit.mcmctree_runner.run_mcmctree` with the same settings.

[IQ2MC preparation](https://iqtree.github.io/doc/Dating) ·
[PAML checkpointing](https://github.com/abacus-gene/paml/wiki/MCMCtree) ·
[R-hat and ESS](https://doi.org/10.1214/20-BA1221)

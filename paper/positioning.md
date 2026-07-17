# Positioning for *Systematic Biology*

Last updated: 2026-07-17

## Article type

Software for Systematics and Evolution

## Central claim, version 0.2

> NWKIT is a composable command-line toolkit for reproducible curation,
> transformation, and analysis of phylogenetic trees.

The earlier wording, "a validated ... environment," is deliberately not used.
NWKIT has a substantial automated test suite and a `validate` command, but that
does not constitute independent or comprehensive validation of every algorithm.
The manuscript will use concrete formulations such as "628 automated tests
passed for version 0.29.0" and will reserve "accurate," "robust," or
"validated" for properties established by external-reference comparisons.

## Reader problem

Phylogenetic studies repeatedly transform trees between inference,
reconciliation, dating, comparative analysis, visualization, and reporting.
Many operations are individually simple but become fragile when implemented as
study-specific scripts with inconsistent assumptions about rooting, support
values, internal labels, taxon names, tree collections, or Newick dialects.

## Proposed contribution

NWKIT exposes 33 functional subcommands through consistent command-line input
and output conventions. Its candidate distinguishing combination is:

1. shell-composable Newick curation and transformation;
2. explicit detection and normalization of structural and format problems;
3. clade-level tree differencing and explicit multi-source composition;
4. topology-aware transfer across exact or partially overlapping tip sets;
5. tabular annotation plus taxonomy- and trait-aware constraints, sampling,
   and monophyly assessment;
6. machine-readable provenance records for every command;
7. tree-collection summaries and consensus support transfer; and
8. categorical ancestral-state reconstruction and stochastic mapping under Mk
   models from the same interface.

This is a hypothesis to test, not yet a novelty claim. The comparison must show
which capabilities are unique, which are shared, and which are better served by
other software.

## Primary comparison set

### Command-line toolkits

- Newick Utilities
- Gotree
- PhyKIT

### Programming libraries

- ETE Toolkit
- DendroPy
- Bio.Phylo
- TreeSwift
- ape, where an R comparison is scientifically informative

The manuscript must not imply that a command-line interface alone is novel.
Newick Utilities, Gotree, and PhyKIT already provide automatable or chainable
tree-processing commands.

## Candidate biological and analytical demonstrations

1. **Tree curation and cross-program handoff**: diagnose ambiguous or malformed
   Newick, normalize it, reconcile taxon sets, root a tree, and transfer node
   metadata without ad hoc scripts.
2. **Tree-collection synthesis**: calculate weighted clade frequencies, build a
   consensus, and transfer support to a fixed reference topology.
3. **Taxonomy- and trait-aware analysis**: create a representative subset or
   topological constraint, assess trait monophyly, and reconstruct categorical
   ancestral states.

At least one demonstration should use a published biological dataset in which
NWKIT already played a documented role. At least one should exercise major
capabilities added in the 0.22--0.28 development cycle.

## Main display-item plan

The journal permits four figures and tables in total. The working allocation is:

1. Figure 1: architecture, command families, and composable data flow;
2. Table 1: capability comparison with related tools;
3. Figure 2: two or three worked biological workflows and outputs; and
4. Figure 3: correctness checks, predefined input cases, and scaling.

All command-level detail and secondary benchmarks belong in supplementary
material.

## Claim boundaries

Do not currently claim that NWKIT is:

- formally or independently validated;
- faster or more memory efficient than related tools;
- the most comprehensive phylogenetic tree toolkit;
- widely used by the community; or
- a replacement for the programming APIs of ETE, DendroPy, Bio.Phylo, or ape.

These statements may be narrowed or revised after comparison, benchmarking,
and a systematic search for published uses.

## Evidence threshold for stronger language

- **Well-tested**: report test count, test scope, supported Python versions, and
  continuous-integration configuration.
- **Correct/accurate**: compare defined outputs with independent implementations
  or analytically known results on documented fixtures.
- **Robust**: predefine malformed and edge-case inputs and report behavior across
  tools without selectively excluding failures.
- **Scalable/fast**: benchmark equivalent operations, versions, hardware,
  wall-clock time, peak memory, and tree sizes using reproducible scripts.
- **Used in research**: cite peer-reviewed studies that explicitly name NWKIT;
  do not equate package downloads with unique users.

## Author decisions required later

- final author list, order, affiliations, and corresponding author;
- which published NWKIT use should anchor the biological example;
- datasets that can be redistributed through Dryad;
- funding, acknowledgments, and conflicts of interest;
- confirmation of software support for at least two years after publication; and
- final disclosure wording for any AI assistance retained in the submitted work.

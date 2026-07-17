# Detailed manuscript outline

Target length: approximately 3,600 words excluding references.

## Working titles

1. **NWKIT: Composable Command-Line Workflows for Phylogenetic Tree Curation and Analysis**
2. **NWKIT: A Command-Line Toolkit for Reproducible Phylogenetic Tree Processing**
3. **Reproducible Curation and Analysis of Phylogenetic Trees with NWKIT**

Working short title: **COMPOSABLE TREE PROCESSING WITH NWKIT**

## One-sentence argument

NWKIT brings routine tree curation, topology-aware transformations,
tree-collection summaries, and taxonomy- or trait-aware analyses into a
consistent command-line interface that can be tested, composed, and reproduced.

## Abstract: 180--220 words

1. Problem: downstream tree processing is recurrent but fragmented and often
   encoded in study-specific scripts.
2. Existing landscape: mature libraries and command-line toolkits exist, so the
   gap must be stated as a specific combination of operations rather than an
   absence of software.
3. Implementation: NWKIT 0.28.0, 33 functional commands, standard stream/file
   interfaces, multi-tree comparison and composition, machine-readable audit
   records, open-source distribution, and seeded stochastic operations.
4. Evidence: capability comparison, independent-output checks, biological
   workflows, and scaling/input-case results.
5. Conclusion: exact benefit supported by those results.

Do not use "validated," "robust," "comprehensive," or "high performance" in
the abstract until the relevant analyses support them.

## Introduction: approximately 550 words

### Paragraph 1: recurring downstream tree operations

- Trees move between inference, reconciliation, dating, comparative analyses,
  visualization, and publication.
- Newick is compact and widespread, but software-specific conventions for
  rooting, support, labels, annotations, and collections create friction.

### Paragraph 2: existing software and the remaining need

- Libraries: ape, Bio.Phylo, DendroPy, ETE, TreeSwift.
- CLI suites: Newick Utilities, Gotree, PhyKIT.
- Acknowledge strengths explicitly.
- Define the candidate gap: a consistent CLI combining defensive Newick
  handling, topology-aware transfer, taxonomy-aware operations, tree-set
  synthesis, sampling, and selected comparative analyses.

### Paragraph 3: NWKIT and paper objectives

- Development since 2019 and major functional expansion in version 0.28.0.
- Objectives: describe design; compare coverage; demonstrate biological
  workflows; evaluate correctness, failure handling, and scale.

## Materials and Methods: approximately 900 words

### Software design and implementation: 220 words

- Python 3.10+, ETE 4 tree representation, dependencies, MIT license.
- Version and commit used in the paper.
- Input/output conventions and command composition.
- Explicit parsing choices for internal labels and quoted names.

### Functional organization: 180 words

Group commands by user task rather than list all 33 in prose:

1. inspect, validate, sanitize, and convert;
2. rename, label, mark, prune, collapse, root, rescale, and transfer;
3. compare and summarize trees or tree collections;
4. use taxonomy or traits for constraints, sampling, and monophyly;
5. perform selected downstream analyses and prepare external-tool inputs; and
6. draw trees and retrieve licensed representative images.

Exact category membership belongs in Figure 1 or supplementary material.

### Software testing and release process: 130 words

- Report 619 tests for the manuscript commit, not a timeless claim.
- Separate unit, regression, CLI, exact-output, round-trip, and documentation
  example tests where possible.
- CI versions, lint, source distribution, wheel build, and installed-wheel smoke
  test.

### Capability comparison: 160 words

- Pin tool versions and evaluation date.
- Compare user-visible tasks, not internal APIs or loosely similar features.
- Distinguish native CLI, bundled script, programming API, and unavailable.
- Have at least one human author review every matrix entry.

### Correctness, input-case, and scaling evaluation: 120 words

- Analytically known small trees.
- Independent implementation comparisons where semantics match.
- Predefined malformed/ambiguous-input corpus.
- Synthetic balanced and imbalanced trees and tree collections over predefined
  sizes; repetitions, hardware, timing, and peak memory.

### Worked biological examples: 90 words

- Select two or three examples only after confirming data redistribution rights.
- Record every command in executable scripts.
- Avoid using a worked example as validation of an algorithm.

## Results: approximately 1,200 words

### A consistent interface spans distinct downstream tree tasks: 250 words

- Figure 1.
- Report functional groups and stream-compatible operations.
- Give one short pipeline, not a command catalogue.

### Capability comparison identifies overlap and distinguishing combinations: 300 words

- Table 1.
- Lead with shared capabilities and acknowledge tasks where alternatives are
  stronger.
- Highlight only audited differences.

### Reproducible workflows connect curation to biological analysis: 350 words

- Figure 2.
- Report input data, transformations, outputs, and biological interpretation for
  each example.
- Include one published-use example and one example of a major recent capability.

### Output checks and scaling define the reliable operating range: 300 words

- Figure 3.
- Report agreements and disagreements, failure behavior, and performance without
  turning the paper into a speed contest.
- State any size, format, platform, or numerical limitations.

## Discussion: approximately 650 words

### Main contribution: 220 words

- Return to the demonstrated combination of composability, explicit input
  handling, and phylogenetically informed operations.
- Explain when a CLI is preferable to a custom script and when it is not.

### Relationship to existing tools: 170 words

- NWKIT complements rather than replaces general libraries.
- Discuss overlaps with Newick Utilities, Gotree, and PhyKIT fairly.

### Limitations: 140 words

- Newick-centered rather than general multi-format data model.
- Some commands depend on external databases or services.
- Breadth creates a documentation and maintenance burden.
- Not every command represents a novel algorithm.

### Availability and future support: 120 words

- Versioned release, GitHub, Bioconda, Zenodo DOI, Dryad data.
- Two-year support commitment after author confirmation.

## Back matter

- Funding: author input is required before submission if funding applies.
- Acknowledgments: author review is required before submission.
- Data Availability: Dryad reviewer URL and Zenodo DOI placeholders until deposit.
- Conflict of Interest: author confirmation is required before submission.

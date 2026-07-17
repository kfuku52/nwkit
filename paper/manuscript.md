---
title: "NWKIT: Composable Command-Line Workflows for Phylogenetic Tree Curation and Analysis"
bibliography: references.bib
link-citations: true
---

**Kenji Fukushima**^1,2^ (ORCID:
<https://orcid.org/0000-0002-2353-9274>)

^1^ Center for Frontier Research, National Institute of Genetics, Mishima,
Japan\
^2^ Genetics Program, Graduate Institute for Advanced Studies, SOKENDAI,
Mishima, Japan

Correspondence: **Kenji Fukushima**, National Institute of Genetics, 1111 Yata,
Mishima, Shizuoka 411-8540, Japan;
telephone: +81-55-981-6751; email: **kenji.fukushima@nig.ac.jp**

Author running head: **FUKUSHIMA**

Title running head: **COMPOSABLE TREE PROCESSING WITH NWKIT**

## Abstract

Modern systematic research combines diverse phylogenetic methods to arrive at
biological interpretations, often generating multiple related trees along the
way. These trees may differ in topology, rooting, branch lengths, support
values, or annotations. Researchers must also inspect and transform trees,
connect tips to taxonomic or trait data, summarize collections, and perform
downstream analyses. Manual edits and study-specific scripts often obscure how
inferred trees become the objects ultimately interpreted. We developed
NWKIT, a command-line toolkit that organizes post-inference tree processing as
composable, pipe-compatible commands. NWKIT reports clade-level differences and
composes a target tree from roots, node labels, support values, branch lengths,
and arbitrary annotations supplied by related trees. It uses exact
descendant-taxon matches when tip sets agree and unique clade projections for
conservative transfer when sampling overlaps only partly. The same interface
supports Newick inspection and normalization, topology-aware transformation,
attachment and aggregation of taxonomic or trait data, synthesis of tree
collections, phylogenetic-diversity sampling, and categorical ancestral-state
reconstruction. Every command can append a machine-readable audit record of
arguments, input and output hashes, interpretation settings, warnings, and
execution status. These conventions make comparisons and transformations
inspectable, executable, and traceable, turning post-inference processing from
an implicit prelude into a reproducible analytical layer of systematic
research.

**Keywords:** ancestral-state reconstruction; command line; Newick;
phylogenetics; reproducibility; tree processing

Phylogenetic inference is rarely the final computational step in a systematic
study. Inferred trees pass to programs for reconciliation, dating, comparative
analysis, visualization, and reporting. Alternative analyses may yield
different topologies or roots, while bootstrap, dating, and annotation add
support, branch lengths, and metadata. Investigators must reconcile labels,
edit taxon sets, reroot topologies, transfer annotations, and summarize tree
collections. These operations determine which tree enters an analysis, yet
manual edits, short scripts, and software-specific conversions can leave the
final tree available without a reconstructable path to it.

The relevant evidence is also rarely confined to one tree. A
maximum-likelihood topology may be accompanied by bootstrap replicates, an
alternative root, a time-calibrated version, or trees annotated by different
programs. Choosing one topology need not imply discarding support, rooting,
branch-length, or annotation evidence carried by the others. Combining these
components, however, requires an explicit rule for deciding which nodes in
nonidentical trees represent the same biological clade. The problem becomes
harder when sampling differs because node positions provide no reliable
correspondence.

Newick is compact and broadly supported, but its conventions vary among
implementations. A numeric internal token may represent support or a node name
[@Czech2017]; rooting, quoted labels, and extended annotations also differ.
Thus, a syntactically readable tree can acquire a different biological
interpretation. Empty or duplicated tips, negative lengths, and mismatched
taxon sets may likewise remain unnoticed until a later step fails or returns an
unintended result.

Mature libraries include ape [@Paradis2004], DendroPy
[@Sukumaran2010; @Moreno2024], and ETE [@HuertaCepas2016]. Newick Utilities
provides automatable filters [@Junier2010], Gotree chainable commands and a Go
API [@Lemoine2021], and PhyKIT alignment, tree, and comparative functions
[@Steenwyk2021]. The gap is therefore not generic tree software, but that
defensive Newick handling, topology-aware transfer, tree-set synthesis,
taxonomy-aware operations, and selected downstream analyses remain distributed
across interfaces with different assumptions.

NWKIT began in 2019 as reusable shell commands for tree transformation. Version
0.28.0 adds preflight reporting, tree differencing, multi-source composition,
provenance, tree-set summaries, taxonomy- and trait-aware operations,
visualization, and categorical ancestral-state reconstruction. We describe its
design, compare task-level capabilities with three command-line toolkits,
demonstrate a biological workflow, check selected outputs independently, test a
problematic-input corpus, and measure scaling. These analyses evaluate NWKIT as
a reproducible layer between inference and interpretation rather than by
command count.

## Materials and Methods

### Software Design and Implementation

Analyses used NWKIT 0.28.0 at commit `9d7353df9c43`. NWKIT requires Python 3.10
or later and uses ETE 4.4.0, Biopython, NumPy, SciPy, pandas, and Matplotlib. It
is MIT-licensed, distributed through Bioconda [@Gruning2018], and documented in
the project wiki.

Individual commands remain small, while tree I/O, label interpretation, clade
mapping, and audit capture are centralized. Commands therefore share stream and
format behavior, whereas biological choices remain explicit at the subcommand
level.

The interface contains 33 functional commands. All read their primary input
from standard input or a path; 31 write their primary result to standard output,
while `draw` and `image` create files. Tree-producing commands can therefore be
piped. Stochastic commands expose an integer seed. Every command accepts an
`--audit` path and appends a JSON Lines record with version, arguments, input
semantics, seeds, external-data settings, warnings, status, and input/output
SHA-256 hashes.

Input handling separates syntax from interpretation. Users may select an ETE
format, request automatic detection, or reject ambiguous unquoted numeric
internal labels. Quoted-name control and shared species-label parsers are also
explicit because guessed label or support conventions can alter an analysis.

### Functional Organization

Commands were grouped into five task families (Fig. 1; Table S9). Comparison
commands report clade differences and rooted Robinson--Foulds distance
[@RobinsonFoulds1981]; tree-set commands calculate clade frequencies and
consensus. Taxonomy and trait commands attach and aggregate data, generate
constraints, diagnose monophyly, and sample phylogenetic diversity [@Faith1992].
Rooting includes minimal ancestor deviation [@Tria2017] and taxonomic reference
trees [@Schoch2020; @Hinchliff2015; @Kumar2022; @APG2016]. Downstream commands
implement Mk ancestral-state reconstruction [@Lewis2001], prepare MCMCtree
inputs, and create visual outputs.

`diff` reports rooted clades or root-independent splits. `compose` takes a
target topology and separate sources for root, names, support, branch lengths,
or NHX properties; `transfer` maps one source. Matching uses descendant-taxon
sets, or informative projections unique in both trees when tip sets overlap
partly. Ambiguous requests remain unchanged or fail under a strict policy;
per-clade TSV reports retain values and reasons.

`annotate` joins table rows to unique tip names and can propagate a column to
internal nodes by uniqueness, mode, count, mean, sum, extrema, or list
aggregation. Explicit policies govern missing tips and unmatched rows, which
can also be reported. Taxonomic or trait data can thus enter the tree stream
without a study-specific join script.

### Software Tests and Release Checks

With ETE 4.4.0 and Biopython 1.87, 612 tests passed. They covered outputs,
errors, regressions, round trips, seeds, threading, partial-taxon ambiguity,
composition, aggregation, provenance hashes, and documentation. Continuous
integration spans Python 3.10--3.14 and checks distributions. An inventory
script extracted interface defaults; a smoke test piped `sanitize`, `rename`,
and `rescale` while requiring three ordered audit records with stdin hashes.

### Task-Level Software Comparison

We compared NWKIT 0.28.0 with Gotree 0.5.2, PhyKIT 2.3.0, and Newick Utilities
1.6 as available on 16 July 2026. Rows represent tasks, not similarly named
functions. “Native” required a documented CLI workflow, “Partial” a narrower or
distributed implementation, and an em dash no identified equivalent. General
libraries were excluded because API and shell workflows differ. Tables S2 and
S3 retain evidence and decision rules; absence does not imply impossibility.

### Independent Output and Input-Case Checks

With seed 20260716, rooted RF distances were compared with DendroPy for 50 tree
pairs, clade frequencies and majority consensus with direct calculations for 20
collections each, and Mk marginal probabilities with exhaustive enumeration
for 100 four-tip trees. Seeded `shuffle` outputs were compared bytewise (Table
S4).

A predefined 12-case corpus tested whether `validate` reported declared
representation, branch, topology, support, rooting, and taxon-set conditions.
It assessed reporting rather than recovery or repair (Table S6).

### Scaling and Cross-Tool Timing

Preflight scaling used balanced and caterpillar trees with 128--32,768 tips;
consensus timing used 10, 50, or 250 rooted 256-tip trees with NWKIT, Gotree,
and PhyKIT. Three runs measured startup-inclusive time and process-tree peak
memory on a 12-core x86-64 Mac with 64 GB RAM. Differing semantics preclude
algorithmic ranking (Tables S7 and S8).

### Worked Biological Example

The worked example used the checksum-pinned CSUBST PEPC tree and C4 definitions
[@Fukushima2023]. It expanded tip traits, checked monophyly, fitted a two-state
equal-rates Mk model, and selected eight C4 leaves maximizing rooted
phylogenetic diversity. `worked_example.py` generates all commands and values;
the PD fraction compares selected and all-candidate root-to-tip branch unions.

## Results

### One Interface Spans Five Tree-Processing Task Families

The 33 commands comprised seven inspection and normalization, 13
transformation, four comparison or tree-set, five taxonomy or trait, and four
downstream operations (Fig. 1). They share parser and label semantics and return
Newick unless another output is requested. Primary trees remain pipeable when
commands use secondary reference trees or tables. Accordingly, the smoke test
removed a singleton, renamed four tips, doubled branch lengths, and appended
three ordered audit records without an intermediate tree file.

![Figure 1](figures/figure1_architecture.png){width=90%}

**Figure 1. Functional organization and shared interfaces in NWKIT 0.28.0.**
Command families describe user tasks and are not claims of algorithmic novelty.
All 33 functional commands read their primary input from standard input
(stdin) by default; 31 emit their primary result to standard output (stdout).
JSONL denotes JSON Lines; TSV denotes tab-separated values.

*Alt text:* A flow diagram passes Newick trees and metadata through shared input
rules to five shaded command-family boxes, then through shared output and audit
rules to Newick, tables, provenance records, and figures.

`diff` separates root, clade, value, and annotation differences; `compose` can
take topology, root, names, support, lengths, and properties from separate
trees. With different tip sets, it uses informative, unique shared-tip
projections. Reports distinguish transferred, unmatched, and ambiguous
requests instead of assigning by node position.

This supports a common multi-tree decision sequence. `diff` first establishes
whether candidates disagree in clades or only in rooting and values. `compose`
then applies each source independently to a selected target. Compatible
requests may succeed when another source is unmatched, whereas strict mode
requires complete transfer. A JSON manifest retains source roles, and the TSV
report records the clade-level basis of the result.

### Comparison Revealed Extensive Overlap and a Distinct Combination

All reviewed toolkits provided routine transformation and visualization, and
three supported standard-stream composition (Table 1). Gotree and PhyKIT
calculated consensus, while PhyKIT covered far more comparative, trait, network,
and phylogenomic analyses than the selected rows. NWKIT is therefore not a
comprehensive replacement.

Differences concerned combinations and interfaces. Among the reviewed CLIs,
NWKIT alone exposed strict ambiguous-Newick rejection and integrated MCMCtree
calibration; it also combined multiple rooting methods, clade differencing,
multi-source composition, and per-command audits. Gotree downloaded and pruned
NCBI taxonomy, whereas NWKIT linked species-labelled tips to NCBI/APG,
OpenTree, or TimeTree. NWKIT and PhyKIT both transferred annotations and
reconstructed categorical states. Thus, NWKIT's position reflects preflight,
shell composition, and co-occurring taxonomy, selection, tree-set, and dating
tasks, not exclusive algorithms.

**Table 1. Task-level capabilities in reviewed command-line toolkits.** Native
denotes a directly documented command-line workflow, Partial a narrower or
distributed implementation, and — that an equivalent was not identified in the
reviewed version. The rows describe operations used in the NWKIT workflows and
are not a general ranking; see Tables S2 and S3 for the complete matrix and
decision rules.

| User task | NWKIT 0.28.0 | Gotree 0.5.2 | PhyKIT 2.3.0 | Newick Utilities 1.6 |
|---|:---:|:---:|:---:|:---:|
| Standard-stream composition | Native | Native | Partial | Native |
| Explicit Newick interpretation and ambiguity rejection | Native | — | — | — |
| Consolidated preflight report | Native | Partial | Partial | Partial |
| Routine tree transformations | Native | Native | Native | Native |
| Rooting strategies | Native | Partial | Partial | Partial |
| Topology-aware annotation transfer | Native | Partial | Native | — |
| Tree-collection summaries | Native | Native | Native | Partial |
| Taxonomy-derived topology | Native | Partial | — | — |
| Trait-aware selection and monophyly | Native | Partial | Partial | — |
| Categorical ancestral-state analysis | Native | Partial | Native | — |
| MCMCtree calibration preparation | Native | — | — | — |
| Tree and trait visualization | Native | Native | Native | Native |

### A PEPC Workflow Connected Diagnosis, Inference, and Sampling

Fukushima and Pollock [-@Fukushima2023] used NWKIT 0.10.0 `shuffle` in 1,000
simulations and `constrain` for NCBI-informed gene-tree constraints. We revisited
their PEPC data with functions added subsequently.

The rooted, binary PEPC tree contained 71 unique tips, including 21 labelled C4.
Neither state was monophyletic; the C4 MRCA contained 48 C3 tips (Fig. 2). The
equal-rates Mk rate was 3.633 per branch-length unit, root P(C4) was 0.108, and
15 edges joined different maximum-a-posteriori states. The last value describes
one marginal reconstruction, not the number of transitions.

Greedy maximum-PD selection retained eight C4 tips, spanning grasses, sedges,
and eudicots, and 75.9% of all-candidate rooted PD. The example tests workflow
reproduction, not PEPC-tree or C4-model accuracy; no tree was edited between
preflight, diagnosis, reconstruction, sampling, and plotting.

![Figure 2](figures/figure2_pepc_workflow.png){width=90%}

**Figure 2. Worked C4-state analysis of the CSUBST phosphoenolpyruvate
carboxylase (PEPC) example.** a) Branch tone gives the NWKIT equal-rates Markov
k-state (Mk) marginal probability of C4; circles and triangles give observed C3
and C4 tip states, respectively. Black outlines identify the eight C4 tips
selected by greedy maximum phylogenetic diversity; labels give abbreviated
species names and sequence accessions. b) Outputs passed between the four
workflow stages. ER, equal rates; PD, phylogenetic diversity.

*Alt text:* A 71-tip rectangular tree has branches that vary along a marginal
C4-probability scale, with observed C3 tips shown as circles and C4 tips as
triangles. Eight C4 triangles are outlined and labelled by species and
accession. Four boxes summarize preflight, polyphyly, Mk reconstruction, and
retention of 75.9% of C4-tip rooted phylogenetic diversity.

### Output Checks and Scaling Delimited the Tested Range

NWKIT matched all 191 predefined output checks (Fig. 3a): 50 rooted RF pairs,
20 clade-frequency collections, 20 consensus collections, 100 Mk calculations
(maximum error 3.33 × 10^-16^), and seeded shuffle reproduction. Preflight
identified all 11 declared problems and passed the valid control (Fig. 3b).

All 30 preflight and 27 consensus runs completed. Median preflight time rose
from 1.49 s at 128 tips to 2.69 s (balanced) and 2.38 s (caterpillar) at 32,768
tips; peak memory was 261 and 267 MB (Fig. 3c). Python startup contributed to
the approximately 1.5-s lower bound.

For 250 256-tip trees, median consensus time was 2.60 s for NWKIT, 0.130 s for
Gotree, and 6.29 s for PhyKIT; peak memory was 158, 21, and 82 MB (Fig. 3d).
Gotree was also fastest at 10 and 50 trees. NWKIT handled the tested sizes, but
startup, runtime, consensus semantics, and input structure preclude a general
performance claim.

![Figure 3](figures/figure3_evaluation.png){width=90%}

**Figure 3. Independent output checks, predefined input cases, and scaling.**
a) Agreement counts for rooted Robinson--Foulds (RF) distance, clade frequency,
consensus, Markov k-state (Mk) marginal probabilities, and seed reproduction.
b) Expected outcomes for the valid control and 11 declared problem conditions.
c) Median wall time and range across three NWKIT preflight runs per tree size
and shape; the horizontal axis is logarithmic. d) Median wall time and range
for three majority-consensus command-line interface (CLI) runs; both axes are
logarithmic.

*Alt text:* Two upper panels show 100% agreement and 12 matched-outcome case
markers. Lower panels show modest growth in NWKIT preflight time to 32,768 tips
and differing consensus-time curves for NWKIT, Gotree, and PhyKIT.

## Discussion

A final tree can be archived even when its derivation cannot be reconstructed.
NWKIT expresses curation and analysis as named, pipe-compatible operations. Its
contribution is not command count, but explicit connections among preflight,
transformation, taxonomic references, tree collections, trait analysis,
sampling, and downstream preparation. Such histories are more inspectable and
reusable than manual edits or disposable scripts.

Tree differencing and composition make this layer more than independent
filters. A selected topology can receive a root, support, lengths, and
annotations from related analyses. Reported exact or projected clade matches
replace positional copying; ambiguous projections are withheld. Composition
reports and audits retain both mapping decisions and the computational path.

This distinction matters biologically. Rooting, support estimation,
time-scaling, and annotation answer different questions and may use related but
nonidentical taxon samples. Composition preserves those origins rather than
presenting the final Newick as if every property arose from one analysis. It
also exposes disagreement: a value lacking a unique clade match remains an
unresolved request rather than silently reassigned evidence.

Gotree was faster in the consensus runs, and PhyKIT includes many analyses
without NWKIT counterparts. ETE, DendroPy, and ape remain preferable for new
algorithms, custom data models, and programmatic control; NWKIT depends on,
rather than replaces, scientific Python.

NWKIT remains Newick-centered. Strict parsing cannot infer authorial intent,
and structural rootedness requires biological confirmation. Projected matching
rejects fewer than two shared descendants and nonunique projections; it cannot
recover information absent from shared sampling. NHX is representationally
limited, some commands use external services, and the Mk models do not replace
broader comparative modeling.

Projected matching is conservative, not an ancestral reconciliation method. It
does not infer where a missing taxon attaches, reconcile gene and species trees,
or adjudicate between roots. Those decisions must be made upstream or expressed
through the selected target.

The evaluation is bounded: the comparison is not exhaustive, checks used
defined small-tree semantics, and scaling used synthetic inputs on one
computer. NWKIT is not universally validated or error-proof. Executable
comparisons, inventories, independent checks, and cross-platform benchmarks
should expand with future releases, maintaining a transparent layer between
inference and specialized analysis.

Likewise, the audit establishes computational provenance, not the scientific
validity of chosen operations. A trace can show exactly how a tree was produced
while leaving its biological assumptions open to review; reproducibility
requires both forms of scrutiny.

## Acknowledgments

OpenAI Codex was used to assist with manuscript drafting, analysis-code
development, document formatting, and consistency checks. The author reviewed
and verified the material and accepts responsibility for the work.

## Funding

This work was supported by the Sofja Kovalevskaja Programme of the Alexander
von Humboldt Foundation, Human Frontier Science Program grant RGY0082/2021,
and JSPS KAKENHI JP23K20050.

## Author Contributions

Kenji Fukushima conceptualized the study, developed NWKIT, performed the
analyses, and wrote and revised the manuscript.

## Conflict of Interest

The author declares no conflict of interest.

## Data Availability

NWKIT source code is available at <https://github.com/kfuku52/nwkit>. A
versioned source and analysis snapshot will be archived in Zenodo. Comparison
records, results, benchmarks, worked-example inputs, and reproduction
instructions will be deposited in Dryad; identifiers will be added before
submission. The PEPC inputs are also available in the CSUBST repository at the
commit and checksums reported above.

## References

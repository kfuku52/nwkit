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

Phylogenetic inference seldom yields a single tree ready for biological
interpretation. Studies instead produce related trees that differ in topology,
rooting, branch lengths, support values, or annotations. Researchers must also
inspect and transform tree files, connect tips to taxonomic or trait data,
summarize tree collections, and perform downstream analyses. These operations
are often dispersed among manual edits and study-specific scripts, obscuring
how inferred trees become the objects ultimately interpreted. We developed
NWKIT, a command-line toolkit that organizes post-inference tree processing as
composable, pipe-compatible commands. NWKIT can transfer roots and map node
labels, support values, and branch lengths between trees according to clades
defined by shared taxa, allowing compatible information from nonidentical
phylogenetic hypotheses to be combined explicitly. The same interface supports
Newick inspection and normalization, topology-aware transformation, taxonomy-
and trait-informed operations, synthesis of tree collections,
phylogenetic-diversity sampling, and categorical ancestral-state
reconstruction. Consistent input, output, and tree-interpretation conventions
allow individual commands to be inspected, substituted, and rerun without
intermediate custom scripts. By making both comparisons among trees and
transformations applied to them executable and traceable, NWKIT turns
post-inference processing from an implicit prelude into a reproducible
analytical layer of systematic research.

**Keywords:** ancestral-state reconstruction; command line; Newick;
phylogenetics; reproducibility; tree processing

Phylogenetic inference is rarely the final computational step in a systematic
study. Although publications often emphasize inference methods or biological
conclusions, substantial practical work occurs between them. Inferred trees
pass to programs for reconciliation, divergence-time estimation, comparative
analysis, visualization, and reporting. Alternative analyses may yield trees
with different topologies or roots, whereas bootstrap, dating, and annotation
steps add support, branch lengths, and metadata to related trees. Investigators
must then reconcile taxon labels, retain or remove leaves, reinterpret support,
reroot topologies, transfer annotations, and summarize bootstrap or posterior
tree collections. Each operation can be simple in isolation, but together they
determine which tree enters a downstream analysis. When these transformations
are recorded only as manual edits, short scripts, and software-specific
conversions, the final tree may remain available while the path that produced
it does not.

Newick is compact and broadly supported, but it leaves several conventions to
implementations. A numeric token following a closing parenthesis, for example,
may be interpreted as support or as an internal-node name, a distinction
handled inconsistently among tree viewers and toolkits [@Czech2017]. Rooting may
be represented only by the degree of the root, quoted labels may be accepted
inconsistently, and extended annotations differ among programs. These
differences are not merely cosmetic: a tree can remain syntactically readable
while acquiring a different biological interpretation. Empty or duplicated tip
labels, a negative branch length, or a mismatched taxon set can likewise pass
unnoticed until a later program fails or, more seriously, produces an
unintended result.

Several mature projects already provide phylogenetic tree representations and
algorithms, including ape in R [@Paradis2004], DendroPy
[@Sukumaran2010; @Moreno2024], and ETE [@HuertaCepas2016] in Python. Other
projects expose tree operations directly in the shell. Newick Utilities introduced automatable
filters for high-throughput Newick processing [@Junier2010]; Gotree provides
chainable commands and a Go application programming interface [@Lemoine2021];
and PhyKIT combines alignment, tree, and comparative functions for
phylogenomics [@Steenwyk2021]. The problem is therefore not a lack of generic
tree objects or shell utilities. It is that defensive Newick handling,
topology-aware transfer, tree-set synthesis, taxonomy-aware operations, and
selected downstream analyses remain distributed across interfaces with
different assumptions.

NWKIT began in 2019 as a collection of commands for making routine tree
transformations explicit and reusable in shell workflows. Version 0.27.0
expands this scope with preflight reporting and conversion, consensus and
clade-frequency calculation, taxonomy- and trait-aware selection,
visualization, and likelihood-based categorical ancestral-state
reconstruction. Here, we describe its design, audit task-level capabilities
against three command-line toolkits, and demonstrate an end-to-end biological
workflow. We further compare selected outputs with independent calculations,
exercise a predefined corpus of problematic inputs, and measure wall time and
peak memory across tree shapes and collection sizes. Together, these analyses
evaluate NWKIT as a reproducible layer between tree inference and biological
interpretation rather than treating command count as evidence of utility.

## Materials and Methods

### Software Design and Implementation

Analyses used NWKIT 0.27.0 at source commit `f71dc345ac83`. NWKIT is implemented
in Python, requires Python 3.10 or later, and uses ETE 4.4.0 for its primary tree
representation. Biopython, NumPy, SciPy, pandas, and Matplotlib support sequence
input, numerical analysis, tabular output, and plotting. The source is
distributed under the MIT License, packaged releases are available through
Bioconda [@Gruning2018], and command-specific documentation is maintained in
the project wiki.

The command-line interface contains 30 functional commands plus a help
pseudo-command. All 30 read their primary input from standard input by default
and also accept an explicit path. Twenty-eight write their primary result to
standard output by default; `draw` and `image` instead create graphical or image
files.
Commands can therefore be joined by shell pipes when one command's tree output
is the next command's input. Stochastic operations in `asr`, `shuffle`, and
`skim` expose an integer seed.

Input handling separates syntax from biological interpretation. Users may
select an ETE Newick parser format, request automatic detection, or use a strict
automatic mode that rejects ambiguous unquoted numeric internal labels. A
separate option controls quoted node names. Shared species-label parsers support
genus--species prefixes, qualifier-aware labels, regular expressions, and
mapping tables. These choices are command options because silently guessing a
species or support convention can alter a downstream analysis.

### Functional Organization

Commands were grouped into five user-task families (Fig. 1; Table S9):
inspection and normalization, transformation, tree-set summaries, taxonomy or
trait operations, and downstream analysis or presentation. The tree-set
commands calculate rooted Robinson--Foulds distance [@RobinsonFoulds1981], clade
frequencies, and consensus trees. Taxonomy and trait commands generate
constraints, diagnose monophyly, and select representatives by phylogenetic
diversity [@Faith1992]. Rooting includes minimal ancestor deviation [@Tria2017]
and NCBI, Open Tree of Life, or TimeTree references [@Schoch2020;
@Hinchliff2015; @Kumar2022], while constraint generation can incorporate APG IV
[@APG2016]. Selected downstream commands implement Mk-model ancestral-state
reconstruction [@Lewis2001], prepare MCMCtree inputs, and create visual outputs.

A root can be transferred from a reference tree, whereas node names, support
values, and branch lengths are mapped by matching descendant-taxon sets.
Unmatched clades remain unchanged or receive a user-specified fill value, so
compatible information can be combined without assuming identical topologies.

### Software Tests and Release Checks

The pytest suite was run in an isolated Python 3.12 environment with ETE 4.4.0
and Biopython 1.87; 597 tests passed at the manuscript commit. Tests covered
exact outputs, error cases, regression, round trips, seeded operations,
threading agreement, and documentation examples. Continuous integration spans
Python 3.10--3.14, performs static and distribution checks, and tests the
installed wheel. An independent inventory script extracted command counts and
defaults, and a process-level smoke test joined `sanitize`, `rename`, and
`rescale` through standard streams.

### Task-Level Software Comparison

We compared NWKIT 0.27.0 with Gotree 0.5.2, PhyKIT 2.3.0, and Newick Utilities
1.6 as available on 16 July 2026, using current documentation and representative
commands. Rows denote user tasks rather than similarly named functions.
“Native” required a documented command-line workflow, “Partial” a narrower or
distributed implementation, and an em dash means no equivalent was identified.
General libraries were excluded because APIs and shell operations impose
different workflow costs. Tables S2 and S3 retain the evidence and decision
rules; absence is not evidence that a task cannot be programmed from lower-level
components.

### Independent Output and Input-Case Checks

All checks used random seed 20260716. We compared rooted RF distances with
DendroPy for 50 tree pairs, clade frequencies and majority consensus with direct
descendant-set calculations for 20 collections each, and Mk marginal
probabilities with exhaustive internal-state enumeration for 100 four-tip
trees. Seeded `shuffle` outputs were compared byte for byte. Table S4 records
the full designs.

A predefined 12-case corpus tested whether `validate` reported declared
representation, branch, topology, support, rooting, and taxon-set conditions.
It assessed reporting rather than recovery or repair (Table S6).

### Scaling and Cross-Tool Timing

Preflight scaling used balanced and caterpillar trees from 128 to 32,768 tips.
Majority-consensus timing used 10, 50, or 250 rooted 256-tip trees with NWKIT,
Gotree, and PhyKIT. Three runs per condition measured startup-inclusive wall
time and process-tree peak memory on a 12-core x86-64 Mac with 64 GB of memory.
Because command semantics differ, timings characterize operating costs rather
than algorithmic rank (Tables S7 and S8).

### Worked Biological Example

The worked example used the CSUBST PEPC tree and C4 photosynthesis foreground
definitions from a published convergence analysis [@Fukushima2023], pinned by
commit and SHA-256 checksum. The workflow expanded trait labels to all tips,
checked monophyly, fitted a two-state equal-rates Mk model, and selected eight C4
leaves maximizing rooted phylogenetic diversity. `worked_example.py` generates
every command and derived value. The rooted-PD fraction divides the selected
root-to-tip branch union by the corresponding union for all C4 candidates.

## Results

### One Interface Spans Five Tree-Processing Task Families

The 30 commands formed five practical groups: seven for inspection and
normalization, 12 for transformation, three for tree-set summaries, four for
taxonomy or trait operations, and four for downstream analysis or presentation
(Fig. 1). The categories share parser and label semantics, so a choice such as
strict automatic Newick interpretation is available across operations rather
than implemented separately in each workflow. Likewise, tree-producing commands
return Newick unless another output is requested. A curation sequence can
therefore be recorded as a shell pipeline rather than as intermediate files and
manual edits. Commands with secondary inputs, such as a reference tree or trait
table, still accept the primary tree through standard input. The three-command
stream smoke test removed a singleton, renamed all four tips, and doubled every
remaining branch length without an intermediate file.

![Figure 1](figures/figure1_architecture.png){width=100%}

**Figure 1. Functional organization and shared interfaces in NWKIT 0.27.0.**
Command families describe user tasks and are not claims of algorithmic novelty.
All 30 functional commands read their primary input from standard input
(stdin) by default; 28 emit their primary result to standard output (stdout).
TSV denotes tab-separated values.

*Alt text:* A flow diagram passes Newick trees and metadata through shared input
rules to five shaded command-family boxes, then through shared output rules to
Newick, tables, and figures.

Tree-to-tree transfer extends this composition beyond a linear pipeline. NWKIT
can take the root from one reference tree and map node labels, support values,
or branch lengths from another Newick tree to clades shared with the target.
The trees must contain the same uniquely named tips, but their topologies need
not be identical; unmatched clades are reported rather than assigned by node
position. Topology, rooting, support, and branch lengths can therefore originate
from different analyses and be combined only where their phylogenetic
hypotheses are compatible.

### Comparison Revealed Extensive Overlap and a Distinct Combination

All four reviewed toolkits provided routine transformation and visualization,
and three supported standard-stream composition directly (Table 1). Gotree and
PhyKIT both calculated consensus trees, and PhyKIT 2.3.0 provided a much broader
set of comparative, trait, network, and phylogenomic analyses than represented
by the selected rows. These overlaps argue against describing NWKIT as a
comprehensive replacement.

The audited differences instead concerned combinations and interface details.
NWKIT alone among the reviewed CLIs exposed an explicit strict mode for
ambiguous Newick interpretation and an integrated MCMCtree-calibration workflow.
Its rooting command combined outgroup and midpoint methods with minimal ancestor
deviation and taxonomy-derived roots. Gotree could download and prune NCBI
taxonomy, whereas NWKIT additionally mapped species-labelled tips to NCBI/APG,
OpenTree, or TimeTree workflows. NWKIT and PhyKIT both transferred annotations
by topology and implemented categorical ancestral-state workflows. The table
therefore positions NWKIT around preflight handling, shell composition, and the
co-occurrence of taxonomy, selection, tree-set, and dating-preparation tasks,
not around exclusive ownership of common algorithms.

**Table 1. Task-level capabilities in reviewed command-line toolkits.** Native
denotes a directly documented command-line workflow, Partial a narrower or
distributed implementation, and — that an equivalent was not identified in the
reviewed version. The rows describe operations used in the NWKIT workflows and
are not a general ranking; see Tables S2 and S3 for the complete matrix and
decision rules.

| User task | NWKIT 0.27.0 | Gotree 0.5.2 | PhyKIT 2.3.0 | Newick Utilities 1.6 |
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

Earlier NWKIT versions had already supported published phylogenetic workflows.
Fukushima and Pollock [-@Fukushima2023] used `shuffle` in NWKIT 0.10.0 to
generate randomized trees in 1,000 simulations and `constrain` to derive NCBI
taxonomy-informed topological constraints for gene-tree inference. The present
worked example revisited the PEPC data from that study with functions added in
later releases.

The pinned PEPC tree contained 71 unique tips and was rooted and binary.
Twenty-one tips were labelled C4. Neither the C4 nor C3 set was monophyletic:
the MRCA of C4-labelled tips contained 48 C3 intruders, consistent with repeated
origins and losses represented in this gene tree (Fig. 2). The fitted equal-rates Mk
model had an estimated transition rate of 3.633 per unit branch length. With an
equal root prior, the marginal probability of C4 at the root was 0.108. Fifteen
edges connected nodes with different maximum-a-posteriori states. This last
count is a descriptive summary of a single marginal reconstruction, not an
estimate of the number of evolutionary transitions.

Greedy maximum-PD selection retained eight of the 21 C4 tips and 75.9% of the
rooted PD represented by all C4 candidates. The selected accessions spanned C4
lineages in grasses, sedges, and eudicots. Importantly, the example did not
establish the accuracy of the PEPC tree or the C4 model. It showed that a
versioned sequence of preflight, trait diagnosis, reconstruction, sampling, and
plotting operations could be rerun without editing the tree between steps.

![Figure 2](figures/figure2_pepc_workflow.png){width=100%}

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

NWKIT matched the comparison result in all 191 predefined output checks (Fig.
3a). Rooted RF distance agreed exactly for 50 tree pairs. Clade sets and
frequencies agreed for all 20 tree collections, as did majority-consensus clades
for another 20 collections. Across 100 Mk checks, the largest absolute
difference from exhaustive enumeration was 2.22 × 10^-16^. Repeated seeded
shuffling produced identical output, whereas the second seed changed the
output. The preflight report identified the declared condition in all 11
problem cases and reported no issue for the valid control (Fig. 3b).

All 30 single-tree preflight runs and all 27 consensus runs completed
successfully. Median preflight time rose
from 1.49 s at 128 tips to 2.69 s for the 32,768-tip balanced tree and 2.38 s for
the caterpillar tree; corresponding peak memory was 261 and 267 MB (Fig. 3c).
The approximately 1.5-s lower bound largely reflected Python process and module
startup on this machine.

For 250 trees of 256 tips, median majority-consensus time was 2.60 s for NWKIT,
0.130 s for Gotree, and 6.29 s for PhyKIT; peak resident memory was 158, 21, and
82 MB, respectively (Fig. 3d). Gotree was also the fastest at 10 and 50 trees.
These results show that NWKIT handled the tested collection sizes, but they do
not support a general performance advantage. Process startup, language runtime,
consensus details, and input structure all contribute to the observed values.

![Figure 3](figures/figure3_evaluation.png){width=100%}

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

The final tree of a phylogenetic study can be archived even when the decisions
that produced it cannot be reconstructed. NWKIT addresses this gap by
expressing tree curation and analysis as named operations with shared,
pipe-compatible conventions. Its contribution is not its command count, but
the way explicit input rules connect preflight reporting, transformations,
taxonomic references, tree collections, trait analysis, sampling, and
preparation for downstream programs. For recurring operations, such a command
history is easier to inspect, test, and reuse than a succession of manual edits
or disposable scripts.

Tree-to-tree transfer illustrates why this layer is more than a collection of
independent filters. A topology selected under one criterion can receive a root
from another tree and support values, branch lengths, or annotations from
related analyses. Transfer is restricted to compatible clades and reports
unmatched nodes rather than copying information by node position. This permits
compatible components of phylogenetic hypotheses to be assembled reproducibly
without treating different topologies as interchangeable.

The comparison also identifies preferable alternatives. Gotree was
substantially faster in the present consensus runs. PhyKIT includes many
alignment, comparative, network, and phylogenomic analyses without NWKIT
counterparts. ETE, DendroPy, and ape remain appropriate for new algorithms,
custom data models, and programmatic control; NWKIT depends on ETE and
scientific Python rather than replacing them.

NWKIT remains Newick-centered rather than a general container for characters,
alignments, and trees. Strict parsing can reject ambiguity but cannot infer
authorial intent, and structural rootedness still requires biological
confirmation. Some commands depend on external services, and the implemented
Mk models do not replace broader comparative modeling.

The evaluation is likewise bounded: the capability table is not an exhaustive
software census, output checks used defined small-tree semantics, and scaling
used synthetic inputs on one computer. These results do not justify describing
NWKIT as universally validated or error-proof. Future releases should preserve
the executable comparison and inventory, expand independent checks, and add
cross-platform benchmarks. NWKIT can thereby remain a transparent layer
between tree inference and specialized downstream analysis while complementing
the libraries and toolkits on which phylogenetic workflows depend.

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

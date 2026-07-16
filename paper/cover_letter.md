# Draft cover letter

Dear Editors,

Please consider my manuscript, “NWKIT: Composable Command-Line Workflows for
Phylogenetic Tree Curation and Analysis,” for publication as a *Software for
Systematics and Evolution* article in *Systematic Biology*.

Phylogenetic trees commonly pass through several programs between inference
and biological interpretation. At these interfaces, small differences in taxon
sets, Newick conventions, rooting, annotations, and metadata can make a
workflow difficult to inspect or reproduce. NWKIT 0.27.0 provides 30
pipe-compatible commands that combine explicit Newick interpretation and
preflight reporting with topology-aware transformations, tree-collection
summaries, taxonomy- and trait-informed operations, categorical
ancestral-state reconstruction, and preparation for downstream programs.

The manuscript does not position NWKIT as a replacement for established
libraries or command-line toolkits. Instead, it documents the particular
combination of operations available through common input and output semantics.
I audit task-level overlap with Gotree, PhyKIT, and Newick Utilities and report
where specialized alternatives are broader or faster. In predefined checks,
NWKIT matched independent calculations for rooted Robinson--Foulds distance,
clade frequency, majority consensus, and Mk marginal probabilities. A public
71-tip PEPC example connects preflight checks, monophyly diagnosis,
ancestral-state reconstruction, and phylogenetic-diversity sampling without
manual tree editing.

I believe the work fits *Systematic Biology* because it addresses a recurrent
reproducibility problem across systematic and evolutionary studies rather than
a single organismal system. The manuscript includes comparative timings, a
real biological example, sample data, accessible figures, and a complete
supplementary audit trail. NWKIT is open source under the MIT License and is
distributed through GitHub and package repositories.

The source snapshot and analysis code have been prepared for Zenodo, and the
input data, raw results, figures, and manifests have been prepared for Dryad.
Active review links will be supplied with the submission.

Thank you for considering this manuscript.

Sincerely,

Kenji Fukushima\
Center for Frontier Research, National Institute of Genetics\
Genetics Program, Graduate Institute for Advanced Studies, SOKENDAI\
1111 Yata, Mishima, Shizuoka 411-8540, Japan\
Telephone: +81-55-981-6751\
<kenji.fukushima@nig.ac.jp>\
ORCID: <https://orcid.org/0000-0002-2353-9274>

## Author confirmation before sending

The signed letter should additionally confirm that the work is not under
consideration elsewhere, that all authors approve the submission, that any
conflicts have been disclosed, and that NWKIT will be supported for at least
two years after publication. The affiliation and contact details above were
compiled from current institutional and published sources and should be checked
before sending.

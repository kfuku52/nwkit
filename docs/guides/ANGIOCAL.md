# AngioCal fossil minimum-age constraints

`nwkit mcmctree --angiocal v1.0` imports the 238-record
[AngioCal v1.0 dataset](https://github.com/eflowerproject/angiocal/tree/v1.0)
and produces a PAML calibration tree. It downloads the original XLS from a
fixed commit and verifies its SHA-256, then reads the matching normalized
records bundled with NWKIT. The bundled TSV is separately checksum-verified
and preserves all imported fields and original worksheet row numbers; reports
and audits still identify the original XLS. Official v1.0 and normalized TSV
inputs need no Excel-reading dependency. Custom XLS input also uses NWKIT's
built-in reader without additional dependencies.
AngioCal v1.1, published with the 2024
angiosperm phylogeny, has a different distribution and is not yet supported.

```sh
nwkit mcmctree -i species.nwk --angiocal v1.0 \
  --species-parser taxonomic --lower-tail-prob 1e-300 \
  --download-dir downloads --report fossils.tsv \
  --add-header yes -o calibrated.tre
```

The tree must be rooted and fully bifurcating, with unique, nonempty tip names.
PAML tip identifiers must be at most 100 ASCII characters without whitespace,
quotes, control characters, Newick delimiters or `#`, and cannot consist solely of digits; use
matching identifiers in the alignment. PAML can interpret digit-only tips as
species numbers instead of preserving their names.
The length limit matches the default PAML species-name buffer. PAML's tree reader does not
apply Newick quoting to species names, so quoting a colon or space cannot make
such a name compatible. NWKIT rejects these inputs rather than silently renaming
tips. `--angiocal` is
exclusive with `--timetree point|ci`, `--posterior`, and manual selection through
`--left-species`, `--right-species`, `--lower-bound`, or `--upper-bound`.
The existing manual and TimeTree modes retain their previous behavior.

## Placement and sampling

A fossil's **minimum age** constrains its assigned divergence; it is not the
estimated divergence age. The dataset's crown/stem placement is preserved.
There are three ways to place a fossil, in this order:

1. `--calibration-map-tsv` supplies two anchors for that fossil ID. Their MRCA
   is the final calibrated node, including for stem fossils. Anchors can be
   original leaf names or species labels parsed by the shared species parser.
   The two anchors must identify different species.
2. An unambiguous internal node label equal to the original `Node calibrated`
   field or `crown CLADE`/`stem CLADE` asserts the final node directly. An
   internal label equal to `CLADE` asserts its crown; a stem fossil is placed
   on that node's parent. For example, `((A_a,B_b)Testaceae,C_c);` identifies
   the crown of Testaceae and its stem divergence. Biological labels are
   interpreted as user assertions, not independently verified phylogenies.
   Tip names are not internal clade assertions, even when they equal `CLADE`.
3. Unlabeled **stem** clades can be resolved with NCBI taxonomy. The sampled
   clade must be monophyletic and have an outside lineage. Taxon unions such
   as `Chloranthus+Sarcandra` use the union of their sampled descendants.
   Ambiguous taxids, unknown tip taxonomy, unavailable clade names and
   nonmonophyletic clades are reported as exclusions.

Crown fossils require explicit anchors or a biological node label. The MRCA
of a few sampled members can be younger than the full clade's crown, so
species names alone are insufficient to place its crown minimum. At least
two distinct species must occur below a crown node. Duplicate tips for a
species are supported only when monophyletic.

For an automatically resolved stem, missing sister lineages can cause its
nearest sampled outside lineage to identify an **older ancestral divergence**.
The fossil minimum remains a conservative lower bound there; the report marks
this as `taxonomy_stem` / `sampled_stem_ancestor`. Use explicit anchors when
the exact stem divergence is required. A clade covering the entire input tree
has no represented stem and is excluded, rather than placed at its crown root.
Automatic stem placement also checks that the candidate node and both of its
branches are compatible with the NCBI clades represented by the sampled tips.
If the rooted tree conflicts with one of those clades, the fossil is excluded
with reason `conflicting_tree_taxonomy`; correct the rooting or topology before
using the tree for dating. Explicit anchors and biological node labels remain
user assertions and are not checked against NCBI taxonomy.

Example mapping (IDs refer to `NFos` in the original table):

```tsv
fossil_id	left_species	right_species
330	Amborella_trichopoda	Ginkgo_biloba
```

This explicitly places the Angiospermae stem fossil at the divergence from
Ginkgo on a suitable rooted tree. The user is responsible for the biological
validity of each explicit assignment. Duplicate or unknown fossil IDs and
missing anchors fail before writing the result.

## Ages, priors and existing constraints

Original ages and report columns use **Ma**. `--time-unit-ma` specifies the
number of Ma in one output unit: the default is `1`; `100` converts 132.1 Ma
to 1.321 output units. Existing input-tree constraints must already use that
output unit. Branch lengths are removed from the PAML calibration tree.

Minimum-only constraints use the existing MCMCtree
`L(age, offset, scale, tail_probability)` syntax. The CLI's
`--lower-offset`, `--lower-scale`, and `--lower-tail-prob` supply its prior
parameters; these parameters are not inferred from fossils. The existing
default lower tail is `0.025`, a soft bound. `1e-300` approximates a hard bound
as in the existing manual mode. The report contains the actual serialized prior.

When several fossils map to one node, their largest minimum age is applied.
All contributing records remain in the report. Quality scores and primary
references are retained without assigning an automatic quality threshold.
Fossil stratigraphic ranges do not create maximum-age constraints. Supply a
justified root maximum in the input tree, for example `U(2.47, 1e-300)` when
the output unit is 100 Ma; NWKIT does not select a root maximum automatically.

Compatible existing upper-only constraints combine with the fossil lower bound
as `B(...)`, retaining the existing upper age and tail without rounding. An
existing lower constraint that already satisfies the minimum is preserved. A weaker
existing lower prior is rejected because silently changing its support would
replace a separately specified prior; remove that calibration before importing.
Upper ages younger than an imported minimum also fail.
Combining an upper-only age equal to the fossil minimum is rejected because
it would create a degenerate bounded prior.
Bounded priors require a positive age range and lower/upper tail probabilities
whose sum is less than one.
Lower priors also require a positive age and representable derived location
and scale; a zero age or numeric overflow is rejected before output.

Existing priors must use fully specified `L(age, offset, scale, tail)`,
`U(age, tail)`, `B(lower, upper, lower_tail, upper_tail)`, or legacy `>`/`<`
bounds. Lowercase `l`/`u`/`b` constructors are normalized to uppercase so that
PAML recognizes them. Shorthand such as `U(30)`, other distributions such as
`G(...)`, and duplication annotations such as `#1` are rejected before import;
expand shorthand using its intended parameters or remove unsupported priors.
They are never silently discarded or replaced by a fossil constraint.

`@age` annotations are rejected in AngioCal mode, including on nodes without
an imported fossil. They are age annotations, not fixed-age fossil priors in
[PAML's calibration parser](https://github.com/abacus-gene/paml/blob/4c7902fe972737ef5e80bb18f159a6e6acace3d3/src/treesub.c#L8653).
Keeping one in place of an imported lower bound would discard the fossil prior.
Remove it or supply an appropriate `L`/`U`/`B` prior. This restriction does not
change the existing manual or TimeTree point-mode interfaces.

Nominal bounds are checked throughout the tree: an ancestor's upper age must
be older than a descendant's lower age to allow positive intervening time.
This is an input-consistency check,
not a calculation of the joint calibration prior or its soft-tail probabilities.
`--min-clade-prop` exclusions are recorded; equal-age constraints on different
nodes are retained rather than silently removed by the legacy cleanup.

## Offline inputs, caches and reports

`--download-dir downloads` stores the pinned XLS under
`downloads/angiocal/v1.0/`. With `auto`, AngioCal uses
`~/.cache/nwkit/angiocal/v1.0/`; NCBI keeps its existing separate ETE cache.
A verified AngioCal cache is reused without a request. A checksum mismatch is
an error and does not silently replace the file.

`--angiocal-file PATH` accepts a local original XLS or a normalized UTF-8 TSV
with these required columns:

```tsv
fossil_id	fossil_taxon	minimum_age_ma	placement	clade
1	Synthetic example	12.5	crown	Testaceae
```

The exact pinned official workbook, modified XLS workbooks and TSV input all
work without `xlrd` or an extra installation. The built-in reader supports
OLE-contained BIFF5/8 workbooks, BIFF4 workbooks and standalone BIFF2/3/4
worksheets. It reads stored cell values, including cached formula results;
it does not recalculate formulas. Encrypted workbooks and XLSX files are
rejected. Date, boolean and error cells retain their types and are rejected
in fossil IDs or ages. The former `xls` installation extra remains an empty
compatibility alias.

Additional TSV columns can retain `node_calibrated`, `safe_minimum_age`,
`age_quality_score`, `node_assignment_score`, `reconciliation_score`,
`relationship_reference`, and `age_reference`. Local files may be subsets;
their actual checksum is reported, without claiming that they are the pinned
official file. IDs must be unique positive integers and ages finite and positive.
Decimal-text identifiers are compared exactly, without binary-float rounding.
An optional `node_calibrated` crown/stem prefix must agree with `placement`.
The original XLS container is detected from its bytes, including files without
an extension and symlinks; source paths and checksums refer to the resolved file.
XLS IDs and ages must use numeric or decimal-text cells. Excel Boolean and date
cells are rejected, as are Excel errors in imported fields, so their internal
codes cannot be mistaken for fossil ages or metadata.

For fully offline placement, use a local source and
`--angiocal-taxonomy no` with biological labels or an explicit mapping.
`--calibration-map-tsv -` accepts stdin; the primary tree and any other input
cannot simultaneously use stdin. The XLS/TSV source and report require file paths.

`--report fossils.tsv` writes one row per fossil with source version, URL/path,
SHA-256, original row and Ma age, quality scores, references, mapping method,
target-tip JSON and its stable clade hash, selected minimum, output unit,
serialized constraint, status and reason. Status is `applied` for the strongest
minimum (including ties), `supporting` for weaker contributing minima, and
`skipped` for exclusions.

The tree and report are staged together; a handled write or commit failure
restores the previous pair. Inputs cannot be overwritten by either output.
The automatic source-cache path is also protected from outputs and `--audit`
before downloading or opening logs. With taxonomy enabled, NCBI cache files are
also protected. Binary cache writes verify their staging-file identity before
writing, as tree/report text writes do. An audit records the source hash even when
the cache is downloaded for the first time during that command.
If no fossil can be placed, the command fails and leaves the previous tree
untouched, while intentionally writing the exclusion report for diagnosis.
Stdout is emitted only after serialization and successful report installation;
already delivered stream bytes cannot be retracted.

Output trees retain prior annotations but omit the biological internal labels
used during placement. Regenerate from the original tree, or retain an explicit
fossil map when reimporting a calibrated output. Use the same `--time-unit-ma`
when reimporting: a PAML tree does not encode its Ma conversion factor.

The synthetic test XLS at `tests/data/angiocal-v1.0-synthetic.xls` contains
invented evidence for parser tests and must not be used for scientific dating.

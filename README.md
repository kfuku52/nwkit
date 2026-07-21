![](logo/logo_nwkit_large.png)

[![Run Tests](https://github.com/kfuku52/nwkit/actions/workflows/tests.yml/badge.svg)](https://github.com/kfuku52/nwkit/actions/workflows/tests.yml)
[![GitHub release](https://img.shields.io/github/v/tag/kfuku52/nwkit?label=release)](https://github.com/kfuku52/nwkit/releases)
[![Bioconda](https://img.shields.io/conda/vn/bioconda/nwkit.svg)](https://anaconda.org/bioconda/nwkit)
[![Python](https://img.shields.io/badge/python-3.10%20%7C%203.11%20%7C%203.12%20%7C%203.13%20%7C%203.14-blue)](https://github.com/kfuku52/nwkit)
[![Platforms](https://img.shields.io/conda/pn/bioconda/nwkit.svg)](https://anaconda.org/bioconda/nwkit)
[![Downloads](https://img.shields.io/conda/dn/bioconda/nwkit.svg)](https://anaconda.org/bioconda/nwkit)
[![License](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)

## Overview
**NWKIT** ([/njuːkit/](https://ipa-reader.com/?text=nju%CB%90kit&voice=Joanna)) is a toolkit for manipulating phylogenetic trees in the [Newick format](https://en.wikipedia.org/wiki/Newick_format).

## Installation

Packaged releases of NWKIT are available from [Bioconda](https://anaconda.org/bioconda/nwkit). The repository can be ahead of the current Bioconda package; compare the release and Bioconda badges above when selecting a version. For users requiring a `conda` installation, please refer to [Miniforge](https://github.com/conda-forge/miniforge) for a lightweight conda environment.

#### Install from Bioconda

```
conda install bioconda::nwkit
```

#### Verify the installation by displaying the available options

```
nwkit -h
```

#### (For advanced users) Install the development version from GitHub

```
pip install git+https://github.com/kfuku52/nwkit
```

NWKIT requires Python 3.10 or newer.

#### Optional dependencies for image post-processing

`nwkit image` can normalize image format, trim margins, and resize/pad output files when the optional image-processing dependencies are installed:

```
pip install "nwkit[image]"
```

## Subcommands
See [Wiki](https://github.com/kfuku52/nwkit/wiki) for usage.

Shared option naming, standard-input rules, TSV schemas, missing-value policy,
and output-column vocabulary are defined in
[CLI and TSV conventions](https://github.com/kfuku52/nwkit/blob/master/CLI_TSV_CONVENTIONS.md).

- [`annotate`](https://github.com/kfuku52/nwkit/wiki/nwkit-annotate): Attaching tip-table values and aggregating them as Newick properties
- [`asr`](https://github.com/kfuku52/nwkit/wiki/nwkit-asr): Inferring categorical ancestral states and imputing missing tip states under an Mk model
- [`constrain`](https://github.com/kfuku52/nwkit/wiki/nwkit-constrain): Generating a species-tree-like Newick file for topological constraint
- [`collapse`](https://github.com/kfuku52/nwkit/wiki/nwkit-collapse): Collapsing internal branches by support and/or branch length
- [`compose`](https://github.com/kfuku52/nwkit/wiki/nwkit-compose): Assembling compatible roots, values, and annotations from multiple trees
- [`cladefreq`](https://github.com/kfuku52/nwkit/wiki/nwkit-cladefreq): Summarizing clade frequencies across a tree collection
- [`consensus`](https://github.com/kfuku52/nwkit/wiki/nwkit-consensus): Generating a consensus tree or transferring consensus support to a reference tree
- [`diff`](https://github.com/kfuku52/nwkit/wiki/nwkit-diff): Reporting interpretable clade, root, value, and annotation differences between trees
- [`dist`](https://github.com/kfuku52/nwkit/wiki/nwkit-dist): Comparing tree topology and branch lengths with multiple distance metrics
- [`draw`](https://github.com/kfuku52/nwkit/wiki/nwkit-draw): Drawing a phylogenetic tree with optional speciation/duplication node markers
- [`drop`](https://github.com/kfuku52/nwkit/wiki/nwkit-drop): Removing node and branch information
- [`image`](https://github.com/kfuku52/nwkit/wiki/nwkit-image): Retrieving representative species images with license-aware filtering
- [`info`](https://github.com/kfuku52/nwkit/wiki/nwkit-info): Printing tree information
- [`intersection`](https://github.com/kfuku52/nwkit/wiki/nwkit-intersection): Dropping non-overlapping leaves/sequences between two trees or between a tree and an alignment
- [`label`](https://github.com/kfuku52/nwkit/wiki/nwkit-label): Adding unique node labels
- [`mark`](https://github.com/kfuku52/nwkit/wiki/nwkit-mark): Adding texts to node labels by identifying the targets with a leaf name regex
- [`mcmctree`](https://github.com/kfuku52/nwkit/wiki/nwkit-mcmctree): Introducing divergence time constraints for PAML's mcmctree
- [`monophyly`](https://github.com/kfuku52/nwkit/wiki/nwkit-monophyly): Assessing whether species or trait-defined groups are monophyletic
- [`nhx2nwk`](https://github.com/kfuku52/nwkit/wiki/nwkit-nhx2nwk): Generating Newick from NHX
- [`nwk2table`](https://github.com/kfuku52/nwkit/wiki/nwkit-nwk2table): Converting a Newick tree into a parent-child table
- [`printlabel`](https://github.com/kfuku52/nwkit/wiki/nwkit-printlabel): Searching and printing node labels
- [`prune`](https://github.com/kfuku52/nwkit/wiki/nwkit-prune): Pruning leaves
- [`rename`](https://github.com/kfuku52/nwkit/wiki/nwkit-rename): Renaming nodes using a TSV mapping or regular expression
- [`rescale`](https://github.com/kfuku52/nwkit/wiki/nwkit-rescale): Rescale branch length with a given factor
- [`root`](https://github.com/kfuku52/nwkit/wiki/nwkit-root): Placing or transferring the tree root
- [`sanitize`](https://github.com/kfuku52/nwkit/wiki/nwkit-sanitize): Eliminating non-standard Newick flavors
- [`sample`](https://github.com/kfuku52/nwkit/wiki/nwkit-sample): Selecting a representative leaf subset by filters, ranks, and sampling method
- [`shuffle`](https://github.com/kfuku52/nwkit/wiki/nwkit-shuffle): Shuffling branches and/or labels
- [`skim`](https://github.com/kfuku52/nwkit/wiki/nwkit-skim): Sampling leaves from clades with shared traits
- [`subtree`](https://github.com/kfuku52/nwkit/wiki/nwkit-subtree): Generating a subtree Newick file
- [`table2nwk`](https://github.com/kfuku52/nwkit/wiki/nwkit-table2nwk): Converting a parent-child table into a Newick tree
- [`transfer`](https://github.com/kfuku52/nwkit/wiki/nwkit-transfer): Transferring information between trees
- [`validate`](https://github.com/kfuku52/nwkit/wiki/nwkit-validate): Validating one or more Newick trees and reporting structural issues

## Tree comparison, composition, and provenance

`diff` reports how two trees differ before information is combined. `compose`
then assembles compatible components from named sources and records every match
or conflict in a TSV report:

```sh
nwkit diff -i topology.nwk -i2 bootstrap.nwk -o differences.tsv
nwkit compose -i topology.nwk \
  --root-source rooted.nwk \
  --support-source bootstrap.nwk \
  --length-source chronogram.nwk \
  --property-source habitat=state@annotated.nhx \
  --report composition.tsv \
  -o combined.nwk
```

`transfer`, `compose`, `diff`, and transferred roots accept
`--taxon-mode intersection`. Reports distinguish `exact_match` from
`projected_match`: the latter is unique only in the trees induced by their
shared tips and does not establish that the original branches are identical.
Ambiguous projections are reported and left unchanged. Projected node names
and NHX properties can be transferred under the default `compatible-only`
policy, whereas projected support values and branch lengths require an explicit
`--allow-projected-values yes`. Strict policy accepts exact matches only.

Matching defaults to rooted descendant clades. For edge-associated values in
trees with different rootings, `--match-basis split` instead uses canonical,
root-independent splits:

```sh
nwkit transfer -i target.nwk -i2 rerooted-source.nwk \
  --support yes --match-basis split --report transfer.tsv \
  -o supported.nwk
```

Root alignment preserves support values, internal names, and NHX properties on
their original unrooted branch splits. A bifurcating root represents one
unrooted edge as two child branches, so `transfer` and `compose` resolve that
representation explicitly. By default, one source annotation is copied to
both target halves, equal source annotations are treated as one value, and
conflicting annotations follow the source half with the same projected
descendant taxa. Branch lengths use the source edge total and preserve the
target root-position ratio; a zero-length target root edge is split equally.
No result depends on traversal order.

Use repeatable `--root-edge-policy TARGET_PROPERTY=POLICY` options to override
the default per target property. `skip` leaves the edge unchanged,
`equal-only` accepts only a unique or equal annotation, `matching-side` follows
projected descendant taxa, `mean`, `min`, and `max` reduce numeric annotations,
and `edge-total` applies the branch-length rule. `auto` is the annotation
default; `edge-total` is the length default. `*` can set a fallback policy.
For example:

```sh
nwkit compose -i topology.nwk \
  --support-source bootstrap.nwk \
  --length-source chronogram.nwk \
  --match-basis split \
  --root-edge-policy support=mean \
  --root-edge-policy length=edge-total \
  --report composition.tsv \
  -o combined.nwk
```

A composition manifest accepts a top-level `root_edge_policies` object, and a
custom `properties` entry can contain `root_edge_policy`. Precedence is
top-level manifest, property entry, then CLI. Reports record the policy,
resolution, candidate counts, and all candidate values. Projected support and
length values still require `--allow-projected-values yes`, and strict mode
still rejects every projected match.

`compose` changes the target rooting only when `--root-source` is supplied.
The `midpoint`, `outgroup`, `mad`, `mv`, `taxonomy`, and `transfer` rooting
methods apply the same split-based annotation preservation. For tree
collections, `validate --require-same-rooting yes` compares the root
bipartition when leaf sets are identical.

Arbitrary NHX properties can be transferred or renamed directly:

```sh
nwkit transfer -i target.nwk -i2 source.nhx \
  --property color --property-map posterior=source_posterior \
  --report transfer.tsv -o annotated.nhx
```

`annotate` joins a TSV keyed by `leaf_name` to tree tips and can aggregate
values onto ancestral clades:

```sh
nwkit annotate -i tree.nwk --table traits.tsv \
  --columns habitat \
  --aggregate habitat:unique:shared_habitat \
  --aggregate body_mass:mean:mean_body_mass \
  -o annotated.nhx
```

Every functional command accepts `--audit PATH`. Each invocation appends one
JSON Lines record containing NWKIT version and arguments, input and output
SHA-256 hashes, input interpretation, random seeds, external-data settings,
warnings, runtime status, and captured command messages:

```sh
nwkit sanitize -i raw.nwk --audit workflow.audit.jsonl |
  nwkit root --method midpoint --audit workflow.audit.jsonl > rooted.nwk
```

## Citation
There is no published paper on NWKIT itself, but we used and cited NWKIT in several papers including [Fukushima & Pollock (2023, Nat Ecol Evol 7: 155-170)](https://www.nature.com/articles/s41559-022-01932-7).

## Development

Install the development and optional image dependencies, then run the same checks used by CI:

```
pip install -e ".[dev,image]"
ruff check nwkit tests
pytest tests/ -q
python -m build
```

See [CHANGELOG.md](CHANGELOG.md) for changes and [RELEASING.md](RELEASING.md) for the release checklist.

# Licensing
This program is MIT-licensed. See [LICENSE](LICENSE) for details.

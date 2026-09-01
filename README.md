![](https://raw.githubusercontent.com/kfuku52/nwkit/master/logo/logo_nwkit_large_02.png)

[![Run Tests](https://github.com/kfuku52/nwkit/actions/workflows/tests.yml/badge.svg)](https://github.com/kfuku52/nwkit/actions/workflows/tests.yml)
[![GitHub release](https://img.shields.io/github/v/tag/kfuku52/nwkit?label=release)](https://github.com/kfuku52/nwkit/releases)
[![Bioconda](https://img.shields.io/conda/vn/bioconda/nwkit.svg)](https://anaconda.org/bioconda/nwkit)
[![Python](https://img.shields.io/badge/python-3.10%20%7C%203.11%20%7C%203.12%20%7C%203.13%20%7C%203.14-blue)](https://github.com/kfuku52/nwkit)
[![Platforms](https://img.shields.io/conda/pn/bioconda/nwkit.svg)](https://anaconda.org/bioconda/nwkit)
[![Downloads](https://img.shields.io/conda/dn/bioconda/nwkit.svg)](https://anaconda.org/bioconda/nwkit)
[![License](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)

## Overview
**NWKIT** ([/njuːkit/](https://ipa-reader.com/?text=nju%CB%90kit&voice=Joanna)) is a command-line toolkit for processing and visualizing phylogenetic trees and for phylogeny-aware comparative analysis.

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

#### Optional SVG rasterization support

Raster-image validation and processing use Pillow, which is installed with
NWKIT. The optional CairoSVG integration also needs the native Cairo library.
Install that library first (for example, `brew install cairo libffi` on macOS
or `sudo apt-get install libcairo2 libffi-dev` on Debian/Ubuntu), then install
NWKIT and its image extra from GitHub:

```
pip install "nwkit[image] @ git+https://github.com/kfuku52/nwkit.git"
```

From an existing source checkout, use `pip install -e ".[image]"`. If NWKIT
was installed from Bioconda, install the native Cairo library and then run
`pip install CairoSVG` in the same environment. These steps enable SVG
normalization, trimming, resizing with `nwkit image`, and rendering SVG tip
images with `nwkit draw`.

## Subcommands
See [Wiki](https://github.com/kfuku52/nwkit/wiki) for usage.

Shared option naming, standard-input rules, TSV schemas, missing-value policy,
and output-column vocabulary are defined in
[CLI and TSV conventions](https://github.com/kfuku52/nwkit/blob/master/CLI_TSV_CONVENTIONS.md).
The reconciled speciation-contrast calculation is derived step by step, with a
minimal worked example, in
[the mathematical guide](https://github.com/kfuku52/nwkit/blob/master/RECONCILED_SPECIATION_CONTRAST_MATH.md).

- [`annotate`](https://github.com/kfuku52/nwkit/wiki/nwkit-annotate): Attaching tip-table values and aggregating them as Newick properties
- [`asr`](https://github.com/kfuku52/nwkit/blob/master/ASR.md): Inferring discrete (Mk) or continuous (Brownian) ancestral traits and imputing missing tips with automatic trait-type detection
- [`constrain`](https://github.com/kfuku52/nwkit/wiki/nwkit-constrain): Generating a species-tree-like Newick file for topological constraint
- [`collapse`](https://github.com/kfuku52/nwkit/wiki/nwkit-collapse): Collapsing internal branches by support and/or branch length
- [`compose`](https://github.com/kfuku52/nwkit/wiki/nwkit-compose): Assembling compatible roots, values, and annotations from multiple trees
- [`cladefreq`](https://github.com/kfuku52/nwkit/wiki/nwkit-cladefreq): Summarizing clade frequencies across a tree collection
- [`consensus`](https://github.com/kfuku52/nwkit/wiki/nwkit-consensus): Generating a consensus tree or transferring consensus support to a reference tree
- [`contrast`](https://github.com/kfuku52/nwkit/wiki/nwkit-contrast): Calculating continuous-trait phylogenetic independent contrasts, with biological/technical replicates, batch adjustment, propagated sampling covariance, and reconciled gene-tree event mappings
- [`diff`](https://github.com/kfuku52/nwkit/wiki/nwkit-diff): Reporting interpretable clade, root, value, and annotation differences between trees
- [`dist`](https://github.com/kfuku52/nwkit/wiki/nwkit-dist): Comparing tree topology and branch lengths with multiple distance metrics
- [`draw`](https://github.com/kfuku52/nwkit/wiki/nwkit-draw): Drawing phylogenetic trees with Cartesian, polar, unrooted, spiral, or fractal geometry, annotation-aware spacing, and auditable layout reports
- [`drop`](https://github.com/kfuku52/nwkit/wiki/nwkit-drop): Removing node and branch information
- [`image`](https://github.com/kfuku52/nwkit/wiki/nwkit-image): Retrieving representative species images with license-aware filtering
- [`info`](https://github.com/kfuku52/nwkit/wiki/nwkit-info): Printing tree information
- [`intersection`](https://github.com/kfuku52/nwkit/wiki/nwkit-intersection): Dropping non-overlapping leaves/sequences between two trees or between a tree and an alignment
- [`label`](https://github.com/kfuku52/nwkit/wiki/nwkit-label): Adding unique node labels
- [`mark`](https://github.com/kfuku52/nwkit/wiki/nwkit-mark): Adding texts to node labels by identifying the targets with a leaf name regex
- [`mcmctree`](https://github.com/kfuku52/nwkit/wiki/nwkit-mcmctree): Preparing PAML MCMCtree calibrations and converting posterior node ages into pipeable dated NHX trees
- [`monophyly`](https://github.com/kfuku52/nwkit/wiki/nwkit-monophyly): Assessing whether species or trait-defined groups are monophyletic
- [`nhx2nwk`](https://github.com/kfuku52/nwkit/wiki/nwkit-nhx2nwk): Generating Newick from NHX
- [`nwk2table`](https://github.com/kfuku52/nwkit/wiki/nwkit-nwk2table): Converting a Newick tree into a parent-child table
- [`regress`](https://github.com/kfuku52/nwkit/wiki/nwkit-regress): Fitting conventional or reconciled Gaussian/multivariate PGLS and categorical, count, zero-inflated, positive, proportion, or censored phylogenetic GLMMs, with partial responses, biological replicates, gene-tree ensembles, latent-predictor measurement error, and automatic shape-parameter estimation
- [`printlabel`](https://github.com/kfuku52/nwkit/wiki/nwkit-printlabel): Searching and printing node labels
- [`prune`](https://github.com/kfuku52/nwkit/wiki/nwkit-prune): Pruning leaves
- [`rename`](https://github.com/kfuku52/nwkit/wiki/nwkit-rename): Renaming nodes using a TSV mapping or regular expression
- [`reconcile`](https://github.com/kfuku52/nwkit/wiki/nwkit-reconcile): Mapping rooted gene-tree nodes and events onto a rooted species tree
- [`rescale`](https://github.com/kfuku52/nwkit/wiki/nwkit-rescale): Rescale branch length with a given factor
- [`root`](https://github.com/kfuku52/nwkit/wiki/nwkit-root): Placing, transferring, or reconciliation-rooting the tree root
- [`rootcompare`](https://github.com/kfuku52/nwkit/wiki/nwkit-rootcompare): Comparing rooting methods in a TSV summary and a branch-marked PDF
- [`sanitize`](https://github.com/kfuku52/nwkit/wiki/nwkit-sanitize): Eliminating non-standard Newick flavors
- [`sample`](https://github.com/kfuku52/nwkit/wiki/nwkit-sample): Selecting a representative leaf subset by filters, ranks, and sampling method
- [`shuffle`](https://github.com/kfuku52/nwkit/wiki/nwkit-shuffle): Shuffling branches and/or labels
- [`skim`](https://github.com/kfuku52/nwkit/wiki/nwkit-skim): Sampling leaves from clades with shared traits
- [`subtree`](https://github.com/kfuku52/nwkit/wiki/nwkit-subtree): Generating a subtree Newick file
- [`table2nwk`](https://github.com/kfuku52/nwkit/wiki/nwkit-table2nwk): Converting a parent-child table into a Newick tree
- [`transfer`](https://github.com/kfuku52/nwkit/wiki/nwkit-transfer): Transferring information between trees
- [`validate`](https://github.com/kfuku52/nwkit/wiki/nwkit-validate): Validating one or more Newick trees and reporting structural issues


## Citation
There is no published paper on NWKIT itself, but we used and cited NWKIT in several papers including [Fukushima & Pollock (2023, Nat Ecol Evol 7: 155-170)](https://www.nature.com/articles/s41559-022-01932-7).

## Development

Create an isolated environment, install the development and optional image
dependencies with the reproducible tool constraints, then use the same check
runner as CI:

```sh
python -m venv .venv
. .venv/bin/activate
python -m pip install -U pip
python -m pip install -c constraints-dev.txt -e ".[dev,image]"
python tools/check.py quick
```

`quick` runs formatting, linting, incremental type checks, and tests excluding
the `slow` group. Use `full` for all tests and the
security, dependency, coverage, and maintainability gates; `dist` for
reproducible package validation; or `release` for the complete pre-release
suite. Focused checks, CI coverage, and benchmarks are described in
[DEVELOPMENT.md](https://github.com/kfuku52/nwkit/blob/master/DEVELOPMENT.md).

See [CHANGELOG.md](https://github.com/kfuku52/nwkit/blob/master/CHANGELOG.md) for changes and
[RELEASING.md](https://github.com/kfuku52/nwkit/blob/master/RELEASING.md) for the release checklist.

# Licensing
This program is MIT-licensed. See
[LICENSE](https://github.com/kfuku52/nwkit/blob/master/LICENSE) for details.

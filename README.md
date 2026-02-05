![](logo/logo_nwkit_large.png)

[![Bioconda](https://img.shields.io/conda/vn/bioconda/nwkit.svg)](https://anaconda.org/bioconda/nwkit)
[![Downloads](https://img.shields.io/conda/dn/bioconda/nwkit.svg)](https://anaconda.org/bioconda/nwkit)
[![Platforms](https://img.shields.io/conda/pn/bioconda/nwkit.svg)](https://anaconda.org/bioconda/nwkit)
[![License](https://img.shields.io/badge/License-BSD_3--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)

## Overview
**NWKIT** ([/njuːkit/](http://ipa-reader.xyz/?text=nju%CB%90kit&voice=Joanna)) is a toolkit for manipulating phylogenetic trees in the [Newick format](https://en.wikipedia.org/wiki/Newick_format). 

## Installation

The latest version of NWKIT is available from [Bioconda](https://anaconda.org/bioconda/nwkit). For users requiring a `conda` installation, please refer to [Miniforge](https://github.com/conda-forge/miniforge) for a lightweight conda environment.

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

## Subcommands
See [Wiki](https://github.com/kfuku52/nwkit/wiki) for usage.

- [`constrain`](https://github.com/kfuku52/nwkit/wiki/nwkit-constrain): Generating a species-tree-like Newick file for topological constraint
- [`dist`](https://github.com/kfuku52/nwkit/wiki/nwkit-dist): Calculating topological distance between two trees
- [`drop`](https://github.com/kfuku52/nwkit/wiki/nwkit-drop): Removing node and branch information
- [`info`](https://github.com/kfuku52/nwkit/wiki/nwkit-info): Printing tree information
- [`intersection`](https://github.com/kfuku52/nwkit/wiki/nwkit-intersection): Dropping non-overlapping leaves/sequences between two trees or between a tree and an alignment
- [`label`](https://github.com/kfuku52/nwkit/wiki/nwkit-label): Adding unique node labels
- [`mark`](https://github.com/kfuku52/nwkit/wiki/nwkit-mark): Adding texts to node labels by identifying the targets with a leaf name regex
- [`mcmctree`](https://github.com/kfuku52/nwkit/wiki/nwkit-mcmctree): Introducing divergence time constraints for PAML's mcmctree
- [`nhx2nwk`](https://github.com/kfuku52/nwkit/wiki/nwkit-nhx2nwk): Generating Newick from NHX
- [`printlabel`](https://github.com/kfuku52/nwkit/wiki/nwkit-printlabel): Searching and printing node labels
- [`prune`](https://github.com/kfuku52/nwkit/wiki/nwkit-prune): Pruning leaves
- [`rescale`](https://github.com/kfuku52/nwkit/wiki/nwkit-rescale): Rescale branch length with a given factor
- [`root`](https://github.com/kfuku52/nwkit/wiki/nwkit-root): Placing or transferring the tree root
- [`sanitize`](https://github.com/kfuku52/nwkit/wiki/nwkit-sanitize): Eliminating non-standard Newick flavors
- [`shuffle`](https://github.com/kfuku52/nwkit/wiki/nwkit-shuffle): Shuffling branches and/or labels
- [`skim`](https://github.com/kfuku52/nwkit/wiki/nwkit-skim): Sampling leaves from clades with shared traits
- [`subtree`](https://github.com/kfuku52/nwkit/wiki/nwkit-subtree): Generating a subtree Newick file
- [`transfer`](https://github.com/kfuku52/nwkit/wiki/nwkit-transfer): Transferring information between trees

## Citation
There is no published paper on NWKIT itself, but we used and cited NWKIT in several papers including [Fukushima & Pollock (2023, Nat Ecol Evol 7: 155-170)](https://www.nature.com/articles/s41559-022-01932-7).

# Licensing
This program is BSD-licensed (3 clause). See [LICENSE](LICENSE) for details.


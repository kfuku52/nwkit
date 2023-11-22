![](logo/logo_nwkit_large.png)

## Overview
**NWKIT** ([/njuːkit/](http://ipa-reader.xyz/?text=nju%CB%90kit&voice=Joanna)) is a toolkit for manipulating phylogenetic trees in the [Newick format](https://en.wikipedia.org/wiki/Newick_format). 

## Installation
```
# Installation with pip
pip install git+https://github.com/kfuku52/nwkit

# This should show complete options
nwkit -h 
```

## Subcommands
See [Wiki](https://github.com/kfuku52/nwkit/wiki) for usage.

- [`constrain`](https://github.com/kfuku52/nwkit/wiki/nwkit-constrain): Generating a species-tree-like Newick file for topological constraint
- [`dist`](https://github.com/kfuku52/nwkit/wiki/nwkit-dist): Calculating topological distance between two trees
- [`drop`](https://github.com/kfuku52/nwkit/wiki/nwkit-drop): Removing node and branch information
- [`intersection`](https://github.com/kfuku52/nwkit/wiki/nwkit-intersection): Dropping non-overlapping leaves/sequences between two trees or between a tree and an alignment
- [`mark`](https://github.com/kfuku52/nwkit/wiki/nwkit-mark): Adding texts to node labels by identifying the targets with a leaf name regex
- [`mcmctree`](https://github.com/kfuku52/nwkit/wiki/nwkit-mcmctree): Introducing divergence time constraints for PAML's mcmctree
- [`nhx2nwk`](https://github.com/kfuku52/nwkit/wiki/nwkit-nhx2nwk): Generating Newick from NHX
- [`printlabel`](https://github.com/kfuku52/nwkit/wiki/nwkit-printlabel): Searching and printing node labels
- [`rescale`](https://github.com/kfuku52/nwkit/wiki/nwkit-rescale): Rescale branch length with a given factor
- [`root`](https://github.com/kfuku52/nwkit/wiki/nwkit-root): Placing or transferring the tree root
- [`sanitize`](https://github.com/kfuku52/nwkit/wiki/nwkit-sanitize): Eliminating non-standard Newick flavors
- [`shuffle`](https://github.com/kfuku52/nwkit/wiki/nwkit-shuffle): Shuffling branches and/or labels
- [`subtree`](https://github.com/kfuku52/nwkit/wiki/nwkit-subtree): Generating a subtree Newick file
- [`transfer`](https://github.com/kfuku52/nwkit/wiki/nwkit-transfer): Transferring information between trees.

## Citation
There is no published paper on NWKIT itself, but we used and cited NWKIT in several papers including [Fukushima & Pollock (2023, Nat Ecol Evol 7: 155-170)](https://www.nature.com/articles/s41559-022-01932-7).

# Licensing
This program is BSD-licensed (3 clause). See [LICENSE](LICENSE) for details.


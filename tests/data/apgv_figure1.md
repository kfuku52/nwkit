# APG V backbone topology reference

Source: [APG V, Figure 1](https://doi.org/10.1111/jse.70096).
The production Newick is our manual transcription, not an author-supplied
tree file. On 2026-09-18, the article's supporting information listed only
the APG V survey DOCX; the linked TU Delft educational poster was a PDF.
The accompanying `apgv_figure1.json` records every parent and its immediate
children independently of the production Newick. The children follow the
figure's top-to-bottom order; the tests ignore sibling order.

The figure has 68 tips (66 orders plus Huaceae and Columelliaceae), 55 internal
nodes, eight three-way nodes and two four-way nodes. Internal identifiers in
the JSON are descriptive test labels, not additional taxonomic assertions.
This is a transcription of the figure's synopsis topology, not an assertion
that unresolved relationships have been resolved by a single analysis.

The complete rooted-clade comparison checks all 123 nodes, including the root
and tips. It rejects missing clades, extra resolutions of polytomies, duplicate
tips and unary nodes. End-to-end tests also check this topology after taxonomic
matching and Newick output, with and without gene collapsing. Their simulated
NCBI lineages place Huaceae within Oxalidales and Columelliaceae within
Bruniales to exercise the family-priority rule while retaining all four tips.

The 2026-09-18 audit corrected four regions in the initial transcription:

| Region | Topology depicted in Figure 1 |
| --- | --- |
| Commelinids | `(Arecales,((Zingiberales,Commelinales),Poales))` |
| Huaceae | `(Celastrales,Huaceae,Malpighiales)` (three-way split) |
| Lamiids after Icacinales | `((Vahliales,Solanales,(Boraginales,Gentianales)),Lamiales)` |
| Campanulids after Asterales | `((Columelliaceae,(Escalloniales,Paracryphiales),Dipsacales),Apiales)` |

The other parent-child relationships were checked against the figure as well.
To rerun the topology and output checks:

```sh
python tools/check.py test -- tests/test_constrain.py -k figure1
```

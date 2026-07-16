# Figure workflow

`paper/analysis/make_figures.py` generates every main figure from the command
inventory and raw result tables. Each figure is written here as PDF and SVG for
vector use and as a 300-ppi PNG for manuscript embedding.

The figures use a colorblind-safe blue/orange/green/purple palette, sans-serif
labels, and lowercase `a)` panel labels. Figure titles remain in the manuscript
legends rather than being duplicated inside the artwork. Selected tips in Figure
2 have italicized abbreviated species names and roman sequence accessions.
Alternative text is maintained in `manuscript.md` and embedded into the
generated Word document by `format_docx.py`.

Display allocation:

1. Figure 1: architecture, command families, and shared input/output semantics.
2. Figure 2: public PEPC worked example.
3. Figure 3: output checks, predefined input cases, and scaling.
4. Table 1: task-level command-line capability comparison.

The four displays collectively meet the target journal's limit of four main
figures and tables.

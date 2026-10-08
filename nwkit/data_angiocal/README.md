# Normalized AngioCal v1.0 records

`v1.0.tsv` contains all 238 records and imported metadata, including original
worksheet row numbers, from `Data2b_CalibrationList.xls` at
[the pinned upstream commit](https://github.com/eflowerproject/angiocal/tree/263e31ce22cbd9df6c7e8644406e859f43f59d70).
NWKIT normalized the imported fields to UTF-8 TSV without changing fossil ages,
identifiers, crown/stem assignments, quality scores or reference text.

Source XLS SHA-256:
`7e6e0b8ba8dfc8b58812a9c38cd2741cb3ca9f534bea096cf4b041279c3f1ee0`.
Normalized TSV SHA-256:
`4f121728b4959c41f9c2fda6117f5ab311177c0b176e6380ad7b853ed86cb5ba`.

The source dataset was compiled by Hervé Sauquet, Santiago Ramírez-Barahona and
Susana Magallón. Its [Zenodo record](https://doi.org/10.5281/zenodo.3828071)
describes it as **other-open**; that dataset license is distinct from NWKIT's MIT
software license. Cite the source dataset and its associated paper when using
these calibrations:

Ramírez-Barahona S, Sauquet H, Magallón S (2020). The delayed and geographically
heterogeneous diversification of flowering plant families. *Nature Ecology &
Evolution* 4, 1232–1238. <https://doi.org/10.1038/s41559-020-1241-3>.

At runtime NWKIT first verifies the original XLS checksum, then reads this
independently verified TSV. Reports and audits still identify the original XLS
path/URL, checksum and worksheet row. An XLS with different bytes uses the
optional `xlrd` reader rather than this version-specific copy.

To reproduce or verify this normalization, install the development or XLS extra
and run from the checkout root:

```sh
python tools/normalize_angiocal.py /path/to/Data2b_CalibrationList.xls --check
```

Omit `--check` to regenerate the TSV from that exact source. The source digest
is checked before parsing; changing the source or TSV requires updating both
digests and comparing all imported records independently.

# Audit and verify a partial topological constraint

Use `nwkit constrain --report matches.tsv` to trace each input label to its
backbone anchor. The NCBI, APG and user modes report the taxonomy query, its
origin, exact name/genus/taxid resolution, anchor and known rank, constrained
or free status, reason, final output label and backbone source. In
`--species-map-tsv` the optional `source` and `anchor_rank` columns record the
evidence for a supplied taxonomy query, e.g. a dated family-to-order mapping.
Supply source versions in this column when constructing your mapping. A free
tip can stay in the sequence alignment. An entirely unmatched input produces
the diagnostic report and fails to generate a tree. An unknown explicitly supplied taxid also fails the command. Ambiguous NCBI name matches remain free. VertLife retains its
dataset-specific report with DOI and checksum.

After IQ-TREE searches with `-g`, verify the result against the list of every
input sequence ID, including those left free.

```sh
nwkit check-constraint --infile ml.treefile --constraint constraint.nwk \
  --species-list alignment-tips.txt --outfile constraint-check.json
```

The command checks exact tip retention independently of the nontrivial unrooted
splits on the constrained tips. Free tips may attach anywhere. A star constraint
has no required split. Missing constrained tips make the split count unassessable.
A nonzero exit and JSON diagnostics identify mismatches. Each violated interval
is a zero-based inclusive slice of `constraint_tip_order` defining one side of
a missing split. This representation keeps even a large set of nested violations
linear in output size. Root placement, branch lengths and support are outside
this test.

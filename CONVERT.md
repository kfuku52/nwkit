# Tree format conversion

`nwkit convert` converts one tree without implicitly renaming its nodes or inferring new
ages or credible intervals. It accepts Newick, NHX, FigTree NEXUS and the
`Species tree for FigTree` block in MCMCtree's main output.

```sh
nwkit convert -i mcmctree.out --from mcmctree-output --to figtree --time-factor 1000 -o FigTree.tre
nwkit convert -i FigTree.tre --to nhx -o dated.nhx
nwkit convert -i dated.nhx --to newick --age-ci drop -o dated.nwk
nwkit label -i dated.nwk --target intnode --prefix s --start 1 --force yes -o labeled.nwk
nwkit validate -i dated.nhx --require-all-lengths yes --fail-on-issue yes
```

Input and output default to stdin/stdout (`-`). `--from` defaults to `auto` and
accepts `newick`, `nhx`, `figtree`, and `mcmctree-output`. `--to` defaults to `nhx`
and accepts `newick`, `nhx`, and `figtree`. FigTree output is a NEXUS document.
`--format` still selects the ETE parser for names versus support; it is not the
container format. `--input-rooted` follows the common NWKIT rooting contract.

`--time-factor` must be finite and positive (default `1`). It multiplies supplied
branch lengths and explicit `age`, `age_mean`, `age_median`, `age_ci_low`, and
`age_ci_high` attributes. FigTree `height`, `height_mean`, and `height_median`
are normalized to the corresponding age attributes. Names, numeric support,
probabilities, and unrelated annotations are not scaled. In particular, a node
named `'0.2,0.3'` remains a name, not an inferred interval.

FigTree `{low,high}` intervals become four NHX attributes: `age_ci_low`,
`age_ci_high`, `age_ci_kind` (`HPD` or `equal-tail`), and `age_ci_level` (for
example `0.95`). Their supplied method and level are retained. `--age-ci drop`
explicitly removes all four interval attributes; its default is `keep`.
Plain Newick output refuses to discard NHX properties implicitly. Use NHX or
FigTree to retain other attributes. Ordinary comments and uninterpreted FigTree
fields are retained as comments; no time semantics are inferred for them.

## Replacing NHX-to-Newick workflows

Use `--properties drop` to explicitly discard all normalized NHX and age
properties, including age intervals. The default `keep` continues to reject
lossy plain-Newick conversion. Ordinary comments, uninterpreted FigTree fields,
and explicit rooting declarations are retained. Invalid annotations are still
rejected before removal.

`--node-label PROPERTY` copies that property's original value into internal
node labels (including the root) before removing or scaling properties. It
replaces existing internal labels; nodes without the property keep their labels.
Tip names never change. Values are quoted as Newick names; numeric values are
labels, not support estimates. Read them with `--format 1` downstream when
necessary. Property copies remain in NHX/FigTree output unless explicitly dropped.
The input-only `--quoted-node-names no` check does not forbid safely quoted
generated labels. `nwkit_rooted` can be copied before its rooting declaration
is canonicalized, including when `--input-rooted` overrides that declaration.

NHX keys `name`, `dist`, and `support` are rejected, even with `--properties drop`:
ETE treats them as overrides of Newick fields, so retaining or removing them can
silently change names, lengths, or support. Put those values in Newick fields or
rename the attributes to non-reserved keys before conversion.

```sh
# Strip NHX properties explicitly.
nwkit convert -i input.nhx --to newick --properties drop -o output.nwk
# Equivalent purpose to nhx2nwk --node-label S, with explicit property removal.
nwkit convert -i input.nhx --to newick --node-label S --properties drop -o labeled.nwk
```

`nhx2nwk` remains available for existing scripts. `convert` retains its stricter
validation and container-selection rules, so it is not a byte-for-byte alias.

Multiple input trees require a one-based `--tree-index`. MCMCtree's known
paired plain and CI-annotated renderings of the same dated tree are recognized
as one result, and the annotated rendering is selected. Its topology/index
rendering is excluded. Distinct dated trees remain ambiguous. An explicit index
counts all original tree statements, including topology/index statements.
NEXUS TRANSLATE tables are rejected; provide trees with direct tip labels.
Malformed, duplicate-tip, non-finite and negative-length inputs fail before
publishing output. Named file output uses NWKIT's output transaction.

`validate --require-all-lengths yes` marks a tree with any missing non-root
branch length as invalid (`missing_branch_length`). The default is `no`.
This does not require ultrametricity: rounded posterior mean branch lengths may
be slightly non-ultrametric. Finite and non-negative length checks always apply.

`label --start INT` sets the first candidate label number (default `0`). Labels
are assigned in level order, skipping names already reserved by other nodes.

# Collect taxon images

`nwkit image` reads a Newick or NHX tree and writes representative images,
`manifest.tsv`, `unmatched.tsv`, and `ATTRIBUTION.md`. The manifest can be
passed to `nwkit draw --tip-image-manifest`.

For a large photo collection:

```sh
nwkit image -i species.nwk --out-dir species-images \
  --style photo --source wikimedia,gbif,inaturalist,ncbi \
  --first-available yes --resume yes --synonyms yes \
  --provider-interval wikimedia=0.5 --output-format jpg --max-edge 1600
```

With `--first-available yes`, NWKIT tries providers in order and stops after
saving `--max-per-species` suitable images. It tries later providers when
earlier media are unreachable or identified as fallback photos. Direct matches
are attempted before synonyms and broader fallbacks. The normal ranked mode
remains the default. Repeat `--provider-interval` to set the minimum spacing
between API calls for each provider across the workers in one run. iNaturalist
uses at least one second between calls.

## Prefer photos that show the organism

For the same taxonomic match, the default ranking places photos identified as
specimen sheets, product packaging, and similar material behind other images.
With `--style photo`, the ranked mode also prioritizes explicit descriptions
of the organism, such as flowering branches and whole plants. NWKIT uses
source metadata for this estimate. Descriptions with little detail keep their
ordinary rank. A flagged photo can still fill a match when better candidates
are unavailable. Look for `organism_photo` or `secondary_photo_fallback` in
the manifest's `selection_reason`. `--style silhouette` keeps its existing
ranking.

## Verify and reuse individual images

Pass `--verified-images-tsv verified.tsv` to pin one image to a specific
leaf. This also works for a photo cropped from a paper or a cultivar that a
provider cannot match reliably. For example:

```tsv
leaf_name	local_path	media_url	source_page_url	attribution	license_url	matched_name
Platanthera_zijinensis	photos/zijinensis.jpg		https://doi.org/10.1038/s41477-022-01127-9	Li et al. (2022), Fig. 1b, cropped	https://creativecommons.org/licenses/by/4.0/	Platanthera zijinensis
```

Each row references an input leaf exactly once and specifies one local file
(relative to the TSV) or an HTTPS media URL. `source_page_url` must identify
the original source. Supply `license_url`, `license_code`, or both. If the
license has a specific version, include its URL. Supply `attribution` for
any licensed image with an attribution requirement. The optional
`matched_name` and `matched_rank` record the identity verified by the
curator. The command applies the same license and media validation and
formatting to these files as to provider results. A pinned image that fails
to load is reported in `unmatched.tsv`.

## Match names below species rank

`--species-map-tsv` overrides a tree label's lookup name. The staged
workflow also preserves explicit qualifiers such as `var.`, `subsp.`, and
`cultivar` in tree labels when no lookup override is provided. If a species
image is used without confirming the infraspecific identity, the manifest
records that limitation. `--synonyms yes` tries up to two names that resolve
unambiguously to the same NCBI Taxonomy ID. A synonym match preserves the
original requested query in `species_name` and records the matching name
in `matched_name`.

`--wikimedia-depicts yes` adds a Commons structured-data search when a
unique Wikidata scientific-name item can be found. A file tagged as depicting
that taxon ranks ahead of a filename or description match within Wikimedia.
The text search remains available for taxa and files lacking those records.

## Resume an interrupted collection

Run with `--resume yes` from the first invocation. Completed taxa are
checkpointed individually. Reuse the same output directory with compatible
selection and image settings to reuse saved files after checking their SHA-256
hashes. New tips can be added to the tree without repeating completed queries.
Missing or damaged files and unfinished taxa are retried. Use
`--refresh-cache yes` to repeat selection for all taxa.
The photo ranking update reselects earlier automatic photo results once, using
fresh provider metadata. Existing pinned results can resume from their
validated checkpoints.

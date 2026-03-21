import csv
import re
from dataclasses import dataclass


DEFAULT_SPECIES_REGEX = r'^([^_]+_[^_]+)(?:_|$)'
DEFAULT_SPECIES_PARSER = 'legacy'
SUPPORTED_SPECIES_PARSERS = ('legacy', 'taxonomic')
TAXONOMIC_SPECIES_SUFFIX_TOKENS = {'cf', 'aff', 'nr'}
GENUS_ONLY_PLACEHOLDER_TOKENS = {'sp', 'sp.', 'spp', 'spp.'}
TAXONOMIC_RANK_ALIASES = {
    'subsp': 'subsp',
    'ssp': 'subsp',
    'subspecies': 'subsp',
    'var': 'var',
    'variety': 'var',
    'forma': 'forma',
    'form': 'forma',
    'f': 'forma',
    'strain': 'strain',
    'substrain': 'substrain',
    'serovar': 'serovar',
    'serotype': 'serotype',
    'serogroup': 'serogroup',
    'pathovar': 'pathovar',
    'pv': 'pathovar',
    'biovar': 'biovar',
    'biotype': 'biotype',
    'chemovar': 'chemovar',
    'morphovar': 'morphovar',
    'cultivar': 'cultivar',
    'cv': 'cultivar',
    'isolate': 'isolate',
    'group': 'group',
    'subgroup': 'subgroup',
    'complex': 'complex',
    'clade': 'clade',
    'lineage': 'lineage',
    'section': 'section',
    'series': 'series',
    'ecotype': 'ecotype',
    'breed': 'breed',
}


@dataclass(frozen=True)
class ParsedSpecies:
    species_label: str = None
    taxonomy_query: str = None


def compile_species_regex(species_regex=None):
    regex_text = DEFAULT_SPECIES_REGEX if (species_regex is None) else str(species_regex).strip()
    if regex_text == '':
        return None
    try:
        return re.compile(regex_text)
    except re.error as exc:
        raise ValueError('--species_regex is not a valid regular expression: {}'.format(exc))


def _normalize_species_label(text):
    if text in [None, '']:
        return None
    normalized = re.sub(r'\s+', '_', str(text).strip())
    normalized = re.sub(r'_+', '_', normalized).strip('_')
    return normalized or None


def _normalize_taxonomy_query(text):
    if text in [None, '']:
        return None
    normalized = str(text).strip().replace('_', ' ')
    normalized = re.sub(r'\s+', ' ', normalized)
    return normalized or None


def _canonical_taxonomic_token(token):
    cleaned = str(token or '').strip()
    lowered = cleaned.lower()
    if lowered in GENUS_ONLY_PLACEHOLDER_TOKENS:
        return 'sp'
    if lowered in TAXONOMIC_SPECIES_SUFFIX_TOKENS:
        return re.sub(r'\.$', '', lowered)
    if lowered in TAXONOMIC_RANK_ALIASES:
        return TAXONOMIC_RANK_ALIASES[lowered]
    return cleaned


def _extract_species_label_with_regex(label, species_pattern):
    if label is None or species_pattern is None:
        return None
    label_text = str(label)
    match = species_pattern.search(label_text)
    if match is None:
        return None
    token = None
    if match.lastindex is not None:
        captured_tokens = list()
        for i in range(1, int(match.lastindex) + 1):
            candidate = match.group(i)
            if candidate is None:
                continue
            candidate = str(candidate).strip()
            if candidate != '':
                captured_tokens.append(candidate)
        if len(captured_tokens) > 0:
            token = '_'.join(captured_tokens)
    if token is None:
        token = str(match.group(0)).strip()
    if token == '':
        return None
    token = re.sub(r'[_\s]+$', '', token)
    return _normalize_species_label(token)


def read_species_map_tsv(species_map_tsv):
    if species_map_tsv in ['', None]:
        return dict()
    with open(species_map_tsv, newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        fieldnames = reader.fieldnames or list()
        if 'leaf_name' not in fieldnames:
            raise ValueError('--species_map_tsv must contain a "leaf_name" column.')
        if ('species_label' not in fieldnames) and ('taxonomy_query' not in fieldnames):
            raise ValueError('--species_map_tsv must contain "species_label" and/or "taxonomy_query" columns.')
        overrides = dict()
        for row in reader:
            leaf_name = str(row.get('leaf_name', '')).strip()
            if leaf_name == '':
                raise ValueError('--species_map_tsv contains an empty "leaf_name" value.')
            if leaf_name in overrides:
                raise ValueError('Duplicate "leaf_name" entries are not supported in --species_map_tsv: {}'.format(leaf_name))
            species_label = _normalize_species_label(row.get('species_label'))
            taxonomy_query = _normalize_taxonomy_query(row.get('taxonomy_query'))
            if (species_label is None) and (taxonomy_query is None):
                raise ValueError(
                    '--species_map_tsv row for "{}" must define "species_label" and/or "taxonomy_query".'.format(
                        leaf_name
                    )
                )
            overrides[leaf_name] = ParsedSpecies(
                species_label=species_label,
                taxonomy_query=taxonomy_query,
            )
    return overrides


class SpeciesParser:
    def __init__(self, preset=None, species_regex=None, species_map_tsv=None):
        self.preset = str(preset or DEFAULT_SPECIES_PARSER).strip().lower()
        if self.preset not in SUPPORTED_SPECIES_PARSERS:
            raise ValueError(
                '--species_parser must be one of: {}'.format(', '.join(SUPPORTED_SPECIES_PARSERS))
            )
        self.species_regex = species_regex
        self.species_pattern = compile_species_regex(species_regex=species_regex)
        self.map_overrides = read_species_map_tsv(species_map_tsv)

    def _parse_legacy(self, label):
        species_label = _extract_species_label_with_regex(label=label, species_pattern=self.species_pattern)
        if species_label is None:
            return ParsedSpecies()
        return ParsedSpecies(
            species_label=species_label,
            taxonomy_query=_normalize_taxonomy_query(species_label),
        )

    def _parse_taxonomic(self, label):
        normalized_label = _normalize_species_label(label)
        if normalized_label is None:
            return ParsedSpecies()
        tokens = [_canonical_taxonomic_token(token) for token in normalized_label.split('_') if token != '']
        if len(tokens) < 2:
            return ParsedSpecies()
        genus_token = tokens[0]
        second_token = tokens[1]
        second_lower = second_token.lower()
        if second_lower in GENUS_ONLY_PLACEHOLDER_TOKENS:
            species_tokens = tokens[:3] if len(tokens) >= 3 else tokens[:2]
            return ParsedSpecies(
                species_label='_'.join(species_tokens),
                taxonomy_query=genus_token,
            )
        if (len(tokens) >= 3) and (second_lower in TAXONOMIC_SPECIES_SUFFIX_TOKENS):
            return ParsedSpecies(
                species_label='_'.join(tokens[:3]),
                taxonomy_query='{} {}'.format(genus_token, tokens[2]),
            )
        if (len(tokens) >= 3) and (tokens[2].lower() in TAXONOMIC_SPECIES_SUFFIX_TOKENS):
            return ParsedSpecies(
                species_label='_'.join(tokens[:3]),
                taxonomy_query='{} {}'.format(genus_token, second_token),
            )
        if (len(tokens) >= 4) and (tokens[2].lower() in set(TAXONOMIC_RANK_ALIASES.values())):
            return ParsedSpecies(
                species_label='_'.join(tokens[:4]),
                taxonomy_query='{} {}'.format(genus_token, second_token),
            )
        return ParsedSpecies(
            species_label='{}_{}'.format(genus_token, second_token),
            taxonomy_query='{} {}'.format(genus_token, second_token),
        )

    def _apply_overrides(self, label, parsed_species):
        override = self.map_overrides.get(str(label), None)
        if override is None:
            return parsed_species
        species_label = override.species_label if override.species_label is not None else parsed_species.species_label
        taxonomy_query = override.taxonomy_query if override.taxonomy_query is not None else parsed_species.taxonomy_query
        if (species_label is None) and (taxonomy_query is not None):
            species_label = _normalize_species_label(taxonomy_query)
        if (taxonomy_query is None) and (species_label is not None):
            taxonomy_query = _normalize_taxonomy_query(species_label)
        return ParsedSpecies(species_label=species_label, taxonomy_query=taxonomy_query)

    def parse(self, label):
        if self.preset == 'legacy':
            parsed_species = self._parse_legacy(label)
        elif self.preset == 'taxonomic':
            parsed_species = self._parse_taxonomic(label)
        else:
            raise ValueError(
                '--species_parser must be one of: {}'.format(', '.join(SUPPORTED_SPECIES_PARSERS))
            )
        return self._apply_overrides(label=label, parsed_species=parsed_species)


def get_species_parser(args=None, species_parser=None, species_regex=None, species_map_tsv=None):
    if (args is not None) and (species_parser is None) and (species_regex is None) and (species_map_tsv is None):
        cache_key = (
            getattr(args, 'species_parser', DEFAULT_SPECIES_PARSER),
            getattr(args, 'species_regex', DEFAULT_SPECIES_REGEX),
            getattr(args, 'species_map_tsv', None),
        )
        cached_key = getattr(args, '_nwkit_species_parser_cache_key', None)
        cached_parser = getattr(args, '_nwkit_species_parser_cache', None)
        if (cached_key == cache_key) and (cached_parser is not None):
            return cached_parser
        parser = SpeciesParser(
            preset=cache_key[0],
            species_regex=cache_key[1],
            species_map_tsv=cache_key[2],
        )
        setattr(args, '_nwkit_species_parser_cache_key', cache_key)
        setattr(args, '_nwkit_species_parser_cache', parser)
        return parser
    return SpeciesParser(
        preset=species_parser if (species_parser is not None) else getattr(args, 'species_parser', DEFAULT_SPECIES_PARSER),
        species_regex=species_regex if (species_regex is not None) else getattr(args, 'species_regex', DEFAULT_SPECIES_REGEX),
        species_map_tsv=species_map_tsv if (species_map_tsv is not None) else getattr(args, 'species_map_tsv', None),
    )


def extract_parsed_species(label, args=None, species_parser=None, species_regex=None, species_map_tsv=None):
    parser = get_species_parser(
        args=args,
        species_parser=species_parser,
        species_regex=species_regex,
        species_map_tsv=species_map_tsv,
    )
    return parser.parse(label)

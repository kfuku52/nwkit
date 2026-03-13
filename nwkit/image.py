import csv
import json
import functools
import glob
import hashlib
import html
import io
import math
import os
import re
import shutil
import sys
import tarfile
import tempfile
import threading
from concurrent.futures import ThreadPoolExecutor, as_completed
from collections import defaultdict
from urllib.parse import urlparse

import requests

from nwkit.util import (
    acquire_exclusive_lock,
    get_ete_ncbitaxa,
    label2sciname,
    read_tree,
    resolve_download_dir,
    validate_unique_named_leaves,
)


PHYLIPIC_API_ROOT = 'https://api.phylopic.org'
INATURALIST_API_ROOT = 'https://api.inaturalist.org/v1'
WIKIMEDIA_API_ROOT = 'https://commons.wikimedia.org/w/api.php'
GBIF_API_ROOT = 'https://api.gbif.org/v1'
EOL_API_ROOT = 'https://eol.org/api'
IDIGBIO_API_ROOT = 'https://search.idigbio.org/v2'
BIOICONS_GITHUB_API_ROOT = 'https://api.github.com/repos/duerrsimon/bioicons/git/trees/main'
BIOICONS_MEDIA_ROOT = 'https://bioicons.com/icons'
OPENVERSE_API_ROOT = 'https://api.openverse.org/v1'
NCBI_NEWTAXDUMP_URL = 'https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz'
REQUEST_TIMEOUT = (10, 60)
DEFAULT_LOOKUP_WORKERS = 4
DEFAULT_DOWNLOAD_WORKERS = 4
IMAGE_QUERY_CACHE_VERSION = 1
LOOKUP_FALLBACK_BUFFER = 2
BIOICONS_CATALOG_CACHE_VERSION = 1
SEMANTIC_MASK_MAX_DIM = 256
SEMANTIC_ALPHA_THRESHOLD = 16
SEMANTIC_DIFF_THRESHOLD_FLOOR = 12
HTTP_HEADERS = {
    'Accept': 'application/json',
    'User-Agent': 'nwkit-image/0.22.9 (+https://github.com/kfuku52/nwkit)',
}
SUPPORTED_SOURCES = ('phylopic', 'bioicons', 'inaturalist', 'wikimedia', 'gbif', 'eol', 'idigbio', 'openverse', 'ncbi')
DEFAULT_SOURCES = {
    'auto': ['phylopic', 'bioicons', 'inaturalist', 'wikimedia', 'gbif', 'eol', 'idigbio', 'openverse', 'ncbi'],
    'photo': ['inaturalist', 'wikimedia', 'gbif', 'eol', 'idigbio', 'openverse', 'ncbi'],
    'silhouette': ['phylopic', 'bioicons', 'wikimedia'],
}
LICENSE_LEVELS = {
    'public-domain': 0,
    'cc-by': 1,
    'cc-by-sa': 2,
    'cc-by-nc': 3,
    'cc-by-nc-sa': 4,
    'mit': 1,
    'bsd': 1,
}
LICENSE_OPENNESS = {
    'public-domain': 70,
    'mit': 65,
    'bsd': 63,
    'cc-by': 60,
    'cc-by-sa': 50,
    'cc-by-nd': 45,
    'cc-by-nc': 40,
    'cc-by-nc-sa': 30,
    'cc-by-nc-nd': 25,
    'unknown': 0,
    'all-rights-reserved': -50,
}
FILENAME_SANITIZE_PATTERN = re.compile(r'[^A-Za-z0-9._-]+')
BIOICONS_DESCRIPTOR_TOKENS = {
    'adult', 'blackeyes', 'blue', 'brown', 'chunky', 'cyan', 'darkgray', 'early',
    'embryo', 'embryoearly', 'embryolate', 'fat', 'gender', 'gray', 'green', 'head',
    'juvenile', 'late', 'new', 'orange', 'pink', 'purple', 'redeyes', 'small',
    'smiling', 'test', 'thin', 'white', 'yellow',
}
BIOICONS_SPECIES_ALIASES = {
    'anopheles gambiae': ('anopheles', 'mosquito'),
    'arabidopsis thaliana': ('arabidopsis_thaliana', 'arabidopsis'),
    'caenorhabditis elegans': ('celegans', 'c_elegans', 'nematode'),
    'danio rerio': ('zebrafish',),
    'drosophila melanogaster': ('drosophila', 'fruit_fly'),
    'escherichia coli': ('e_coli', 'coli', 'bacteria'),
    'macaca mulatta': ('rhesus_monkey', 'macaque', 'monkey'),
    'mus musculus': ('mouse',),
    'rattus norvegicus': ('rat',),
    'saccharomyces cerevisiae': ('budding_yeast', 'yeast'),
    'schizosaccharomyces pombe': ('fission_yeast', 'pombe'),
    'xenopus laevis': ('xenopus_laevis', 'xenopus'),
}
MIME_TYPE_TO_EXTENSION = {
    'image/gif': '.gif',
    'image/jpeg': '.jpg',
    'image/jpg': '.jpg',
    'image/png': '.png',
    'image/svg+xml': '.svg',
    'image/tiff': '.tif',
    'image/webp': '.webp',
}
PILLOW_INSTALL_HINT = 'Install optional image-processing dependencies with: pip install "nwkit[image]"'
RASTER_OUTPUT_EXTENSIONS = {
    '.gif': 'GIF',
    '.jpg': 'JPEG',
    '.jpeg': 'JPEG',
    '.png': 'PNG',
    '.tif': 'TIFF',
    '.tiff': 'TIFF',
    '.webp': 'WEBP',
}
BIOICONS_CATALOG_MEMORY_CACHE = dict()
BIOICONS_CATALOG_MEMORY_CACHE_LOCK = threading.Lock()


def _stderr(message):
    sys.stderr.write(str(message).rstrip() + '\n')


def _close_ncbi_db(ncbi):
    if ncbi is None:
        return
    if hasattr(ncbi, 'close') and callable(ncbi.close):
        ncbi.close()
        return
    db = getattr(ncbi, 'db', None)
    if db is not None:
        try:
            db.close()
        except Exception:
            pass


def normalize_species_name(name):
    if name is None:
        return None
    normalized = str(name).strip().replace('_', ' ')
    normalized = re.sub(r'\s+', ' ', normalized)
    return normalized if normalized != '' else None


def sanitize_filename_component(value):
    normalized = FILENAME_SANITIZE_PATTERN.sub('_', str(value).strip())
    normalized = normalized.strip('._')
    return normalized or 'item'


def normalize_license_code(raw_code=None, raw_url=None, attribution=None):
    if raw_code is not None:
        code = str(raw_code).strip().lower()
        if code in ('', 'none', 'null', 'nan'):
            code = None
        elif code in ('cc0', 'cc-0', 'pd', 'pdm', 'public-domain', 'public domain'):
            return 'public-domain'
        elif code in ('mit', 'mit license'):
            return 'mit'
        elif code in ('bsd', 'bsd license', 'bsd-2-clause', 'bsd-3-clause'):
            return 'bsd'
        elif code in ('by', 'cc-by', 'cc_by'):
            return 'cc-by'
        elif code in ('by-sa', 'cc-by-sa', 'cc_by_sa'):
            return 'cc-by-sa'
        elif code in ('by-nc', 'cc-by-nc', 'cc_by_nc'):
            return 'cc-by-nc'
        elif code in ('by-nc-sa', 'cc-by-nc-sa', 'cc_by_nc_sa'):
            return 'cc-by-nc-sa'
        elif code in ('by-nd', 'cc-by-nd', 'cc_by_nd'):
            return 'cc-by-nd'
        elif code in ('by-nc-nd', 'cc-by-nc-nd', 'cc_by_nc_nd'):
            return 'cc-by-nc-nd'
        elif code in ('all-rights-reserved', 'arr', 'all rights reserved'):
            return 'all-rights-reserved'
        elif ('public domain' in code) or ('cc0' in code) or ('pdm' == code):
            return 'public-domain'
        elif 'mit' in code:
            return 'mit'
        elif 'bsd' in code:
            return 'bsd'
        elif ('cc by-nc-nd' in code) or ('cc-by-nc-nd' in code):
            return 'cc-by-nc-nd'
        elif ('cc by-nd' in code) or ('cc-by-nd' in code):
            return 'cc-by-nd'
        elif ('cc by-nc-sa' in code) or ('cc-by-nc-sa' in code):
            return 'cc-by-nc-sa'
        elif ('cc by-nc' in code) or ('cc-by-nc' in code):
            return 'cc-by-nc'
        elif ('cc by-sa' in code) or ('cc-by-sa' in code):
            return 'cc-by-sa'
        elif ('cc by' in code) or ('cc-by' in code):
            return 'cc-by'

    if raw_url:
        url = str(raw_url).strip().lower()
        if any(token in url for token in ('publicdomain/zero', 'publicdomain/mark', '/zero/1.0', '/mark/1.0')):
            return 'public-domain'
        if '/licenses/by-nc-nd/' in url:
            return 'cc-by-nc-nd'
        if '/licenses/by-nd/' in url:
            return 'cc-by-nd'
        if '/licenses/by-nc-sa/' in url:
            return 'cc-by-nc-sa'
        if '/licenses/by-nc/' in url:
            return 'cc-by-nc'
        if '/licenses/by-sa/' in url:
            return 'cc-by-sa'
        if '/licenses/by/' in url:
            return 'cc-by'

    if attribution and ('all rights reserved' in str(attribution).lower()):
        return 'all-rights-reserved'
    return 'unknown'


def canonical_license_url(license_code):
    mapping = {
        'public-domain': 'https://creativecommons.org/publicdomain/zero/1.0/',
        'mit': 'https://opensource.org/licenses/MIT',
        'bsd': 'https://opensource.org/licenses/BSD-3-Clause',
        'cc-by': 'https://creativecommons.org/licenses/by/4.0/',
        'cc-by-sa': 'https://creativecommons.org/licenses/by-sa/4.0/',
        'cc-by-nc': 'https://creativecommons.org/licenses/by-nc/4.0/',
        'cc-by-nc-sa': 'https://creativecommons.org/licenses/by-nc-sa/4.0/',
        'cc-by-nd': 'https://creativecommons.org/licenses/by-nd/4.0/',
        'cc-by-nc-nd': 'https://creativecommons.org/licenses/by-nc-nd/4.0/',
    }
    return mapping.get(license_code, '')


def license_allowed(license_code, license_max='any', allow_nd=False):
    if license_code in ('unknown', 'all-rights-reserved', None):
        return False

    if license_code.endswith('-nd') and (not allow_nd):
        return False

    if license_max == 'any':
        return True

    if license_code in ('mit', 'bsd'):
        return LICENSE_LEVELS[license_code] <= LICENSE_LEVELS[license_max]

    if license_code in LICENSE_LEVELS:
        return LICENSE_LEVELS[license_code] <= LICENSE_LEVELS[license_max]

    if license_code == 'cc-by-nd':
        return allow_nd and (LICENSE_LEVELS['cc-by'] <= LICENSE_LEVELS[license_max])

    if license_code == 'cc-by-nc-nd':
        return allow_nd and (LICENSE_LEVELS['cc-by-nc'] <= LICENSE_LEVELS[license_max])

    return False


def parse_sources(style, source_arg):
    if source_arg in (None, ''):
        return list(DEFAULT_SOURCES[style])
    sources = [s.strip().lower() for s in str(source_arg).split(',') if s.strip() != '']
    if len(sources) == 0:
        return list(DEFAULT_SOURCES[style])
    unsupported = [s for s in sources if s not in SUPPORTED_SOURCES]
    if unsupported:
        raise ValueError(
            'Unsupported --source value(s) in the current implementation: {}. '
            'Supported sources are: {}.'.format(
                ', '.join(sorted(set(unsupported))),
                ', '.join(SUPPORTED_SOURCES),
            )
        )
    return sources


def read_name_tsv(path):
    with open(path, newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        if reader.fieldnames is None:
            raise ValueError('--name_tsv is empty.')
        required = {'leaf_name', 'species_name'}
        missing = required.difference(reader.fieldnames)
        if missing:
            raise ValueError('--name_tsv must contain "leaf_name" and "species_name" columns.')
        mapping = dict()
        for row in reader:
            leaf_name = str(row.get('leaf_name', '')).strip()
            species_name = normalize_species_name(row.get('species_name'))
            if leaf_name == '':
                raise ValueError('--name_tsv contains an empty "leaf_name" value.')
            if leaf_name in mapping:
                raise ValueError('Duplicate values in the "leaf_name" column of --name_tsv are not supported.')
            if species_name is None:
                raise ValueError('--name_tsv contains an empty "species_name" value.')
            mapping[leaf_name] = species_name
    if len(mapping) == 0:
        raise ValueError('--name_tsv is empty.')
    return mapping


def extract_species_mapping(tree, name_mapping=None):
    validate_unique_named_leaves(tree, '--infile', context=' for nwkit image')
    name_mapping = name_mapping or dict()
    leaf_names = set(tree.leaf_names())
    unknown_names = set(name_mapping.keys()).difference(leaf_names)
    if unknown_names:
        raise ValueError(
            '--name_tsv contains leaf names not found in --infile: {}'.format(
                ', '.join(sorted(unknown_names)),
            )
        )

    leaf_to_species = dict()
    unmatched_rows = list()
    for leaf in tree.leaves():
        species_name = name_mapping.get(leaf.name)
        if species_name is None:
            species_name = normalize_species_name(label2sciname(leaf.name, out_delim=' '))
        if species_name is None:
            unmatched_rows.append({
                'leaf_name': leaf.name,
                'species_name': '',
                'reason': 'unparsable leaf label',
                'details': "Expected the 'GENUS_SPECIES[_...]' convention or a matching --name_tsv entry.",
            })
            continue
        leaf_to_species[leaf.name] = species_name
    return leaf_to_species, unmatched_rows


def get_taxonomic_queries(species_name, fallback_rank='none', ncbi=None):
    queries = [('species', species_name)]
    parts = species_name.split()
    genus_name = parts[0] if len(parts) >= 1 else None
    family_name = None

    if fallback_rank in ('genus', 'family') and genus_name is not None:
        queries.append(('genus', genus_name))

    if fallback_rank == 'family' and ncbi is not None:
        lookup_names = [species_name]
        if genus_name is not None and genus_name != species_name:
            lookup_names.append(genus_name)
        name_to_taxid = ncbi.get_name_translator(lookup_names)
        taxid = None
        for lookup_name in lookup_names:
            taxids = name_to_taxid.get(lookup_name)
            if taxids:
                taxid = int(taxids[0])
                break
        if taxid is not None:
            lineage = ncbi.get_lineage(taxid)
            rank_by_taxid = ncbi.get_rank(lineage)
            name_by_taxid = ncbi.get_taxid_translator(lineage)
            for lineage_taxid in lineage:
                if rank_by_taxid.get(lineage_taxid) == 'family':
                    family_name = name_by_taxid.get(lineage_taxid)
                    break
    if family_name:
        queries.append(('family', family_name))

    deduped = list()
    seen = set()
    for matched_rank, query_name in queries:
        key = (matched_rank, normalize_species_name(query_name))
        if key in seen or key[1] is None:
            continue
        seen.add(key)
        deduped.append((matched_rank, key[1]))
    return deduped


def parse_size(text):
    if not text:
        return None, None
    parts = str(text).lower().split('x')
    if len(parts) != 2:
        return None, None
    try:
        width = int(float(parts[0]))
        height = int(float(parts[1]))
    except ValueError:
        return None, None
    return width, height


def strip_html_markup(value):
    if value in (None, ''):
        return ''
    text = re.sub(r'<[^>]+>', ' ', str(value))
    return re.sub(r'\s+', ' ', html.unescape(text)).strip()


def tokenize_search_terms(value):
    return re.findall(r'[a-z0-9]+', str(value or '').lower())


def canonicalize_search_terms(value):
    return ''.join(tokenize_search_terms(value))


def bioicons_display_author(author_slug):
    display = str(author_slug or '').replace('_', ' ').replace('-', ' ')
    return re.sub(r'\s+', ' ', display).strip()


def bioicons_match_quality(icon_name, query_name, matched_rank):
    icon_tokens = tokenize_search_terms(icon_name)
    query_tokens = tokenize_search_terms(query_name)
    if not icon_tokens:
        return 0

    best = 0
    if query_tokens:
        if icon_tokens == query_tokens:
            best = 90 if matched_rank == 'species' else 80
        elif icon_tokens[:len(query_tokens)] == query_tokens:
            best = max(best, 70 if matched_rank == 'species' else 60)
        elif (len(query_tokens) > 1) and all(token in icon_tokens for token in query_tokens):
            best = max(best, 55)

    scientific_aliases = BIOICONS_SPECIES_ALIASES.get(str(query_name or '').strip().lower(), ())
    for alias in scientific_aliases:
        alias_tokens = tokenize_search_terms(alias)
        if not alias_tokens:
            continue
        if icon_tokens == alias_tokens:
            best = max(best, 80)
        elif icon_tokens[:len(alias_tokens)] == alias_tokens:
            best = max(best, 60)
        elif all(token in icon_tokens for token in alias_tokens):
            best = max(best, 45)

    if best == 0:
        return 0

    descriptor_penalty = 0
    for token in icon_tokens[1:]:
        if token in BIOICONS_DESCRIPTOR_TOKENS:
            descriptor_penalty += 8
    return max(1, best - min(descriptor_penalty, 30))


def bioicons_index_keys_for_name(name):
    tokens = tokenize_search_terms(name)
    if not tokens:
        return set()
    keys = {
        canonicalize_search_terms(name),
        canonicalize_search_terms(' '.join(tokens)),
    }
    non_descriptor_tokens = [token for token in tokens if token not in BIOICONS_DESCRIPTOR_TOKENS]
    if non_descriptor_tokens:
        keys.add(canonicalize_search_terms(' '.join(non_descriptor_tokens)))
        keys.add(non_descriptor_tokens[0])
    keys.add(tokens[0])
    return {key for key in keys if key not in ('', None)}


def wikimedia_page_mentions_query(page, query_name):
    image_info = (page.get('imageinfo') or [{}])[0]
    metadata = image_info.get('extmetadata') or {}
    title_text = str(page.get('title', '')).replace('File:', ' ')
    object_name = strip_html_markup(metadata.get('ObjectName', {}).get('value', ''))
    description = strip_html_markup(metadata.get('ImageDescription', {}).get('value', ''))
    combined_text = ' '.join([title_text, object_name, description]).lower()
    return normalize_species_name(query_name).lower() in combined_text


def search_text_mentions_query(text_fragments, query_name):
    query_tokens = tokenize_search_terms(query_name)
    if not query_tokens:
        return False
    combined_tokens = tokenize_search_terms(' '.join([str(fragment or '') for fragment in text_fragments]))
    if not combined_tokens:
        return False
    combined_set = set(combined_tokens)
    return all(token in combined_set for token in query_tokens)


def infer_extension(url, default_ext='.bin'):
    path = urlparse(url).path
    _, ext = os.path.splitext(path)
    return ext.lower() if ext != '' else default_ext


def replace_extension(path, ext):
    root, _ = os.path.splitext(path)
    return root + ext


def infer_extension_from_content_type(content_type, default_ext='.bin'):
    normalized = str(content_type or '').split(';', 1)[0].strip().lower()
    return MIME_TYPE_TO_EXTENSION.get(normalized, default_ext)


def infer_extension_from_bytes_prefix(prefix, default_ext='.bin'):
    if prefix.startswith(b'\xff\xd8\xff'):
        return '.jpg'
    if prefix.startswith(b'\x89PNG\r\n\x1a\n'):
        return '.png'
    if prefix.startswith((b'GIF87a', b'GIF89a')):
        return '.gif'
    if prefix.startswith((b'II*\x00', b'MM\x00*')):
        return '.tif'
    if len(prefix) >= 12 and prefix[:4] == b'RIFF' and prefix[8:12] == b'WEBP':
        return '.webp'
    if prefix.lstrip().startswith(b'<svg') or prefix.lstrip().startswith(b'<?xml'):
        return '.svg'
    return default_ext


def infer_extension_from_response(response, media_url, first_chunk=b'', default_ext='.bin'):
    response_headers = getattr(response, 'headers', {}) or {}
    response_url = getattr(response, 'url', media_url)
    inferred_ext = infer_extension_from_content_type(response_headers.get('Content-Type'), default_ext='')
    if inferred_ext != '':
        return inferred_ext
    inferred_ext = infer_extension(response_url, default_ext='')
    if inferred_ext != '':
        return inferred_ext
    inferred_ext = infer_extension_from_bytes_prefix(first_chunk, default_ext='')
    if inferred_ext != '':
        return inferred_ext
    return infer_extension(media_url, default_ext=default_ext)


def find_existing_media_variant(path):
    if os.path.exists(path):
        return path
    base, _ = os.path.splitext(path)
    candidates = sorted(
        candidate for candidate in glob.glob(base + '.*')
        if (not candidate.endswith('.tmp')) and os.path.isfile(candidate)
    )
    return candidates[0] if candidates else None


def normalize_existing_media_path(path, default_ext='.bin'):
    if path is None:
        return None
    existing_path = find_existing_media_variant(path)
    if existing_path is None:
        return path
    _, ext = os.path.splitext(existing_path)
    if ext.lower() != default_ext:
        return existing_path
    with open(existing_path, 'rb') as handle:
        prefix = handle.read(256)
    inferred_ext = infer_extension_from_bytes_prefix(prefix, default_ext=default_ext)
    if inferred_ext == default_ext:
        return existing_path
    normalized_path = replace_extension(existing_path, inferred_ext)
    if normalized_path == existing_path:
        return existing_path
    os.replace(existing_path, normalized_path)
    return normalized_path


def image_postprocessing_requested(args):
    return any((
        getattr(args, 'output_format', 'original') != 'original',
        getattr(args, 'max_edge', None) not in (None, 0),
        getattr(args, 'canvas', 'none') != 'none',
        getattr(args, 'trim', 'off') != 'off',
        getattr(args, 'trim_shape', 'bbox') != 'bbox',
    ))


def load_pillow_modules():
    try:
        from PIL import Image, ImageChops, ImageOps
    except ImportError as exc:
        raise RuntimeError(
            'Image post-processing requires the optional Pillow dependency. {}'.format(PILLOW_INSTALL_HINT)
        ) from exc
    return Image, ImageChops, ImageOps


def load_cairosvg_module():
    try:
        import cairosvg
    except ImportError as exc:
        raise RuntimeError(
            'SVG image post-processing requires the optional CairoSVG dependency. {}'.format(PILLOW_INSTALL_HINT)
        ) from exc
    return cairosvg


def get_resampling_filter(Image):
    resampling = getattr(Image, 'Resampling', Image)
    return getattr(resampling, 'LANCZOS')


def get_nearest_resampling_filter(Image):
    resampling = getattr(Image, 'Resampling', Image)
    return getattr(resampling, 'NEAREST')


def output_extension_for_path(source_path, args, rasterized=False):
    source_ext = os.path.splitext(source_path)[1].lower()
    output_format = getattr(args, 'output_format', 'original')
    if output_format == 'original':
        if rasterized and source_ext == '.svg':
            return '.png'
        return source_ext or '.bin'
    if output_format == 'png':
        return '.png'
    if output_format == 'jpg':
        return '.jpg'
    raise ValueError('Unsupported --output_format: {}'.format(output_format))


def estimate_background_color_from_border(rgb_image):
    width, height = rgb_image.size
    if width <= 0 or height <= 0:
        return (255, 255, 255)
    if width == 1 and height == 1:
        return rgb_image.getpixel((0, 0))
    samples = []
    x_step = max(1, width // 64)
    y_step = max(1, height // 64)
    for x in range(0, width, x_step):
        samples.append(rgb_image.getpixel((x, 0)))
        if height > 1:
            samples.append(rgb_image.getpixel((x, height - 1)))
    for y in range(y_step, max(y_step, height - 1), y_step):
        samples.append(rgb_image.getpixel((0, y)))
        if width > 1:
            samples.append(rgb_image.getpixel((width - 1, y)))
    if not samples:
        return rgb_image.getpixel((0, 0))
    buckets = defaultdict(lambda: [0, 0, 0, 0])
    for red, green, blue in samples:
        key = (red // 16, green // 16, blue // 16)
        buckets[key][0] += red
        buckets[key][1] += green
        buckets[key][2] += blue
        buckets[key][3] += 1
    best_bucket = max(buckets.values(), key=lambda item: item[3])
    return tuple(int(round(best_bucket[index] / best_bucket[3])) for index in range(3))


def otsu_threshold_from_histogram(histogram):
    total = sum(histogram)
    if total <= 0:
        return SEMANTIC_DIFF_THRESHOLD_FLOOR
    total_weight = sum(index * count for index, count in enumerate(histogram))
    background_weight = 0
    background_sum = 0
    max_variance = -1
    threshold = SEMANTIC_DIFF_THRESHOLD_FLOOR
    for index, count in enumerate(histogram):
        background_weight += count
        if background_weight == 0:
            continue
        foreground_weight = total - background_weight
        if foreground_weight == 0:
            break
        background_sum += index * count
        background_mean = background_sum / background_weight
        foreground_mean = (total_weight - background_sum) / foreground_weight
        between_variance = background_weight * foreground_weight * ((background_mean - foreground_mean) ** 2)
        if between_variance > max_variance:
            max_variance = between_variance
            threshold = index
    return max(SEMANTIC_DIFF_THRESHOLD_FLOOR, min(96, int(threshold)))


def largest_component_bbox(mask, Image):
    bbox = mask.getbbox()
    if bbox is None:
        return None
    cropped_mask = mask.crop(bbox)
    source_width, source_height = cropped_mask.size
    working_mask = cropped_mask
    if max(source_width, source_height) > SEMANTIC_MASK_MAX_DIM:
        if source_width >= source_height:
            working_width = SEMANTIC_MASK_MAX_DIM
            working_height = max(1, int(round(source_height * SEMANTIC_MASK_MAX_DIM / float(source_width))))
        else:
            working_height = SEMANTIC_MASK_MAX_DIM
            working_width = max(1, int(round(source_width * SEMANTIC_MASK_MAX_DIM / float(source_height))))
        working_mask = cropped_mask.resize((working_width, working_height), get_nearest_resampling_filter(Image))
    working_width, working_height = working_mask.size
    pixels = working_mask.load()
    visited = bytearray(working_width * working_height)
    best_area = 0
    best_bbox = None
    for y in range(working_height):
        row_offset = y * working_width
        for x in range(working_width):
            index = row_offset + x
            if visited[index]:
                continue
            visited[index] = 1
            if pixels[x, y] == 0:
                continue
            stack = [(x, y)]
            area = 0
            min_x = max_x = x
            min_y = max_y = y
            while stack:
                current_x, current_y = stack.pop()
                area += 1
                if current_x < min_x:
                    min_x = current_x
                elif current_x > max_x:
                    max_x = current_x
                if current_y < min_y:
                    min_y = current_y
                elif current_y > max_y:
                    max_y = current_y
                if current_x > 0:
                    neighbor_index = current_y * working_width + (current_x - 1)
                    if not visited[neighbor_index]:
                        visited[neighbor_index] = 1
                        if pixels[current_x - 1, current_y] != 0:
                            stack.append((current_x - 1, current_y))
                if current_x + 1 < working_width:
                    neighbor_index = current_y * working_width + (current_x + 1)
                    if not visited[neighbor_index]:
                        visited[neighbor_index] = 1
                        if pixels[current_x + 1, current_y] != 0:
                            stack.append((current_x + 1, current_y))
                if current_y > 0:
                    neighbor_index = (current_y - 1) * working_width + current_x
                    if not visited[neighbor_index]:
                        visited[neighbor_index] = 1
                        if pixels[current_x, current_y - 1] != 0:
                            stack.append((current_x, current_y - 1))
                if current_y + 1 < working_height:
                    neighbor_index = (current_y + 1) * working_width + current_x
                    if not visited[neighbor_index]:
                        visited[neighbor_index] = 1
                        if pixels[current_x, current_y + 1] != 0:
                            stack.append((current_x, current_y + 1))
            if area > best_area:
                best_area = area
                best_bbox = (min_x, min_y, max_x + 1, max_y + 1)
    if best_bbox is None or best_area < 4:
        return bbox
    scale_x = source_width / float(working_width)
    scale_y = source_height / float(working_height)
    left = bbox[0] + int(math.floor(best_bbox[0] * scale_x))
    top = bbox[1] + int(math.floor(best_bbox[1] * scale_y))
    right = bbox[0] + int(math.ceil(best_bbox[2] * scale_x))
    bottom = bbox[1] + int(math.ceil(best_bbox[3] * scale_y))
    right = min(bbox[0] + source_width, max(left + 1, right))
    bottom = min(bbox[1] + source_height, max(top + 1, bottom))
    return (left, top, right, bottom)


def semantic_foreground_bbox(image, Image, ImageChops):
    rgba = image.convert('RGBA')
    alpha = rgba.getchannel('A')
    alpha_bbox = alpha.getbbox()
    alpha_min, alpha_max = alpha.getextrema()
    if alpha_bbox and alpha_max > 0 and alpha_min < 255:
        alpha_mask = alpha.point(lambda value: 255 if value > SEMANTIC_ALPHA_THRESHOLD else 0)
        return largest_component_bbox(alpha_mask, Image) or alpha_bbox
    rgb = image.convert('RGB')
    background = estimate_background_color_from_border(rgb)
    background_image = Image.new('RGB', rgb.size, background)
    diff = ImageChops.difference(rgb, background_image)
    red, green, blue = diff.split()
    max_diff = ImageChops.lighter(ImageChops.lighter(red, green), blue)
    threshold = otsu_threshold_from_histogram(max_diff.histogram())
    semantic_mask = max_diff.point(lambda value: 255 if value > threshold else 0)
    bbox = largest_component_bbox(semantic_mask, Image)
    if bbox is not None:
        return bbox
    return max_diff.getbbox()


def trim_image(image, trim_mode, trim_shape, background_mode, Image, ImageChops):
    if trim_mode == 'off':
        trimmed = image
    elif trim_mode == 'transparent':
        rgba = image.convert('RGBA')
        bbox = rgba.getchannel('A').getbbox()
        trimmed = rgba.crop(bbox) if bbox else rgba
    elif trim_mode == 'white':
        rgb = image.convert('RGB')
        background = rgb.getpixel((0, 0)) if rgb.size[0] > 0 and rgb.size[1] > 0 else (255, 255, 255)
        diff = ImageChops.difference(rgb, Image.new('RGB', rgb.size, background))
        bbox = diff.getbbox()
        trimmed = image.crop(bbox) if bbox else image
    elif trim_mode == 'semantic':
        bbox = semantic_foreground_bbox(image, Image, ImageChops)
        trimmed = image.crop(bbox) if bbox else image
    else:
        raise ValueError('Unsupported --trim mode: {}'.format(trim_mode))
    if trim_shape == 'bbox':
        return trimmed
    if trim_shape == 'square':
        width, height = trimmed.size
        side = min(width, height)
        left = max(0, (width - side) // 2)
        top = max(0, (height - side) // 2)
        return trimmed.crop((left, top, left + side, top + side))
    raise ValueError('Unsupported --trim_shape: {}'.format(trim_shape))


def resize_image_max_edge(image, max_edge, Image):
    if max_edge in (None, 0):
        return image
    width, height = image.size
    if max(width, height) <= int(max_edge):
        return image
    resized = image.copy()
    resized.thumbnail((int(max_edge), int(max_edge)), get_resampling_filter(Image))
    return resized


def square_canvas_image(image, background_mode, Image):
    if background_mode == 'transparent':
        canvas = Image.new('RGBA', (max(image.size), max(image.size)), (0, 0, 0, 0))
        overlay = image.convert('RGBA')
        offset = ((canvas.size[0] - overlay.size[0]) // 2, (canvas.size[1] - overlay.size[1]) // 2)
        canvas.paste(overlay, offset, overlay)
        return canvas
    canvas = Image.new('RGBA', (max(image.size), max(image.size)), (255, 255, 255, 255))
    overlay = image.convert('RGBA')
    offset = ((canvas.size[0] - overlay.size[0]) // 2, (canvas.size[1] - overlay.size[1]) // 2)
    canvas.paste(overlay, offset, overlay)
    return canvas


def save_processed_raster_image(image, destination_path):
    ensure_directory(os.path.dirname(destination_path))
    dest_ext = os.path.splitext(destination_path)[1].lower()
    image_format = RASTER_OUTPUT_EXTENSIONS.get(dest_ext)
    if image_format is None:
        raise RuntimeError('Unsupported output extension for processed image: {}'.format(dest_ext or '(none)'))
    image_to_save = image
    save_kwargs = dict()
    if image_format == 'JPEG':
        if image.mode not in ('RGB', 'L'):
            Image, _, _ = load_pillow_modules()
            background = Image.new('RGBA', image.size, (255, 255, 255, 255))
            background.alpha_composite(image.convert('RGBA'))
            image_to_save = background.convert('RGB')
        save_kwargs.update({'quality': 95})
    elif image_format == 'PNG':
        if image.mode not in ('RGB', 'RGBA', 'L', 'LA'):
            image_to_save = image.convert('RGBA')
    tmp_path = destination_path + '.tmp'
    image_to_save.save(tmp_path, format=image_format, **save_kwargs)
    os.replace(tmp_path, destination_path)


def rasterize_svg_to_image(source_path):
    cairosvg = load_cairosvg_module()
    Image, _, _ = load_pillow_modules()
    png_bytes = cairosvg.svg2png(url=source_path)
    with Image.open(io.BytesIO(png_bytes)) as image:
        rasterized = image.convert('RGBA')
        rasterized.load()
    return rasterized


def process_raster_image(source_path, args):
    Image, ImageChops, ImageOps = load_pillow_modules()
    with Image.open(source_path) as source_image:
        image = ImageOps.exif_transpose(source_image)
        image.load()
    image = trim_image(
        image,
        getattr(args, 'trim', 'off'),
        getattr(args, 'trim_shape', 'bbox'),
        getattr(args, 'background', 'white'),
        Image,
        ImageChops,
    )
    image = resize_image_max_edge(image, getattr(args, 'max_edge', None), Image)
    if getattr(args, 'canvas', 'none') == 'square':
        image = square_canvas_image(image, getattr(args, 'background', 'white'), Image)
    destination_ext = output_extension_for_path(source_path, args=args, rasterized=False)
    if destination_ext in ('.jpg', '.jpeg') and getattr(args, 'background', 'white') == 'transparent':
        raise RuntimeError('Transparent background is incompatible with JPEG output.')
    destination_path = replace_extension(source_path, destination_ext)
    save_processed_raster_image(image, destination_path)
    if destination_path != source_path and os.path.exists(source_path):
        os.remove(source_path)
    return destination_path


def process_svg_image(source_path, args):
    rasterized = rasterize_svg_to_image(source_path)
    Image, ImageChops, _ = load_pillow_modules()
    image = trim_image(
        rasterized,
        getattr(args, 'trim', 'off'),
        getattr(args, 'trim_shape', 'bbox'),
        getattr(args, 'background', 'white'),
        Image,
        ImageChops,
    )
    image = resize_image_max_edge(image, getattr(args, 'max_edge', None), Image)
    if getattr(args, 'canvas', 'none') == 'square':
        image = square_canvas_image(image, getattr(args, 'background', 'white'), Image)
    destination_ext = output_extension_for_path(source_path, args=args, rasterized=True)
    if destination_ext in ('.jpg', '.jpeg') and getattr(args, 'background', 'white') == 'transparent':
        raise RuntimeError('Transparent background is incompatible with JPEG output.')
    destination_path = replace_extension(source_path, destination_ext)
    save_processed_raster_image(image, destination_path)
    if destination_path != source_path and os.path.exists(source_path):
        os.remove(source_path)
    return destination_path


def postprocess_media_file(source_path, args):
    normalized_path = normalize_existing_media_path(source_path)
    if not image_postprocessing_requested(args):
        return normalized_path
    source_ext = os.path.splitext(normalized_path)[1].lower()
    if source_ext == '.svg':
        return process_svg_image(normalized_path, args=args)
    return process_raster_image(normalized_path, args=args)


def resolve_image_cache_dir(args=None):
    download_dir = resolve_download_dir(args)
    if download_dir is None:
        return None
    return os.path.join(download_dir, 'nwkit', 'image-cache')


def resolve_image_query_cache_dir(args=None):
    download_dir = resolve_download_dir(args)
    if download_dir is not None:
        return os.path.join(download_dir, 'nwkit', 'image-query-cache')
    out_dir = getattr(args, 'out_dir', None) if args is not None else None
    if out_dir not in (None, ''):
        return os.path.join(os.path.realpath(out_dir), '.nwkit-cache', 'image-query-cache')
    return os.path.join(tempfile.gettempdir(), 'nwkit', 'image-query-cache')


def resolve_bioicons_catalog_cache_path(args=None):
    return os.path.join(
        resolve_image_query_cache_dir(args),
        'bioicons-catalog-v{}.json'.format(BIOICONS_CATALOG_CACHE_VERSION),
    )


def build_media_cache_path(cache_dir, media_url, provider, provider_record_id):
    ext = infer_extension(media_url, default_ext='.bin')
    digest = hashlib.sha256(str(media_url).encode('utf-8')).hexdigest()[:16]
    filename = '{}__{}__{}{}'.format(
        sanitize_filename_component(provider),
        sanitize_filename_component(provider_record_id),
        digest,
        ext,
    )
    return os.path.join(cache_dir, filename)


def build_query_cache_path(cache_dir, provider, species_name, fallback_rank):
    digest = hashlib.sha256(
        '{}\0{}\0{}'.format(provider, fallback_rank, species_name).encode('utf-8')
    ).hexdigest()[:16]
    filename = '{}__{}__{}__{}.json'.format(
        sanitize_filename_component(provider),
        sanitize_filename_component(str(species_name).replace(' ', '_')),
        sanitize_filename_component(fallback_rank),
        digest,
    )
    return os.path.join(cache_dir, filename)


def load_cached_provider_candidates(cache_path):
    try:
        with open(cache_path) as handle:
            payload = json.load(handle)
    except (FileNotFoundError, json.JSONDecodeError, OSError):
        return None
    if payload.get('version') != IMAGE_QUERY_CACHE_VERSION:
        return None
    candidates = payload.get('candidates')
    if not isinstance(candidates, list):
        return None
    return candidates


def write_cached_provider_candidates(cache_path, candidates):
    ensure_directory(os.path.dirname(cache_path))
    payload = {
        'version': IMAGE_QUERY_CACHE_VERSION,
        'candidates': candidates,
    }
    tmp_path = cache_path + '.tmp'
    with open(tmp_path, 'w') as handle:
        json.dump(payload, handle, sort_keys=True)
    os.replace(tmp_path, cache_path)


def load_cached_bioicons_catalog(cache_path):
    try:
        with open(cache_path) as handle:
            payload = json.load(handle)
    except (FileNotFoundError, json.JSONDecodeError, OSError):
        return None
    if payload.get('version') != BIOICONS_CATALOG_CACHE_VERSION:
        return None
    catalog = payload.get('catalog')
    return catalog if isinstance(catalog, list) else None


def write_cached_bioicons_catalog(cache_path, catalog):
    ensure_directory(os.path.dirname(cache_path))
    payload = {
        'version': BIOICONS_CATALOG_CACHE_VERSION,
        'catalog': catalog,
    }
    tmp_path = cache_path + '.tmp'
    with open(tmp_path, 'w') as handle:
        json.dump(payload, handle, sort_keys=True)
    os.replace(tmp_path, cache_path)


def get_style_priority(candidate, style='auto'):
    asset_type = candidate.get('asset_type')
    if style == 'silhouette':
        return 2 if asset_type == 'silhouette' else 0
    if style == 'photo':
        return 2 if asset_type == 'photo' else 0
    return 1 if asset_type in ('photo', 'silhouette') else 0


def get_provider_quality_bonus(candidate):
    return int(candidate.get('provider_quality', 0) or 0)


def dedupe_sorted_candidates(candidates):
    deduped = list()
    seen_urls = set()
    for candidate in sorted(candidates, key=lambda item: item['score'], reverse=True):
        media_url = candidate.get('media_url')
        if media_url in seen_urls:
            continue
        seen_urls.add(media_url)
        deduped.append(candidate)
    return deduped


def allowed_candidates_from_scored_candidates(candidates, args):
    return [
        candidate for candidate in dedupe_sorted_candidates(candidates)
        if license_allowed(
            candidate['license_code'],
            license_max=args.license_max,
            allow_nd=args.allow_nd,
        )
    ]


def should_stop_after_provider(candidates, args):
    max_per_species = int(getattr(args, 'max_per_species', 1) or 1)
    allowed_candidates = allowed_candidates_from_scored_candidates(candidates, args=args)
    if len(allowed_candidates) < (max_per_species + LOOKUP_FALLBACK_BUFFER):
        return False
    top_candidates = allowed_candidates[:max_per_species]
    return all(candidate.get('matched_rank') == 'species' for candidate in top_candidates)


def resolve_provider_fetch_limit(args=None, minimum=10, maximum=30, extra_buffer=2):
    try:
        max_per_species = int(getattr(args, 'max_per_species', 1) or 1)
    except (TypeError, ValueError):
        max_per_species = 1
    desired = max_per_species + LOOKUP_FALLBACK_BUFFER + int(extra_buffer)
    return max(int(minimum), min(int(maximum), desired))


def resolve_ncbi_taxonomy_image_cache_dir(args=None):
    download_dir = resolve_download_dir(args)
    if download_dir is not None:
        return os.path.join(download_dir, 'ncbi-taxonomy-images')
    out_dir = getattr(args, 'out_dir', None) if args is not None else None
    if out_dir not in (None, ''):
        return os.path.join(os.path.realpath(out_dir), '.nwkit-cache', 'ncbi-taxonomy-images')
    return os.path.join(tempfile.gettempdir(), 'nwkit', 'ncbi-taxonomy-images')


def _download_to_path(session, url, destination_path):
    ensure_directory(os.path.dirname(destination_path))
    tmp_path = destination_path + '.tmp'
    response = session.get(url, stream=True, timeout=REQUEST_TIMEOUT, headers=HTTP_HEADERS)
    response.raise_for_status()
    with open(tmp_path, 'wb') as handle:
        for chunk in response.iter_content(chunk_size=1024 * 256):
            if chunk:
                handle.write(chunk)
    os.replace(tmp_path, destination_path)


def _extract_ncbi_images_table(tarball_path, destination_path):
    ensure_directory(os.path.dirname(destination_path))
    tmp_path = destination_path + '.tmp'
    with tarfile.open(tarball_path, 'r:gz') as handle:
        member = handle.getmember('images.dmp')
        extracted = handle.extractfile(member)
        if extracted is None:
            raise FileNotFoundError('images.dmp not found in {}'.format(tarball_path))
        with open(tmp_path, 'wb') as out_handle:
            shutil.copyfileobj(extracted, out_handle)
    os.replace(tmp_path, destination_path)


def ensure_ncbi_images_table(args, session):
    cache_dir = resolve_ncbi_taxonomy_image_cache_dir(args)
    images_path = os.path.join(cache_dir, 'images.dmp')
    if os.path.exists(images_path) and os.path.getsize(images_path) > 0:
        return images_path
    lock_path = os.path.join(cache_dir, '.ncbi_images.lock')
    tarball_path = os.path.join(cache_dir, 'new_taxdump.tar.gz')
    with acquire_exclusive_lock(lock_path=lock_path, lock_label='NCBI taxonomy images'):
        if os.path.exists(images_path) and os.path.getsize(images_path) > 0:
            return images_path
        if (not os.path.exists(tarball_path)) or (os.path.getsize(tarball_path) == 0):
            _download_to_path(session=session, url=NCBI_NEWTAXDUMP_URL, destination_path=tarball_path)
        _extract_ncbi_images_table(tarball_path=tarball_path, destination_path=images_path)
        try:
            os.remove(tarball_path)
        except FileNotFoundError:
            pass
    return images_path


def parse_ncbi_images_dmp_line(line):
    parts = [part.strip() for part in str(line).rstrip('\n').split('|')]
    if parts and parts[-1] == '':
        parts = parts[:-1]
    if len(parts) < 8:
        raise ValueError('Unexpected images.dmp row: {}'.format(line))
    license_text = parts[3]
    license_url_match = re.search(r'\((https?://[^)]+)\)', license_text)
    license_url = license_url_match.group(1) if license_url_match else ''
    license_code_text = re.sub(r'\s*\(https?://[^)]+\)\s*', '', license_text).strip()
    taxids = [int(token) for token in parts[7].split() if token.strip() != '']
    return {
        'record_id': parts[0],
        'title': re.sub(r'^image:', '', parts[1]),
        'image_url': parts[2],
        'license_text': license_text,
        'license_code_text': license_code_text,
        'license_url': license_url,
        'attribution': parts[4],
        'source_name': parts[5],
        'taxids': taxids,
    }


@functools.lru_cache(maxsize=16)
def load_ncbi_images_index(images_path):
    index = defaultdict(list)
    with open(images_path, 'r') as handle:
        for raw_line in handle:
            stripped = raw_line.strip()
            if stripped == '':
                continue
            record = parse_ncbi_images_dmp_line(raw_line)
            if len(record['taxids']) == 0:
                continue
            for taxid in record['taxids']:
                index[taxid].append(record)
    return {taxid: records for taxid, records in index.items()}


def candidate_score(candidate, provider_index=0, style='auto'):
    rank_priority = {'species': 3, 'genus': 2, 'family': 1}.get(candidate['matched_rank'], 0)
    provider_priority = max(1, 100 - int(provider_index))
    license_priority = max(0, LICENSE_OPENNESS.get(candidate['license_code'], 0))
    style_priority = get_style_priority(candidate, style=style)
    quality_priority = min(max(candidate.get('width') or 0, candidate.get('height') or 0), 10000)
    provider_quality = max(0, min(get_provider_quality_bonus(candidate), 99))
    score = (
        rank_priority * 10**12
        + provider_priority * 10**10
        + license_priority * 10**8
        + style_priority * 10**7
        + quality_priority * 10**2
        + provider_quality
    )
    return int(score)


class ProviderError(RuntimeError):
    pass


class LazyNCBITaxa:
    def __init__(self, args):
        self.args = args
        self._ncbi = None

    def _get_ncbi(self):
        if self._ncbi is None:
            self._ncbi = get_ete_ncbitaxa(args=self.args)
        return self._ncbi

    def get_name_translator(self, names):
        return self._get_ncbi().get_name_translator(names)

    def get_lineage(self, taxid):
        return self._get_ncbi().get_lineage(taxid)

    def get_rank(self, lineage):
        return self._get_ncbi().get_rank(lineage)

    def get_taxid_translator(self, lineage):
        return self._get_ncbi().get_taxid_translator(lineage)

    def close(self):
        if self._ncbi is None:
            return
        _close_ncbi_db(self._ncbi)
        self._ncbi = None


class PhylopicProvider:
    provider_name = 'phylopic'

    def __init__(self, session, ncbi):
        self.session = session
        self.ncbi = ncbi

    def _get_taxid(self, query_name):
        name_to_taxid = self.ncbi.get_name_translator([query_name])
        taxids = name_to_taxid.get(query_name, [])
        return int(taxids[0]) if taxids else None

    def _resolve_node(self, taxid):
        response = self.session.get(
            f'{PHYLIPIC_API_ROOT}/resolve/ncbi.nlm.nih.gov/taxid/{taxid}',
            timeout=REQUEST_TIMEOUT,
            headers=HTTP_HEADERS,
        )
        if response.status_code == 404:
            return None
        response.raise_for_status()
        return response.json()

    def _fetch_node_images(self, node_uuid, build):
        response = self.session.get(
            f'{PHYLIPIC_API_ROOT}/images',
            params={
                'build': build,
                'filter_node': node_uuid,
                'page': 0,
                'embed_items': 'true',
            },
            timeout=REQUEST_TIMEOUT,
            headers=HTTP_HEADERS,
        )
        response.raise_for_status()
        payload = response.json()
        return payload.get('_embedded', {}).get('items', [])

    def fetch_candidates(self, species_name, fallback_rank='none'):
        candidates = list()
        for matched_rank, query_name in get_taxonomic_queries(species_name, fallback_rank=fallback_rank, ncbi=self.ncbi):
            taxid = self._get_taxid(query_name)
            if taxid is None:
                continue
            node = self._resolve_node(taxid)
            if node is None:
                continue
            build = node.get('build')
            node_uuid = node.get('uuid')
            if build is None or node_uuid is None:
                continue
            for image_item in self._fetch_node_images(node_uuid=node_uuid, build=build):
                candidate = self._candidate_from_image(
                    image_item=image_item,
                    matched_name=query_name,
                    matched_rank=matched_rank,
                )
                if candidate is not None:
                    candidates.append(candidate)
        return candidates

    def _candidate_from_image(self, image_item, matched_name, matched_rank):
        links = image_item.get('_links', {})
        vector_file = links.get('vectorFile')
        source_file = links.get('sourceFile')
        raster_files = links.get('rasterFiles', [])
        selected_file = vector_file or source_file
        if selected_file is None and raster_files:
            selected_file = max(
                raster_files,
                key=lambda item: max(parse_size(item.get('sizes'))[0] or 0, parse_size(item.get('sizes'))[1] or 0),
            )
        if selected_file is None:
            return None

        license_url = links.get('license', {}).get('href', '')
        license_code = normalize_license_code(raw_url=license_url, attribution=image_item.get('attribution'))
        width, height = parse_size(selected_file.get('sizes'))
        self_link = links.get('self', {}).get('href', '')

        return {
            'provider': self.provider_name,
            'provider_record_id': image_item.get('uuid', ''),
            'matched_name': matched_name,
            'matched_rank': matched_rank,
            'license_code': license_code,
            'license_url': license_url,
            'attribution': image_item.get('attribution', ''),
            'source_page_url': f'{PHYLIPIC_API_ROOT}{self_link}' if self_link else '',
            'media_url': selected_file.get('href', ''),
            'width': width,
            'height': height,
            'asset_type': 'silhouette',
            'provider_quality': 30 if vector_file is not None else 10,
        }


class BioiconsProvider:
    provider_name = 'bioicons'

    def __init__(self, session, ncbi=None, args=None):
        self.session = session
        self.ncbi = ncbi
        self.args = args
        self._catalog = None
        self._catalog_index = None

    def _load_catalog(self):
        if self._catalog is not None:
            return self._catalog
        cache_path = resolve_bioicons_catalog_cache_path(args=self.args)
        with BIOICONS_CATALOG_MEMORY_CACHE_LOCK:
            cached_catalog = BIOICONS_CATALOG_MEMORY_CACHE.get(cache_path)
        if cached_catalog is not None:
            self._catalog = cached_catalog
            return self._catalog

        catalog = load_cached_bioicons_catalog(cache_path)
        if catalog is None:
            lock_path = cache_path + '.lock'
            with acquire_exclusive_lock(lock_path=lock_path, lock_label='Bioicons catalog'):
                catalog = load_cached_bioicons_catalog(cache_path)
                if catalog is None:
                    response = self.session.get(
                        BIOICONS_GITHUB_API_ROOT,
                        params={'recursive': 1},
                        timeout=REQUEST_TIMEOUT,
                        headers=HTTP_HEADERS,
                    )
                    response.raise_for_status()
                    payload = response.json()
                    catalog = list()
                    for item in payload.get('tree', []):
                        path = item.get('path', '')
                        if (not path.startswith('static/icons/')) or (not path.endswith('.svg')):
                            continue
                        relative_path = path[len('static/icons/'):]
                        parts = relative_path.split('/')
                        if len(parts) != 4:
                            continue
                        license_slug, category, author_slug, filename = parts
                        name, _ = os.path.splitext(filename)
                        catalog.append({
                            'license_slug': license_slug,
                            'category': category,
                            'author_slug': author_slug,
                            'name': name,
                            'relative_path': relative_path,
                        })
                    write_cached_bioicons_catalog(cache_path, catalog)
        with BIOICONS_CATALOG_MEMORY_CACHE_LOCK:
            BIOICONS_CATALOG_MEMORY_CACHE[cache_path] = catalog
        self._catalog = catalog
        return self._catalog

    def _load_catalog_index(self):
        if self._catalog_index is not None:
            return self._catalog_index
        index = defaultdict(list)
        for icon in self._load_catalog():
            for key in bioicons_index_keys_for_name(icon['name']):
                index[key].append(icon)
        self._catalog_index = dict(index)
        return self._catalog_index

    def fetch_candidates(self, species_name, fallback_rank='none'):
        candidates = list()
        catalog = self._load_catalog()
        catalog_index = self._load_catalog_index()
        seen_paths = set()
        for matched_rank, query_name in get_taxonomic_queries(species_name, fallback_rank=fallback_rank, ncbi=self.ncbi):
            query_keys = bioicons_index_keys_for_name(query_name)
            for alias in BIOICONS_SPECIES_ALIASES.get(str(query_name or '').strip().lower(), ()):
                query_keys.update(bioicons_index_keys_for_name(alias))
            candidate_icons = list()
            seen_icon_paths = set()
            for key in query_keys:
                for icon in catalog_index.get(key, []):
                    relative_path = icon['relative_path']
                    if relative_path in seen_icon_paths:
                        continue
                    seen_icon_paths.add(relative_path)
                    candidate_icons.append(icon)
            if len(candidate_icons) == 0:
                candidate_icons = catalog
            for icon in candidate_icons:
                quality = bioicons_match_quality(
                    icon_name=icon['name'],
                    query_name=query_name,
                    matched_rank=matched_rank,
                )
                if quality <= 0:
                    continue
                relative_path = icon['relative_path']
                if relative_path in seen_paths:
                    continue
                seen_paths.add(relative_path)
                license_code = normalize_license_code(raw_code=icon['license_slug'])
                media_url = '{}/{}'.format(BIOICONS_MEDIA_ROOT.rstrip('/'), relative_path)
                candidates.append({
                    'provider': self.provider_name,
                    'provider_record_id': relative_path,
                    'matched_name': query_name,
                    'matched_rank': matched_rank,
                    'license_code': license_code,
                    'license_url': canonical_license_url(license_code),
                    'attribution': bioicons_display_author(icon['author_slug']),
                    'source_page_url': media_url,
                    'media_url': media_url,
                    'width': None,
                    'height': None,
                    'asset_type': 'silhouette',
                    'provider_quality': quality,
                })
        return candidates


class INaturalistProvider:
    provider_name = 'inaturalist'

    def __init__(self, session, ncbi=None, args=None):
        self.session = session
        self.ncbi = ncbi
        self.result_limit = resolve_provider_fetch_limit(args=args, minimum=10, maximum=30)

    def _find_taxon(self, query_name, expected_rank):
        response = self.session.get(
            f'{INATURALIST_API_ROOT}/taxa',
            params={'q': query_name, 'per_page': self.result_limit},
            timeout=REQUEST_TIMEOUT,
            headers=HTTP_HEADERS,
        )
        response.raise_for_status()
        results = response.json().get('results', [])
        expected_name = query_name.lower()
        matches = [
            item for item in results
            if str(item.get('name', '')).lower() == expected_name and item.get('rank') == expected_rank
        ]
        if matches:
            return matches[0]
        return None

    def _fetch_observations(self, taxon_id, per_page):
        response = self.session.get(
            f'{INATURALIST_API_ROOT}/observations',
            params={
                'taxon_id': taxon_id,
                'photos': 'true',
                'per_page': per_page,
                'order': 'desc',
                'order_by': 'votes',
            },
            timeout=REQUEST_TIMEOUT,
            headers=HTTP_HEADERS,
        )
        response.raise_for_status()
        return response.json().get('results', [])

    def fetch_candidates(self, species_name, fallback_rank='none'):
        candidates = list()
        per_page = self.result_limit
        for matched_rank, query_name in get_taxonomic_queries(species_name, fallback_rank=fallback_rank, ncbi=self.ncbi):
            taxon = self._find_taxon(query_name=query_name, expected_rank=matched_rank)
            if taxon is None:
                continue
            observations = self._fetch_observations(taxon_id=taxon['id'], per_page=per_page)
            for observation in observations:
                for photo in observation.get('photos', []):
                    candidate = self._candidate_from_photo(
                        photo=photo,
                        observation=observation,
                        taxon=taxon,
                        matched_rank=matched_rank,
                    )
                    if candidate is not None:
                        candidates.append(candidate)
        return candidates

    def _candidate_from_photo(self, photo, observation, taxon, matched_rank):
        photo_url = photo.get('url')
        if not photo_url:
            return None
        original_url = str(photo_url).replace('/square.', '/original.')
        dimensions = photo.get('original_dimensions', {}) or {}
        license_code = normalize_license_code(raw_code=photo.get('license_code'), attribution=photo.get('attribution'))
        quality_grade_bonus = {
            'research': 30,
            'needs_id': 10,
            'casual': 0,
        }.get(str(observation.get('quality_grade', '')).lower(), 0)
        popularity_bonus = min(int(observation.get('cached_votes_total') or observation.get('faves_count') or 0), 20)
        return {
            'provider': self.provider_name,
            'provider_record_id': str(photo.get('id', '')),
            'matched_name': taxon.get('name', ''),
            'matched_rank': matched_rank,
            'license_code': license_code,
            'license_url': canonical_license_url(license_code),
            'attribution': photo.get('attribution', ''),
            'source_page_url': observation.get('uri', ''),
            'media_url': original_url,
            'width': dimensions.get('width'),
            'height': dimensions.get('height'),
            'asset_type': 'photo',
            'provider_quality': quality_grade_bonus + popularity_bonus,
        }


class EOLProvider:
    provider_name = 'eol'

    def __init__(self, session, ncbi=None, args=None):
        self.session = session
        self.ncbi = ncbi
        self.result_limit = resolve_provider_fetch_limit(args=args, minimum=10, maximum=30)
        self.page_limit = resolve_provider_fetch_limit(args=args, minimum=3, maximum=5, extra_buffer=0)

    def _search_pages(self, query_name):
        response = self.session.get(
            '{}/search/1.0.json'.format(EOL_API_ROOT),
            params={
                'q': query_name,
                'page': 1,
                'exact': 'true',
            },
            timeout=REQUEST_TIMEOUT,
            headers=HTTP_HEADERS,
        )
        response.raise_for_status()
        return response.json().get('results', [])

    def _fetch_page(self, page_id):
        response = self.session.get(
            '{}/pages/1.0/{}.json'.format(EOL_API_ROOT, page_id),
            params={
                'details': 'true',
                'licenses': 'all',
                'images_per_page': self.result_limit,
                'images_page': 1,
                'videos_per_page': 0,
                'sounds_per_page': 0,
                'maps_per_page': 0,
                'texts_per_page': 0,
                'common_names': 'false',
                'synonyms': 'false',
                'references': 'false',
                'taxonomy': 'false',
            },
            timeout=REQUEST_TIMEOUT,
            headers=HTTP_HEADERS,
        )
        response.raise_for_status()
        return response.json()

    def fetch_candidates(self, species_name, fallback_rank='none'):
        candidates = list()
        for matched_rank, query_name in get_taxonomic_queries(species_name, fallback_rank=fallback_rank, ncbi=self.ncbi):
            normalized_query = normalize_species_name(query_name).lower()
            seen_page_ids = set()
            exact_matches = list()
            for result in self._search_pages(query_name=query_name):
                page_id = result.get('id')
                title = normalize_species_name(result.get('title', ''))
                if page_id in seen_page_ids:
                    continue
                seen_page_ids.add(page_id)
                if title and title.lower() == normalized_query:
                    exact_matches.append(result)
            for result in exact_matches[:self.page_limit]:
                payload = self._fetch_page(page_id=result['id'])
                taxon_concept = payload.get('taxonConcept', {})
                for data_object in taxon_concept.get('dataObjects') or []:
                    candidate = self._candidate_from_data_object(
                        data_object=data_object,
                        matched_name=query_name,
                        matched_rank=matched_rank,
                    )
                    if candidate is not None:
                        candidates.append(candidate)
        return candidates

    def _candidate_from_data_object(self, data_object, matched_name, matched_rank):
        data_type = str(data_object.get('dataType', ''))
        medium_type = str(data_object.get('mediumType', '')).lower()
        if (medium_type != 'image') and ('StillImage' not in data_type):
            return None

        media_url = data_object.get('mediaURL') or data_object.get('eolMediaURL')
        if not media_url:
            return None

        license_value = data_object.get('license', '')
        creator_names = [
            strip_html_markup(agent.get('full_name', ''))
            for agent in data_object.get('agents', []) or []
            if str(agent.get('role', '')).lower() in ('creator', 'photographer', 'illustrator')
        ]
        rights_holder = strip_html_markup(data_object.get('rightsHolder', ''))
        attribution = rights_holder or ', '.join([name for name in creator_names if name != ''])
        license_code = normalize_license_code(
            raw_code=license_value if '://' not in str(license_value) else None,
            raw_url=license_value if '://' in str(license_value) else None,
            attribution=attribution,
        )
        vetted_bonus = {
            'trusted': 20,
            'unreviewed': 5,
            'untrusted': 0,
        }.get(str(data_object.get('vettedStatus', '')).lower(), 0)
        try:
            rating_bonus = int(float(data_object.get('dataRating') or 0) * 4)
        except (TypeError, ValueError):
            rating_bonus = 0
        return {
            'provider': self.provider_name,
            'provider_record_id': str(data_object.get('dataObjectVersionID') or data_object.get('identifier', '')),
            'matched_name': matched_name,
            'matched_rank': matched_rank,
            'license_code': license_code,
            'license_url': license_value if '://' in str(license_value) else canonical_license_url(license_code),
            'attribution': attribution,
            'source_page_url': data_object.get('source') or '',
            'media_url': media_url,
            'width': data_object.get('width'),
            'height': data_object.get('height'),
            'asset_type': 'photo',
            'provider_quality': vetted_bonus + rating_bonus,
        }


class IDigBioProvider:
    provider_name = 'idigbio'

    def __init__(self, session, ncbi=None, args=None):
        self.session = session
        self.ncbi = ncbi
        self.result_limit = resolve_provider_fetch_limit(args=args, minimum=10, maximum=30)

    def _build_record_query(self, query_name, matched_rank):
        normalized = normalize_species_name(query_name)
        if matched_rank == 'species':
            return {'scientificname': normalized}
        if matched_rank == 'genus':
            return {'genus': normalized}
        if matched_rank == 'family':
            return {'family': normalized}
        raise ValueError('Unsupported matched rank for iDigBio: {}'.format(matched_rank))

    def _search_media(self, query_name, matched_rank, limit=None):
        body = {
            'rq': self._build_record_query(query_name=query_name, matched_rank=matched_rank),
            'limit': int(self.result_limit if limit is None else limit),
            'offset': 0,
        }
        response = self.session.post(
            '{}/search/media'.format(IDIGBIO_API_ROOT),
            json=body,
            timeout=REQUEST_TIMEOUT,
            headers=HTTP_HEADERS,
        )
        response.raise_for_status()
        return response.json().get('items', [])

    def fetch_candidates(self, species_name, fallback_rank='none'):
        candidates = list()
        for matched_rank, query_name in get_taxonomic_queries(species_name, fallback_rank=fallback_rank, ncbi=self.ncbi):
            for item in self._search_media(query_name=query_name, matched_rank=matched_rank):
                candidate = self._candidate_from_item(
                    item=item,
                    matched_name=query_name,
                    matched_rank=matched_rank,
                )
                if candidate is not None:
                    candidates.append(candidate)
        return candidates

    def _candidate_from_item(self, item, matched_name, matched_rank):
        data = item.get('data') or {}
        index_terms = item.get('indexTerms') or {}
        media_url = data.get('ac:accessURI') or data.get('dcterms:identifier') or index_terms.get('accessuri')
        if not media_url:
            return None

        media_kind = str(index_terms.get('mediatype') or index_terms.get('type') or data.get('dc:type') or '').lower()
        if ('image' not in media_kind) and ('stillimage' not in media_kind):
            return None

        scientific_name = normalize_species_name(data.get('dwc:scientificName') or matched_name)

        license_url = data.get('xmpRights:UsageTerms') or data.get('xmpRights:WebStatement') or ''
        license_code = normalize_license_code(
            raw_code=data.get('dcterms:rights') or data.get('dc:rights') or '',
            raw_url=license_url,
            attribution=data.get('xmpRights:Owner') or data.get('ac:providerLiteral') or '',
        )
        try:
            width = int(data.get('exif:PixelXDimension') or index_terms.get('xpixels') or 0) or None
        except (TypeError, ValueError):
            width = None
        try:
            height = int(data.get('exif:PixelYDimension') or index_terms.get('ypixels') or 0) or None
        except (TypeError, ValueError):
            height = None
        try:
            dqs_bonus = int(float(index_terms.get('dqs') or 0) * 20)
        except (TypeError, ValueError):
            dqs_bonus = 0
        provider_bonus = 10 if data.get('ac:providerLiteral') else 0
        exact_name_bonus = 10 if scientific_name and (scientific_name.lower() == normalize_species_name(matched_name).lower()) else 0
        return {
            'provider': self.provider_name,
            'provider_record_id': str(item.get('uuid', '')),
            'matched_name': scientific_name or matched_name,
            'matched_rank': matched_rank,
            'license_code': license_code,
            'license_url': license_url if '://' in str(license_url) else canonical_license_url(license_code),
            'attribution': data.get('xmpRights:Owner') or data.get('ac:providerLiteral') or '',
            'source_page_url': data.get('dcterms:identifier') or data.get('xmpRights:WebStatement') or media_url,
            'media_url': media_url,
            'width': width,
            'height': height,
            'asset_type': 'photo',
            'provider_quality': dqs_bonus + provider_bonus + exact_name_bonus,
        }


class OpenverseProvider:
    provider_name = 'openverse'

    def __init__(self, session, ncbi=None, args=None):
        self.session = session
        self.ncbi = ncbi
        self.result_limit = resolve_provider_fetch_limit(args=args, minimum=10, maximum=30)

    def _search_images(self, query_name, page_size=None):
        response = self.session.get(
            '{}/images/'.format(OPENVERSE_API_ROOT),
            params={
                'q': query_name,
                'page_size': int(self.result_limit if page_size is None else page_size),
            },
            timeout=REQUEST_TIMEOUT,
            headers=HTTP_HEADERS,
        )
        response.raise_for_status()
        return response.json().get('results', [])

    def fetch_candidates(self, species_name, fallback_rank='none'):
        candidates = list()
        for matched_rank, query_name in get_taxonomic_queries(species_name, fallback_rank=fallback_rank, ncbi=self.ncbi):
            for result in self._search_images(query_name=query_name):
                candidate = self._candidate_from_result(
                    result=result,
                    matched_name=query_name,
                    matched_rank=matched_rank,
                )
                if candidate is not None:
                    candidates.append(candidate)
        return candidates

    def _candidate_from_result(self, result, matched_name, matched_rank):
        media_url = result.get('url')
        if not media_url:
            return None

        text_fragments = [
            result.get('title', ''),
            ' '.join([tag.get('name', '') for tag in result.get('tags', []) or []]),
        ]
        if not search_text_mentions_query(text_fragments, matched_name):
            return None

        license_code = normalize_license_code(
            raw_code=result.get('license'),
            raw_url=result.get('license_url'),
            attribution=result.get('attribution') or result.get('creator'),
        )
        matched_fields = [str(field).lower() for field in result.get('fields_matched', []) or []]
        provider_quality = 0
        if 'title' in matched_fields:
            provider_quality += 20
        if 'tags.name' in matched_fields:
            provider_quality += 10
        if search_text_mentions_query([result.get('title', '')], matched_name):
            provider_quality += 20
        return {
            'provider': self.provider_name,
            'provider_record_id': str(result.get('id', '')),
            'matched_name': matched_name,
            'matched_rank': matched_rank,
            'license_code': license_code,
            'license_url': result.get('license_url') or canonical_license_url(license_code),
            'attribution': result.get('creator', '') or '',
            'source_page_url': result.get('foreign_landing_url') or result.get('detail_url') or '',
            'media_url': media_url,
            'width': result.get('width'),
            'height': result.get('height'),
            'asset_type': 'photo',
            'provider_quality': provider_quality,
        }


class WikimediaProvider:
    provider_name = 'wikimedia'

    def __init__(self, session, ncbi=None, args=None):
        self.session = session
        self.ncbi = ncbi
        self.result_limit = resolve_provider_fetch_limit(args=args, minimum=5, maximum=10)

    def _search_pages(self, query_name, limit=10):
        response = self.session.get(
            WIKIMEDIA_API_ROOT,
            params={
                'action': 'query',
                'generator': 'search',
                'gsrsearch': '"{}" filetype:bitmap'.format(query_name),
                'gsrnamespace': 6,
                'gsrlimit': self.result_limit if limit == 10 else limit,
                'prop': 'imageinfo',
                'iiprop': 'url|size|extmetadata',
                'iiurlwidth': 1600,
                'format': 'json',
            },
            timeout=REQUEST_TIMEOUT,
            headers=HTTP_HEADERS,
        )
        response.raise_for_status()
        payload = response.json()
        pages = list(payload.get('query', {}).get('pages', {}).values())
        return sorted(pages, key=lambda page: int(page.get('index', 10**9)))

    def fetch_candidates(self, species_name, fallback_rank='none'):
        candidates = list()
        for matched_rank, query_name in get_taxonomic_queries(species_name, fallback_rank=fallback_rank, ncbi=self.ncbi):
            for page in self._search_pages(query_name=query_name):
                candidate = self._candidate_from_page(page=page, matched_name=query_name, matched_rank=matched_rank)
                if candidate is not None:
                    candidates.append(candidate)
        return candidates

    def _candidate_from_page(self, page, matched_name, matched_rank):
        if not wikimedia_page_mentions_query(page=page, query_name=matched_name):
            return None
        image_info = (page.get('imageinfo') or [{}])[0]
        media_url = image_info.get('url')
        if not media_url:
            return None
        metadata = image_info.get('extmetadata') or {}
        license_short_name = strip_html_markup(metadata.get('LicenseShortName', {}).get('value', ''))
        license_url = strip_html_markup(metadata.get('LicenseUrl', {}).get('value', ''))
        attribution = strip_html_markup(metadata.get('Artist', {}).get('value', ''))
        width = image_info.get('width') or image_info.get('thumbwidth')
        height = image_info.get('height') or image_info.get('thumbheight')
        assessments = str(metadata.get('Assessments', {}).get('value', '')).lower()
        provider_quality = 0
        if 'featured' in assessments:
            provider_quality += 30
        if 'quality' in assessments:
            provider_quality += 20
        return {
            'provider': self.provider_name,
            'provider_record_id': str(page.get('pageid', page.get('title', ''))),
            'matched_name': matched_name,
            'matched_rank': matched_rank,
            'license_code': normalize_license_code(
                raw_code=license_short_name,
                raw_url=license_url,
                attribution=attribution,
            ),
            'license_url': license_url,
            'attribution': attribution,
            'source_page_url': image_info.get('descriptionurl', ''),
            'media_url': media_url,
            'width': width,
            'height': height,
            'asset_type': 'photo',
            'provider_quality': provider_quality,
        }


class GBIFProvider:
    provider_name = 'gbif'

    def __init__(self, session, ncbi=None, args=None):
        self.session = session
        self.ncbi = ncbi
        self.result_limit = resolve_provider_fetch_limit(args=args, minimum=10, maximum=30)

    def _match_taxon(self, query_name, matched_rank):
        response = self.session.get(
            f'{GBIF_API_ROOT}/species/match',
            params={
                'name': query_name,
                'rank': matched_rank,
                'strict': 'true',
                'verbose': 'true',
            },
            timeout=REQUEST_TIMEOUT,
            headers=HTTP_HEADERS,
        )
        response.raise_for_status()
        payload = response.json()
        if str(payload.get('matchType', '')).upper() not in ('EXACT', 'HIGHERRANK'):
            return None
        usage_key = payload.get('usageKey')
        if usage_key is None:
            key_name = '{}Key'.format(matched_rank)
            usage_key = payload.get(key_name)
        if usage_key is None:
            return None
        return {
            'usage_key': int(usage_key),
            'matched_name': payload.get('canonicalName') or query_name,
        }

    def _fetch_occurrences(self, usage_key, limit=None):
        response = self.session.get(
            f'{GBIF_API_ROOT}/occurrence/search',
            params={
                'taxon_key': usage_key,
                'media_type': 'StillImage',
                'limit': self.result_limit if limit is None else limit,
            },
            timeout=REQUEST_TIMEOUT,
            headers=HTTP_HEADERS,
        )
        response.raise_for_status()
        return response.json().get('results', [])

    def fetch_candidates(self, species_name, fallback_rank='none'):
        candidates = list()
        for matched_rank, query_name in get_taxonomic_queries(species_name, fallback_rank=fallback_rank, ncbi=self.ncbi):
            match = self._match_taxon(query_name=query_name, matched_rank=matched_rank)
            if match is None:
                continue
            for occurrence in self._fetch_occurrences(usage_key=match['usage_key']):
                for media in occurrence.get('media', []) or []:
                    candidate = self._candidate_from_media(
                        occurrence=occurrence,
                        media=media,
                        matched_name=match['matched_name'],
                        matched_rank=matched_rank,
                    )
                    if candidate is not None:
                        candidates.append(candidate)
        return candidates

    def _candidate_from_media(self, occurrence, media, matched_name, matched_rank):
        media_url = media.get('identifier')
        if not media_url:
            return None
        license_url = media.get('license') or occurrence.get('license') or ''
        attribution = media.get('creator') or media.get('rightsHolder') or occurrence.get('rightsHolder') or ''
        provider_record_id = media.get('references') or occurrence.get('key') or media_url
        provider_quality = 0
        if media.get('publisher'):
            provider_quality += 5
        if media.get('references'):
            provider_quality += 5
        return {
            'provider': self.provider_name,
            'provider_record_id': str(provider_record_id),
            'matched_name': matched_name,
            'matched_rank': matched_rank,
            'license_code': normalize_license_code(
                raw_url=license_url,
                attribution=attribution,
            ),
            'license_url': license_url,
            'attribution': attribution,
            'source_page_url': media.get('references') or occurrence.get('references') or '',
            'media_url': media_url,
            'width': None,
            'height': None,
            'asset_type': 'photo',
            'provider_quality': provider_quality,
        }


class NCBIProvider:
    provider_name = 'ncbi'

    def __init__(self, session, ncbi, args):
        self.session = session
        self.ncbi = ncbi
        self.args = args
        self.images_path = None
        self.images_index = None

    def _ensure_images_index(self):
        if self.images_index is not None:
            return
        self.images_path = ensure_ncbi_images_table(args=self.args, session=self.session)
        self.images_index = load_ncbi_images_index(self.images_path)

    def _get_taxid(self, query_name):
        name_to_taxid = self.ncbi.get_name_translator([query_name])
        taxids = name_to_taxid.get(query_name, [])
        return int(taxids[0]) if taxids else None

    def fetch_candidates(self, species_name, fallback_rank='none'):
        self._ensure_images_index()
        candidates = list()
        for matched_rank, query_name in get_taxonomic_queries(species_name, fallback_rank=fallback_rank, ncbi=self.ncbi):
            taxid = self._get_taxid(query_name)
            if taxid is None:
                continue
            for record in self.images_index.get(taxid, []):
                candidates.append(self._candidate_from_record(record=record, matched_name=query_name, matched_rank=matched_rank))
        return candidates

    def _candidate_from_record(self, record, matched_name, matched_rank):
        source_page_url = record['image_url']
        attribution_bits = [record['attribution'], record['source_name']]
        attribution = ', '.join([bit for bit in attribution_bits if bit not in ('', None)])
        provider_quality = 0
        if record['attribution'] not in ('', None):
            provider_quality += 5
        if record['source_name'] not in ('', None):
            provider_quality += 5
        return {
            'provider': self.provider_name,
            'provider_record_id': record['record_id'],
            'matched_name': matched_name,
            'matched_rank': matched_rank,
            'license_code': normalize_license_code(
                raw_code=record['license_code_text'],
                raw_url=record['license_url'],
                attribution=attribution,
            ),
            'license_url': record['license_url'],
            'attribution': attribution,
            'source_page_url': source_page_url,
            'media_url': record['image_url'],
            'width': None,
            'height': None,
            'asset_type': 'photo',
            'provider_quality': provider_quality,
        }


def build_providers(args, sources, session=None):
    session = session or requests.Session()
    ncbi = None
    if ('phylopic' in sources) or ('ncbi' in sources) or (args.fallback_rank == 'family'):
        ncbi = LazyNCBITaxa(args=args)
    providers = dict()
    for source in sources:
        if source == 'phylopic':
            if ncbi is None:
                raise ValueError('PhyloPic lookups require an initialized taxonomy database.')
            providers[source] = PhylopicProvider(session=session, ncbi=ncbi)
        elif source == 'bioicons':
            providers[source] = BioiconsProvider(session=session, ncbi=ncbi, args=args)
        elif source == 'inaturalist':
            providers[source] = INaturalistProvider(session=session, ncbi=ncbi, args=args)
        elif source == 'wikimedia':
            providers[source] = WikimediaProvider(session=session, ncbi=ncbi, args=args)
        elif source == 'gbif':
            providers[source] = GBIFProvider(session=session, ncbi=ncbi, args=args)
        elif source == 'eol':
            providers[source] = EOLProvider(session=session, ncbi=ncbi, args=args)
        elif source == 'idigbio':
            providers[source] = IDigBioProvider(session=session, ncbi=ncbi, args=args)
        elif source == 'openverse':
            providers[source] = OpenverseProvider(session=session, ncbi=ncbi, args=args)
        elif source == 'ncbi':
            if ncbi is None:
                raise ValueError('NCBI lookups require an initialized taxonomy database.')
            providers[source] = NCBIProvider(session=session, ncbi=ncbi, args=args)
    return session, ncbi, providers


def build_download_session():
    return requests.Session()


def collect_candidates_for_species(species_name, args, sources, providers):
    provider_errors = list()
    candidates = list()
    query_cache_dir = resolve_image_query_cache_dir(args)
    ensure_directory(query_cache_dir)
    for provider_index, source in enumerate(sources):
        if source == 'ncbi':
            has_allowed_before_ncbi = any(
                license_allowed(
                    candidate.get('license_code'),
                    license_max=args.license_max,
                    allow_nd=args.allow_nd,
                )
                for candidate in candidates
            )
            if has_allowed_before_ncbi:
                continue
        provider = providers[source]
        cache_path = build_query_cache_path(
            cache_dir=query_cache_dir,
            provider=source,
            species_name=species_name,
            fallback_rank=args.fallback_rank,
        )
        try:
            provider_candidates = load_cached_provider_candidates(cache_path)
            if provider_candidates is None:
                provider_candidates = provider.fetch_candidates(species_name, fallback_rank=args.fallback_rank)
                write_cached_provider_candidates(cache_path, provider_candidates)
        except requests.RequestException as exc:
            message = '{} lookup failed for {}: {}'.format(source, species_name, exc)
            _stderr(message)
            provider_errors.append(message)
            continue
        except Exception as exc:
            message = '{} lookup failed for {}: {}'.format(source, species_name, exc)
            _stderr(message)
            provider_errors.append(message)
            continue
        for candidate in provider_candidates:
            candidate = dict(candidate)
            candidate['score'] = candidate_score(candidate, provider_index=provider_index, style=args.style)
            candidates.append(candidate)
        if should_stop_after_provider(candidates, args=args):
            break
    return dedupe_sorted_candidates(candidates), provider_errors


def ensure_directory(path):
    if path in ('', None):
        return
    os.makedirs(path, exist_ok=True)


def download_media(session, media_url, destination_path, cache_path=None):
    ensure_directory(os.path.dirname(destination_path))
    destination_path = normalize_existing_media_path(destination_path)
    if os.path.exists(destination_path) and os.path.getsize(destination_path) > 0:
        return {
            'status': 'cached',
            'destination_path': destination_path,
            'cache_path': normalize_existing_media_path(cache_path) if cache_path is not None else None,
        }
    if cache_path is not None:
        ensure_directory(os.path.dirname(cache_path))
        cache_path = normalize_existing_media_path(cache_path)
        if os.path.exists(cache_path) and os.path.getsize(cache_path) > 0:
            _, cache_ext = os.path.splitext(cache_path)
            resolved_destination_path = replace_extension(destination_path, cache_ext or infer_extension(media_url, default_ext='.bin'))
            shutil.copyfile(cache_path, resolved_destination_path)
            return {
                'status': 'cached',
                'destination_path': resolved_destination_path,
                'cache_path': cache_path,
            }
    target_path = cache_path if cache_path is not None else destination_path
    tmp_path = target_path + '.tmp'
    response = None
    try:
        response = session.get(media_url, stream=True, timeout=REQUEST_TIMEOUT, headers=HTTP_HEADERS)
        response.raise_for_status()
        chunk_iter = response.iter_content(chunk_size=1024 * 64)
        first_chunk = b''
        for chunk in chunk_iter:
            if chunk:
                first_chunk = chunk
                break
        resolved_ext = infer_extension_from_response(
            response=response,
            media_url=media_url,
            first_chunk=first_chunk,
            default_ext=infer_extension(media_url, default_ext='.bin'),
        )
        resolved_target_path = replace_extension(target_path, resolved_ext)
        resolved_destination_path = replace_extension(destination_path, resolved_ext)
        tmp_path = resolved_target_path + '.tmp'
        with open(tmp_path, 'wb') as handle:
            if first_chunk:
                handle.write(first_chunk)
            for chunk in chunk_iter:
                if chunk:
                    handle.write(chunk)
        os.replace(tmp_path, resolved_target_path)
        if cache_path is not None:
            shutil.copyfile(resolved_target_path, resolved_destination_path)
            resolved_cache_path = resolved_target_path
        else:
            resolved_cache_path = None
        return {
            'status': 'downloaded',
            'destination_path': resolved_destination_path,
            'cache_path': resolved_cache_path,
        }
    except Exception:
        if os.path.exists(tmp_path):
            os.remove(tmp_path)
        raise
    finally:
        if response is not None and hasattr(response, 'close'):
            response.close()


def default_output_paths(args):
    out_dir = os.path.realpath(args.out_dir)
    images_dir = os.path.join(out_dir, 'images')
    manifest_path = os.path.realpath(args.manifest) if args.manifest else os.path.join(out_dir, 'manifest.tsv')
    attribution_path = os.path.realpath(args.attribution) if args.attribution else os.path.join(out_dir, 'ATTRIBUTION.md')
    unmatched_path = os.path.join(out_dir, 'unmatched.tsv')
    return out_dir, images_dir, manifest_path, attribution_path, unmatched_path


def write_tsv(path, rows, fieldnames):
    ensure_directory(os.path.dirname(path))
    with open(path, 'w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter='\t', extrasaction='ignore')
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def write_attribution_markdown(path, selected_assets):
    ensure_directory(os.path.dirname(path))
    grouped = defaultdict(list)
    for asset in selected_assets:
        grouped[asset['local_path']].append(asset)

    lines = ['# Attribution', '']
    for local_path in sorted(grouped.keys()):
        assets = grouped[local_path]
        representative = assets[0]
        species_names = sorted({asset['species_name'] for asset in assets})
        lines.append('## {}'.format(', '.join(species_names)))
        lines.append('')
        lines.append('Provider: {}'.format(representative['provider']))
        lines.append('Matched taxon: {} ({})'.format(representative['matched_name'], representative['matched_rank']))
        lines.append('Creator / attribution: {}'.format(representative['attribution'] or ''))
        lines.append('License: {}'.format(representative['license_code']))
        if representative['license_url']:
            lines.append('License URL: {}'.format(representative['license_url']))
        if representative['source_page_url']:
            lines.append('Source: {}'.format(representative['source_page_url']))
        lines.append('Local file: {}'.format(local_path))
        lines.append('')

    with open(path, 'w') as handle:
        handle.write('\n'.join(lines).rstrip() + '\n')


def validate_args(args):
    if args.out_dir in (None, ''):
        raise ValueError('--out_dir must be specified.')
    if int(args.max_per_species) <= 0:
        raise ValueError('--max_per_species must be > 0.')
    max_edge = getattr(args, 'max_edge', None)
    if max_edge not in (None, 0) and int(max_edge) <= 0:
        raise ValueError('--max_edge must be > 0 when specified.')
    if (getattr(args, 'output_format', 'original') == 'jpg') and (getattr(args, 'background', 'white') == 'transparent'):
        raise ValueError('--background transparent is incompatible with --output_format jpg.')


def _close_provider_bundle(session, ncbi):
    _close_ncbi_db(ncbi)
    if session is not None:
        session.close()


class ThreadLocalProviderPool:
    def __init__(self, args, sources):
        self.args = args
        self.sources = sources
        self._local = threading.local()
        self._bundles = list()
        self._lock = threading.Lock()

    def _get_bundle(self):
        bundle = getattr(self._local, 'bundle', None)
        if bundle is not None:
            return bundle
        session, ncbi, providers = build_providers(args=self.args, sources=self.sources)
        bundle = (session, ncbi, providers)
        self._local.bundle = bundle
        with self._lock:
            self._bundles.append(bundle)
        return bundle

    def lookup_species(self, species_name):
        _, _, providers = self._get_bundle()
        return collect_candidates_for_species(
            species_name=species_name,
            args=self.args,
            sources=self.sources,
            providers=providers,
        )

    def close(self):
        with self._lock:
            bundles = list(self._bundles)
            self._bundles = list()
        for session, ncbi, _ in bundles:
            _close_provider_bundle(session=session, ncbi=ncbi)


def resolve_lookup_worker_count(args, sources, species_count):
    env_value = os.environ.get('NWKIT_IMAGE_LOOKUP_WORKERS')
    if env_value not in (None, ''):
        try:
            worker_count = int(env_value)
        except ValueError:
            worker_count = DEFAULT_LOOKUP_WORKERS
        else:
            if worker_count <= 0:
                worker_count = DEFAULT_LOOKUP_WORKERS
        return min(worker_count, species_count)
    uses_taxonomy_db = ('phylopic' in sources) or ('ncbi' in sources) or (getattr(args, 'fallback_rank', 'none') == 'family')
    default_workers = 2 if uses_taxonomy_db else max(DEFAULT_LOOKUP_WORKERS, 8)
    return min(default_workers, species_count)


def resolve_download_worker_count(species_count):
    env_value = os.environ.get('NWKIT_IMAGE_DOWNLOAD_WORKERS')
    if env_value not in (None, ''):
        try:
            worker_count = int(env_value)
        except ValueError:
            worker_count = DEFAULT_DOWNLOAD_WORKERS
        else:
            if worker_count <= 0:
                worker_count = DEFAULT_DOWNLOAD_WORKERS
        return min(worker_count, species_count)
    return min(max(DEFAULT_DOWNLOAD_WORKERS, 8), species_count)


def collect_candidates_for_species_map(species_names, args, sources):
    species_names = list(species_names)
    if len(species_names) == 0:
        return dict()
    max_workers = resolve_lookup_worker_count(args=args, sources=sources, species_count=len(species_names))
    if max_workers <= 1:
        session = None
        ncbi = None
        providers = None
        try:
            session, ncbi, providers = build_providers(args=args, sources=sources)
            return {
                species_name: collect_candidates_for_species(
                    species_name=species_name,
                    args=args,
                    sources=sources,
                    providers=providers,
                )
                for species_name in species_names
            }
        finally:
            _close_provider_bundle(session=session, ncbi=ncbi)

    lookup_pool = ThreadLocalProviderPool(args=args, sources=sources)
    future_to_species = dict()
    results = dict()
    try:
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            for species_name in species_names:
                future_to_species[executor.submit(lookup_pool.lookup_species, species_name)] = species_name
            for future in as_completed(future_to_species):
                species_name = future_to_species[future]
                results[species_name] = future.result()
    finally:
        lookup_pool.close()
    return results


class SharedMediaMaterializer:
    def __init__(self, args, out_dir, images_dir, shared_cache_dir, media_plan, session=None, session_factory=None):
        self.session = session
        self.args = args
        self.out_dir = out_dir
        self.images_dir = images_dir
        self.shared_cache_dir = shared_cache_dir
        self.media_plan = media_plan
        self.session_factory = session_factory
        self._session_local = threading.local()
        self._session_lock = threading.Lock()
        self._managed_sessions = list()
        self._lock = threading.Lock()
        self._tasks = dict()

    def _get_download_session(self):
        if self.session_factory is None:
            return self.session
        session = getattr(self._session_local, 'session', None)
        if session is None:
            session = self.session_factory()
            self._session_local.session = session
            with self._session_lock:
                self._managed_sessions.append(session)
        return session

    def close(self):
        with self._session_lock:
            sessions = list(self._managed_sessions)
            self._managed_sessions = list()
        for session in sessions:
            try:
                session.close()
            except Exception:
                pass
        if self.session is not None:
            try:
                self.session.close()
            except Exception:
                pass

    def _download_candidate(self, media_url):
        plan_entry = self.media_plan.get(media_url)
        if plan_entry is None:
            raise KeyError('Missing media plan for {}'.format(media_url))
        species_name = plan_entry['species_name']
        candidate = plan_entry['candidate']
        ext = infer_extension(media_url, default_ext='.bin')
        filename = '{}__{}__{}{}'.format(
            sanitize_filename_component(species_name.replace(' ', '_')),
            candidate['provider'],
            sanitize_filename_component(candidate['provider_record_id']),
            ext,
        )
        absolute_local_path = os.path.join(self.images_dir, filename)
        cache_path = None
        if self.shared_cache_dir is not None:
            cache_path = build_media_cache_path(
                cache_dir=self.shared_cache_dir,
                media_url=media_url,
                provider=candidate['provider'],
                provider_record_id=candidate['provider_record_id'],
            )
        download_result = download_media(
            session=self._get_download_session(),
            media_url=media_url,
            destination_path=absolute_local_path,
            cache_path=cache_path,
        )
        if isinstance(download_result, dict):
            download_status = download_result['status']
            absolute_local_path = download_result.get('destination_path', absolute_local_path)
        else:
            download_status = download_result
        absolute_local_path = postprocess_media_file(absolute_local_path, args=self.args)
        return {
            'local_path': os.path.relpath(absolute_local_path, self.out_dir),
            'download_status': download_status,
        }

    def materialize(self, species_name, candidate):
        media_url = candidate['media_url']
        plan_entry = self.media_plan.get(media_url, {'species_name': species_name, 'candidate': candidate})
        owner = False
        with self._lock:
            task = self._tasks.get(media_url)
            if task is None:
                task = {
                    'event': threading.Event(),
                    'result': None,
                    'error': None,
                }
                self._tasks[media_url] = task
                owner = True
        if owner:
            try:
                task['result'] = self._download_candidate(media_url=media_url)
            except Exception as exc:
                task['error'] = exc
                with self._lock:
                    if self._tasks.get(media_url) is task:
                        self._tasks.pop(media_url, None)
            finally:
                task['event'].set()
        else:
            task['event'].wait()
        if task['error'] is not None:
            raise task['error']
        result = dict(task['result'])
        if species_name != plan_entry['species_name']:
            result['download_status'] = 'reused'
        return result


def build_media_plan(species_names, lookup_results, args):
    media_plan = dict()
    for species_name in species_names:
        candidates, _ = lookup_results[species_name]
        for candidate in allowed_candidates_from_scored_candidates(candidates, args=args):
            media_url = candidate['media_url']
            if media_url not in media_plan:
                media_plan[media_url] = {
                    'species_name': species_name,
                    'candidate': candidate,
                }
    return media_plan


def process_species_assets(
    species_name,
    leaf_names,
    candidates,
    provider_errors,
    args,
    materializer,
):
    allowed_candidates = allowed_candidates_from_scored_candidates(candidates, args=args)
    manifest_rows = list()
    selected_assets = list()
    unmatched_rows = list()
    if not allowed_candidates:
        reason = 'no exact or fallback match found'
        if candidates and (not allowed_candidates):
            reason = 'only disallowed-license candidates found'
        elif provider_errors and (not candidates):
            reason = 'provider_error'
        details = '; '.join(provider_errors) if provider_errors else ''
        for leaf_name in leaf_names:
            unmatched_rows.append({
                'leaf_name': leaf_name,
                'species_name': species_name,
                'reason': reason,
                'details': details,
            })
        return manifest_rows, selected_assets, unmatched_rows

    selected_count = 0
    download_errors = list()
    for candidate in allowed_candidates:
        if selected_count >= int(args.max_per_species):
            break
        try:
            materialized = materializer.materialize(species_name=species_name, candidate=candidate)
        except (requests.RequestException, OSError) as exc:
            message = '{} download failed for {}: {}'.format(
                candidate['provider'],
                species_name,
                exc,
            )
            _stderr(message)
            download_errors.append(message)
            continue

        for leaf_name in leaf_names:
            manifest_row = {
                'leaf_name': leaf_name,
                'species_name': species_name,
                'provider': candidate['provider'],
                'provider_record_id': candidate['provider_record_id'],
                'matched_name': candidate['matched_name'],
                'matched_rank': candidate['matched_rank'],
                'license_code': candidate['license_code'],
                'license_url': candidate['license_url'],
                'attribution': candidate['attribution'],
                'source_page_url': candidate['source_page_url'],
                'media_url': candidate['media_url'],
                'local_path': materialized['local_path'],
                'score': '{:.1f}'.format(candidate['score']),
                'status': materialized['download_status'],
            }
            manifest_rows.append(manifest_row)
            selected_assets.append(manifest_row)
        selected_count += 1

    if selected_count == 0:
        unmatched_details = '; '.join(download_errors) if download_errors else ''
        for leaf_name in leaf_names:
            unmatched_rows.append({
                'leaf_name': leaf_name,
                'species_name': species_name,
                'reason': 'download_error',
                'details': unmatched_details,
            })
    return manifest_rows, selected_assets, unmatched_rows


def image_main(args):
    validate_args(args)
    sources = parse_sources(style=args.style, source_arg=args.source)
    out_dir, images_dir, manifest_path, attribution_path, unmatched_path = default_output_paths(args)
    shared_cache_dir = resolve_image_cache_dir(args)
    ensure_directory(out_dir)
    ensure_directory(images_dir)
    ensure_directory(shared_cache_dir)

    name_mapping = read_name_tsv(args.name_tsv) if args.name_tsv else None
    tree = read_tree(args.infile, args.format, args.quoted_node_names, quiet=True)
    leaf_to_species, unmatched_rows = extract_species_mapping(tree, name_mapping=name_mapping)
    leaf_names_by_species = defaultdict(list)
    for leaf_name, species_name in leaf_to_species.items():
        leaf_names_by_species[species_name].append(leaf_name)

    manifest_rows = list()
    selected_assets = list()
    species_names = sorted(leaf_names_by_species.keys())
    lookup_results = collect_candidates_for_species_map(species_names=species_names, args=args, sources=sources)
    materializer = None
    try:
        media_plan = build_media_plan(species_names=species_names, lookup_results=lookup_results, args=args)
        materializer = SharedMediaMaterializer(
            args=args,
            out_dir=out_dir,
            images_dir=images_dir,
            shared_cache_dir=shared_cache_dir,
            media_plan=media_plan,
            session_factory=build_download_session,
        )
        download_results = dict()
        max_download_workers = resolve_download_worker_count(species_count=len(species_names))
        if max_download_workers <= 1:
            for species_name in species_names:
                candidates, provider_errors = lookup_results[species_name]
                download_results[species_name] = process_species_assets(
                    species_name=species_name,
                    leaf_names=leaf_names_by_species[species_name],
                    candidates=candidates,
                    provider_errors=provider_errors,
                    args=args,
                    materializer=materializer,
                )
        else:
            future_to_species = dict()
            with ThreadPoolExecutor(max_workers=max_download_workers) as executor:
                for species_name in species_names:
                    candidates, provider_errors = lookup_results[species_name]
                    future_to_species[
                        executor.submit(
                            process_species_assets,
                            species_name,
                            leaf_names_by_species[species_name],
                            candidates,
                            provider_errors,
                            args,
                            materializer,
                        )
                    ] = species_name
                for future in as_completed(future_to_species):
                    download_results[future_to_species[future]] = future.result()

        for species_name in species_names:
            species_manifest_rows, species_selected_assets, species_unmatched_rows = download_results[species_name]
            manifest_rows.extend(species_manifest_rows)
            selected_assets.extend(species_selected_assets)
            unmatched_rows.extend(species_unmatched_rows)
    finally:
        if materializer is not None:
            materializer.close()

    manifest_fieldnames = [
        'leaf_name',
        'species_name',
        'provider',
        'provider_record_id',
        'matched_name',
        'matched_rank',
        'license_code',
        'license_url',
        'attribution',
        'source_page_url',
        'media_url',
        'local_path',
        'score',
        'status',
    ]
    unmatched_fieldnames = ['leaf_name', 'species_name', 'reason', 'details']

    write_tsv(manifest_path, manifest_rows, manifest_fieldnames)
    write_tsv(unmatched_path, unmatched_rows, unmatched_fieldnames)
    write_attribution_markdown(attribution_path, selected_assets)

    if unmatched_rows and args.fail_on_missing:
        raise SystemExit(1)

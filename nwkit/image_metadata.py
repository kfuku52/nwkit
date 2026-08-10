"""Pure metadata, license, search-term, and media-type helpers for image lookup."""

import hashlib
import html
import os
import re
from urllib.parse import urlparse


def _hostname_matches_allowlist(hostname, allowed_hosts):
    normalized = str(hostname or "").lower().rstrip(".")
    return any(
        normalized == allowed or normalized.endswith("." + allowed)
        for allowed in allowed_hosts
    )


LICENSE_LEVELS = {
    "public-domain": 0,
    "cc-by": 1,
    "cc-by-sa": 2,
    "cc-by-nc": 3,
    "cc-by-nc-sa": 4,
    "mit": 1,
    "bsd": 1,
}

LICENSE_OPENNESS = {
    "public-domain": 70,
    "mit": 65,
    "bsd": 63,
    "cc-by": 60,
    "cc-by-sa": 50,
    "cc-by-nd": 45,
    "cc-by-nc": 40,
    "cc-by-nc-sa": 30,
    "cc-by-nc-nd": 25,
    "unknown": 0,
    "all-rights-reserved": -50,
}

FILENAME_SANITIZE_PATTERN = re.compile(r"[^A-Za-z0-9._-]+")

MAX_FILENAME_COMPONENT_LENGTH = 64

FILENAME_COMPONENT_HASH_LENGTH = 12

BIOICONS_DESCRIPTOR_TOKENS = {
    "adult",
    "blackeyes",
    "blue",
    "brown",
    "chunky",
    "cyan",
    "darkgray",
    "early",
    "embryo",
    "embryoearly",
    "embryolate",
    "fat",
    "gender",
    "gray",
    "green",
    "head",
    "juvenile",
    "late",
    "new",
    "orange",
    "pink",
    "purple",
    "redeyes",
    "small",
    "smiling",
    "test",
    "thin",
    "white",
    "yellow",
}

BIOICONS_SPECIES_ALIASES = {
    "anopheles gambiae": ("anopheles", "mosquito"),
    "arabidopsis thaliana": ("arabidopsis_thaliana", "arabidopsis"),
    "caenorhabditis elegans": ("celegans", "c_elegans", "nematode"),
    "danio rerio": ("zebrafish",),
    "drosophila melanogaster": ("drosophila", "fruit_fly"),
    "escherichia coli": ("e_coli", "coli", "bacteria"),
    "macaca mulatta": ("rhesus_monkey", "macaque", "monkey"),
    "mus musculus": ("mouse",),
    "rattus norvegicus": ("rat",),
    "saccharomyces cerevisiae": ("budding_yeast", "yeast"),
    "schizosaccharomyces pombe": ("fission_yeast", "pombe"),
    "xenopus laevis": ("xenopus_laevis", "xenopus"),
}

MIME_TYPE_TO_EXTENSION = {
    "image/gif": ".gif",
    "image/jpeg": ".jpg",
    "image/jpg": ".jpg",
    "image/png": ".png",
    "image/svg+xml": ".svg",
    "image/tiff": ".tif",
    "image/webp": ".webp",
}

RASTER_OUTPUT_EXTENSIONS = {
    ".gif": "GIF",
    ".jpg": "JPEG",
    ".jpeg": "JPEG",
    ".png": "PNG",
    ".tif": "TIFF",
    ".tiff": "TIFF",
    ".webp": "WEBP",
}

SAFE_IMAGE_EXTENSIONS = frozenset(RASTER_OUTPUT_EXTENSIONS) | frozenset((".svg",))


def normalize_species_name(name):
    if name is None:
        return None
    normalized = str(name).strip().replace("_", " ")
    normalized = re.sub(r"\s+", " ", normalized)
    return normalized if normalized != "" else None


def sanitize_filename_component(value):
    raw_value = str(value).strip()
    normalized = FILENAME_SANITIZE_PATTERN.sub("_", raw_value)
    normalized = normalized.strip("._")
    normalized = normalized or "item"
    if (
        len(normalized) > MAX_FILENAME_COMPONENT_LENGTH
        or len(raw_value.encode("utf-8")) > MAX_FILENAME_COMPONENT_LENGTH
    ):
        digest = hashlib.sha256(raw_value.encode("utf-8")).hexdigest()[
            :FILENAME_COMPONENT_HASH_LENGTH
        ]
        prefix_length = (
            MAX_FILENAME_COMPONENT_LENGTH - FILENAME_COMPONENT_HASH_LENGTH - 1
        )
        prefix = normalized[:prefix_length].rstrip("._-") or "item"
        normalized = "{}-{}".format(prefix, digest)
    return normalized


def normalize_license_code(raw_code=None, raw_url=None, attribution=None):
    code = None
    if raw_code is not None:
        code = str(raw_code).strip().lower()
        if code in ("", "none", "null", "nan"):
            code = None

    attribution_text = str(attribution or "").strip().lower()
    combined_text = " ".join(value for value in (code, attribution_text) if value)
    if "all rights reserved" in combined_text:
        return "all-rights-reserved"
    if re.search(
        r"(?:"
        r"\b(?:not|non)[ -]+public[ -]+domain\b|"
        r"\bno\s+(?:redistribution|reuse|reproduction)\b|"
        r"\blimited(?:\s+[a-z]+){0,3}\s+use\b|"
        r"\b(?:personal|educational|editorial|research|"
        r"noncommercial|non-commercial)\s+use\s+only\b"
        r")",
        combined_text,
    ):
        return "unknown"

    if code is not None:
        if code in ("cc0", "cc-0", "pd", "pdm", "public-domain", "public domain"):
            return "public-domain"
        elif code in ("mit", "mit license"):
            return "mit"
        elif code in ("bsd", "bsd license", "bsd-2-clause", "bsd-3-clause"):
            return "bsd"
        elif code in ("by", "cc-by", "cc_by"):
            return "cc-by"
        elif code in ("by-sa", "cc-by-sa", "cc_by_sa"):
            return "cc-by-sa"
        elif code in ("by-nc", "cc-by-nc", "cc_by_nc"):
            return "cc-by-nc"
        elif code in ("by-nc-sa", "cc-by-nc-sa", "cc_by_nc_sa"):
            return "cc-by-nc-sa"
        elif code in ("by-nd", "cc-by-nd", "cc_by_nd"):
            return "cc-by-nd"
        elif code in ("by-nc-nd", "cc-by-nc-nd", "cc_by_nc_nd"):
            return "cc-by-nc-nd"
        elif code in ("all-rights-reserved", "arr", "all rights reserved"):
            return "all-rights-reserved"
        elif re.fullmatch(r"cc[- ]?0(?:\s+\d+(?:\.\d+)*)?", code) or re.fullmatch(
            r"public[ -]+domain(?:[ -]+mark)?(?:\s+\d+(?:\.\d+)*)?",
            code,
        ):
            return "public-domain"
        elif re.fullmatch(
            r"mit(?:\s+software)?\s+licen[cs]e",
            code,
        ):
            return "mit"
        elif re.fullmatch(
            r"bsd(?:[- ](?:2|3)[- ]clause)?(?:\s+licen[cs]e)?",
            code,
        ):
            return "bsd"
        elif re.search(r"\bcc[ -]by[ -]nc[ -]nd\b", code):
            return "cc-by-nc-nd"
        elif re.search(r"\bcc[ -]by[ -]nd\b", code):
            return "cc-by-nd"
        elif re.search(r"\bcc[ -]by[ -]nc[ -]sa\b", code):
            return "cc-by-nc-sa"
        elif re.search(r"\bcc[ -]by[ -]nc\b", code):
            return "cc-by-nc"
        elif re.search(r"\bcc[ -]by[ -]sa\b", code):
            return "cc-by-sa"
        elif re.search(r"\bcc[ -]by\b", code):
            return "cc-by"

    if raw_url:
        parsed_url = urlparse(str(raw_url).strip())
        hostname = str(parsed_url.hostname or "").lower().rstrip(".")
        path = parsed_url.path.lower()
        try:
            port = parsed_url.port
        except ValueError:
            port = -1
        has_canonical_origin = (
            parsed_url.scheme.lower() in ("http", "https")
            and port
            in (
                None,
                80 if parsed_url.scheme.lower() == "http" else 443,
            )
            and _hostname_matches_allowlist(
                hostname,
                ("creativecommons.org",),
            )
        )
        if has_canonical_origin:
            if re.fullmatch(
                r"/publicdomain/(?:zero|mark)/1\.0/?",
                path,
            ):
                return "public-domain"
            match = re.fullmatch(
                r"/licenses/"
                r"(by-nc-nd|by-nd|by-nc-sa|by-nc|by-sa|by)"
                r"/\d+(?:\.\d+)*/?",
                path,
            )
            if match is not None:
                return "cc-{}".format(match.group(1))
    return "unknown"


def canonical_license_url(license_code):
    mapping = {
        "public-domain": "https://creativecommons.org/publicdomain/zero/1.0/",
        "mit": "https://opensource.org/licenses/MIT",
        "bsd": "https://opensource.org/licenses/BSD-3-Clause",
        "cc-by": "https://creativecommons.org/licenses/by/4.0/",
        "cc-by-sa": "https://creativecommons.org/licenses/by-sa/4.0/",
        "cc-by-nc": "https://creativecommons.org/licenses/by-nc/4.0/",
        "cc-by-nc-sa": "https://creativecommons.org/licenses/by-nc-sa/4.0/",
        "cc-by-nd": "https://creativecommons.org/licenses/by-nd/4.0/",
        "cc-by-nc-nd": "https://creativecommons.org/licenses/by-nc-nd/4.0/",
    }
    return mapping.get(license_code, "")


def license_allowed(license_code, license_max="any", allow_nd=False):
    if license_code in ("unknown", "all-rights-reserved", None):
        return False

    if license_code.endswith("-nd") and (not allow_nd):
        return False

    if license_max == "any":
        return True

    if license_code in ("mit", "bsd"):
        return LICENSE_LEVELS[license_code] <= LICENSE_LEVELS[license_max]

    if license_code in LICENSE_LEVELS:
        return LICENSE_LEVELS[license_code] <= LICENSE_LEVELS[license_max]

    if license_code == "cc-by-nd":
        return allow_nd and (LICENSE_LEVELS["cc-by"] <= LICENSE_LEVELS[license_max])

    if license_code == "cc-by-nc-nd":
        return allow_nd and (LICENSE_LEVELS["cc-by-nc"] <= LICENSE_LEVELS[license_max])

    return False


def parse_size(text):
    if not text:
        return None, None
    parts = str(text).lower().split("x")
    if len(parts) != 2:
        return None, None
    try:
        width = int(float(parts[0]))
        height = int(float(parts[1]))
    except ValueError:
        return None, None
    return width, height


def strip_html_markup(value):
    if value in (None, ""):
        return ""
    text = re.sub(r"<[^>]+>", " ", str(value))
    return re.sub(r"\s+", " ", html.unescape(text)).strip()


def tokenize_search_terms(value):
    return re.findall(r"[a-z0-9]+", str(value or "").lower())


def canonicalize_search_terms(value):
    return "".join(tokenize_search_terms(value))


def bioicons_display_author(author_slug):
    display = str(author_slug or "").replace("_", " ").replace("-", " ")
    return re.sub(r"\s+", " ", display).strip()


def bioicons_match_quality(icon_name, query_name, matched_rank):
    icon_tokens = tokenize_search_terms(icon_name)
    query_tokens = tokenize_search_terms(query_name)
    if not icon_tokens:
        return 0

    best = 0
    if query_tokens:
        if icon_tokens == query_tokens:
            best = 90 if matched_rank == "species" else 80
        elif icon_tokens[: len(query_tokens)] == query_tokens:
            best = max(best, 70 if matched_rank == "species" else 60)
        elif (len(query_tokens) > 1) and all(
            token in icon_tokens for token in query_tokens
        ):
            best = max(best, 55)

    scientific_aliases = BIOICONS_SPECIES_ALIASES.get(
        str(query_name or "").strip().lower(), ()
    )
    for alias in scientific_aliases:
        alias_tokens = tokenize_search_terms(alias)
        if not alias_tokens:
            continue
        if icon_tokens == alias_tokens:
            best = max(best, 80)
        elif icon_tokens[: len(alias_tokens)] == alias_tokens:
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
        canonicalize_search_terms(" ".join(tokens)),
    }
    non_descriptor_tokens = [
        token for token in tokens if token not in BIOICONS_DESCRIPTOR_TOKENS
    ]
    if non_descriptor_tokens:
        keys.add(canonicalize_search_terms(" ".join(non_descriptor_tokens)))
        keys.add(non_descriptor_tokens[0])
    keys.add(tokens[0])
    return {key for key in keys if key not in ("", None)}


def wikimedia_page_mentions_query(page, query_name):
    image_info = (page.get("imageinfo") or [{}])[0]
    metadata = image_info.get("extmetadata") or {}
    title_text = str(page.get("title", "")).replace("File:", " ")
    object_name = strip_html_markup(metadata.get("ObjectName", {}).get("value", ""))
    description = strip_html_markup(
        metadata.get("ImageDescription", {}).get("value", "")
    )
    combined_text = " ".join([title_text, object_name, description]).lower()
    return normalize_species_name(query_name).lower() in combined_text


def classify_wikimedia_asset(page):
    image_info = (page.get("imageinfo") or [{}])[0]
    metadata = image_info.get("extmetadata") or {}
    descriptive_text = " ".join(
        [
            str(page.get("title", "")),
            strip_html_markup(metadata.get("ObjectName", {}).get("value", "")),
            strip_html_markup(metadata.get("ImageDescription", {}).get("value", "")),
            strip_html_markup(metadata.get("Categories", {}).get("value", "")),
        ]
    ).lower()
    silhouette_terms = (
        "silhouette",
        "outline",
        "pictogram",
        "black shape",
        "shadow profile",
    )
    illustration_terms = (
        "anatomical plate",
        "chart",
        "cladogram",
        "diagram",
        "drawing",
        "figure",
        "graph",
        "icon",
        "illustration",
        "infographic",
        "logo",
        "map",
        "micrograph montage",
        "poster",
        "schematic",
        "sequence alignment",
        "taxonomic plate",
    )
    if any(term in descriptive_text for term in silhouette_terms):
        return "silhouette"
    if any(term in descriptive_text for term in illustration_terms):
        return "illustration"
    mime_type = str(image_info.get("mime") or "").lower()
    if (
        mime_type == "image/gif"
        or infer_extension(image_info.get("url", ""), default_ext="") == ".gif"
    ):
        return "illustration"
    return "photo"


def search_text_mentions_query(text_fragments, query_name):
    query_tokens = tokenize_search_terms(query_name)
    if not query_tokens:
        return False
    combined_tokens = tokenize_search_terms(
        " ".join([str(fragment or "") for fragment in text_fragments])
    )
    if not combined_tokens:
        return False
    combined_set = set(combined_tokens)
    return all(token in combined_set for token in query_tokens)


def infer_extension(url, default_ext=".bin"):
    try:
        path = urlparse(str(url or "")).path
    except (TypeError, ValueError):
        return default_ext
    if any(ord(character) < 32 or ord(character) == 127 for character in path):
        return default_ext
    _, ext = os.path.splitext(path)
    normalized_ext = ext.lower()
    if normalized_ext in SAFE_IMAGE_EXTENSIONS:
        return normalized_ext
    return default_ext


def replace_extension(path, ext):
    root, _ = os.path.splitext(path)
    return root + ext


def infer_extension_from_content_type(content_type, default_ext=".bin"):
    normalized = str(content_type or "").split(";", 1)[0].strip().lower()
    return MIME_TYPE_TO_EXTENSION.get(normalized, default_ext)


def infer_extension_from_bytes_prefix(prefix, default_ext=".bin"):
    if prefix.startswith(b"\xff\xd8\xff"):
        return ".jpg"
    if prefix.startswith(b"\x89PNG\r\n\x1a\n"):
        return ".png"
    if prefix.startswith((b"GIF87a", b"GIF89a")):
        return ".gif"
    if prefix.startswith((b"II*\x00", b"MM\x00*")):
        return ".tif"
    if len(prefix) >= 12 and prefix[:4] == b"RIFF" and prefix[8:12] == b"WEBP":
        return ".webp"
    normalized_prefix = prefix.lstrip(b"\xef\xbb\xbf\x00\t\r\n ")
    if re.search(rb"<svg(?:\s|>)", normalized_prefix[:4096], flags=re.IGNORECASE):
        return ".svg"
    return default_ext


def infer_extension_from_response(
    response, media_url, first_chunk=b"", default_ext=".bin"
):
    response_headers = getattr(response, "headers", {}) or {}
    response_url = getattr(response, "url", media_url)
    inferred_ext = infer_extension_from_bytes_prefix(first_chunk, default_ext="")
    if inferred_ext != "":
        return inferred_ext
    inferred_ext = infer_extension_from_content_type(
        response_headers.get("Content-Type"), default_ext=""
    )
    if inferred_ext != "":
        return inferred_ext
    inferred_ext = infer_extension(response_url, default_ext="")
    if inferred_ext != "":
        return inferred_ext
    return infer_extension(media_url, default_ext=default_ext)

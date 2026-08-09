import csv
import json
import glob
import hashlib
import html
import ipaddress
import io
import math
import os
import re
import secrets
import shutil
import socket
import sqlite3
import stat
import sys
import tarfile
import tempfile
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from collections import defaultdict
from urllib.parse import quote, urljoin, urlparse

import requests
from defusedxml import ElementTree
from defusedxml.common import DefusedXmlException
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

from nwkit import __version__
from nwkit.species_parser import DEFAULT_SPECIES_REGEX
from nwkit.util import (
    acquire_exclusive_lock,
    extract_taxonomy_query,
    get_ete_ncbitaxa,
    read_tree,
    resolve_download_dir,
    validate_unique_named_leaves,
    validate_distinct_output_paths,
    warn_cleanup_failure,
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
NCBI_LEGACY_IMAGE_HOSTS = frozenset(('ncbi.nlm.nih.gov', 'www.ncbi.nlm.nih.gov'))
NCBI_LEGACY_IMAGE_PATH_PREFIX = '/Taxonomy/taxi/images/'
REQUEST_TIMEOUT = (10, 60)
DEFAULT_LOOKUP_WORKERS = 4
DEFAULT_DOWNLOAD_WORKERS = 4
IMAGE_QUERY_CACHE_VERSION = 3
LOOKUP_FALLBACK_BUFFER = 2
BIOICONS_CATALOG_CACHE_VERSION = 2
DEFAULT_QUERY_CACHE_MAX_AGE_HOURS = 168.0
DEFAULT_MAX_DOWNLOAD_BYTES = 100 * 1024 * 1024
DEFAULT_MAX_PROVIDER_RESPONSE_BYTES = 8 * 1024 * 1024
DEFAULT_MAX_DECODED_PIXELS = 40_000_000
DEFAULT_MAX_IMAGE_DIMENSION = 32_768
DEFAULT_MAX_SVG_BYTES = 10 * 1024 * 1024
NCBI_NEWTAXDUMP_MAX_BYTES = 512 * 1024 * 1024
NCBI_IMAGES_TABLE_MAX_BYTES = 512 * 1024 * 1024
MAX_SAFE_REDIRECTS = 5
REDIRECT_SAFE_HEADERS = frozenset((
    'accept',
    'accept-encoding',
    'accept-language',
    'range',
    'user-agent',
))
SEMANTIC_MASK_MAX_DIM = 256
SEMANTIC_ALPHA_THRESHOLD = 16
SEMANTIC_DIFF_THRESHOLD_FLOOR = 12
REQUEST_RETRY_TOTAL = 4
REQUEST_RETRY_BACKOFF_FACTOR = 1.0
REQUEST_RETRY_STATUS_CODES = (429, 500, 502, 503, 504)
HTTP_SESSION_POOL_SIZE = 16
HTTP_USER_AGENT = 'nwkit-image/{} (+https://github.com/kfuku52/nwkit)'.format(
    __version__
)
API_HTTP_HEADERS = {
    'Accept': 'application/json',
    'User-Agent': HTTP_USER_AGENT,
}
MEDIA_HTTP_HEADERS = {
    'Accept': 'image/*, application/octet-stream;q=0.9, */*;q=0.1',
    'User-Agent': HTTP_USER_AGENT,
}
BINARY_HTTP_HEADERS = {
    'Accept': '*/*',
    'User-Agent': HTTP_USER_AGENT,
}
# Backward-compatible name for callers that customize provider API requests.
HTTP_HEADERS = API_HTTP_HEADERS
SUPPORTED_SOURCES = ('phylopic', 'bioicons', 'inaturalist', 'wikimedia', 'gbif', 'eol', 'idigbio', 'openverse', 'ncbi')
PROVIDER_API_HOSTS = {
    'phylopic': ('api.phylopic.org',),
    'bioicons': ('api.github.com',),
    'inaturalist': ('api.inaturalist.org',),
    'wikimedia': ('commons.wikimedia.org',),
    'gbif': ('api.gbif.org',),
    'eol': ('eol.org',),
    'idigbio': ('search.idigbio.org',),
    'openverse': ('api.openverse.org',),
    'ncbi': ('ftp.ncbi.nlm.nih.gov',),
}
PROVIDER_MEDIA_HOSTS = {
    'phylopic': ('phylopic.org',),
    'bioicons': ('bioicons.com',),
    'inaturalist': (
        'inaturalist.org',
        'inaturalist-open-data.s3.amazonaws.com',
        'inaturalist-open-data.s3.us-west-2.amazonaws.com',
    ),
    'wikimedia': ('wikimedia.org', 'wikimediausercontent.org'),
    # Aggregators below intentionally index publisher-hosted media. Their
    # effective allowlist is any DNS-resolved public host over HTTPS.
    'gbif': None,
    'eol': None,
    'idigbio': None,
    'openverse': None,
    'ncbi': None,
}
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
MAX_FILENAME_COMPONENT_LENGTH = 64
FILENAME_COMPONENT_HASH_LENGTH = 12
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
PILLOW_INSTALL_HINT = (
    'Pillow is a required NWKIT dependency; reinstall with: pip install --upgrade nwkit'
)
CAIROSVG_INSTALL_HINT = (
    'Install CairoSVG and the native Cairo library; from a source checkout run: '
    'pip install -e ".[image]". See the NWKIT README for platform-specific steps.'
)
RASTER_OUTPUT_EXTENSIONS = {
    '.gif': 'GIF',
    '.jpg': 'JPEG',
    '.jpeg': 'JPEG',
    '.png': 'PNG',
    '.tif': 'TIFF',
    '.tiff': 'TIFF',
    '.webp': 'WEBP',
}
SAFE_IMAGE_EXTENSIONS = frozenset(RASTER_OUTPUT_EXTENSIONS) | frozenset(('.svg',))
BIOICONS_CATALOG_MEMORY_CACHE: dict[str, dict[str, object]] = {}
BIOICONS_CATALOG_MEMORY_CACHE_LOCK = threading.Lock()


class MediaDownloadError(RuntimeError):
    pass


def _stderr(message):
    sys.stderr.write(str(message).rstrip() + '\n')


def build_retry_config():
    retry_kwargs = dict(
        total=REQUEST_RETRY_TOTAL,
        connect=REQUEST_RETRY_TOTAL,
        read=REQUEST_RETRY_TOTAL,
        status=REQUEST_RETRY_TOTAL,
        backoff_factor=REQUEST_RETRY_BACKOFF_FACTOR,
        status_forcelist=REQUEST_RETRY_STATUS_CODES,
        respect_retry_after_header=True,
        raise_on_status=False,
    )
    try:
        return Retry(allowed_methods=frozenset(['GET']), **retry_kwargs)
    except TypeError:
        return Retry(method_whitelist=frozenset(['GET']), **retry_kwargs)


def build_http_session():
    session = requests.Session()
    retry = build_retry_config()
    adapter = HTTPAdapter(
        max_retries=retry,
        pool_connections=HTTP_SESSION_POOL_SIZE,
        pool_maxsize=HTTP_SESSION_POOL_SIZE,
    )
    session.mount('https://', adapter)
    session.mount('http://', adapter)
    return session


def _hostname_matches_allowlist(hostname, allowed_hosts):
    hostname = str(hostname or '').lower().rstrip('.')
    return any(
        hostname == allowed_host or hostname.endswith('.' + allowed_host)
        for allowed_host in (str(item).lower().rstrip('.') for item in allowed_hosts or ())
    )


def _validate_public_ip_address(value, url):
    try:
        address = ipaddress.ip_address(str(value).split('%', 1)[0])
    except ValueError as exc:
        raise MediaDownloadError("URL host has an invalid IP address: {!r}.".format(url)) from exc
    if not address.is_global:
        raise MediaDownloadError(
            "URL resolves to a non-public address ({}) and was refused: {!r}.".format(address, url)
        )


def resolve_hostname_addresses(hostname):
    """Resolve a host for SSRF checks; kept separate so callers can test deterministically."""
    return {
        entry[4][0]
        for entry in socket.getaddrinfo(hostname, 443, type=socket.SOCK_STREAM)
    }


def _resolve_public_addresses(hostname, url):
    try:
        addresses = resolve_hostname_addresses(hostname)
    except OSError as exc:
        raise MediaDownloadError(
            "Could not resolve external image URL host {!r}: {}.".format(hostname, exc)
        ) from exc
    if not addresses:
        raise MediaDownloadError(
            "External image URL host did not resolve: {!r}.".format(hostname)
        )
    for address in addresses:
        _validate_public_ip_address(address, url)
    return tuple(sorted(
        addresses,
        key=lambda address: (
            ipaddress.ip_address(str(address).split('%', 1)[0]).version,
            int(ipaddress.ip_address(str(address).split('%', 1)[0])),
        ),
    ))


def validate_external_url(url, allowed_hosts=None, resolve_dns=False):
    parsed = urlparse(str(url or ''))
    if parsed.scheme.lower() != 'https':
        raise MediaDownloadError("External image URLs must use HTTPS: {!r}.".format(url))
    if parsed.username is not None or parsed.password is not None:
        raise MediaDownloadError("External image URLs must not contain credentials: {!r}.".format(url))
    hostname = parsed.hostname
    if hostname in (None, ''):
        raise MediaDownloadError("External image URL has no hostname: {!r}.".format(url))
    try:
        port = parsed.port
    except ValueError as exc:
        raise MediaDownloadError("External image URL has an invalid port: {!r}.".format(url)) from exc
    if port not in (None, 443):
        raise MediaDownloadError("External image URL must use HTTPS port 443: {!r}.".format(url))
    if allowed_hosts and not _hostname_matches_allowlist(hostname, allowed_hosts):
        raise MediaDownloadError(
            "URL host {!r} is not allowed for this provider.".format(hostname)
        )
    try:
        literal_address = ipaddress.ip_address(hostname.split('%', 1)[0])
    except ValueError:
        literal_address = None
    if literal_address is not None:
        _validate_public_ip_address(literal_address, url)
    elif resolve_dns:
        _resolve_public_addresses(hostname, url)
    return str(url)


def normalize_ncbi_image_url(url):
    """Upgrade the legacy NCBI taxonomy-image endpoint without relaxing HTTPS."""
    parsed = urlparse(str(url or ''))
    hostname = str(parsed.hostname or '').lower().rstrip('.')
    if (
        parsed.scheme.lower() == 'http'
        and hostname in NCBI_LEGACY_IMAGE_HOSTS
        and parsed.username is None
        and parsed.password is None
        and parsed.port in (None, 80)
        and parsed.path.startswith(NCBI_LEGACY_IMAGE_PATH_PREFIX)
    ):
        return parsed._replace(scheme='https', netloc=hostname).geturl()
    return str(url)


def _safe_cross_origin_redirect_headers(headers):
    return {
        key: value
        for key, value in (headers or {}).items()
        if str(key).lower() in REDIRECT_SAFE_HEADERS
    }


def _redirect_changes_origin(source_url, destination_url):
    def origin(url):
        parsed = urlparse(str(url))
        try:
            port = parsed.port
        except ValueError:
            return None
        scheme = parsed.scheme.lower()
        if port is None:
            port = 443 if scheme == 'https' else 80 if scheme == 'http' else None
        return (
            scheme,
            str(parsed.hostname or '').lower().rstrip('.'),
            port,
        )

    return origin(source_url) != origin(destination_url)


class _PinnedHTTPSAdapter(HTTPAdapter):
    """Connect to one validated address while retaining the URL host for TLS."""

    def __init__(self, address, tls_hostname, **kwargs):
        self._nwkit_address = str(address)
        self._nwkit_tls_hostname = str(tls_hostname)
        super().__init__(**kwargs)

    def _connection_from_pinned_host(self, pool_kwargs=None):
        pool_kwargs = dict(pool_kwargs or {})
        pool_kwargs.update({
            'assert_hostname': self._nwkit_tls_hostname,
            'server_hostname': self._nwkit_tls_hostname,
        })
        return self.poolmanager.connection_from_host(
            host=self._nwkit_address,
            port=443,
            scheme='https',
            pool_kwargs=pool_kwargs,
        )

    def get_connection_with_tls_context(self, request, verify, proxies=None, cert=None):
        build_pool_key = getattr(self, 'build_connection_pool_key_attributes', None)
        pool_kwargs = {}
        if callable(build_pool_key):
            _, pool_kwargs = build_pool_key(request, verify, cert)
        return self._connection_from_pinned_host(pool_kwargs=pool_kwargs)

    def get_connection(self, url, proxies=None):
        # Compatibility with Requests versions predating
        # get_connection_with_tls_context().
        return self._connection_from_pinned_host()


def _request_pinned_address(
    session,
    method,
    url,
    address,
    suppress_sensitive_headers=False,
    **kwargs
):
    parsed = urlparse(url)
    hostname = str(parsed.hostname or '')
    host_header = '[{}]'.format(hostname) if ':' in hostname else hostname
    configured_proxies = dict(getattr(session, 'proxies', None) or {})
    if getattr(session, 'trust_env', False):
        configured_proxies.update(requests.utils.get_environ_proxies(url))
    configured_proxies.update(kwargs.get('proxies', None) or {})
    if any(configured_proxies.values()):
        raise MediaDownloadError(
            'Proxy-based external image requests are disabled because the '
            'validated destination address cannot be pinned safely.'
        )
    original_adapter = session.get_adapter(url)
    adapter = _PinnedHTTPSAdapter(
        address=address,
        tls_hostname=hostname,
        max_retries=getattr(original_adapter, 'max_retries', 0),
        pool_connections=getattr(original_adapter, '_pool_connections', 1),
        pool_maxsize=getattr(original_adapter, '_pool_maxsize', 1),
        pool_block=getattr(original_adapter, '_pool_block', False),
    )
    pinned_session = requests.Session()
    pinned_session.trust_env = False
    session_headers = session.headers
    if suppress_sensitive_headers:
        session_headers = _safe_cross_origin_redirect_headers(session_headers)
    pinned_session.headers.update(session_headers)
    if not suppress_sensitive_headers:
        pinned_session.cookies.update(session.cookies)
        pinned_session.params.update(session.params)
    pinned_session.auth = None if suppress_sensitive_headers else session.auth
    pinned_session.cert = None if suppress_sensitive_headers else session.cert
    pinned_session.verify = session.verify
    pinned_session.hooks = {
        event: list(hooks)
        for event, hooks in session.hooks.items()
    }
    pinned_session.mount('https://', adapter)
    request_headers = dict(kwargs.pop('headers', None) or {})
    if suppress_sensitive_headers:
        request_headers = _safe_cross_origin_redirect_headers(request_headers)
    request_headers['Host'] = host_header
    kwargs['headers'] = request_headers
    kwargs.pop('proxies', None)
    response = None
    try:
        response = getattr(pinned_session, method)(url, **kwargs)
        session.cookies.update(pinned_session.cookies)
    except Exception:
        if response is not None and hasattr(response, 'close'):
            response.close()
        pinned_session.close()
        raise

    original_close = response.close

    def close_response_and_session():
        try:
            original_close()
        finally:
            pinned_session.close()

    response.close = close_response_and_session
    return response


def safe_external_request(session, method, url, allowed_hosts=None, **kwargs):
    """Make an external request while validating every redirect target."""
    method = str(method).lower()
    current_url = str(url)
    resolve_dns = isinstance(session, requests.Session)
    suppress_sensitive_headers = False
    for redirect_count in range(MAX_SAFE_REDIRECTS + 1):
        validate_external_url(
            current_url,
            allowed_hosts=allowed_hosts,
            resolve_dns=False,
        )
        parsed = urlparse(current_url)
        addresses = None
        if resolve_dns:
            literal_address = None
            try:
                literal_address = ipaddress.ip_address(
                    str(parsed.hostname or '').split('%', 1)[0]
                )
            except ValueError:
                pass
            if literal_address is not None:
                _validate_public_ip_address(literal_address, current_url)
                addresses = (str(literal_address),)
            else:
                addresses = _resolve_public_addresses(parsed.hostname, current_url)
        request_method = getattr(session, method)
        if resolve_dns:
            response = None
            last_error = None
            for address in addresses:
                try:
                    response = _request_pinned_address(
                        session,
                        method,
                        current_url,
                        address,
                        suppress_sensitive_headers=suppress_sensitive_headers,
                        allow_redirects=False,
                        **kwargs
                    )
                    break
                except requests.RequestException as exc:
                    last_error = exc
            if response is None:
                raise last_error
        else:
            try:
                response = request_method(current_url, allow_redirects=False, **kwargs)
            except TypeError as exc:
                if 'allow_redirects' not in str(exc):
                    raise
                raise MediaDownloadError(
                    'External request sessions must support allow_redirects=False.'
                ) from exc
        status_code = int(getattr(response, 'status_code', 0) or 0)
        location = (getattr(response, 'headers', {}) or {}).get('Location')
        if status_code not in (301, 302, 303, 307, 308) or not location:
            response_url = getattr(response, 'url', current_url) or current_url
            validate_external_url(
                response_url,
                allowed_hosts=allowed_hosts,
                resolve_dns=resolve_dns and str(response_url) != current_url,
            )
            return response
        if hasattr(response, 'close'):
            response.close()
        if redirect_count >= MAX_SAFE_REDIRECTS:
            raise requests.TooManyRedirects(
                'External request exceeded {} redirects.'.format(MAX_SAFE_REDIRECTS)
            )
        next_url = urljoin(current_url, str(location))
        if _redirect_changes_origin(current_url, next_url):
            suppress_sensitive_headers = True
            if kwargs.get('headers') is not None:
                kwargs['headers'] = _safe_cross_origin_redirect_headers(kwargs['headers'])
            kwargs.pop('auth', None)
            kwargs.pop('cert', None)
            kwargs.pop('cookies', None)
        current_url = next_url
        kwargs.pop('params', None)
        if status_code == 303 or (status_code in (301, 302) and method == 'post'):
            method = 'get'
            kwargs.pop('json', None)
            kwargs.pop('data', None)
    raise requests.TooManyRedirects(
        'External request exceeded {} redirects.'.format(MAX_SAFE_REDIRECTS)
    )


def bounded_response_json(response, max_bytes=DEFAULT_MAX_PROVIDER_RESPONSE_BYTES):
    max_bytes = int(max_bytes)
    content_length = response_content_length(response)
    if content_length is not None and content_length > max_bytes:
        raise MediaDownloadError(
            'Provider response declares {} bytes, exceeding {}.'.format(content_length, max_bytes)
        )
    if hasattr(response, 'iter_content'):
        chunks = list()
        total = 0
        for chunk in response.iter_content(chunk_size=64 * 1024):
            if not chunk:
                continue
            total += len(chunk)
            if total > max_bytes:
                raise MediaDownloadError('Provider response exceeded {} bytes.'.format(max_bytes))
            chunks.append(chunk)
        try:
            return json.loads(b''.join(chunks).decode('utf-8'))
        except (UnicodeDecodeError, json.JSONDecodeError) as exc:
            raise MediaDownloadError('Provider returned invalid JSON: {}.'.format(exc)) from exc
    payload = response.json()
    try:
        encoded_size = len(json.dumps(payload, ensure_ascii=False).encode('utf-8'))
    except (TypeError, ValueError) as exc:
        raise MediaDownloadError('Provider returned invalid JSON: {}.'.format(exc)) from exc
    if encoded_size > max_bytes:
        raise MediaDownloadError('Provider response exceeded {} bytes.'.format(max_bytes))
    return payload


def provider_json_request(session, provider, method, url, **kwargs):
    kwargs.setdefault('stream', True)
    response = safe_external_request(
        session,
        method,
        url,
        allowed_hosts=PROVIDER_API_HOSTS.get(provider, ()),
        **kwargs
    )
    try:
        response.raise_for_status()
        return bounded_response_json(response)
    finally:
        if hasattr(response, 'close'):
            response.close()


def _close_ncbi_db(ncbi):
    if ncbi is None:
        return
    if hasattr(ncbi, 'close') and callable(ncbi.close):
        try:
            ncbi.close()
        except Exception as exc:
            warn_cleanup_failure('NCBI taxonomy database handle', exc)
        return
    db = getattr(ncbi, 'db', None)
    if db is not None:
        try:
            db.close()
        except Exception as exc:
            warn_cleanup_failure('NCBI taxonomy database handle', exc)


def normalize_species_name(name):
    if name is None:
        return None
    normalized = str(name).strip().replace('_', ' ')
    normalized = re.sub(r'\s+', ' ', normalized)
    return normalized if normalized != '' else None


def sanitize_filename_component(value):
    raw_value = str(value).strip()
    normalized = FILENAME_SANITIZE_PATTERN.sub('_', raw_value)
    normalized = normalized.strip('._')
    normalized = normalized or 'item'
    if (
        len(normalized) > MAX_FILENAME_COMPONENT_LENGTH
        or len(raw_value.encode('utf-8')) > MAX_FILENAME_COMPONENT_LENGTH
    ):
        digest = hashlib.sha256(raw_value.encode('utf-8')).hexdigest()[
            :FILENAME_COMPONENT_HASH_LENGTH
        ]
        prefix_length = (
            MAX_FILENAME_COMPONENT_LENGTH
            - FILENAME_COMPONENT_HASH_LENGTH
            - 1
        )
        prefix = normalized[:prefix_length].rstrip('._-') or 'item'
        normalized = '{}-{}'.format(prefix, digest)
    return normalized


def normalize_license_code(raw_code=None, raw_url=None, attribution=None):
    code = None
    if raw_code is not None:
        code = str(raw_code).strip().lower()
        if code in ('', 'none', 'null', 'nan'):
            code = None

    attribution_text = str(attribution or '').strip().lower()
    combined_text = ' '.join(
        value for value in (code, attribution_text) if value
    )
    if 'all rights reserved' in combined_text:
        return 'all-rights-reserved'
    if re.search(
        r'(?:'
        r'\b(?:not|non)[ -]+public[ -]+domain\b|'
        r'\bno\s+(?:redistribution|reuse|reproduction)\b|'
        r'\blimited(?:\s+[a-z]+){0,3}\s+use\b|'
        r'\b(?:personal|educational|editorial|research|'
        r'noncommercial|non-commercial)\s+use\s+only\b'
        r')',
        combined_text,
    ):
        return 'unknown'

    if code is not None:
        if code in ('cc0', 'cc-0', 'pd', 'pdm', 'public-domain', 'public domain'):
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
        elif (
            re.fullmatch(r'cc[- ]?0(?:\s+\d+(?:\.\d+)*)?', code)
            or re.fullmatch(
                r'public[ -]+domain(?:[ -]+mark)?(?:\s+\d+(?:\.\d+)*)?',
                code,
            )
        ):
            return 'public-domain'
        elif re.fullmatch(
            r'mit(?:\s+software)?\s+licen[cs]e',
            code,
        ):
            return 'mit'
        elif re.fullmatch(
            r'bsd(?:[- ](?:2|3)[- ]clause)?(?:\s+licen[cs]e)?',
            code,
        ):
            return 'bsd'
        elif re.search(r'\bcc[ -]by[ -]nc[ -]nd\b', code):
            return 'cc-by-nc-nd'
        elif re.search(r'\bcc[ -]by[ -]nd\b', code):
            return 'cc-by-nd'
        elif re.search(r'\bcc[ -]by[ -]nc[ -]sa\b', code):
            return 'cc-by-nc-sa'
        elif re.search(r'\bcc[ -]by[ -]nc\b', code):
            return 'cc-by-nc'
        elif re.search(r'\bcc[ -]by[ -]sa\b', code):
            return 'cc-by-sa'
        elif re.search(r'\bcc[ -]by\b', code):
            return 'cc-by'

    if raw_url:
        parsed_url = urlparse(str(raw_url).strip())
        hostname = str(parsed_url.hostname or '').lower().rstrip('.')
        path = parsed_url.path.lower()
        try:
            port = parsed_url.port
        except ValueError:
            port = -1
        has_canonical_origin = (
            parsed_url.scheme.lower() in ('http', 'https')
            and port in (
                None,
                80 if parsed_url.scheme.lower() == 'http' else 443,
            )
            and _hostname_matches_allowlist(
                hostname,
                ('creativecommons.org',),
            )
        )
        if has_canonical_origin:
            if re.fullmatch(
                r'/publicdomain/(?:zero|mark)/1\.0/?',
                path,
            ):
                return 'public-domain'
            match = re.fullmatch(
                r'/licenses/'
                r'(by-nc-nd|by-nd|by-nc-sa|by-nc|by-sa|by)'
                r'/\d+(?:\.\d+)*/?',
                path,
            )
            if match is not None:
                return 'cc-{}'.format(match.group(1))
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
    handle = sys.stdin if path == '-' else open(path, newline='')
    try:
        reader = csv.DictReader(handle, delimiter='\t')
        if reader.fieldnames is None:
            raise ValueError('--species-name-tsv is empty.')
        required = {'leaf_name', 'species_name'}
        missing = required.difference(reader.fieldnames)
        if missing:
            raise ValueError('--species-name-tsv must contain "leaf_name" and "species_name" columns.')
        mapping = dict()
        for row in reader:
            raw_leaf_name = row.get('leaf_name')
            leaf_name = '' if raw_leaf_name is None else str(raw_leaf_name)
            species_name = normalize_species_name(row.get('species_name'))
            if leaf_name.strip() == '':
                raise ValueError('--species-name-tsv contains an empty "leaf_name" value.')
            if leaf_name in mapping:
                raise ValueError('Duplicate values in the "leaf_name" column of --species-name-tsv are not supported.')
            if species_name is None:
                raise ValueError('--species-name-tsv contains an empty "species_name" value.')
            mapping[leaf_name] = species_name
    finally:
        if handle is not sys.stdin:
            handle.close()
    if len(mapping) == 0:
        raise ValueError('--species-name-tsv is empty.')
    return mapping


def extract_species_mapping(tree, name_mapping=None, species_regex=DEFAULT_SPECIES_REGEX, args=None):
    validate_unique_named_leaves(tree, '--infile', context=' for nwkit image')
    name_mapping = name_mapping or dict()
    leaf_names = set(tree.leaf_names())
    unknown_names = set(name_mapping.keys()).difference(leaf_names)
    if unknown_names:
        raise ValueError(
            '--species-name-tsv contains leaf names not found in --infile: {}'.format(
                ', '.join(sorted(unknown_names)),
            )
        )

    leaf_to_species = dict()
    unmatched_rows = list()
    for leaf in tree.leaves():
        species_name = name_mapping.get(leaf.name)
        if species_name is None:
            species_name = normalize_species_name(extract_taxonomy_query(
                leaf.name,
                args=args,
                species_regex=species_regex if args is None else None,
                out_delim=' ',
            ))
        if species_name is None:
            unmatched_rows.append({
                'leaf_name': leaf.name,
                'species_name': '',
                'reason': 'unparsable leaf label',
                'details': "Expected the configured species parser or a matching --species-name-tsv entry.",
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


def classify_wikimedia_asset(page):
    image_info = (page.get('imageinfo') or [{}])[0]
    metadata = image_info.get('extmetadata') or {}
    descriptive_text = ' '.join([
        str(page.get('title', '')),
        strip_html_markup(metadata.get('ObjectName', {}).get('value', '')),
        strip_html_markup(metadata.get('ImageDescription', {}).get('value', '')),
        strip_html_markup(metadata.get('Categories', {}).get('value', '')),
    ]).lower()
    silhouette_terms = (
        'silhouette', 'outline', 'pictogram', 'black shape', 'shadow profile',
    )
    illustration_terms = (
        'anatomical plate', 'chart', 'cladogram', 'diagram', 'drawing', 'figure',
        'graph', 'icon', 'illustration', 'infographic', 'logo', 'map', 'micrograph montage',
        'poster', 'schematic', 'sequence alignment', 'taxonomic plate',
    )
    if any(term in descriptive_text for term in silhouette_terms):
        return 'silhouette'
    if any(term in descriptive_text for term in illustration_terms):
        return 'illustration'
    mime_type = str(image_info.get('mime') or '').lower()
    if mime_type == 'image/gif' or infer_extension(image_info.get('url', ''), default_ext='') == '.gif':
        return 'illustration'
    return 'photo'


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
    try:
        path = urlparse(str(url or '')).path
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
    normalized_prefix = prefix.lstrip(b'\xef\xbb\xbf\x00\t\r\n ')
    if re.search(br'<svg(?:\s|>)', normalized_prefix[:4096], flags=re.IGNORECASE):
        return '.svg'
    return default_ext


def infer_extension_from_response(response, media_url, first_chunk=b'', default_ext='.bin'):
    response_headers = getattr(response, 'headers', {}) or {}
    response_url = getattr(response, 'url', media_url)
    inferred_ext = infer_extension_from_bytes_prefix(first_chunk, default_ext='')
    if inferred_ext != '':
        return inferred_ext
    inferred_ext = infer_extension_from_content_type(response_headers.get('Content-Type'), default_ext='')
    if inferred_ext != '':
        return inferred_ext
    inferred_ext = infer_extension(response_url, default_ext='')
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
            'Image validation and post-processing require Pillow. {}'.format(PILLOW_INSTALL_HINT)
        ) from exc
    return Image, ImageChops, ImageOps


def load_cairosvg_module():
    try:
        import cairosvg
    except (ImportError, OSError) as exc:
        raise RuntimeError(
            'SVG image post-processing requires CairoSVG and its native Cairo '
            'library. {}'.format(
                CAIROSVG_INSTALL_HINT
            )
        ) from exc
    return cairosvg


def validate_image_dimensions(
    width,
    height,
    max_pixels=DEFAULT_MAX_DECODED_PIXELS,
    max_dimension=DEFAULT_MAX_IMAGE_DIMENSION,
    label='Image',
):
    try:
        width = int(math.ceil(float(width)))
        height = int(math.ceil(float(height)))
    except (TypeError, ValueError, OverflowError) as exc:
        raise MediaDownloadError('{} has invalid dimensions.'.format(label)) from exc
    if width <= 0 or height <= 0:
        raise MediaDownloadError('{} has non-positive dimensions {}x{}.'.format(label, width, height))
    if width > int(max_dimension) or height > int(max_dimension):
        raise MediaDownloadError(
            '{} dimensions {}x{} exceed the maximum edge {}.'.format(
                label, width, height, int(max_dimension)
            )
        )
    if width * height > int(max_pixels):
        raise MediaDownloadError(
            '{} dimensions {}x{} exceed the maximum decoded pixel count {}.'.format(
                label, width, height, int(max_pixels)
            )
        )
    return width, height


def _svg_length_to_pixels(value):
    match = re.fullmatch(
        r'\s*([+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)\s*(px|pt|pc|in|cm|mm|q)?\s*',
        str(value or ''),
        flags=re.IGNORECASE,
    )
    if match is None:
        return None
    number = float(match.group(1))
    unit = (match.group(2) or 'px').lower()
    factors = {
        'px': 1.0,
        'pt': 96.0 / 72.0,
        'pc': 16.0,
        'in': 96.0,
        'cm': 96.0 / 2.54,
        'mm': 96.0 / 25.4,
        'q': 96.0 / 101.6,
    }
    return number * factors[unit]


def _svg_dimensions_from_root(root, enforce_limits=True):
    if root.tag.rsplit('}', 1)[-1].lower() != 'svg':
        raise MediaDownloadError('SVG image root element must be <svg>.')
    view_box = str(root.attrib.get('viewBox') or root.attrib.get('viewbox') or '').replace(',', ' ').split()
    view_width = view_height = None
    if len(view_box) == 4:
        try:
            view_width = float(view_box[2])
            view_height = float(view_box[3])
        except (TypeError, ValueError):
            view_width = view_height = None
    width = _svg_length_to_pixels(root.attrib.get('width'))
    height = _svg_length_to_pixels(root.attrib.get('height'))
    if width is None:
        width = view_width if view_width is not None else 300.0
    if height is None:
        height = view_height if view_height is not None else 150.0
    if not math.isfinite(width) or not math.isfinite(height) or width <= 0 or height <= 0:
        raise MediaDownloadError('SVG image has invalid or non-positive dimensions.')
    if enforce_limits:
        return validate_image_dimensions(width, height, label='SVG image')
    return int(math.ceil(width)), int(math.ceil(height))


def inspect_svg_dimensions(source_path, enforce_limits=True):
    try:
        size = os.path.getsize(source_path)
    except OSError as exc:
        raise MediaDownloadError('SVG image is unavailable: {}.'.format(exc)) from exc
    if size <= 0 or size > DEFAULT_MAX_SVG_BYTES:
        raise MediaDownloadError(
            'SVG image size {} bytes is outside the supported range (1-{}).'.format(
                size, DEFAULT_MAX_SVG_BYTES
            )
        )
    try:
        root = ElementTree.parse(source_path).getroot()
    except (ElementTree.ParseError, DefusedXmlException, OSError) as exc:
        raise MediaDownloadError('SVG image is not well-formed XML: {}.'.format(exc)) from exc
    return _svg_dimensions_from_root(root, enforce_limits=enforce_limits)


def validate_safe_svg(
    source_path,
    allow_oversized_dimensions=False,
    rewrite_sanitized=True,
    return_sanitized_bytes=False,
):
    """Reject active or externally loaded SVG content while retaining safe SVG output."""
    try:
        size = os.path.getsize(source_path)
        with open(source_path, 'rb') as handle:
            raw = handle.read(DEFAULT_MAX_SVG_BYTES + 1)
    except OSError as exc:
        raise MediaDownloadError('SVG image is unavailable: {}.'.format(exc)) from exc
    if size <= 0 or size > DEFAULT_MAX_SVG_BYTES or len(raw) > DEFAULT_MAX_SVG_BYTES:
        raise MediaDownloadError(
            'SVG image size {} bytes is outside the supported range (1-{}).'.format(
                size, DEFAULT_MAX_SVG_BYTES
            )
        )
    if raw.startswith((b'\xff\xfe', b'\xfe\xff', b'\x00\x00\xfe\xff')) or b'\x00' in raw:
        raise MediaDownloadError(
            'SVG image uses an unsupported UTF-16/UTF-32 or NUL-containing encoding.'
        )
    if re.search(br'<!\s*ENTITY\b', raw, flags=re.IGNORECASE):
        raise MediaDownloadError('SVG image contains a forbidden entity declaration.')
    if re.search(br'<!\s*DOCTYPE\b[^[]*\[', raw, flags=re.IGNORECASE):
        raise MediaDownloadError('SVG image contains a forbidden internal document type subset.')
    sanitized_raw = re.sub(
        br'<!\s*DOCTYPE\b(?:(?:"[^"]*")|(?:\'[^\']*\')|[^>])*>\s*',
        b'',
        raw,
        flags=re.IGNORECASE,
    )
    sanitized_raw = re.sub(
        br'<\?(?!xml(?:\s|$))[\s\S]*?\?>\s*',
        b'',
        sanitized_raw,
        flags=re.IGNORECASE,
    )
    if re.search(br'<!\s*DOCTYPE\b', sanitized_raw, flags=re.IGNORECASE):
        raise MediaDownloadError('SVG image contains an unsupported document type declaration.')
    try:
        root = ElementTree.fromstring(sanitized_raw)
    except (ElementTree.ParseError, DefusedXmlException) as exc:
        raise MediaDownloadError('SVG image is not well-formed XML: {}.'.format(exc)) from exc
    if root.tag.rsplit('}', 1)[-1].lower() != 'svg':
        raise MediaDownloadError('SVG image root element must be <svg>.')
    forbidden_tags = {
        'a', 'animate', 'animatemotion', 'animatetransform', 'audio', 'embed',
        'foreignobject', 'iframe', 'image', 'object', 'script', 'set', 'video',
    }
    css_url = re.compile(r'url\s*\(\s*([^)]+?)\s*\)', flags=re.IGNORECASE)

    def css_has_unsafe_reference(value):
        value = str(value or '')
        # CSS escapes can conceal tokens such as url(), @import, or a URL
        # scheme from text-level checks. Reject them conservatively rather
        # than attempting to implement a full CSS tokenizer.
        if '\\' in value:
            return True
        if re.search(r'(?:@import\b|expression\s*\(|(?:java|vb)script\s*:)', value, flags=re.IGNORECASE):
            return True
        # Colons separate ordinary CSS declarations (for example,
        # ``fill:red``), so only treat protocol-relative text and explicit
        # resource-bearing url() functions as URL references here.
        if '//' in value:
            return True
        for match in css_url.finditer(value):
            target = match.group(1).strip().strip('"\'').strip()
            if target and not target.startswith('#'):
                return True
        return False

    for element_count, element in enumerate(root.iter(), start=1):
        if element_count > 100_000:
            raise MediaDownloadError('SVG image contains too many XML elements.')
        local_tag = str(element.tag).rsplit('}', 1)[-1].lower()
        if local_tag in forbidden_tags:
            raise MediaDownloadError('SVG image contains forbidden <{}> content.'.format(local_tag))
        if local_tag == 'style' and css_has_unsafe_reference(element.text):
            raise MediaDownloadError('SVG image contains an external or executable CSS reference.')
        for raw_name, raw_value in element.attrib.items():
            local_name = str(raw_name).rsplit('}', 1)[-1].lower()
            value = str(raw_value or '').strip()
            if local_name.startswith('on'):
                raise MediaDownloadError('SVG image contains an event-handler attribute.')
            if local_name == 'base':
                raise MediaDownloadError('SVG image contains a forbidden base URL.')
            if local_name in ('href', 'src') and value not in ('',) and not value.startswith('#'):
                raise MediaDownloadError('SVG image contains an external resource reference.')
            if css_has_unsafe_reference(value):
                raise MediaDownloadError('SVG image contains an external or executable CSS reference.')
            if re.match(r'\s*(?:javascript|vbscript)\s*:', value, flags=re.IGNORECASE):
                raise MediaDownloadError('SVG image contains an executable URL.')
    _svg_dimensions_from_root(
        root,
        enforce_limits=not allow_oversized_dimensions,
    )
    if sanitized_raw != raw and rewrite_sanitized:
        source_mode = stat.S_IMODE(os.stat(source_path).st_mode)
        tmp_path = make_temporary_sibling_path(source_path)
        try:
            with open(tmp_path, 'wb') as handle:
                handle.write(sanitized_raw)
            os.chmod(tmp_path, source_mode)
            os.replace(tmp_path, source_path)
        except Exception:
            remove_temporary_path(tmp_path, 'temporary sanitized SVG')
            raise
    if return_sanitized_bytes:
        return sanitized_raw
    return source_path


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
    raise ValueError('Unsupported --output-format: {}'.format(output_format))


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
    raise ValueError('Unsupported --trim-shape: {}'.format(trim_shape))


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
    output_mode = output_file_mode(destination_path)
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
    tmp_path = make_temporary_sibling_path(destination_path)
    try:
        image_to_save.save(tmp_path, format=image_format, **save_kwargs)
        os.chmod(tmp_path, output_mode)
        os.replace(tmp_path, destination_path)
    except Exception:
        remove_temporary_path(tmp_path, 'temporary processed image')
        raise


def rasterize_svg_to_image(source_path, max_edge=None):
    allow_oversized = max_edge not in (None, 0)
    sanitized_svg = validate_safe_svg(
        source_path,
        allow_oversized_dimensions=allow_oversized,
        rewrite_sanitized=False,
        return_sanitized_bytes=True,
    )
    natural_width, natural_height = _svg_dimensions_from_root(
        ElementTree.fromstring(sanitized_svg),
        enforce_limits=not allow_oversized,
    )
    scale = 1.0
    if max_edge not in (None, 0):
        scale = min(1.0, float(max_edge) / float(max(natural_width, natural_height)))
    output_width = max(1, int(round(natural_width * scale)))
    output_height = max(1, int(round(natural_height * scale)))
    validate_image_dimensions(output_width, output_height, label='Rasterized SVG image')
    cairosvg = load_cairosvg_module()
    Image, _, _ = load_pillow_modules()
    png_bytes = cairosvg.svg2png(
        bytestring=sanitized_svg,
        output_width=output_width,
        output_height=output_height,
    )
    with Image.open(io.BytesIO(png_bytes)) as image:
        validate_image_dimensions(*image.size, label='Rasterized SVG image')
        rasterized = image.convert('RGBA')
        rasterized.load()
    return rasterized


def process_raster_image(source_path, args):
    Image, ImageChops, ImageOps = load_pillow_modules()
    with Image.open(source_path) as source_image:
        validate_image_dimensions(*source_image.size, label='Raster image')
        max_edge = getattr(args, 'max_edge', None)
        can_downsample_before_decode = (
            max_edge not in (None, 0)
            and getattr(args, 'trim', 'off') == 'off'
            and getattr(args, 'trim_shape', 'bbox') == 'bbox'
        )
        if can_downsample_before_decode:
            source_image.draft(source_image.mode, (int(max_edge), int(max_edge)))
        if can_downsample_before_decode:
            source_image.thumbnail(
                (int(max_edge), int(max_edge)),
                get_resampling_filter(Image),
            )
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


def validate_safe_raster(source_path):
    """Validate raster dimensions and file integrity without retaining pixels."""
    Image, _, _ = load_pillow_modules()
    try:
        with Image.open(source_path) as image:
            validate_image_dimensions(*image.size, label='Raster image')
            image.verify()
        with Image.open(source_path) as image:
            validate_image_dimensions(*image.size, label='Raster image')
            image.load()
    except MediaDownloadError:
        raise
    except Exception as exc:
        raise MediaDownloadError(
            'Raster image is truncated, corrupt, or unsupported: {}.'.format(exc)
        ) from exc
    return source_path


def process_svg_image(source_path, args):
    rasterized = rasterize_svg_to_image(source_path, max_edge=getattr(args, 'max_edge', None))
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
    source_ext = os.path.splitext(normalized_path)[1].lower()
    processing_requested = image_postprocessing_requested(args)
    if source_ext == '.svg':
        validate_safe_svg(
            normalized_path,
            allow_oversized_dimensions=(
                processing_requested
                and getattr(args, 'max_edge', None) not in (None, 0)
            ),
        )
    elif not processing_requested:
        validate_safe_raster(normalized_path)
    if not processing_requested:
        return normalized_path
    if source_ext == '.svg':
        return process_svg_image(normalized_path, args=args)
    return process_raster_image(normalized_path, args=args)


def resolve_image_cache_dir(args=None):
    download_dir = resolve_download_dir(args)
    if download_dir is not None:
        return os.path.join(download_dir, 'nwkit', 'image-cache')
    out_dir = getattr(args, 'out_dir', None) if args is not None else None
    if out_dir in (None, ''):
        return None
    return os.path.join(os.path.realpath(out_dir), '.nwkit-cache', 'image-cache')


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


def build_local_media_filename(species_name, candidate, media_url):
    ext = infer_extension(media_url, default_ext='.bin')
    digest = hashlib.sha256(str(media_url).encode('utf-8')).hexdigest()[:12]
    return '{}__{}__{}__{}{}'.format(
        sanitize_filename_component(str(species_name).replace(' ', '_')),
        sanitize_filename_component(candidate['provider']),
        sanitize_filename_component(candidate['provider_record_id']),
        digest,
        ext,
    )


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


def cache_max_age_seconds_from_args(args):
    max_age_hours = getattr(args, 'query_cache_max_age_hours', DEFAULT_QUERY_CACHE_MAX_AGE_HOURS)
    if max_age_hours in (None, 0, 0.0):
        return None
    return float(max_age_hours) * 60 * 60


def cache_timestamp_is_fresh(cached_at, max_age_seconds):
    if max_age_seconds is None:
        return True
    try:
        age_seconds = max(0.0, time.time() - float(cached_at))
    except (TypeError, ValueError):
        return False
    return age_seconds <= float(max_age_seconds)


def load_cached_provider_candidates(
    cache_path,
    minimum_fetch_limit=None,
    max_age_seconds=None,
    refresh=False,
):
    if refresh:
        return None
    try:
        with open(cache_path) as handle:
            payload = json.load(handle)
    except (FileNotFoundError, json.JSONDecodeError, OSError):
        return None
    if not isinstance(payload, dict):
        return None
    if payload.get('version') != IMAGE_QUERY_CACHE_VERSION:
        return None
    if not cache_timestamp_is_fresh(payload.get('cached_at'), max_age_seconds):
        return None
    if minimum_fetch_limit is not None:
        try:
            if int(payload.get('fetch_limit', 0)) < int(minimum_fetch_limit):
                return None
        except (TypeError, ValueError):
            return None
    candidates = payload.get('candidates')
    if not isinstance(candidates, list):
        return None
    try:
        for candidate in candidates:
            validate_candidate_media_url(candidate, resolve_dns=False)
    except (MediaDownloadError, TypeError, ValueError):
        return None
    return candidates


def validate_candidate_media_url(candidate, resolve_dns=False):
    if not isinstance(candidate, dict):
        raise MediaDownloadError('Image provider candidate must be an object.')
    provider = str(candidate.get('provider') or '')
    if provider not in SUPPORTED_SOURCES:
        raise MediaDownloadError('Image provider candidate has an unknown provider {!r}.'.format(provider))
    media_url = candidate.get('media_url')
    allowed_hosts = PROVIDER_MEDIA_HOSTS.get(provider)
    return validate_external_url(
        media_url,
        allowed_hosts=allowed_hosts,
        resolve_dns=resolve_dns,
    )


def normalize_provider_candidate(candidate, expected_provider=None):
    if not isinstance(candidate, dict):
        raise ValueError('Image provider candidate must be an object.')
    candidate = dict(candidate)
    provider = candidate.get('provider')
    if not isinstance(provider, str) or provider not in SUPPORTED_SOURCES:
        raise ValueError('Image provider candidate has an invalid provider.')
    if expected_provider is not None and provider != expected_provider:
        raise ValueError(
            'Image provider candidate source does not match its cache/provider.'
        )
    for field in (
        'provider_record_id',
        'matched_name',
        'matched_rank',
        'license_code',
        'media_url',
    ):
        value = candidate.get(field)
        if not isinstance(value, (str, int, float)) or isinstance(value, bool):
            raise ValueError(
                "Image provider candidate has an invalid '{}' value.".format(
                    field
                )
            )
        candidate[field] = str(value)
    for field in ('license_url', 'attribution', 'source_page_url', 'asset_type'):
        value = candidate.get(field)
        candidate[field] = (
            str(value)
            if isinstance(value, (str, int, float)) and not isinstance(value, bool)
            else ''
        )
    candidate['width'] = _finite_nonnegative_number(candidate.get('width'))
    candidate['height'] = _finite_nonnegative_number(candidate.get('height'))
    candidate['provider_quality'] = int(
        _finite_nonnegative_number(candidate.get('provider_quality'))
    )
    for field in ('is_primary', 'is_vector'):
        value = candidate.get(field)
        if value is None:
            continue
        if isinstance(value, bool):
            candidate[field] = value
            continue
        if isinstance(value, str) and value.strip().lower() in ('true', 'false'):
            candidate[field] = value.strip().lower() == 'true'
            continue
        raise ValueError(
            "Image provider candidate has an invalid '{}' value; "
            "expected true or false.".format(field)
        )
    validate_candidate_media_url(candidate, resolve_dns=False)
    return candidate


def make_temporary_sibling_path(destination_path):
    ensure_directory(os.path.dirname(destination_path))
    file_descriptor, tmp_path = tempfile.mkstemp(
        prefix=os.path.basename(destination_path) + '.',
        suffix='.tmp',
        dir=os.path.dirname(destination_path) or '.',
    )
    os.close(file_descriptor)
    return tmp_path


def remove_temporary_path(path, label):
    if path in (None, ''):
        return
    try:
        os.remove(path)
    except FileNotFoundError:
        pass
    except Exception as exc:
        warn_cleanup_failure('{} {}'.format(label, path), exc)


def write_cached_provider_candidates(cache_path, candidates, fetch_limit=None):
    ensure_directory(os.path.dirname(cache_path))
    payload = {
        'version': IMAGE_QUERY_CACHE_VERSION,
        'cached_at': time.time(),
        'fetch_limit': int(fetch_limit or 0),
        'candidates': candidates,
    }
    tmp_path = make_temporary_sibling_path(cache_path)
    try:
        with open(tmp_path, 'w') as handle:
            json.dump(payload, handle, sort_keys=True)
        os.replace(tmp_path, cache_path)
    except Exception:
        remove_temporary_path(tmp_path, 'temporary provider cache file')
        raise


def load_cached_bioicons_catalog(cache_path, max_age_seconds=None, refresh=False):
    if refresh:
        return None
    try:
        with open(cache_path) as handle:
            payload = json.load(handle)
    except (FileNotFoundError, json.JSONDecodeError, OSError):
        return None
    if not isinstance(payload, dict):
        return None
    if payload.get('version') != BIOICONS_CATALOG_CACHE_VERSION:
        return None
    if not cache_timestamp_is_fresh(payload.get('cached_at'), max_age_seconds):
        return None
    catalog = payload.get('catalog')
    return catalog if isinstance(catalog, list) else None


def write_cached_bioicons_catalog(cache_path, catalog):
    ensure_directory(os.path.dirname(cache_path))
    payload = {
        'version': BIOICONS_CATALOG_CACHE_VERSION,
        'cached_at': time.time(),
        'catalog': catalog,
    }
    tmp_path = make_temporary_sibling_path(cache_path)
    try:
        with open(tmp_path, 'w') as handle:
            json.dump(payload, handle, sort_keys=True)
        os.replace(tmp_path, cache_path)
    except Exception:
        remove_temporary_path(tmp_path, 'temporary Bioicons catalog cache file')
        raise


def get_style_priority(candidate, style='auto'):
    asset_type = candidate.get('asset_type')
    if style == 'silhouette':
        return 2 if asset_type == 'silhouette' else 0
    if style == 'photo':
        return 2 if asset_type == 'photo' else 0
    return 1 if asset_type in ('photo', 'silhouette') else 0


def _finite_nonnegative_number(value, default=0.0):
    try:
        number = float(value)
    except (TypeError, ValueError):
        return float(default)
    if not math.isfinite(number) or number < 0.0:
        return float(default)
    return number


def get_provider_quality_bonus(candidate):
    return int(_finite_nonnegative_number(candidate.get('provider_quality')))


def candidate_is_vector(candidate):
    explicit_value = candidate.get('is_vector')
    if explicit_value is not None:
        return bool(explicit_value)
    media_url = str(candidate.get('media_url') or '')
    return urlparse(media_url).path.lower().endswith('.svg')


def get_aspect_fit_bonus(candidate):
    width = _finite_nonnegative_number(candidate.get('width'))
    height = _finite_nonnegative_number(candidate.get('height'))
    if width <= 0 or height <= 0:
        return 0
    aspect_ratio = max(width, height) / min(width, height)
    if aspect_ratio <= 1.5:
        return 100
    if aspect_ratio <= 2:
        return 90
    if aspect_ratio <= 3:
        return 75
    if aspect_ratio <= 4:
        return 60
    if aspect_ratio <= 6:
        return 40
    if aspect_ratio <= 10:
        return 20
    return 0


def describe_candidate_selection(candidate):
    rank_reason = {
        'species': 'exact_species_match',
        'genus': 'genus_fallback',
        'family': 'family_fallback',
    }.get(candidate.get('matched_rank'), 'taxon_match')
    if candidate.get('provider') == 'phylopic':
        choice_reason = (
            'phylopic_primary_image'
            if candidate.get('is_primary')
            else 'phylopic_ranked_fallback'
        )
    else:
        choice_reason = 'provider_ranked_candidate'
    return '{};{}'.format(rank_reason, choice_reason)


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


def _file_md5(path):
    # NCBI publishes MD5 for transfer-integrity checking.
    digest = hashlib.md5(
        usedforsecurity=False
    )
    with open(path, 'rb') as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def _fetch_ncbi_archive_md5(session):
    checksum_url = NCBI_NEWTAXDUMP_URL + '.md5'
    response = safe_external_request(
        session,
        'get',
        checksum_url,
        allowed_hosts=PROVIDER_API_HOSTS['ncbi'],
        stream=True,
        timeout=REQUEST_TIMEOUT,
        headers=BINARY_HTTP_HEADERS,
    )
    try:
        response.raise_for_status()
        content_length = response_content_length(response)
        if content_length is not None and content_length > 4096:
            raise MediaDownloadError('NCBI checksum response is unexpectedly large.')
        content = bytearray()
        for chunk in response.iter_content(chunk_size=4096):
            if chunk:
                content.extend(chunk)
                if len(content) > 4096:
                    raise MediaDownloadError('NCBI checksum response is unexpectedly large.')
        match = re.search(rb'\b([0-9a-fA-F]{32})\b', bytes(content))
        if match is None:
            raise MediaDownloadError('NCBI checksum response did not contain an MD5 digest.')
        return match.group(1).decode('ascii').lower()
    finally:
        if hasattr(response, 'close'):
            response.close()


def _download_to_path(
    session,
    url,
    destination_path,
    max_bytes=NCBI_NEWTAXDUMP_MAX_BYTES,
    expected_md5=None,
):
    ensure_directory(os.path.dirname(destination_path))
    tmp_path = make_temporary_sibling_path(destination_path)
    response = None
    try:
        response = safe_external_request(
            session,
            'get',
            url,
            allowed_hosts=PROVIDER_API_HOSTS['ncbi'],
            stream=True,
            timeout=REQUEST_TIMEOUT,
            headers=BINARY_HTTP_HEADERS,
        )
        response.raise_for_status()
        content_length = response_content_length(response)
        if content_length is not None and content_length > int(max_bytes):
            raise MediaDownloadError(
                'NCBI archive declares {} bytes, exceeding {}.'.format(content_length, int(max_bytes))
            )
        downloaded_bytes = 0
        # This is a transfer-integrity check against NCBI's published digest.
        digest = hashlib.md5(
            usedforsecurity=False
        )
        with open(tmp_path, 'wb') as handle:
            for chunk in response.iter_content(chunk_size=1024 * 256):
                if chunk:
                    downloaded_bytes += len(chunk)
                    if downloaded_bytes > int(max_bytes):
                        raise MediaDownloadError(
                            'NCBI archive exceeded {} bytes.'.format(int(max_bytes))
                        )
                    digest.update(chunk)
                    handle.write(chunk)
        if downloaded_bytes == 0:
            raise MediaDownloadError('Downloaded NCBI archive is empty.')
        if expected_md5 and digest.hexdigest().lower() != str(expected_md5).lower():
            raise MediaDownloadError('Downloaded NCBI archive checksum did not match the published digest.')
        os.replace(tmp_path, destination_path)
    except Exception:
        remove_temporary_path(tmp_path, 'temporary NCBI archive download')
        raise
    finally:
        if response is not None and hasattr(response, 'close'):
            try:
                response.close()
            except Exception as exc:
                warn_cleanup_failure('NCBI archive HTTP response', exc)


def _extract_ncbi_images_table(tarball_path, destination_path):
    ensure_directory(os.path.dirname(destination_path))
    tmp_path = make_temporary_sibling_path(destination_path)
    try:
        with tarfile.open(tarball_path, 'r:gz') as handle:
            member = handle.getmember('images.dmp')
            if not member.isfile():
                raise MediaDownloadError('NCBI archive images.dmp entry is not a regular file.')
            if member.size <= 0 or member.size > NCBI_IMAGES_TABLE_MAX_BYTES:
                raise MediaDownloadError(
                    'NCBI images.dmp size {} is outside the supported range.'.format(member.size)
                )
            extracted = handle.extractfile(member)
            if extracted is None:
                raise FileNotFoundError('images.dmp not found in {}'.format(tarball_path))
            try:
                with open(tmp_path, 'wb') as out_handle:
                    copied = 0
                    while True:
                        chunk = extracted.read(1024 * 1024)
                        if not chunk:
                            break
                        copied += len(chunk)
                        if copied > NCBI_IMAGES_TABLE_MAX_BYTES:
                            raise MediaDownloadError('Extracted NCBI images.dmp exceeded its size limit.')
                        out_handle.write(chunk)
                if copied != member.size:
                    raise MediaDownloadError(
                        'Extracted NCBI images.dmp size did not match the archive metadata.'
                    )
            finally:
                extracted.close()
        os.replace(tmp_path, destination_path)
    except Exception:
        remove_temporary_path(tmp_path, 'temporary extracted NCBI images table')
        raise


def _ncbi_images_table_is_valid(path):
    try:
        size = os.path.getsize(path)
        if size <= 0 or size > NCBI_IMAGES_TABLE_MAX_BYTES:
            return False
        with open(path, 'r') as handle:
            for line in handle:
                if line.strip():
                    parse_ncbi_images_dmp_line(line)
                    return True
    except (OSError, UnicodeError, ValueError):
        return False
    return False


def _ncbi_images_database_is_valid(path):
    try:
        if os.path.getsize(path) <= 0:
            return False
        uri = 'file:{}?mode=ro'.format(
            quote(os.path.abspath(path), safe='/')
        )
        connection = sqlite3.connect(uri, uri=True)
        try:
            quick_check = connection.execute('PRAGMA quick_check').fetchone()
            if not quick_check or quick_check[0] != 'ok':
                return False
            columns = {
                row[1]
                for row in connection.execute('PRAGMA table_info(images)').fetchall()
            }
            if not {'taxid', 'record_id', 'image_url'}.issubset(columns):
                return False
            return connection.execute(
                'SELECT 1 FROM images LIMIT 1'
            ).fetchone() is not None
        finally:
            connection.close()
    except (OSError, sqlite3.DatabaseError):
        return False


def ensure_ncbi_images_table(args, session):
    cache_dir = resolve_ncbi_taxonomy_image_cache_dir(args)
    images_path = os.path.join(cache_dir, 'images.dmp')
    if _ncbi_images_table_is_valid(images_path):
        return images_path
    lock_path = os.path.join(cache_dir, '.ncbi_images.lock')
    tarball_path = os.path.join(cache_dir, 'new_taxdump.tar.gz')
    with acquire_exclusive_lock(lock_path=lock_path, lock_label='NCBI taxonomy images'):
        if _ncbi_images_table_is_valid(images_path):
            return images_path
        remove_temporary_path(images_path, 'invalid cached NCBI images table')
        try:
            expected_md5 = _fetch_ncbi_archive_md5(session)
        except (OSError, requests.RequestException, MediaDownloadError) as exc:
            expected_md5 = None
            _stderr('Warning: could not verify the published NCBI archive checksum: {}'.format(exc))
        archive_is_valid = False
        if os.path.isfile(tarball_path):
            try:
                archive_is_valid = (
                    0 < os.path.getsize(tarball_path) <= NCBI_NEWTAXDUMP_MAX_BYTES
                    and (expected_md5 is None or _file_md5(tarball_path) == expected_md5)
                )
                if archive_is_valid:
                    with tarfile.open(tarball_path, 'r:gz') as archive:
                        member = archive.getmember('images.dmp')
                        archive_is_valid = (
                            member.isfile()
                            and 0 < member.size <= NCBI_IMAGES_TABLE_MAX_BYTES
                        )
            except (OSError, KeyError, tarfile.TarError):
                archive_is_valid = False
        downloaded_archive = False
        if not archive_is_valid:
            remove_temporary_path(tarball_path, 'invalid cached NCBI archive')
            _download_to_path(
                session=session,
                url=NCBI_NEWTAXDUMP_URL,
                destination_path=tarball_path,
                max_bytes=NCBI_NEWTAXDUMP_MAX_BYTES,
                expected_md5=expected_md5,
            )
            downloaded_archive = True
        try:
            _extract_ncbi_images_table(tarball_path=tarball_path, destination_path=images_path)
            if not _ncbi_images_table_is_valid(images_path):
                raise MediaDownloadError('Extracted NCBI images.dmp failed validation.')
        except Exception:
            remove_temporary_path(tarball_path, 'invalid NCBI archive')
            remove_temporary_path(images_path, 'invalid extracted NCBI images table')
            if downloaded_archive:
                raise
            _download_to_path(
                session=session,
                url=NCBI_NEWTAXDUMP_URL,
                destination_path=tarball_path,
                max_bytes=NCBI_NEWTAXDUMP_MAX_BYTES,
                expected_md5=expected_md5,
            )
            try:
                _extract_ncbi_images_table(
                    tarball_path=tarball_path,
                    destination_path=images_path,
                )
                if not _ncbi_images_table_is_valid(images_path):
                    raise MediaDownloadError('Extracted NCBI images.dmp failed validation.')
            except Exception:
                remove_temporary_path(tarball_path, 'invalid NCBI archive')
                remove_temporary_path(images_path, 'invalid extracted NCBI images table')
                raise
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


def build_ncbi_images_database(images_path, database_path):
    ensure_directory(os.path.dirname(database_path))
    tmp_path = make_temporary_sibling_path(database_path)
    connection = None
    try:
        connection = sqlite3.connect(tmp_path)
        connection.execute('PRAGMA journal_mode=OFF')
        connection.execute('PRAGMA synchronous=OFF')
        connection.execute(
            'CREATE TABLE images ('
            'taxid INTEGER NOT NULL, record_id TEXT NOT NULL, title TEXT, image_url TEXT NOT NULL, '
            'license_code_text TEXT, license_url TEXT, attribution TEXT, source_name TEXT)'
        )
        insert_sql = (
            'INSERT INTO images '
            '(taxid, record_id, title, image_url, license_code_text, license_url, attribution, source_name) '
            'VALUES (?, ?, ?, ?, ?, ?, ?, ?)'
        )
        pending_rows = list()
        with open(images_path, 'r') as handle:
            for raw_line in handle:
                if raw_line.strip() == '':
                    continue
                record = parse_ncbi_images_dmp_line(raw_line)
                for taxid in record['taxids']:
                    pending_rows.append((
                        taxid,
                        record['record_id'],
                        record['title'],
                        record['image_url'],
                        record['license_code_text'],
                        record['license_url'],
                        record['attribution'],
                        record['source_name'],
                    ))
                if len(pending_rows) >= 1000:
                    connection.executemany(insert_sql, pending_rows)
                    pending_rows = list()
        if pending_rows:
            connection.executemany(insert_sql, pending_rows)
        connection.execute('CREATE INDEX images_taxid_idx ON images (taxid)')
        connection.commit()
        connection.close()
        connection = None
        os.replace(tmp_path, database_path)
    except Exception:
        if connection is not None:
            connection.close()
        remove_temporary_path(tmp_path, 'temporary NCBI images database')
        raise


def ensure_ncbi_images_database(args, session):
    cache_dir = resolve_ncbi_taxonomy_image_cache_dir(args)
    database_path = os.path.join(cache_dir, 'images.sqlite3')
    if _ncbi_images_database_is_valid(database_path):
        return database_path
    images_path = ensure_ncbi_images_table(args=args, session=session)
    lock_path = os.path.join(cache_dir, '.ncbi_images_database.lock')
    with acquire_exclusive_lock(lock_path=lock_path, lock_label='NCBI taxonomy images database'):
        if _ncbi_images_database_is_valid(database_path):
            return database_path
        remove_temporary_path(database_path, 'invalid cached NCBI images database')
        for attempt in range(2):
            try:
                build_ncbi_images_database(
                    images_path=images_path,
                    database_path=database_path,
                )
            except (UnicodeError, ValueError):
                remove_temporary_path(
                    database_path,
                    'invalid rebuilt NCBI images database',
                )
                if attempt > 0:
                    raise
            else:
                if _ncbi_images_database_is_valid(database_path):
                    return database_path
                remove_temporary_path(
                    database_path,
                    'invalid rebuilt NCBI images database',
                )
                if attempt > 0:
                    raise MediaDownloadError(
                        'Rebuilt NCBI images database failed validation.'
                    )
            remove_temporary_path(
                images_path,
                'invalid cached NCBI images table',
            )
            images_path = ensure_ncbi_images_table(
                args=args,
                session=session,
            )
    return database_path


def load_ncbi_images_records(database_path, taxid):
    connection = sqlite3.connect(database_path)
    try:
        rows = connection.execute(
            'SELECT record_id, title, image_url, license_code_text, license_url, attribution, source_name '
            'FROM images WHERE taxid = ? ORDER BY rowid',
            (int(taxid),),
        ).fetchall()
    finally:
        connection.close()
    return [
        {
            'record_id': row[0],
            'title': row[1],
            'image_url': row[2],
            'license_code_text': row[3],
            'license_url': row[4],
            'attribution': row[5],
            'source_name': row[6],
        }
        for row in rows
    ]


def candidate_score(candidate, provider_index=0, style='auto'):
    matched_rank = candidate.get('matched_rank')
    if not isinstance(matched_rank, str):
        matched_rank = ''
    license_code = candidate.get('license_code')
    if not isinstance(license_code, str):
        license_code = ''
    rank_priority = {'species': 3, 'genus': 2, 'family': 1}.get(matched_rank, 0)
    provider_priority = max(
        1,
        100 - int(_finite_nonnegative_number(provider_index)),
    )
    primary_priority = 1 if candidate.get('is_primary') else 0
    license_priority = max(0, LICENSE_OPENNESS.get(license_code, 0))
    style_priority = get_style_priority(candidate, style=style)
    vector_priority = 1 if candidate_is_vector(candidate) else 0
    aspect_priority = get_aspect_fit_bonus(candidate)
    quality_priority = min(
        max(
            _finite_nonnegative_number(candidate.get('width')),
            _finite_nonnegative_number(candidate.get('height')),
        ),
        10000,
    )
    provider_quality = max(0, min(get_provider_quality_bonus(candidate), 99))
    score = (
        rank_priority * 10**15
        + provider_priority * 10**12
        + primary_priority * 10**11
        + license_priority * 10**9
        + style_priority * 10**8
        + vector_priority * 10**7
        + aspect_priority * 10**4
        + quality_priority * 9
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
        response = safe_external_request(
            self.session,
            'get',
            f'{PHYLIPIC_API_ROOT}/resolve/ncbi.nlm.nih.gov/taxid/{taxid}',
            allowed_hosts=PROVIDER_API_HOSTS[self.provider_name],
            stream=True,
            timeout=REQUEST_TIMEOUT,
            headers=HTTP_HEADERS,
        )
        try:
            if response.status_code == 404:
                return None
            response.raise_for_status()
            return bounded_response_json(response)
        finally:
            if hasattr(response, 'close'):
                response.close()

    def _fetch_node_images(self, node_uuid, build):
        payload = provider_json_request(
            self.session,
            self.provider_name,
            'get',
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
        return payload.get('_embedded', {}).get('items', [])

    @staticmethod
    def _primary_image_link(node):
        primary_image = node.get('_links', {}).get('primaryImage')
        return primary_image if isinstance(primary_image, dict) else {}

    @classmethod
    def _primary_image_uuid(cls, node):
        href = cls._primary_image_link(node).get('href', '')
        path_parts = [part for part in urlparse(href).path.split('/') if part]
        if len(path_parts) >= 2 and path_parts[-2] == 'images':
            return path_parts[-1]
        return ''

    def _fetch_linked_image(self, href):
        if str(href).startswith(('http://', 'https://')):
            url = href
        else:
            url = '{}{}'.format(PHYLIPIC_API_ROOT, href)
        return provider_json_request(
            self.session,
            self.provider_name,
            'get',
            url,
            timeout=REQUEST_TIMEOUT,
            headers=HTTP_HEADERS,
        )

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
            image_items = self._fetch_node_images(node_uuid=node_uuid, build=build)
            primary_uuid = self._primary_image_uuid(node)
            if primary_uuid and not any(item.get('uuid') == primary_uuid for item in image_items):
                primary_href = self._primary_image_link(node).get('href')
                try:
                    primary_item = self._fetch_linked_image(primary_href)
                except requests.RequestException as exc:
                    _stderr(
                        'PhyloPic primary image lookup failed for {}: {}; '
                        'using ranked fallback candidates.'.format(query_name, exc)
                    )
                else:
                    image_items = [primary_item] + image_items
            for image_item in image_items:
                candidate = self._candidate_from_image(
                    image_item=image_item,
                    matched_name=query_name,
                    matched_rank=matched_rank,
                    is_primary=bool(primary_uuid and image_item.get('uuid') == primary_uuid),
                )
                if candidate is not None:
                    candidates.append(candidate)
        return candidates

    def _candidate_from_image(self, image_item, matched_name, matched_rank, is_primary=False):
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
        selected_type = str(selected_file.get('type') or '').lower()
        selected_href = str(selected_file.get('href') or '')
        is_vector = (
            vector_file is not None
            or selected_type == 'image/svg+xml'
            or urlparse(selected_href).path.lower().endswith('.svg')
        )

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
            'is_primary': bool(is_primary),
            'is_vector': is_vector,
            'provider_quality': 30 if is_vector else 10,
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
        refresh_cache = bool(getattr(self.args, 'refresh_cache', False))
        refresh_token = self.args if refresh_cache else None
        max_age_seconds = cache_max_age_seconds_from_args(self.args)

        def load_memory_catalog():
            with BIOICONS_CATALOG_MEMORY_CACHE_LOCK:
                memory_entry = BIOICONS_CATALOG_MEMORY_CACHE.get(cache_path)
            if memory_entry is None:
                return None
            if refresh_cache:
                if memory_entry.get('refresh_token') is refresh_token:
                    return memory_entry['catalog']
                return None
            if cache_timestamp_is_fresh(memory_entry.get('cached_at'), max_age_seconds):
                return memory_entry['catalog']
            return None

        def remember_catalog(catalog):
            with BIOICONS_CATALOG_MEMORY_CACHE_LOCK:
                BIOICONS_CATALOG_MEMORY_CACHE[cache_path] = {
                    'cached_at': time.time(),
                    'catalog': catalog,
                    'refresh_token': refresh_token,
                }

        catalog = load_memory_catalog()
        if catalog is not None:
            self._catalog = catalog
            return self._catalog

        catalog = load_cached_bioicons_catalog(
            cache_path,
            max_age_seconds=max_age_seconds,
            refresh=refresh_cache,
        )
        if catalog is None:
            lock_path = cache_path + '.lock'
            with acquire_exclusive_lock(lock_path=lock_path, lock_label='Bioicons catalog'):
                catalog = load_memory_catalog()
                if catalog is None:
                    catalog = load_cached_bioicons_catalog(
                        cache_path,
                        max_age_seconds=max_age_seconds,
                        refresh=refresh_cache,
                    )
                    if catalog is None:
                        payload = provider_json_request(
                            self.session,
                            self.provider_name,
                            'get',
                            BIOICONS_GITHUB_API_ROOT,
                            params={'recursive': 1},
                            timeout=REQUEST_TIMEOUT,
                            headers=HTTP_HEADERS,
                        )
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
                remember_catalog(catalog)
        else:
            remember_catalog(catalog)
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
        payload = provider_json_request(
            self.session,
            self.provider_name,
            'get',
            f'{INATURALIST_API_ROOT}/taxa',
            params={'q': query_name, 'per_page': self.result_limit},
            timeout=REQUEST_TIMEOUT,
            headers=HTTP_HEADERS,
        )
        results = payload.get('results', [])
        expected_name = query_name.lower()
        matches = [
            item for item in results
            if str(item.get('name', '')).lower() == expected_name and item.get('rank') == expected_rank
        ]
        if matches:
            return matches[0]
        return None

    def _fetch_observations(self, taxon_id, per_page):
        payload = provider_json_request(
            self.session,
            self.provider_name,
            'get',
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
        return payload.get('results', [])

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
        payload = provider_json_request(
            self.session,
            self.provider_name,
            'get',
            '{}/search/1.0.json'.format(EOL_API_ROOT),
            params={
                'q': query_name,
                'page': 1,
                'exact': 'true',
            },
            timeout=REQUEST_TIMEOUT,
            headers=HTTP_HEADERS,
        )
        return payload.get('results', [])

    def _fetch_page(self, page_id):
        return provider_json_request(
            self.session,
            self.provider_name,
            'get',
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
        payload = provider_json_request(
            self.session,
            self.provider_name,
            'post',
            '{}/search/media'.format(IDIGBIO_API_ROOT),
            json=body,
            timeout=REQUEST_TIMEOUT,
            headers=HTTP_HEADERS,
        )
        return payload.get('items', [])

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
        payload = provider_json_request(
            self.session,
            self.provider_name,
            'get',
            '{}/images/'.format(OPENVERSE_API_ROOT),
            params={
                'q': query_name,
                'page_size': int(self.result_limit if page_size is None else page_size),
            },
            timeout=REQUEST_TIMEOUT,
            headers=HTTP_HEADERS,
        )
        return payload.get('results', [])

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
        payload = provider_json_request(
            self.session,
            self.provider_name,
            'get',
            WIKIMEDIA_API_ROOT,
            params={
                'action': 'query',
                'generator': 'search',
                'gsrsearch': '"{}" filetype:bitmap'.format(query_name),
                'gsrnamespace': 6,
                'gsrlimit': self.result_limit if limit == 10 else limit,
                'prop': 'imageinfo',
                'iiprop': 'url|size|mime|mediatype|extmetadata',
                'iiurlwidth': 1600,
                'format': 'json',
            },
            timeout=REQUEST_TIMEOUT,
            headers=HTTP_HEADERS,
        )
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
        asset_type = classify_wikimedia_asset(page)
        if asset_type == 'illustration':
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
            'asset_type': asset_type,
            'provider_quality': provider_quality,
        }


class GBIFProvider:
    provider_name = 'gbif'

    def __init__(self, session, ncbi=None, args=None):
        self.session = session
        self.ncbi = ncbi
        self.result_limit = resolve_provider_fetch_limit(args=args, minimum=10, maximum=30)

    def _match_taxon(self, query_name, matched_rank):
        payload = provider_json_request(
            self.session,
            self.provider_name,
            'get',
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
        payload = provider_json_request(
            self.session,
            self.provider_name,
            'get',
            f'{GBIF_API_ROOT}/occurrence/search',
            params={
                'taxon_key': usage_key,
                'media_type': 'StillImage',
                'limit': self.result_limit if limit is None else limit,
            },
            timeout=REQUEST_TIMEOUT,
            headers=HTTP_HEADERS,
        )
        return payload.get('results', [])

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
        self.images_database_path = None

    def _ensure_images_database(self):
        if self.images_database_path is not None:
            return
        self.images_database_path = ensure_ncbi_images_database(args=self.args, session=self.session)

    def _get_taxid(self, query_name):
        name_to_taxid = self.ncbi.get_name_translator([query_name])
        taxids = name_to_taxid.get(query_name, [])
        return int(taxids[0]) if taxids else None

    def fetch_candidates(self, species_name, fallback_rank='none'):
        self._ensure_images_database()
        candidates = list()
        for matched_rank, query_name in get_taxonomic_queries(species_name, fallback_rank=fallback_rank, ncbi=self.ncbi):
            taxid = self._get_taxid(query_name)
            if taxid is None:
                continue
            for record in load_ncbi_images_records(self.images_database_path, taxid):
                candidates.append(self._candidate_from_record(record=record, matched_name=query_name, matched_rank=matched_rank))
        return candidates

    def _candidate_from_record(self, record, matched_name, matched_rank):
        media_url = normalize_ncbi_image_url(record['image_url'])
        source_page_url = media_url
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
            'media_url': media_url,
            'width': None,
            'height': None,
            'asset_type': 'photo',
            'provider_quality': provider_quality,
        }


def build_providers(args, sources, session=None):
    session = session or build_http_session()
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
    return build_http_session()


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
        fetch_limit = int(getattr(provider, 'result_limit', resolve_provider_fetch_limit(args=args)))
        cache_path = build_query_cache_path(
            cache_dir=query_cache_dir,
            provider=source,
            species_name=species_name,
            fallback_rank=args.fallback_rank,
        )
        fetched_candidates = False
        try:
            provider_candidates = load_cached_provider_candidates(
                cache_path,
                minimum_fetch_limit=fetch_limit,
                max_age_seconds=cache_max_age_seconds_from_args(args),
                refresh=bool(getattr(args, 'refresh_cache', False)),
            )
            if provider_candidates is None:
                provider_candidates = provider.fetch_candidates(species_name, fallback_rank=args.fallback_rank)
                fetched_candidates = True
            if not isinstance(provider_candidates, list):
                raise ValueError('Image provider candidates must be a list.')
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
        normalized_candidates = list()
        for candidate_index, candidate in enumerate(provider_candidates, start=1):
            try:
                candidate = normalize_provider_candidate(
                    candidate,
                    expected_provider=source,
                )
                candidate['score'] = candidate_score(
                    candidate,
                    provider_index=provider_index,
                    style=args.style,
                )
            except Exception as exc:
                message = (
                    '{} candidate {} was skipped for {}: {}'.format(
                        source,
                        candidate_index,
                        species_name,
                        exc,
                    )
                )
                _stderr(message)
                provider_errors.append(message)
                continue
            normalized_candidates.append(candidate)
            candidates.append(candidate)
        if fetched_candidates:
            try:
                write_cached_provider_candidates(
                    cache_path,
                    [
                        {
                            key: value
                            for key, value in candidate.items()
                            if key != 'score'
                        }
                        for candidate in normalized_candidates
                    ],
                    fetch_limit=fetch_limit,
                )
            except Exception as exc:
                message = '{} cache write failed for {}: {}'.format(
                    source,
                    species_name,
                    exc,
                )
                _stderr(message)
                provider_errors.append(message)
        if should_stop_after_provider(candidates, args=args):
            break
    return dedupe_sorted_candidates(candidates), provider_errors


def ensure_directory(path):
    if path in ('', None):
        return
    os.makedirs(path, exist_ok=True)


def validate_media_file(path, max_download_bytes=DEFAULT_MAX_DOWNLOAD_BYTES):
    try:
        size = os.path.getsize(path)
    except OSError as exc:
        raise MediaDownloadError('Downloaded media is unavailable: {}'.format(exc)) from exc
    if size <= 0:
        raise MediaDownloadError('Downloaded media is empty.')
    if size > int(max_download_bytes):
        raise MediaDownloadError(
            'Downloaded media is {} bytes, exceeding --max-download-bytes {}.'.format(
                size,
                int(max_download_bytes),
            )
        )
    with open(path, 'rb') as handle:
        prefix = handle.read(4096)
    extension = infer_extension_from_bytes_prefix(prefix, default_ext='')
    if extension == '':
        raise MediaDownloadError('Downloaded response is not a recognized JPEG, PNG, GIF, TIFF, WebP, or SVG image.')
    return extension


def normalize_valid_media_path(path, max_download_bytes=DEFAULT_MAX_DOWNLOAD_BYTES):
    extension = validate_media_file(path, max_download_bytes=max_download_bytes)
    resolved_path = replace_extension(path, extension)
    if resolved_path != path:
        os.replace(path, resolved_path)
    return resolved_path


def response_content_length(response):
    headers = getattr(response, 'headers', {}) or {}
    raw_length = headers.get('Content-Length')
    if raw_length in (None, ''):
        return None
    try:
        return int(raw_length)
    except (TypeError, ValueError):
        return None


def response_content_type_is_plausible_image(response):
    headers = getattr(response, 'headers', {}) or {}
    content_type = str(headers.get('Content-Type') or '').split(';', 1)[0].strip().lower()
    return (
        content_type == ''
        or content_type.startswith('image/')
        or content_type in ('application/octet-stream', 'binary/octet-stream')
    )


def atomic_copyfile(source_path, destination_path):
    ensure_directory(os.path.dirname(destination_path))
    output_mode = output_file_mode(destination_path)
    tmp_path = make_temporary_sibling_path(destination_path)
    try:
        shutil.copyfile(source_path, tmp_path)
        os.chmod(tmp_path, output_mode)
        os.replace(tmp_path, destination_path)
    except Exception:
        remove_temporary_path(tmp_path, 'temporary copied media file')
        raise


def output_file_mode(destination_path):
    try:
        return stat.S_IMODE(os.stat(destination_path).st_mode)
    except FileNotFoundError:
        pass
    directory = os.path.dirname(os.path.abspath(destination_path))
    ensure_directory(directory)
    for _ in range(100):
        probe_path = os.path.join(
            directory,
            '.nwkit-mode-probe-{}'.format(secrets.token_hex(16)),
        )
        try:
            fd = os.open(
                probe_path,
                os.O_CREAT | os.O_EXCL | os.O_WRONLY,
                0o666,
            )
        except FileExistsError:
            continue
        try:
            return stat.S_IMODE(os.fstat(fd).st_mode)
        finally:
            os.close(fd)
            os.remove(probe_path)
    raise FileExistsError('Could not allocate a temporary output-mode probe.')


def download_media(
    session,
    media_url,
    destination_path,
    cache_path=None,
    max_download_bytes=DEFAULT_MAX_DOWNLOAD_BYTES,
    provider=None,
    reuse_destination=True,
):
    max_download_bytes = int(max_download_bytes)
    if provider not in (None, '') and provider not in SUPPORTED_SOURCES:
        raise MediaDownloadError('Unknown media provider {!r}.'.format(provider))
    allowed_hosts = PROVIDER_MEDIA_HOSTS.get(provider)
    validate_external_url(media_url, allowed_hosts=allowed_hosts, resolve_dns=False)
    ensure_directory(os.path.dirname(destination_path))
    if reuse_destination:
        destination_path = normalize_existing_media_path(destination_path)
        if os.path.exists(destination_path) and os.path.getsize(destination_path) > 0:
            try:
                destination_path = normalize_valid_media_path(
                    destination_path,
                    max_download_bytes=max_download_bytes,
                )
            except MediaDownloadError:
                remove_temporary_path(destination_path, 'invalid existing media file')
            else:
                return {
                    'status': 'cached',
                    'destination_path': destination_path,
                    'cache_path': normalize_existing_media_path(cache_path) if cache_path is not None else None,
                }
    if cache_path is not None:
        ensure_directory(os.path.dirname(cache_path))
        cache_path = normalize_existing_media_path(cache_path)
        if os.path.exists(cache_path) and os.path.getsize(cache_path) > 0:
            try:
                cache_path = normalize_valid_media_path(
                    cache_path,
                    max_download_bytes=max_download_bytes,
                )
            except MediaDownloadError:
                remove_temporary_path(cache_path, 'invalid cached media file')
            else:
                _, cache_ext = os.path.splitext(cache_path)
                resolved_destination_path = replace_extension(destination_path, cache_ext)
                atomic_copyfile(cache_path, resolved_destination_path)
                return {
                    'status': 'cached',
                    'destination_path': resolved_destination_path,
                    'cache_path': cache_path,
                }
    target_path = cache_path if cache_path is not None else destination_path
    tmp_path = make_temporary_sibling_path(target_path)
    response = None
    try:
        response = safe_external_request(
            session,
            'get',
            media_url,
            allowed_hosts=allowed_hosts,
            stream=True,
            timeout=REQUEST_TIMEOUT,
            headers=MEDIA_HTTP_HEADERS,
        )
        response.raise_for_status()
        if not response_content_type_is_plausible_image(response):
            content_type = (getattr(response, 'headers', {}) or {}).get('Content-Type', '')
            raise MediaDownloadError('Media URL returned non-image Content-Type {!r}.'.format(content_type))
        content_length = response_content_length(response)
        if content_length is not None and content_length > max_download_bytes:
            raise MediaDownloadError(
                'Media response declares {} bytes, exceeding --max-download-bytes {}.'.format(
                    content_length,
                    max_download_bytes,
                )
            )
        downloaded_bytes = 0
        with open(tmp_path, 'wb') as handle:
            for chunk in response.iter_content(chunk_size=1024 * 64):
                if chunk:
                    downloaded_bytes += len(chunk)
                    if downloaded_bytes > max_download_bytes:
                        raise MediaDownloadError(
                            'Media response exceeded --max-download-bytes {}.'.format(max_download_bytes)
                        )
                    handle.write(chunk)
        resolved_ext = validate_media_file(tmp_path, max_download_bytes=max_download_bytes)
        resolved_target_path = replace_extension(target_path, resolved_ext)
        resolved_destination_path = replace_extension(destination_path, resolved_ext)
        if cache_path is None:
            os.chmod(tmp_path, output_file_mode(resolved_target_path))
        os.replace(tmp_path, resolved_target_path)
        if cache_path is not None:
            atomic_copyfile(resolved_target_path, resolved_destination_path)
            resolved_cache_path = resolved_target_path
        else:
            resolved_cache_path = None
        return {
            'status': 'downloaded',
            'destination_path': resolved_destination_path,
            'cache_path': resolved_cache_path,
        }
    except Exception:
        remove_temporary_path(tmp_path, 'temporary media download file')
        raise
    finally:
        if response is not None and hasattr(response, 'close'):
            try:
                response.close()
            except Exception as exc:
                warn_cleanup_failure('HTTP response', exc)


def default_output_paths(args):
    out_dir = os.path.realpath(args.out_dir)
    images_dir = os.path.join(out_dir, 'images')
    manifest_out = getattr(args, 'manifest_out', None) or getattr(args, 'manifest', None)
    attribution_out = getattr(args, 'attribution_out', None) or getattr(args, 'attribution', None)
    manifest_path = os.path.realpath(manifest_out) if manifest_out else os.path.join(out_dir, 'manifest.tsv')
    attribution_path = os.path.realpath(attribution_out) if attribution_out else os.path.join(out_dir, 'ATTRIBUTION.md')
    unmatched_path = os.path.join(out_dir, 'unmatched.tsv')
    return out_dir, images_dir, manifest_path, attribution_path, unmatched_path


def write_tsv(path, rows, fieldnames):
    ensure_directory(os.path.dirname(path))
    with open(path, 'w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter='\t', extrasaction='ignore')
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def rebase_local_paths_for_output(rows, out_dir, output_path):
    output_dir = os.path.dirname(os.path.realpath(output_path))
    rebased_rows = list()
    for row in rows:
        rebased_row = dict(row)
        local_path = row.get('local_path')
        if local_path not in (None, ''):
            absolute_path = (
                os.path.abspath(local_path)
                if os.path.isabs(local_path)
                else os.path.abspath(os.path.join(out_dir, local_path))
            )
            rebased_row['local_path'] = os.path.relpath(absolute_path, output_dir)
        rebased_rows.append(rebased_row)
    return rebased_rows


def _safe_markdown_text(value, escape_markdown=True):
    text = ''.join(
        ' ' if (ord(character) < 32 or ord(character) == 127) else character
        for character in str(value or '')
    )
    text = ' '.join(text.split())
    text = html.escape(text, quote=True)
    if escape_markdown:
        text = re.sub(r'([\\`*_\[\]\(\)!|])', r'\\\1', text)
    return text


def write_attribution_markdown(path, selected_assets):
    ensure_directory(os.path.dirname(path))
    grouped = defaultdict(list)
    for asset in selected_assets:
        grouped[asset['local_path']].append(asset)

    lines = ['# Attribution', '']
    for local_path in sorted(grouped.keys()):
        assets = grouped[local_path]
        species_names = sorted({asset['species_name'] for asset in assets})
        unique_records = dict()
        for asset in assets:
            identity = (
                asset['provider'],
                asset['matched_name'],
                asset['matched_rank'],
                asset['attribution'] or '',
                asset['license_code'],
                asset['license_url'] or '',
                asset['source_page_url'] or '',
            )
            unique_records.setdefault(identity, asset)
        attribution_records = [unique_records[key] for key in sorted(unique_records)]
        lines.append(
            '## {}'.format(
                ', '.join(_safe_markdown_text(name) for name in species_names)
            )
        )
        lines.append('')
        lines.append(
            'Local file: {}'.format(
                _safe_markdown_text(local_path, escape_markdown=False)
            )
        )
        lines.append('')
        for index, record in enumerate(attribution_records, start=1):
            if len(attribution_records) > 1:
                lines.append('### Attribution record {}'.format(index))
                lines.append('')
            lines.append('Provider: {}'.format(_safe_markdown_text(record['provider'])))
            lines.append(
                'Matched taxon: {} ({})'.format(
                    _safe_markdown_text(record['matched_name']),
                    _safe_markdown_text(record['matched_rank']),
                )
            )
            lines.append(
                'Creator / attribution: {}'.format(
                    _safe_markdown_text(record['attribution'] or '')
                )
            )
            lines.append('License: {}'.format(_safe_markdown_text(record['license_code'])))
            if record['license_url']:
                lines.append(
                    'License URL: {}'.format(
                        _safe_markdown_text(record['license_url'])
                    )
                )
            if record['source_page_url']:
                lines.append(
                    'Source: {}'.format(
                        _safe_markdown_text(record['source_page_url'])
                    )
                )
            lines.append('')

    with open(path, 'w') as handle:
        handle.write('\n'.join(lines).rstrip() + '\n')


def validate_args(args):
    if args.out_dir in (None, ''):
        raise ValueError('--out-dir must be specified.')
    auxiliary_outputs = {
        '--manifest-out': getattr(args, 'manifest_out', None) or getattr(args, 'manifest', None),
        '--attribution-out': getattr(args, 'attribution_out', None) or getattr(args, 'attribution', None),
    }
    stdout_auxiliary_outputs = [
        option_name for option_name, path in auxiliary_outputs.items() if path == '-'
    ]
    if stdout_auxiliary_outputs:
        raise ValueError(
            'Auxiliary outputs require file paths, not STDOUT: {}'.format(
                ', '.join(stdout_auxiliary_outputs)
            )
        )
    if int(args.max_per_species) <= 0:
        raise ValueError('--max-per-species must be > 0.')
    if int(getattr(args, 'max_download_bytes', DEFAULT_MAX_DOWNLOAD_BYTES)) <= 0:
        raise ValueError('--max-download-bytes must be > 0.')
    if float(getattr(args, 'query_cache_max_age_hours', DEFAULT_QUERY_CACHE_MAX_AGE_HOURS)) < 0:
        raise ValueError('--query-cache-max-age-hours must be >= 0.')
    max_edge = getattr(args, 'max_edge', None)
    if max_edge not in (None, 0) and int(max_edge) <= 0:
        raise ValueError('--max-edge must be > 0 when specified.')
    if (getattr(args, 'output_format', 'original') == 'jpg') and (getattr(args, 'background', 'white') == 'transparent'):
        raise ValueError('--background transparent is incompatible with --output-format jpg.')


def _close_provider_bundle(session, ncbi):
    _close_ncbi_db(ncbi)
    if session is not None:
        session.close()


def _collect_candidates_for_species_batch(species_names, args, sources):
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
    default_workers = 2 if uses_taxonomy_db else DEFAULT_LOOKUP_WORKERS
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
    return min(DEFAULT_DOWNLOAD_WORKERS, species_count)


def collect_candidates_for_species_map(species_names, args, sources):
    species_names = list(species_names)
    if len(species_names) == 0:
        return dict()
    max_workers = resolve_lookup_worker_count(args=args, sources=sources, species_count=len(species_names))
    if max_workers <= 1:
        return _collect_candidates_for_species_batch(
            species_names=species_names,
            args=args,
            sources=sources,
        )

    species_batches = [species_names[index::max_workers] for index in range(max_workers)]
    future_to_batch = dict()
    results = dict()
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        for species_batch in species_batches:
            future = executor.submit(
                _collect_candidates_for_species_batch,
                species_batch,
                args,
                sources,
            )
            future_to_batch[future] = species_batch
        for future in as_completed(future_to_batch):
            results.update(future.result())
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
            except Exception as exc:
                warn_cleanup_failure('HTTP session', exc)
        if self.session is not None:
            try:
                self.session.close()
            except Exception as exc:
                warn_cleanup_failure('HTTP session', exc)

    def _download_candidate(self, media_url):
        plan_entry = self.media_plan.get(media_url)
        if plan_entry is None:
            raise KeyError('Missing media plan for {}'.format(media_url))
        species_name = plan_entry['species_name']
        candidate = plan_entry['candidate']
        filename = build_local_media_filename(
            species_name=species_name,
            candidate=candidate,
            media_url=media_url,
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
            max_download_bytes=getattr(self.args, 'max_download_bytes', DEFAULT_MAX_DOWNLOAD_BYTES),
            provider=candidate.get('provider'),
            reuse_destination=False,
        )
        resolved_cache_path = None
        if isinstance(download_result, dict):
            download_status = download_result['status']
            absolute_local_path = download_result.get('destination_path', absolute_local_path)
            resolved_cache_path = download_result.get('cache_path')
        else:
            download_status = download_result
        try:
            absolute_local_path = postprocess_media_file(absolute_local_path, args=self.args)
        except (MediaDownloadError, OSError, ValueError):
            remove_temporary_path(absolute_local_path, 'invalid materialized media file')
            if resolved_cache_path and os.path.realpath(resolved_cache_path) != os.path.realpath(absolute_local_path):
                remove_temporary_path(resolved_cache_path, 'invalid cached media file')
            raise
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
        except (
            requests.RequestException,
            OSError,
            MediaDownloadError,
            ValueError,
        ) as exc:
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
                'is_primary': 'yes' if candidate.get('is_primary') else 'no',
                'is_vector': 'yes' if candidate_is_vector(candidate) else 'no',
                'selection_reason': describe_candidate_selection(candidate),
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
    validate_distinct_output_paths([
        ('--out-dir', out_dir),
        ('images directory', images_dir),
        ('--manifest-out', manifest_path),
        ('--attribution-out', attribution_path),
        ('unmatched output', unmatched_path),
    ])
    shared_cache_dir = resolve_image_cache_dir(args)
    ensure_directory(out_dir)
    ensure_directory(images_dir)
    ensure_directory(shared_cache_dir)

    species_name_tsv = getattr(args, 'species_name_tsv', None) or getattr(args, 'name_tsv', None)
    name_mapping = read_name_tsv(species_name_tsv) if species_name_tsv else None
    tree = read_tree(args.infile, args.format, args.quoted_node_names, quiet=True)
    leaf_to_species, unmatched_rows = extract_species_mapping(
        tree,
        name_mapping=name_mapping,
        args=args,
    )
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
        'is_primary',
        'is_vector',
        'selection_reason',
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

    manifest_output_rows = rebase_local_paths_for_output(
        manifest_rows,
        out_dir=out_dir,
        output_path=manifest_path,
    )
    attribution_output_assets = rebase_local_paths_for_output(
        selected_assets,
        out_dir=out_dir,
        output_path=attribution_path,
    )
    write_tsv(manifest_path, manifest_output_rows, manifest_fieldnames)
    write_tsv(unmatched_path, unmatched_rows, unmatched_fieldnames)
    write_attribution_markdown(attribution_path, attribution_output_assets)

    if unmatched_rows and args.fail_on_missing:
        raise ValueError('{} tree tip(s) could not be resolved to an image.'.format(len(unmatched_rows)))

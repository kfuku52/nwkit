import csv
import io
import os
import sqlite3
import stat
import struct
import tarfile
import threading
import time
import zlib
from argparse import Namespace
from concurrent.futures import ThreadPoolExecutor

import pytest
import requests

from nwkit import image as image_module
from nwkit.image import (
    BioiconsProvider,
    EOLProvider,
    IDigBioProvider,
    MediaDownloadError,
    NCBIProvider,
    OpenverseProvider,
    PhylopicProvider,
    allowed_candidates_from_scored_candidates,
    build_download_session,
    build_local_media_filename,
    build_retry_config,
    build_providers,
    candidate_score,
    classify_wikimedia_asset,
    collect_candidates_for_species,
    collect_candidates_for_species_map,
    download_media,
    extract_species_mapping,
    get_aspect_fit_bonus,
    get_provider_quality_bonus,
    get_style_priority,
    image_main,
    postprocess_media_file,
    license_allowed,
    normalize_license_code,
    parse_ncbi_images_dmp_line,
    parse_sources,
    resolve_download_worker_count,
    resolve_image_cache_dir,
    resolve_lookup_worker_count,
    resolve_ncbi_taxonomy_image_cache_dir,
    resolve_provider_fetch_limit,
    wikimedia_page_mentions_query,
)


def make_image_args(**kwargs):
    defaults = {
        'infile': '-',
        'format': 'auto',
        'quoted_node_names': True,
        'download_dir': 'auto',
        'out_dir': None,
        'style': 'auto',
        'source': None,
        'license_max': 'cc-by-nc-sa',
        'allow_nd': False,
        'fallback_rank': 'none',
        'max_per_species': 1,
        'max_download_bytes': 104857600,
        'query_cache_max_age_hours': 168.0,
        'refresh_cache': False,
        'species_name_tsv': None,
        'manifest_out': None,
        'attribution_out': None,
        'fail_on_missing': False,
        'output_format': 'original',
        'max_edge': None,
        'canvas': 'none',
        'background': 'white',
        'trim': 'off',
        'trim_shape': 'bbox',
        'species_regex': r'^([^_]+_[^_]+)(?:_|$)',
    }
    defaults.update(kwargs)
    return Namespace(**defaults)


class DummySession:
    def close(self):
        return None


class JSONResponse:
    def __init__(self, payload, status_code=200, url='https://example.org'):
        self.payload = payload
        self.status_code = status_code
        self.url = url
        self.headers = {'Content-Type': 'application/json'}

    def json(self):
        return self.payload

    def raise_for_status(self):
        if self.status_code >= 400:
            raise requests.HTTPError('HTTP {}'.format(self.status_code))

    def close(self):
        return None


class DummyProvider:
    def __init__(self, candidates_by_species):
        self.candidates_by_species = candidates_by_species

    def fetch_candidates(self, species_name, fallback_rank='none'):
        return [dict(candidate) for candidate in self.candidates_by_species.get(species_name, [])]


def read_tsv(path):
    with open(path, newline='') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def write_valid_test_media(path):
    if str(path).lower().endswith('.svg'):
        with open(path, 'wb') as handle:
            handle.write(b'<svg xmlns="http://www.w3.org/2000/svg"></svg>')
        return
    from PIL import Image

    extension = os.path.splitext(str(path))[1].lower()
    image_format = {
        '.gif': 'GIF',
        '.jpg': 'JPEG',
        '.jpeg': 'JPEG',
        '.png': 'PNG',
        '.tif': 'TIFF',
        '.tiff': 'TIFF',
        '.webp': 'WEBP',
    }.get(extension, 'PNG')
    Image.new('RGB', (2, 2), 'white').save(path, format=image_format)


class TestLicenseHelpers:
    def test_normalize_license_code_from_url(self):
        assert normalize_license_code(raw_url='https://creativecommons.org/licenses/by/4.0/') == 'cc-by'
        assert normalize_license_code(raw_url='https://creativecommons.org/publicdomain/zero/1.0/') == 'public-domain'

    def test_normalize_license_code_from_code_and_attribution(self):
        assert normalize_license_code(raw_code='cc-by-nc-sa') == 'cc-by-nc-sa'
        assert normalize_license_code(raw_code='CC BY-SA 3.0') == 'cc-by-sa'
        assert normalize_license_code(raw_code='Public domain') == 'public-domain'
        assert normalize_license_code(raw_code='CC0 1.0') == 'public-domain'
        assert normalize_license_code(raw_code='MIT') == 'mit'
        assert normalize_license_code(raw_code='BSD-3-Clause') == 'bsd'
        assert normalize_license_code(raw_code='by-sa') == 'cc-by-sa'
        assert normalize_license_code(raw_code='pdm') == 'public-domain'
        assert normalize_license_code(raw_code=None, attribution='(c) Someone, all rights reserved') == 'all-rights-reserved'

    @pytest.mark.parametrize(
        'raw_code',
        [
            'limited use only',
            'no redistribution; limited educational use',
            'submitted by user',
            'not public domain',
        ],
    )
    def test_restrictive_or_unrelated_metadata_is_not_open_license(
        self,
        raw_code,
    ):
        license_code = normalize_license_code(raw_code=raw_code)

        assert license_code == 'unknown'
        assert not license_allowed(license_code, license_max='any')

    @pytest.mark.parametrize(
        ('raw_code', 'expected'),
        [
            ('MIT License', 'mit'),
            ('BSD-2-Clause', 'bsd'),
            ('Public Domain Mark 1.0', 'public-domain'),
        ],
    )
    def test_explicit_open_license_labels_remain_allowed(
        self,
        raw_code,
        expected,
    ):
        license_code = normalize_license_code(raw_code=raw_code)

        assert license_code == expected
        assert license_allowed(license_code, license_max='any')

    @pytest.mark.parametrize(
        'raw_url',
        [
            'https://evil.example/licenses/by/4.0/',
            'https://example.com/?next=/publicdomain/zero/1.0',
            'https://notcreativecommons.org/licenses/by-nc/4.0/',
            'javascript://creativecommons.org/licenses/by/4.0/',
            'https://creativecommons.org:444/licenses/by/4.0/',
        ],
    )
    def test_license_url_requires_an_official_host(self, raw_url):
        assert normalize_license_code(raw_url=raw_url) == 'unknown'

    def test_license_url_accepts_official_subdomains_and_ignores_query(self):
        assert normalize_license_code(
            raw_url=(
                'https://licenses.creativecommons.org/licenses/by/4.0/'
                '?next=/licenses/by-nc/4.0/'
            )
        ) == 'cc-by'

    def test_license_allowed_respects_nd_and_ceiling(self):
        assert license_allowed('cc-by', license_max='cc-by', allow_nd=False) is True
        assert license_allowed('cc-by-nc', license_max='cc-by-sa', allow_nd=False) is False
        assert license_allowed('cc-by-nd', license_max='cc-by', allow_nd=False) is False
        assert license_allowed('cc-by-nd', license_max='cc-by', allow_nd=True) is True
        assert license_allowed('mit', license_max='cc-by', allow_nd=False) is True
        assert license_allowed('bsd', license_max='cc-by-sa', allow_nd=False) is True
        assert license_allowed('mit', license_max='public-domain', allow_nd=False) is False
        assert license_allowed('all-rights-reserved', license_max='any', allow_nd=True) is False


class TestFetchLimits:
    def test_resolve_provider_fetch_limit_scales_with_max_per_species(self):
        assert resolve_provider_fetch_limit(make_image_args(max_per_species=1)) == 10
        assert resolve_provider_fetch_limit(make_image_args(max_per_species=8)) == 12
        assert resolve_provider_fetch_limit(make_image_args(max_per_species=99)) == 30

    def test_candidate_score_prioritizes_source_order_before_license_and_quality(self):
        first_source = {
            'matched_rank': 'species',
            'license_code': 'cc-by',
            'asset_type': 'photo',
            'width': 1200,
            'height': 900,
            'provider_quality': 0,
        }
        later_source = {
            'matched_rank': 'species',
            'license_code': 'public-domain',
            'asset_type': 'photo',
            'width': 8000,
            'height': 8000,
            'provider_quality': 99,
        }

        assert candidate_score(first_source, provider_index=0, style='photo') > candidate_score(
            later_source, provider_index=1, style='photo'
        )

    def test_candidate_score_uses_provider_quality_only_after_quality(self):
        lower_quality = {
            'matched_rank': 'species',
            'license_code': 'cc-by',
            'asset_type': 'photo',
            'width': 1000,
            'height': 1000,
            'provider_quality': 90,
        }
        higher_quality = {
            'matched_rank': 'species',
            'license_code': 'cc-by',
            'asset_type': 'photo',
            'width': 2000,
            'height': 2000,
            'provider_quality': 0,
        }

        assert candidate_score(higher_quality, provider_index=0, style='photo') > candidate_score(
            lower_quality, provider_index=0, style='photo'
        )

    def test_candidate_score_prefers_primary_then_vector_and_drawable_aspect(self):
        base = {
            'matched_rank': 'species',
            'license_code': 'cc-by',
            'asset_type': 'silhouette',
            'width': 1000,
            'height': 1000,
            'provider_quality': 0,
            'media_url': 'https://example.org/image.png',
        }
        primary = dict(base, is_primary=True)
        vector = dict(base, is_vector=True, media_url='https://example.org/image.svg')
        raster = dict(base, is_vector=False)
        extreme_vector = dict(vector, width=10000, height=100)
        balanced_vector = dict(vector, width=1000, height=800)

        assert candidate_score(primary, provider_index=0, style='silhouette') > candidate_score(
            vector, provider_index=0, style='silhouette'
        )
        assert candidate_score(vector, provider_index=0, style='silhouette') > candidate_score(
            raster, provider_index=0, style='silhouette'
        )
        assert candidate_score(balanced_vector, provider_index=0, style='silhouette') > candidate_score(
            extreme_vector, provider_index=0, style='silhouette'
        )
        assert get_aspect_fit_bonus(balanced_vector) == 100
        assert get_aspect_fit_bonus(extreme_vector) == 0

    def test_candidate_score_keeps_exact_match_ahead_of_fallback_primary(self):
        exact = {
            'matched_rank': 'species',
            'license_code': 'cc-by',
            'asset_type': 'silhouette',
            'width': 500,
            'height': 500,
            'is_primary': False,
        }
        genus_primary = dict(exact, matched_rank='genus', is_primary=True)

        assert candidate_score(exact, provider_index=0, style='silhouette') > candidate_score(
            genus_primary, provider_index=0, style='silhouette'
        )

    def test_disallowed_primary_falls_back_to_allowed_ranked_candidate(self):
        primary = {
            'matched_rank': 'species',
            'license_code': 'cc-by-nc-sa',
            'asset_type': 'silhouette',
            'width': 1000,
            'height': 1000,
            'is_primary': True,
            'media_url': 'https://example.org/primary.svg',
        }
        fallback = dict(
            primary,
            license_code='public-domain',
            is_primary=False,
            media_url='https://example.org/fallback.svg',
        )
        for candidate in (primary, fallback):
            candidate['score'] = candidate_score(
                candidate,
                provider_index=0,
                style='silhouette',
            )

        allowed = allowed_candidates_from_scored_candidates(
            [primary, fallback],
            args=make_image_args(license_max='public-domain'),
        )

        assert allowed == [fallback]

    def test_style_and_provider_quality_helpers(self):
        candidate = {'asset_type': 'silhouette', 'provider_quality': 12}
        assert get_style_priority(candidate, style='silhouette') == 2
        assert get_style_priority(candidate, style='photo') == 0
        assert get_provider_quality_bonus(candidate) == 12

    def test_candidate_score_normalizes_numeric_strings_and_nonfinite_values(
        self,
    ):
        numeric_strings = {
            'matched_rank': 'species',
            'license_code': 'cc-by',
            'width': '1024',
            'height': '512',
            'provider_quality': '1e3',
        }
        malformed_numbers = dict(
            numeric_strings,
            width=float('nan'),
            height=float('inf'),
            provider_quality='not-a-number',
        )

        assert isinstance(candidate_score(numeric_strings), int)
        assert isinstance(candidate_score(malformed_numbers), int)
        assert get_aspect_fit_bonus(malformed_numbers) == 0

    def test_candidate_normalization_parses_only_explicit_boolean_values(self):
        candidate = {
            'provider': 'phylopic',
            'provider_record_id': 'candidate-1',
            'matched_name': 'Species alpha',
            'matched_rank': 'species',
            'license_code': 'cc-by',
            'media_url': 'https://images.phylopic.org/candidate.svg',
            'is_primary': 'false',
            'is_vector': 'false',
        }

        normalized = image_module.normalize_provider_candidate(
            candidate,
            expected_provider='phylopic',
        )

        assert normalized['is_primary'] is False
        assert normalized['is_vector'] is False
        assert image_module.candidate_is_vector(normalized) is False

        with pytest.raises(ValueError, match="invalid 'is_vector' value"):
            image_module.normalize_provider_candidate(
                dict(candidate, is_vector='yes'),
                expected_provider='phylopic',
            )

    def test_candidate_collection_skips_malformed_records_individually(
        self,
        tmp_path,
    ):
        base = {
            'provider': 'phylopic',
            'provider_record_id': 'valid',
            'matched_name': 'Species alpha',
            'matched_rank': 'species',
            'license_code': 'cc-by',
            'license_url': 'https://creativecommons.org/licenses/by/4.0/',
            'attribution': 'Artist',
            'source_page_url': 'https://phylopic.org/images/valid',
            'media_url': 'https://images.phylopic.org/valid.svg',
            'asset_type': 'silhouette',
            'width': '1024',
            'height': '512',
        }
        candidates = [
            dict(base, provider_record_id='bad-rank', matched_rank={'bad': 1}),
            dict(base, provider_record_id='bad-license', license_code=['cc-by']),
            base,
        ]
        args = make_image_args(
            out_dir=str(tmp_path / 'out'),
            source='phylopic',
            style='silhouette',
        )

        collected, errors = collect_candidates_for_species(
            'Species alpha',
            args=args,
            sources=['phylopic'],
            providers={
                'phylopic': DummyProvider({'Species alpha': candidates}),
            },
        )

        assert [candidate['provider_record_id'] for candidate in collected] == [
            'valid'
        ]
        assert collected[0]['width'] == 1024.0
        assert len(errors) == 2


class TestPhylopicProvider:
    @staticmethod
    def _image_item(uuid, license_url, media_url, sizes='1000x800', vector=True):
        selected_file = {
            'href': media_url,
            'sizes': sizes,
            'type': 'image/svg+xml' if vector else 'image/png',
        }
        links = {
            'license': {'href': license_url},
            'self': {'href': '/images/{}'.format(uuid)},
            'sourceFile': selected_file,
        }
        if vector:
            links['vectorFile'] = selected_file
        return {
            '_links': links,
            'attribution': 'Example Artist',
            'uuid': uuid,
        }

    def test_fetch_candidates_marks_and_fetches_linked_primary_image(self, monkeypatch):
        primary_uuid = 'primary-uuid'
        fallback_item = self._image_item(
            uuid='fallback-uuid',
            license_url='https://creativecommons.org/publicdomain/zero/1.0/',
            media_url='https://images.example.org/fallback.svg',
        )
        primary_item = self._image_item(
            uuid=primary_uuid,
            license_url='https://creativecommons.org/licenses/by/4.0/',
            media_url='https://images.example.org/primary.png',
            vector=False,
        )

        class FakeNCBI:
            def get_name_translator(self, names):
                return {'Apis mellifera': [7460]}

        provider = PhylopicProvider(session=DummySession(), ncbi=FakeNCBI())
        monkeypatch.setattr(
            provider,
            '_resolve_node',
            lambda taxid: {
                '_links': {
                    'primaryImage': {
                        'href': '/images/{}?build=545'.format(primary_uuid),
                    },
                },
                'build': 545,
                'uuid': 'node-uuid',
            },
        )
        monkeypatch.setattr(provider, '_fetch_node_images', lambda node_uuid, build: [fallback_item])
        linked_hrefs = list()

        def fake_fetch_linked_image(href):
            linked_hrefs.append(href)
            return primary_item

        monkeypatch.setattr(provider, '_fetch_linked_image', fake_fetch_linked_image)

        candidates = provider.fetch_candidates('Apis mellifera')
        candidates_by_uuid = {candidate['provider_record_id']: candidate for candidate in candidates}

        assert linked_hrefs == ['/images/{}?build=545'.format(primary_uuid)]
        assert candidates_by_uuid[primary_uuid]['is_primary'] is True
        assert candidates_by_uuid[primary_uuid]['is_vector'] is False
        assert candidates_by_uuid['fallback-uuid']['is_primary'] is False
        assert candidates_by_uuid['fallback-uuid']['is_vector'] is True
        assert candidate_score(
            candidates_by_uuid[primary_uuid], provider_index=0, style='silhouette'
        ) > candidate_score(
            candidates_by_uuid['fallback-uuid'], provider_index=0, style='silhouette'
        )


class TestSourceParsing:
    def test_parse_sources_uses_style_defaults(self):
        assert parse_sources('auto', None) == [
            'phylopic', 'bioicons', 'inaturalist', 'wikimedia', 'gbif', 'eol', 'idigbio', 'openverse', 'ncbi'
        ]
        assert parse_sources('photo', None) == ['inaturalist', 'wikimedia', 'gbif', 'eol', 'idigbio', 'openverse', 'ncbi']
        assert parse_sources('silhouette', None) == ['phylopic', 'bioicons', 'wikimedia']

    def test_parse_sources_accepts_all_implemented_sources(self):
        assert parse_sources('auto', 'phylopic,bioicons,wikimedia,gbif,inaturalist,eol,idigbio,openverse,ncbi') == [
            'phylopic', 'bioicons', 'wikimedia', 'gbif', 'inaturalist', 'eol', 'idigbio', 'openverse', 'ncbi'
        ]

    def test_parse_sources_rejects_unimplemented_sources(self):
        with pytest.raises(ValueError, match='Unsupported --source'):
            parse_sources('auto', 'phylopic,example')


class TestWikimediaHelpers:
    def test_wikimedia_page_mentions_query_filters_irrelevant_results(self):
        relevant_page = {
            'title': 'File:Lion (Panthera leo) old male Chobe.jpg',
            'imageinfo': [{
                'extmetadata': {
                    'ObjectName': {'value': 'Lion (Panthera leo) old male Chobe'},
                    'ImageDescription': {'value': 'Lion (<i>Panthera leo</i>), male, Chobe National Park, Botswana'},
                }
            }],
        }
        irrelevant_page = {
            'title': 'File:La Bohémienne endormie.jpg',
            'imageinfo': [{
                'extmetadata': {
                    'ObjectName': {'value': 'La Bohémienne endormie'},
                    'ImageDescription': {'value': ''},
                }
            }],
        }
        assert wikimedia_page_mentions_query(relevant_page, 'Panthera leo') is True
        assert wikimedia_page_mentions_query(irrelevant_page, 'Panthera leo') is False

    def test_classify_wikimedia_asset_rejects_research_figures_and_gifs(self):
        page = {
            'title': 'File:Arabidopsis thaliana expression figure.gif',
            'imageinfo': [{
                'url': 'https://upload.wikimedia.org/expression.gif',
                'mime': 'image/gif',
                'extmetadata': {
                    'ImageDescription': {'value': 'Gene-expression graph for Arabidopsis thaliana'},
                },
            }],
        }

        assert classify_wikimedia_asset(page) == 'illustration'

    def test_classify_wikimedia_asset_preserves_silhouettes(self):
        page = {
            'title': 'File:Panthera leo silhouette.svg',
            'imageinfo': [{'url': 'https://upload.wikimedia.org/lion.svg', 'mime': 'image/svg+xml'}],
        }

        assert classify_wikimedia_asset(page) == 'silhouette'


class TestNCBIHelpers:
    def test_parse_ncbi_images_dmp_line(self):
        line = '64373\t|\timage:Alternaria brassicae\t|\thttp://www.ncbi.nlm.nih.gov/Taxonomy/taxi/images/12\t|\tCC BY-NC (https://creativecommons.org/licenses/by-nc/4.0/)\t|\tgen_ok\t|\tiNaturalist\t|\t\t|\t29911\t|\n'
        record = parse_ncbi_images_dmp_line(line)

        assert record['record_id'] == '64373'
        assert record['title'] == 'Alternaria brassicae'
        assert record['image_url'] == 'http://www.ncbi.nlm.nih.gov/Taxonomy/taxi/images/12'
        assert record['license_code_text'] == 'CC BY-NC'
        assert record['license_url'] == 'https://creativecommons.org/licenses/by-nc/4.0/'
        assert record['attribution'] == 'gen_ok'
        assert record['source_name'] == 'iNaturalist'
        assert record['taxids'] == [29911]

    def test_parse_ncbi_images_dmp_line_supports_multiple_taxids(self):
        line = '64698\t|\timage:Abudefduf saxatilis\t|\thttp://www.ncbi.nlm.nih.gov/Taxonomy/taxi/images/31217\t|\tCC BY-SA 3.0 (https://creativecommons.org/licenses/by-sa/3.0/)\t|\tCralize\t|\tWikimedia Commons\t|\t\t|\t50731 1567454\t|\n'
        record = parse_ncbi_images_dmp_line(line)
        assert record['taxids'] == [50731, 1567454]

    def test_build_download_session_configures_retries_for_rate_limits(self):
        session = build_download_session()
        adapter = session.get_adapter('https://example.org/')
        retries = adapter.max_retries
        assert retries.total == 4
        assert 429 in retries.status_forcelist

    def test_build_retry_config_limits_retries_to_get_requests(self):
        retries = build_retry_config()
        allowed_methods = getattr(retries, 'allowed_methods', None)
        if allowed_methods is None:
            allowed_methods = getattr(retries, 'method_whitelist', None)
        assert allowed_methods == frozenset(['GET'])

    def test_resolve_lookup_worker_count_defaults_to_four_without_taxonomy(self):
        workers = resolve_lookup_worker_count(
            args=make_image_args(),
            sources=['wikimedia', 'openverse'],
            species_count=20,
        )
        assert workers == 4

    def test_resolve_lookup_worker_count_defaults_to_two_with_taxonomy(self):
        workers = resolve_lookup_worker_count(
            args=make_image_args(),
            sources=['phylopic', 'wikimedia'],
            species_count=20,
        )
        assert workers == 2

    def test_resolve_download_worker_count_defaults_to_four(self):
        assert resolve_download_worker_count(species_count=20) == 4

    def test_parallel_lookup_closes_taxonomy_handles_in_their_worker_threads(self, monkeypatch, tmp_path):
        import threading

        barrier = threading.Barrier(2)
        close_threads = list()

        class ThreadBoundNCBI:
            def __init__(self):
                self.created_thread = threading.get_ident()

            def close(self):
                close_threads.append((self.created_thread, threading.get_ident()))

        class BarrierProvider:
            def fetch_candidates(self, species_name, fallback_rank='none'):
                barrier.wait(timeout=5)
                return []

        def fake_build_providers(args, sources, session=None):
            return DummySession(), ThreadBoundNCBI(), {'phylopic': BarrierProvider()}

        monkeypatch.setattr('nwkit.image.build_providers', fake_build_providers)
        args = make_image_args(
            out_dir=str(tmp_path / 'out'),
            download_dir=str(tmp_path / 'cache'),
            source='phylopic',
        )

        results = collect_candidates_for_species_map(
            species_names=['Species alpha', 'Species beta'],
            args=args,
            sources=['phylopic'],
        )

        assert set(results) == {'Species alpha', 'Species beta'}
        assert len(close_threads) == 2
        assert len({created_thread for created_thread, _ in close_threads}) == 2
        assert all(created_thread == closed_thread for created_thread, closed_thread in close_threads)

    def test_resolve_ncbi_taxonomy_image_cache_dir(self, tmp_path):
        shared_args = make_image_args(download_dir=str(tmp_path / 'shared'), out_dir=str(tmp_path / 'out'))
        auto_args = make_image_args(download_dir='auto', out_dir=str(tmp_path / 'out'))

        assert resolve_ncbi_taxonomy_image_cache_dir(shared_args) == str(tmp_path / 'shared' / 'ncbi-taxonomy-images')
        assert resolve_ncbi_taxonomy_image_cache_dir(auto_args) == str(tmp_path / 'out' / '.nwkit-cache' / 'ncbi-taxonomy-images')

    def test_ncbi_provider_fetches_candidates_from_images_table(self, monkeypatch, tmp_path):
        images_path = tmp_path / 'images.dmp'
        images_path.write_text(
            '64365\t|\timage:Cyanophora paradoxa\t|\thttp://www.ncbi.nlm.nih.gov/Taxonomy/taxi/images/4\t|\tCC BY-SA 3.0 (https://creativecommons.org/licenses/by-sa/3.0/)\t|\tWolfgang Bettighofer\t|\tWikimedia Commons\t|\t\t|\t2762\t|\n'
        )

        class FakeNCBI:
            def get_name_translator(self, names):
                mapping = {
                    'Cyanophora paradoxa': [2762],
                }
                return {name: mapping[name] for name in names if name in mapping}

        monkeypatch.setattr('nwkit.image.ensure_ncbi_images_table', lambda args, session: str(images_path))

        args = make_image_args(
            out_dir=str(tmp_path / 'out'),
            download_dir=str(tmp_path / 'cache'),
        )
        provider = NCBIProvider(session=DummySession(), ncbi=FakeNCBI(), args=args)
        candidates, errors = collect_candidates_for_species(
            'Cyanophora paradoxa',
            args=args,
            sources=['ncbi'],
            providers={'ncbi': provider},
        )

        assert errors == []
        assert len(candidates) == 1
        candidate = candidates[0]
        assert candidate['provider'] == 'ncbi'
        assert candidate['provider_record_id'] == '64365'
        assert candidate['matched_name'] == 'Cyanophora paradoxa'
        assert candidate['matched_rank'] == 'species'
        assert candidate['license_code'] == 'cc-by-sa'
        assert candidate['license_url'] == 'https://creativecommons.org/licenses/by-sa/3.0/'
        assert candidate['attribution'] == 'Wolfgang Bettighofer, Wikimedia Commons'
        assert candidate['media_url'] == 'https://www.ncbi.nlm.nih.gov/Taxonomy/taxi/images/4'

    def test_ncbi_provider_does_not_upgrade_arbitrary_http_media(self):
        provider = NCBIProvider(
            session=DummySession(),
            ncbi=None,
            args=make_image_args(out_dir='unused'),
        )
        candidate = provider._candidate_from_record(
            record={
                'record_id': '1',
                'image_url': 'http://example.org/Taxonomy/taxi/images/4',
                'attribution': '',
                'source_name': '',
                'license_code_text': '',
                'license_url': '',
            },
            matched_name='Example species',
            matched_rank='species',
        )

        with pytest.raises(MediaDownloadError, match='HTTPS'):
            image_module.validate_candidate_media_url(candidate)

    def test_build_providers_does_not_initialize_ncbi_eagerly(self, monkeypatch, tmp_path):
        call_counter = {'count': 0}

        def fake_get_ete_ncbitaxa(args=None):
            call_counter['count'] += 1
            raise AssertionError('NCBI taxonomy should not initialize during provider construction')

        monkeypatch.setattr('nwkit.image.get_ete_ncbitaxa', fake_get_ete_ncbitaxa)

        session, ncbi, providers = build_providers(
            args=make_image_args(out_dir=str(tmp_path / 'out')),
            sources=['ncbi'],
            session=DummySession(),
        )

        assert session is not None
        assert ncbi is not None
        assert 'ncbi' in providers
        assert call_counter['count'] == 0

    def test_collect_candidates_skips_ncbi_when_earlier_provider_has_allowed_candidate(self):
        call_counter = {'ncbi': 0}

        class CountingProvider(DummyProvider):
            def fetch_candidates(self, species_name, fallback_rank='none'):
                call_counter['ncbi'] += 1
                return super().fetch_candidates(species_name, fallback_rank=fallback_rank)

        providers = {
            'inaturalist': DummyProvider({
                'Apis mellifera': [{
                    'provider': 'inaturalist',
                    'provider_record_id': 'inat-1',
                    'matched_name': 'Apis mellifera',
                    'matched_rank': 'species',
                    'license_code': 'cc-by',
                    'license_url': 'https://creativecommons.org/licenses/by/4.0/',
                    'attribution': 'A',
                    'source_page_url': 'https://example.org/obs',
                    'media_url': 'https://static.inaturalist.org/photo.jpg',
                    'width': 1000,
                    'height': 900,
                    'asset_type': 'photo',
                }],
            }),
            'ncbi': CountingProvider({
                'Apis mellifera': [{
                    'provider': 'ncbi',
                    'provider_record_id': 'ncbi-1',
                    'matched_name': 'Apis mellifera',
                    'matched_rank': 'species',
                    'license_code': 'public-domain',
                    'license_url': 'https://creativecommons.org/publicdomain/zero/1.0/',
                    'attribution': 'NCBI',
                    'source_page_url': 'https://www.ncbi.nlm.nih.gov/Taxonomy/taxi/images/1',
                    'media_url': 'https://www.ncbi.nlm.nih.gov/Taxonomy/taxi/images/1',
                    'width': None,
                    'height': None,
                    'asset_type': 'photo',
                }],
            }),
        }

        candidates, provider_errors = collect_candidates_for_species(
            species_name='Apis mellifera',
            args=make_image_args(style='photo', source='inaturalist,ncbi'),
            sources=['inaturalist', 'ncbi'],
            providers=providers,
        )

        assert provider_errors == []
        assert len(candidates) == 1
        assert candidates[0]['provider'] == 'inaturalist'
        assert call_counter['ncbi'] == 0

    def test_collect_candidates_stops_after_exact_species_allowed_candidate(self, tmp_path):
        call_counter = {'openverse': 0}

        class CountingProvider(DummyProvider):
            def fetch_candidates(self, species_name, fallback_rank='none'):
                call_counter['openverse'] += 1
                return super().fetch_candidates(species_name, fallback_rank=fallback_rank)

        providers = {
            'inaturalist': DummyProvider({
                'Apis mellifera': [
                    {
                        'provider': 'inaturalist',
                        'provider_record_id': 'inat-early-1',
                        'matched_name': 'Apis mellifera',
                        'matched_rank': 'species',
                        'license_code': 'cc-by',
                        'license_url': 'https://creativecommons.org/licenses/by/4.0/',
                        'attribution': 'A1',
                        'source_page_url': 'https://example.org/obs/1',
                        'media_url': 'https://static.inaturalist.org/photo-1.jpg',
                        'width': 1000,
                        'height': 900,
                        'asset_type': 'photo',
                    },
                    {
                        'provider': 'inaturalist',
                        'provider_record_id': 'inat-early-2',
                        'matched_name': 'Apis mellifera',
                        'matched_rank': 'species',
                        'license_code': 'cc-by',
                        'license_url': 'https://creativecommons.org/licenses/by/4.0/',
                        'attribution': 'A2',
                        'source_page_url': 'https://example.org/obs/2',
                        'media_url': 'https://static.inaturalist.org/photo-2.jpg',
                        'width': 1200,
                        'height': 950,
                        'asset_type': 'photo',
                    },
                    {
                        'provider': 'inaturalist',
                        'provider_record_id': 'inat-early-3',
                        'matched_name': 'Apis mellifera',
                        'matched_rank': 'species',
                        'license_code': 'cc-by',
                        'license_url': 'https://creativecommons.org/licenses/by/4.0/',
                        'attribution': 'A3',
                        'source_page_url': 'https://example.org/obs/3',
                        'media_url': 'https://static.inaturalist.org/photo-3.jpg',
                        'width': 1100,
                        'height': 920,
                        'asset_type': 'photo',
                    },
                ],
            }),
            'openverse': CountingProvider({
                'Apis mellifera': [{
                    'provider': 'openverse',
                    'provider_record_id': 'ov-1',
                    'matched_name': 'Apis mellifera',
                    'matched_rank': 'species',
                    'license_code': 'public-domain',
                    'license_url': 'https://creativecommons.org/publicdomain/zero/1.0/',
                    'attribution': 'B',
                    'source_page_url': 'https://example.org/ov/1',
                    'media_url': 'https://static.inaturalist.org/photo-2.jpg',
                    'width': 3000,
                    'height': 2000,
                    'asset_type': 'photo',
                }],
            }),
        }

        candidates, provider_errors = collect_candidates_for_species(
            species_name='Apis mellifera',
            args=make_image_args(out_dir=str(tmp_path / 'out'), style='photo', source='inaturalist,openverse'),
            sources=['inaturalist', 'openverse'],
            providers=providers,
        )

        assert provider_errors == []
        assert len(candidates) == 3
        assert candidates[0]['provider'] == 'inaturalist'
        assert call_counter['openverse'] == 0

    def test_collect_candidates_reuses_query_cache(self, tmp_path):
        call_counter = {'count': 0}

        class CountingProvider(DummyProvider):
            def fetch_candidates(self, species_name, fallback_rank='none'):
                call_counter['count'] += 1
                return super().fetch_candidates(species_name, fallback_rank=fallback_rank)

        providers = {
            'eol': CountingProvider({
                'Apis mellifera': [{
                    'provider': 'eol',
                    'provider_record_id': 'eol-1',
                    'matched_name': 'Apis mellifera',
                    'matched_rank': 'species',
                    'license_code': 'cc-by',
                    'license_url': 'https://creativecommons.org/licenses/by/4.0/',
                    'attribution': 'A',
                    'source_page_url': 'https://example.org/eol/1',
                    'media_url': 'https://example.org/eol-photo.jpg',
                    'width': 1200,
                    'height': 800,
                    'asset_type': 'photo',
                }],
            }),
        }
        args = make_image_args(out_dir=str(tmp_path / 'out'), style='photo', source='eol')

        first_candidates, first_errors = collect_candidates_for_species(
            species_name='Apis mellifera',
            args=args,
            sources=['eol'],
            providers=providers,
        )
        second_candidates, second_errors = collect_candidates_for_species(
            species_name='Apis mellifera',
            args=args,
            sources=['eol'],
            providers=providers,
        )

        assert first_errors == []
        assert second_errors == []
        assert call_counter['count'] == 1
        assert [candidate['provider_record_id'] for candidate in first_candidates] == ['eol-1']
        assert [candidate['provider_record_id'] for candidate in second_candidates] == ['eol-1']

    def test_query_cache_refetches_when_requested_capacity_increases(self, tmp_path):
        call_counter = {'count': 0}

        class CapacityProvider(DummyProvider):
            result_limit = 10

            def fetch_candidates(self, species_name, fallback_rank='none'):
                call_counter['count'] += 1
                return []

        provider = CapacityProvider({})
        providers = {'eol': provider}
        base_args = make_image_args(out_dir=str(tmp_path / 'out'), source='eol')
        collect_candidates_for_species('Apis mellifera', base_args, ['eol'], providers)

        provider.result_limit = 24
        expanded_args = make_image_args(out_dir=str(tmp_path / 'out'), source='eol', max_per_species=20)
        collect_candidates_for_species('Apis mellifera', expanded_args, ['eol'], providers)

        assert call_counter['count'] == 2

    def test_query_cache_honors_expiration_and_explicit_refresh(self, monkeypatch, tmp_path):
        call_counter = {'count': 0}
        current_time = {'value': 100.0}

        class CountingProvider(DummyProvider):
            result_limit = 10

            def fetch_candidates(self, species_name, fallback_rank='none'):
                call_counter['count'] += 1
                return []

        monkeypatch.setattr('nwkit.image.time.time', lambda: current_time['value'])
        providers = {'eol': CountingProvider({})}
        args = make_image_args(
            out_dir=str(tmp_path / 'out'),
            source='eol',
            query_cache_max_age_hours=1.0,
        )
        collect_candidates_for_species('Apis mellifera', args, ['eol'], providers)
        current_time['value'] = 3801.0
        collect_candidates_for_species('Apis mellifera', args, ['eol'], providers)
        collect_candidates_for_species(
            'Apis mellifera',
            make_image_args(out_dir=str(tmp_path / 'out'), source='eol', refresh_cache=True),
            ['eol'],
            providers,
        )

        assert call_counter['count'] == 3


class TestBioiconsProvider:
    def test_bioicons_provider_fetches_matching_svg_candidates(self, tmp_path):
        class RoutingSession:
            def get(self, url, params=None, timeout=None, headers=None, stream=None, allow_redirects=None):
                assert url.endswith('/git/trees/main')
                return JSONResponse({
                    'tree': [
                        {'path': 'static/icons/cc-0/Animals/Ben-Murrell/Mouse.svg'},
                        {'path': 'static/icons/cc-by-4.0/Animals/DBCLS/Xenopus_laevis.svg'},
                        {'path': 'static/icons/cc-by-3.0/Animals/Servier/rat-adult.svg'},
                        {'path': 'static/icons/categories.json'},
                    ],
                }, url=url)

            def close(self):
                return None

        provider = BioiconsProvider(session=RoutingSession(), args=make_image_args(out_dir=str(tmp_path / 'out')))
        candidates = provider.fetch_candidates('Mus musculus', fallback_rank='none')

        assert len(candidates) == 1
        candidate = candidates[0]
        assert candidate['provider'] == 'bioicons'
        assert candidate['provider_record_id'] == 'cc-0/Animals/Ben-Murrell/Mouse.svg'
        assert candidate['matched_name'] == 'Mus musculus'
        assert candidate['matched_rank'] == 'species'
        assert candidate['license_code'] == 'public-domain'
        assert candidate['license_url'] == 'https://creativecommons.org/publicdomain/zero/1.0/'
        assert candidate['attribution'] == 'Ben Murrell'
        assert candidate['media_url'] == 'https://bioicons.com/icons/cc-0/Animals/Ben-Murrell/Mouse.svg'
        assert candidate['asset_type'] == 'silhouette'
        assert candidate['provider_quality'] > 0

    @pytest.mark.parametrize('parallel', [False, True], ids=['sequential', 'parallel'])
    def test_bioicons_refresh_is_coalesced_across_provider_instances(
        self,
        tmp_path,
        parallel,
    ):
        request_count = {'value': 0}
        request_count_lock = threading.Lock()
        args = make_image_args(
            out_dir=str(tmp_path / 'out'),
            refresh_cache=True,
        )

        class CountingSession:
            def get(
                self,
                url,
                params=None,
                timeout=None,
                headers=None,
                stream=None,
                allow_redirects=None,
            ):
                with request_count_lock:
                    request_count['value'] += 1
                time.sleep(0.05)
                return JSONResponse({
                    'tree': [{
                        'path': 'static/icons/cc-0/Animals/Ben-Murrell/Mouse.svg',
                    }],
                }, url=url)

        providers = [
            BioiconsProvider(session=CountingSession(), args=args)
            for _ in range(4)
        ]
        if parallel:
            with ThreadPoolExecutor(max_workers=len(providers)) as executor:
                catalogs = list(executor.map(
                    lambda provider: provider._load_catalog(),
                    providers,
                ))
        else:
            catalogs = [provider._load_catalog() for provider in providers]

        assert request_count['value'] == 1
        assert all(catalog is catalogs[0] for catalog in catalogs)
        assert catalogs[0][0]['relative_path'].endswith('/Mouse.svg')


class TestEOLProvider:
    def test_eol_provider_fetches_candidates_from_page_media(self):
        class RoutingSession:
            def get(self, url, params=None, timeout=None, headers=None, stream=None, allow_redirects=None):
                if url.endswith('/search/1.0.json'):
                    return JSONResponse({
                        'results': [{
                            'id': 491995,
                            'title': 'Hapalochlaena lunulata',
                        }],
                    }, url=url)
                if url.endswith('/pages/1.0/491995.json'):
                    assert params['images_per_page'] == 10
                    return JSONResponse({
                        'taxonConcept': {
                            'dataObjects': [{
                                'identifier': 'EOL-media-509-demo',
                                'dataObjectVersionID': 28895676,
                                'dataType': 'http://purl.org/dc/dcmitype/StillImage',
                                'mediumType': 'image',
                                'dataRating': '2.5',
                                'vettedStatus': 'Trusted',
                                'license': 'http://creativecommons.org/licenses/by/3.0/',
                                'rightsHolder': 'Elias Levy',
                                'source': 'https://commons.wikimedia.org/wiki/File:Blue-Ringed_Octopus.jpg',
                                'mediaURL': 'https://upload.wikimedia.org/blue-ringed.jpg',
                                'agents': [{
                                    'full_name': 'Elias Levy',
                                    'role': 'creator',
                                }],
                            }],
                        },
                    }, url=url)
                raise AssertionError('Unexpected URL: {}'.format(url))

            def close(self):
                return None

        provider = EOLProvider(session=RoutingSession(), args=make_image_args(max_per_species=1))
        candidates = provider.fetch_candidates('Hapalochlaena lunulata', fallback_rank='none')

        assert len(candidates) == 1
        candidate = candidates[0]
        assert candidate['provider'] == 'eol'
        assert candidate['provider_record_id'] == '28895676'
        assert candidate['matched_name'] == 'Hapalochlaena lunulata'
        assert candidate['matched_rank'] == 'species'
        assert candidate['license_code'] == 'cc-by'
        assert candidate['license_url'] == 'http://creativecommons.org/licenses/by/3.0/'
        assert candidate['attribution'] == 'Elias Levy'
        assert candidate['source_page_url'] == 'https://commons.wikimedia.org/wiki/File:Blue-Ringed_Octopus.jpg'
        assert candidate['media_url'] == 'https://upload.wikimedia.org/blue-ringed.jpg'
        assert candidate['asset_type'] == 'photo'
        assert candidate['provider_quality'] > 0


class TestIDigBioProvider:
    def test_idigbio_provider_fetches_candidates_from_media_search(self):
        class RoutingSession:
            def post(self, url, json=None, timeout=None, headers=None, stream=None, allow_redirects=None):
                assert url.endswith('/search/media')
                assert json['rq'] == {'scientificname': 'Panthera leo'}
                assert json['limit'] == 10
                return JSONResponse({
                    'items': [{
                        'uuid': 'idigbio-1',
                        'data': {
                            'ac:accessURI': 'https://collections.example.org/panthera_leo.jpg',
                            'dcterms:identifier': 'https://collections.example.org/object/1',
                            'xmpRights:UsageTerms': 'https://creativecommons.org/publicdomain/zero/1.0/',
                            'dcterms:rights': 'CC0',
                            'xmpRights:Owner': 'Museum Example',
                            'xmpRights:WebStatement': 'https://collections.example.org/license',
                            'dwc:scientificName': 'Panthera leo',
                            'exif:PixelXDimension': '2048',
                            'exif:PixelYDimension': '1024',
                        },
                        'indexTerms': {
                            'mediatype': 'images',
                            'dqs': 0.8,
                        },
                    }],
                }, url=url)

            def close(self):
                return None

        provider = IDigBioProvider(session=RoutingSession(), args=make_image_args(max_per_species=1))
        candidates = provider.fetch_candidates('Panthera leo', fallback_rank='none')

        assert len(candidates) == 1
        candidate = candidates[0]
        assert candidate['provider'] == 'idigbio'
        assert candidate['provider_record_id'] == 'idigbio-1'
        assert candidate['matched_name'] == 'Panthera leo'
        assert candidate['matched_rank'] == 'species'
        assert candidate['license_code'] == 'public-domain'
        assert candidate['license_url'] == 'https://creativecommons.org/publicdomain/zero/1.0/'
        assert candidate['attribution'] == 'Museum Example'
        assert candidate['source_page_url'] == 'https://collections.example.org/object/1'
        assert candidate['media_url'] == 'https://collections.example.org/panthera_leo.jpg'
        assert candidate['width'] == 2048
        assert candidate['height'] == 1024
        assert candidate['asset_type'] == 'photo'
        assert candidate['provider_quality'] > 0


class TestOpenverseProvider:
    def test_openverse_provider_fetches_relevant_candidates(self):
        class RoutingSession:
            def get(self, url, params=None, timeout=None, headers=None, stream=None, allow_redirects=None):
                assert url.endswith('/images/')
                assert params['q'] == 'Danio rerio'
                assert params['page_size'] == 10
                return JSONResponse({
                    'results': [
                        {
                            'id': 'openverse-1',
                            'title': 'Danio rerio (Peix Zebra/Zebrafish)',
                            'url': 'https://cdn.example.org/danio.jpg',
                            'creator': 'berarma',
                            'license': 'by-sa',
                            'license_url': 'https://creativecommons.org/licenses/by-sa/2.0/',
                            'foreign_landing_url': 'https://www.flickr.com/photos/example/1',
                            'width': 1024,
                            'height': 683,
                            'tags': [],
                            'fields_matched': ['title'],
                        },
                        {
                            'id': 'openverse-2',
                            'title': 'Zebra fish in aquarium',
                            'url': 'https://cdn.example.org/irrelevant.jpg',
                            'creator': 'someone',
                            'license': 'by',
                            'license_url': 'https://creativecommons.org/licenses/by/4.0/',
                            'foreign_landing_url': 'https://www.flickr.com/photos/example/2',
                            'width': 1200,
                            'height': 800,
                            'tags': [],
                            'fields_matched': ['title'],
                        },
                    ],
                }, url=url)

            def close(self):
                return None

        provider = OpenverseProvider(session=RoutingSession(), args=make_image_args(max_per_species=1))
        candidates = provider.fetch_candidates('Danio rerio', fallback_rank='none')

        assert len(candidates) == 1
        candidate = candidates[0]
        assert candidate['provider'] == 'openverse'
        assert candidate['provider_record_id'] == 'openverse-1'
        assert candidate['matched_name'] == 'Danio rerio'
        assert candidate['matched_rank'] == 'species'
        assert candidate['license_code'] == 'cc-by-sa'
        assert candidate['license_url'] == 'https://creativecommons.org/licenses/by-sa/2.0/'
        assert candidate['attribution'] == 'berarma'
        assert candidate['source_page_url'] == 'https://www.flickr.com/photos/example/1'
        assert candidate['media_url'] == 'https://cdn.example.org/danio.jpg'
        assert candidate['asset_type'] == 'photo'
        assert candidate['provider_quality'] > 0


class TestImageMain:
    def test_resolve_image_cache_dir_uses_download_dir(self, tmp_path):
        args = make_image_args(download_dir=str(tmp_path / 'cache'))
        assert resolve_image_cache_dir(args) == str(tmp_path / 'cache' / 'nwkit' / 'image-cache')

    def test_extract_species_mapping_reports_unparsable_labels(self, tmp_path):
        tree_path = tmp_path / 'tree.nwk'
        tree_path.write_text('(Homo_sapiens_A,BadLabel);')

        from nwkit.util import read_tree

        tree = read_tree(str(tree_path), format='auto', quoted_node_names=True, quiet=True)
        leaf_to_species, unmatched_rows = extract_species_mapping(tree)

        assert leaf_to_species == {'Homo_sapiens_A': 'Homo sapiens'}
        assert unmatched_rows == [{
            'leaf_name': 'BadLabel',
            'species_name': '',
            'reason': 'unparsable leaf label',
            'details': 'Expected the configured species parser or a matching --species-name-tsv entry.',
        }]

    def test_extract_species_mapping_accepts_custom_species_regex(self, tmp_path):
        tree_path = tmp_path / 'tree.nwk'
        tree_path.write_text('(Homo.sapiens|A,Mus.musculus|B);')

        from nwkit.util import read_tree

        tree = read_tree(str(tree_path), format='auto', quoted_node_names=True, quiet=True)
        leaf_to_species, unmatched_rows = extract_species_mapping(
            tree,
            species_regex=r'^([A-Za-z]+)\.([A-Za-z]+)\|',
        )

        assert leaf_to_species == {
            'Homo.sapiens|A': 'Homo sapiens',
            'Mus.musculus|B': 'Mus musculus',
        }
        assert unmatched_rows == []

    def test_image_main_writes_manifest_attribution_and_unmatched(self, monkeypatch, tmp_path):
        tree_path = tmp_path / 'tree.nwk'
        tree_path.write_text('((Homo_sapiens_A,Homo_sapiens_B),Panthera_leo_C,BadLabel);')
        out_dir = tmp_path / 'out'

        phylopic_candidates = {
            'Homo sapiens': [{
                'provider': 'phylopic',
                'provider_record_id': 'phy-1',
                'matched_name': 'Homo sapiens',
                'matched_rank': 'species',
                'license_code': 'cc-by',
                'license_url': 'https://creativecommons.org/licenses/by/4.0/',
                'attribution': 'PhyloPic Artist',
                'source_page_url': 'https://api.phylopic.org/images/phy-1',
                'media_url': 'https://images.phylopic.org/homo.svg',
                'width': 1200,
                'height': 800,
                'asset_type': 'silhouette',
                'is_primary': True,
                'is_vector': True,
            }],
        }
        inaturalist_candidates = {
            'Panthera leo': [{
                'provider': 'inaturalist',
                'provider_record_id': 'inat-1',
                'matched_name': 'Panthera leo',
                'matched_rank': 'species',
                'license_code': 'cc-by-nc-sa',
                'license_url': 'https://creativecommons.org/licenses/by-nc-sa/4.0/',
                'attribution': '(c) Example Photographer, some rights reserved',
                'source_page_url': 'https://www.inaturalist.org/observations/1',
                'media_url': 'https://static.inaturalist.org/photos/1/original.jpg',
                'width': 1600,
                'height': 1200,
                'asset_type': 'photo',
            }],
        }

        def fake_build_providers(args, sources, session=None):
            providers = {
                'phylopic': DummyProvider(phylopic_candidates),
                'inaturalist': DummyProvider(inaturalist_candidates),
            }
            return DummySession(), None, providers

        def fake_download_media(session, media_url, destination_path, cache_path=None, max_download_bytes=None, provider=None, **kwargs):
            write_valid_test_media(destination_path)
            return 'downloaded'

        monkeypatch.setattr('nwkit.image.build_providers', fake_build_providers)
        monkeypatch.setattr('nwkit.image.download_media', fake_download_media)

        args = make_image_args(
            infile=str(tree_path),
            out_dir=str(out_dir),
            source='phylopic,inaturalist',
        )
        image_main(args)

        manifest_rows = read_tsv(out_dir / 'manifest.tsv')
        unmatched_rows = read_tsv(out_dir / 'unmatched.tsv')
        attribution_text = (out_dir / 'ATTRIBUTION.md').read_text()

        assert len(manifest_rows) == 3
        assert {row['leaf_name'] for row in manifest_rows} == {'Homo_sapiens_A', 'Homo_sapiens_B', 'Panthera_leo_C'}
        assert {row['species_name'] for row in manifest_rows} == {'Homo sapiens', 'Panthera leo'}
        assert {row['status'] for row in manifest_rows} == {'downloaded'}
        assert any(row['local_path'].endswith('.svg') for row in manifest_rows)
        assert any(row['local_path'].endswith('.jpg') for row in manifest_rows)
        homo_row = next(row for row in manifest_rows if row['species_name'] == 'Homo sapiens')
        assert homo_row['is_primary'] == 'yes'
        assert homo_row['is_vector'] == 'yes'
        assert homo_row['selection_reason'] == 'exact_species_match;phylopic_primary_image'

        assert unmatched_rows == [{
            'leaf_name': 'BadLabel',
            'species_name': '',
            'reason': 'unparsable leaf label',
            'details': 'Expected the configured species parser or a matching --species-name-tsv entry.',
        }]
        assert 'Homo sapiens' in attribution_text
        assert 'Panthera leo' in attribution_text

    def test_attribution_preserves_distinct_records_for_shared_media(self, tmp_path):
        from nwkit.image import write_attribution_markdown

        path = tmp_path / 'ATTRIBUTION.md'
        shared_path = 'images/shared.jpg'
        rows = [
            {
                'local_path': shared_path,
                'species_name': 'Species alpha',
                'provider': 'provider-a',
                'matched_name': 'Species alpha',
                'matched_rank': 'species',
                'attribution': 'Alice',
                'license_code': 'cc-by',
                'license_url': 'https://example.org/license-a',
                'source_page_url': 'https://example.org/source-a',
            },
            {
                'local_path': shared_path,
                'species_name': 'Species beta',
                'provider': 'provider-b',
                'matched_name': 'Species beta',
                'matched_rank': 'species',
                'attribution': 'Bob',
                'license_code': 'cc-by-sa',
                'license_url': 'https://example.org/license-b',
                'source_page_url': 'https://example.org/source-b',
            },
        ]

        write_attribution_markdown(str(path), rows)

        text = path.read_text()
        assert text.count('### Attribution record') == 2
        assert 'Creator / attribution: Alice' in text
        assert 'Creator / attribution: Bob' in text
        assert 'License: cc-by\n' in text
        assert 'License: cc-by-sa\n' in text

    def test_attribution_escapes_provider_controlled_markdown(self, tmp_path):
        from nwkit.image import write_attribution_markdown

        path = tmp_path / 'ATTRIBUTION.md'
        rows = [{
            'local_path': 'images/example.jpg',
            'species_name': 'Species alpha\n## Forged species',
            'provider': 'provider',
            'matched_name': '<b>Species alpha</b>',
            'matched_rank': 'species',
            'attribution': '</p>\n## Fake license\n![track](https://evil.example/pixel)',
            'license_code': 'cc-by',
            'license_url': 'https://example.org/license',
            'source_page_url': 'https://example.org/source',
        }]

        write_attribution_markdown(str(path), rows)

        text = path.read_text()
        assert '\n## Forged species' not in text
        assert '\n## Fake license' not in text
        assert '<b>' not in text
        assert '</p>' not in text
        assert '![track](' not in text
        assert '&lt;b&gt;Species alpha&lt;/b&gt;' in text


class TestMediaFilenames:
    @pytest.mark.parametrize(
        ('url', 'expected'),
        [
            ('https://example.org/image.JPG?size=large', '.jpg'),
            ('https://example.org/image.jpeg', '.jpeg'),
            ('https://example.org/image.tiff', '.tiff'),
            ('https://example.org/image.svg#preview', '.svg'),
            ('https://example.org/image.php', '.bin'),
            ('https://example.org/image.tar.gz', '.bin'),
            ('https://example.org/image.\0jpg', '.bin'),
            ('https://example.org/image\0.jpg', '.bin'),
            ('https://example.org/image.' + ('x' * 4096), '.bin'),
        ],
    )
    def test_infer_extension_only_accepts_known_image_suffixes(
            self, url, expected):
        assert image_module.infer_extension(url) == expected

    def test_local_filename_distinguishes_urls_with_the_same_provider_record_id(self):
        candidate = {
            'provider': 'gbif',
            'provider_record_id': '12345',
        }

        first = build_local_media_filename('Species alpha', candidate, 'https://example.org/front.jpg')
        second = build_local_media_filename('Species alpha', candidate, 'https://example.org/back.jpg')

        assert first != second
        assert first.endswith('.jpg')
        assert second.endswith('.jpg')

    def test_long_filename_components_are_bounded_and_collision_resistant(self):
        long_species = 'Species_' + ('a' * 500)
        shared_prefix = 'record-' + ('x' * 500)
        first_candidate = {
            'provider': 'provider-' + ('p' * 500),
            'provider_record_id': shared_prefix + '-first',
        }
        second_candidate = {
            **first_candidate,
            'provider_record_id': shared_prefix + '-second',
        }

        first_local = build_local_media_filename(
            long_species,
            first_candidate,
            'https://example.org/image.jpg',
        )
        second_local = build_local_media_filename(
            long_species,
            second_candidate,
            'https://example.org/image.jpg',
        )
        cache_name = os.path.basename(image_module.build_media_cache_path(
            '/tmp/cache',
            'https://example.org/image.jpg',
            first_candidate['provider'],
            first_candidate['provider_record_id'],
        ))
        query_name = os.path.basename(image_module.build_query_cache_path(
            '/tmp/cache',
            first_candidate['provider'],
            long_species,
            'fallback-' + ('f' * 500),
        ))

        assert first_local != second_local
        assert image_module.sanitize_filename_component('Apis mellifera') == 'Apis_mellifera'
        for filename in (first_local, second_local, cache_name, query_name):
            assert len(os.fsencode(filename)) < 255

    def test_image_main_uses_name_tsv_override_and_strict_mode(self, monkeypatch, tmp_path):
        tree_path = tmp_path / 'tree.nwk'
        tree_path.write_text('(Sample_1,Unknown);')
        name_tsv = tmp_path / 'names.tsv'
        name_tsv.write_text('leaf_name\tspecies_name\nSample_1\tApis mellifera\n')
        out_dir = tmp_path / 'out'

        phylopic_candidates = {
            'Apis mellifera': [{
                'provider': 'phylopic',
                'provider_record_id': 'phy-2',
                'matched_name': 'Apis mellifera',
                'matched_rank': 'species',
                'license_code': 'public-domain',
                'license_url': 'https://creativecommons.org/publicdomain/zero/1.0/',
                'attribution': 'PhyloPic Artist',
                'source_page_url': 'https://api.phylopic.org/images/phy-2',
                'media_url': 'https://images.phylopic.org/apis.svg',
                'width': 1000,
                'height': 900,
                'asset_type': 'silhouette',
            }],
        }

        def fake_build_providers(args, sources, session=None):
            providers = {
                'phylopic': DummyProvider(phylopic_candidates),
            }
            return DummySession(), None, providers

        def fake_download_media(session, media_url, destination_path, cache_path=None, max_download_bytes=None, provider=None, **kwargs):
            write_valid_test_media(destination_path)
            return 'downloaded'

        monkeypatch.setattr('nwkit.image.build_providers', fake_build_providers)
        monkeypatch.setattr('nwkit.image.download_media', fake_download_media)

        args = make_image_args(
            infile=str(tree_path),
            out_dir=str(out_dir),
            species_name_tsv=str(name_tsv),
            source='phylopic',
            fail_on_missing=True,
        )
        with pytest.raises(ValueError, match='could not be resolved'):
            image_main(args)

        manifest_rows = read_tsv(out_dir / 'manifest.tsv')
        unmatched_rows = read_tsv(out_dir / 'unmatched.tsv')

        assert manifest_rows[0]['leaf_name'] == 'Sample_1'
        assert manifest_rows[0]['species_name'] == 'Apis mellifera'
        assert unmatched_rows == [{
            'leaf_name': 'Unknown',
            'species_name': '',
            'reason': 'unparsable leaf label',
            'details': 'Expected the configured species parser or a matching --species-name-tsv entry.',
        }]

    def test_image_main_reports_filtered_by_license(self, monkeypatch, tmp_path):
        tree_path = tmp_path / 'tree.nwk'
        tree_path.write_text('(Panthera_leo_A);')
        out_dir = tmp_path / 'out'

        inaturalist_candidates = {
            'Panthera leo': [{
                'provider': 'inaturalist',
                'provider_record_id': 'inat-2',
                'matched_name': 'Panthera leo',
                'matched_rank': 'species',
                'license_code': 'all-rights-reserved',
                'license_url': '',
                'attribution': '(c) Example Photographer, all rights reserved',
                'source_page_url': 'https://www.inaturalist.org/observations/2',
                'media_url': 'https://static.inaturalist.org/photos/2/original.jpg',
                'width': 1600,
                'height': 1200,
                'asset_type': 'photo',
            }],
        }

        def fake_build_providers(args, sources, session=None):
            providers = {
                'inaturalist': DummyProvider(inaturalist_candidates),
            }
            return DummySession(), None, providers

        monkeypatch.setattr('nwkit.image.build_providers', fake_build_providers)

        args = make_image_args(
            infile=str(tree_path),
            out_dir=str(out_dir),
            style='photo',
            source='inaturalist',
        )
        image_main(args)

        manifest_rows = read_tsv(out_dir / 'manifest.tsv')
        unmatched_rows = read_tsv(out_dir / 'unmatched.tsv')

        assert manifest_rows == []
        assert unmatched_rows == [{
            'leaf_name': 'Panthera_leo_A',
            'species_name': 'Panthera leo',
            'reason': 'only disallowed-license candidates found',
            'details': '',
        }]

    def test_image_main_reuses_shared_download_cache_across_runs(self, monkeypatch, tmp_path):
        tree_path = tmp_path / 'tree.nwk'
        tree_path.write_text('(Apis_mellifera_A);')
        shared_download_dir = tmp_path / 'shared-cache'
        out_dir1 = tmp_path / 'out1'
        out_dir2 = tmp_path / 'out2'

        phylopic_candidates = {
            'Apis mellifera': [{
                'provider': 'phylopic',
                'provider_record_id': 'phy-cache',
                'matched_name': 'Apis mellifera',
                'matched_rank': 'species',
                'license_code': 'cc-by',
                'license_url': 'https://creativecommons.org/licenses/by/4.0/',
                'attribution': 'PhyloPic Artist',
                'source_page_url': 'https://api.phylopic.org/images/phy-cache',
                'media_url': 'https://images.phylopic.org/apis-cache.svg',
                'width': 1000,
                'height': 900,
                'asset_type': 'silhouette',
            }],
        }
        call_counter = {'count': 0}

        class DownloadOnlyResponse:
            def raise_for_status(self):
                return None

            def iter_content(self, chunk_size=65536):
                yield b'<svg xmlns="http://www.w3.org/2000/svg"></svg>'

        class CountingSession:
            def get(self, media_url, stream=True, timeout=None, headers=None, allow_redirects=False):
                call_counter['count'] += 1
                return DownloadOnlyResponse()

            def close(self):
                return None

        def fake_build_providers(args, sources, session=None):
            providers = {
                'phylopic': DummyProvider(phylopic_candidates),
            }
            return CountingSession(), None, providers

        monkeypatch.setattr('nwkit.image.build_providers', fake_build_providers)
        monkeypatch.setattr('nwkit.image.build_download_session', CountingSession)

        args1 = make_image_args(
            infile=str(tree_path),
            out_dir=str(out_dir1),
            source='phylopic',
            download_dir=str(shared_download_dir),
        )
        image_main(args1)

        args2 = make_image_args(
            infile=str(tree_path),
            out_dir=str(out_dir2),
            source='phylopic',
            download_dir=str(shared_download_dir),
        )
        image_main(args2)

        manifest_rows_1 = read_tsv(out_dir1 / 'manifest.tsv')
        manifest_rows_2 = read_tsv(out_dir2 / 'manifest.tsv')

        assert call_counter['count'] == 1
        assert manifest_rows_1[0]['status'] == 'downloaded'
        assert manifest_rows_2[0]['status'] == 'cached'
        assert (shared_download_dir / 'nwkit' / 'image-cache').is_dir()

    @pytest.mark.parametrize(
        ('use_shared_cache', 'expected_downloads'),
        [(False, 1), (True, 1)],
        ids=['reuse-local-raw-cache', 'reuse-shared-raw-cache'],
    )
    def test_image_main_reprocesses_raw_media_when_options_change(
        self,
        monkeypatch,
        tmp_path,
        use_shared_cache,
        expected_downloads,
    ):
        Image = pytest.importorskip('PIL.Image')
        tree_path = tmp_path / 'tree.nwk'
        tree_path.write_text('(Apis_mellifera_A);')
        out_dir = tmp_path / 'out'
        shared_download_dir = tmp_path / 'shared-cache'
        image_bytes = io.BytesIO()
        Image.new('RGB', (80, 40), 'black').save(image_bytes, format='PNG')
        raw_png = image_bytes.getvalue()
        call_counter = {'count': 0}
        candidates = {
            'Apis mellifera': [{
                'provider': 'phylopic',
                'provider_record_id': 'phy-options',
                'matched_name': 'Apis mellifera',
                'matched_rank': 'species',
                'license_code': 'cc-by',
                'license_url': 'https://creativecommons.org/licenses/by/4.0/',
                'attribution': 'PhyloPic Artist',
                'source_page_url': 'https://api.phylopic.org/images/phy-options',
                'media_url': 'https://images.phylopic.org/apis-options.png',
                'width': 80,
                'height': 40,
                'asset_type': 'silhouette',
            }],
        }

        class ImageResponse:
            status_code = 200
            headers = {'Content-Type': 'image/png'}
            url = 'https://images.phylopic.org/apis-options.png'

            def raise_for_status(self):
                return None

            def iter_content(self, chunk_size=65536):
                yield raw_png

            def close(self):
                return None

        class CountingSession:
            def get(self, url, **kwargs):
                call_counter['count'] += 1
                return ImageResponse()

            def close(self):
                return None

        def fake_build_providers(args, sources, session=None):
            return DummySession(), None, {
                'phylopic': DummyProvider(candidates),
            }

        monkeypatch.setattr('nwkit.image.build_providers', fake_build_providers)
        monkeypatch.setattr('nwkit.image.build_download_session', CountingSession)
        download_dir = str(shared_download_dir) if use_shared_cache else 'auto'

        image_main(make_image_args(
            infile=str(tree_path),
            out_dir=str(out_dir),
            source='phylopic',
            download_dir=download_dir,
            output_format='png',
            max_edge=20,
        ))
        first_manifest = read_tsv(out_dir / 'manifest.tsv')
        first_path = out_dir / first_manifest[0]['local_path']
        with Image.open(first_path) as first_image:
            assert first_image.size == (20, 10)

        image_main(make_image_args(
            infile=str(tree_path),
            out_dir=str(out_dir),
            source='phylopic',
            download_dir=download_dir,
            output_format='png',
            max_edge=40,
        ))
        second_manifest = read_tsv(out_dir / 'manifest.tsv')
        second_path = out_dir / second_manifest[0]['local_path']
        with Image.open(second_path) as second_image:
            assert second_image.size == (40, 20)

        assert call_counter['count'] == expected_downloads

    def test_custom_manifest_and_attribution_paths_feed_draw_with_ranked_rows(
        self,
        monkeypatch,
        tmp_path,
    ):
        tree_path = tmp_path / 'tree.nwk'
        tree_path.write_text('(Apis_mellifera_A);')
        out_dir = tmp_path / 'image-output'
        metadata_dir = tmp_path / 'metadata'
        manifest_path = metadata_dir / 'ranked-images.tsv'
        attribution_path = metadata_dir / 'ATTRIBUTION.md'
        candidates = {
            'Apis mellifera': [
                {
                    'provider': 'phylopic',
                    'provider_record_id': 'phy-first',
                    'matched_name': 'Apis mellifera',
                    'matched_rank': 'species',
                    'license_code': 'cc-by',
                    'license_url': 'https://creativecommons.org/licenses/by/4.0/',
                    'attribution': 'First Artist',
                    'source_page_url': 'https://api.phylopic.org/images/phy-first',
                    'media_url': 'https://images.phylopic.org/apis-first.png',
                    'width': 80,
                    'height': 40,
                    'asset_type': 'silhouette',
                    'provider_quality': 20,
                },
                {
                    'provider': 'phylopic',
                    'provider_record_id': 'phy-second',
                    'matched_name': 'Apis mellifera',
                    'matched_rank': 'species',
                    'license_code': 'cc-by',
                    'license_url': 'https://creativecommons.org/licenses/by/4.0/',
                    'attribution': 'Second Artist',
                    'source_page_url': 'https://api.phylopic.org/images/phy-second',
                    'media_url': 'https://images.phylopic.org/apis-second.png',
                    'width': 80,
                    'height': 40,
                    'asset_type': 'silhouette',
                    'provider_quality': 10,
                },
            ],
        }

        def fake_build_providers(args, sources, session=None):
            return DummySession(), None, {
                'phylopic': DummyProvider(candidates),
            }

        def fake_download_media(
            session,
            media_url,
            destination_path,
            **kwargs
        ):
            write_valid_test_media(destination_path)
            return 'downloaded'

        monkeypatch.setattr('nwkit.image.build_providers', fake_build_providers)
        monkeypatch.setattr('nwkit.image.download_media', fake_download_media)

        image_main(make_image_args(
            infile=str(tree_path),
            out_dir=str(out_dir),
            source='phylopic',
            max_per_species=2,
            manifest_out=str(manifest_path),
            attribution_out=str(attribution_path),
        ))

        manifest_rows = read_tsv(manifest_path)
        assert len(manifest_rows) == 2
        assert [row['provider_record_id'] for row in manifest_rows] == [
            'phy-first',
            'phy-second',
        ]
        for row in manifest_rows:
            assert os.path.isfile(metadata_dir / row['local_path'])
        local_file_lines = [
            line.split(': ', 1)[1]
            for line in attribution_path.read_text().splitlines()
            if line.startswith('Local file: ')
        ]
        assert len(local_file_lines) == 2
        assert all(os.path.isfile(metadata_dir / path) for path in local_file_lines)

        from nwkit.cli import main as cli_main

        draw_path = tmp_path / 'ranked-images.svg'
        cli_main([
            'draw',
            '-i',
            str(tree_path),
            '--species-overlap-node-plot',
            'no',
            '--tip-image-manifest',
            str(manifest_path),
            '-o',
            str(draw_path),
        ])

        assert draw_path.read_text(encoding='utf-8').count('<image') == 1

    def test_image_main_does_not_rebuild_providers_for_download_stage(self, monkeypatch, tmp_path):
        tree_path = tmp_path / 'tree.nwk'
        tree_path.write_text('(Apis_mellifera_A);')
        out_dir = tmp_path / 'out'
        call_counter = {'count': 0}

        phylopic_candidates = {
            'Apis mellifera': [{
                'provider': 'phylopic',
                'provider_record_id': 'phy-once',
                'matched_name': 'Apis mellifera',
                'matched_rank': 'species',
                'license_code': 'cc-by',
                'license_url': 'https://creativecommons.org/licenses/by/4.0/',
                'attribution': 'PhyloPic Artist',
                'source_page_url': 'https://api.phylopic.org/images/phy-once',
                'media_url': 'https://images.phylopic.org/apis-once.svg',
                'width': 1000,
                'height': 900,
                'asset_type': 'silhouette',
            }],
        }

        def fake_build_providers(args, sources, session=None):
            call_counter['count'] += 1
            providers = {
                'phylopic': DummyProvider(phylopic_candidates),
            }
            return DummySession(), None, providers

        def fake_download_media(session, media_url, destination_path, cache_path=None, max_download_bytes=None, provider=None, **kwargs):
            write_valid_test_media(destination_path)
            return 'downloaded'

        monkeypatch.setattr('nwkit.image.build_providers', fake_build_providers)
        monkeypatch.setattr('nwkit.image.download_media', fake_download_media)

        image_main(
            make_image_args(
                infile=str(tree_path),
                out_dir=str(out_dir),
                source='phylopic',
            )
        )

        assert call_counter['count'] == 1

    def test_image_main_reuses_same_media_across_species_within_run(self, monkeypatch, tmp_path):
        tree_path = tmp_path / 'tree.nwk'
        tree_path.write_text('(Felis_catus_A,Panthera_leo_B);')
        out_dir = tmp_path / 'out'

        shared_candidates = {
            'Felis catus': [{
                'provider': 'wikimedia',
                'provider_record_id': 'wiki-shared-cat',
                'matched_name': 'Felis catus',
                'matched_rank': 'species',
                'license_code': 'cc-by',
                'license_url': 'https://creativecommons.org/licenses/by/4.0/',
                'attribution': 'Shared Photographer',
                'source_page_url': 'https://commons.wikimedia.org/wiki/File:SharedCatLion.jpg',
                'media_url': 'https://upload.wikimedia.org/shared-cat-lion.jpg',
                'width': 1500,
                'height': 1000,
                'asset_type': 'photo',
            }],
            'Panthera leo': [{
                'provider': 'wikimedia',
                'provider_record_id': 'wiki-shared-lion',
                'matched_name': 'Panthera leo',
                'matched_rank': 'species',
                'license_code': 'cc-by',
                'license_url': 'https://creativecommons.org/licenses/by/4.0/',
                'attribution': 'Shared Photographer',
                'source_page_url': 'https://commons.wikimedia.org/wiki/File:SharedCatLion.jpg',
                'media_url': 'https://upload.wikimedia.org/shared-cat-lion.jpg',
                'width': 1500,
                'height': 1000,
                'asset_type': 'photo',
            }],
        }
        call_counter = {'count': 0}

        def fake_build_providers(args, sources, session=None):
            providers = {
                'wikimedia': DummyProvider(shared_candidates),
            }
            return DummySession(), None, providers

        def fake_download_media(session, media_url, destination_path, cache_path=None, max_download_bytes=None, provider=None, **kwargs):
            call_counter['count'] += 1
            write_valid_test_media(destination_path)
            return 'downloaded'

        monkeypatch.setattr('nwkit.image.build_providers', fake_build_providers)
        monkeypatch.setattr('nwkit.image.download_media', fake_download_media)

        image_main(
            make_image_args(
                infile=str(tree_path),
                out_dir=str(out_dir),
                style='photo',
                source='wikimedia',
            )
        )

        manifest_rows = read_tsv(out_dir / 'manifest.tsv')
        rows_by_species = {row['species_name']: row for row in manifest_rows}

        assert call_counter['count'] == 1
        assert rows_by_species['Felis catus']['status'] == 'downloaded'
        assert rows_by_species['Panthera leo']['status'] == 'reused'
        assert rows_by_species['Felis catus']['local_path'] == rows_by_species['Panthera leo']['local_path']

    def test_image_main_records_resolved_extension_from_download(self, monkeypatch, tmp_path):
        tree_path = tmp_path / 'tree.nwk'
        tree_path.write_text('(Cyanophora_paradoxa_A);')
        out_dir = tmp_path / 'out'

        ncbi_candidates = {
            'Cyanophora paradoxa': [{
                'provider': 'ncbi',
                'provider_record_id': '64365',
                'matched_name': 'Cyanophora paradoxa',
                'matched_rank': 'species',
                'license_code': 'cc-by-sa',
                'license_url': 'https://creativecommons.org/licenses/by-sa/3.0/',
                'attribution': 'NCBI',
                'source_page_url': 'https://www.ncbi.nlm.nih.gov/Taxonomy/taxi/images/4',
                'media_url': 'https://www.ncbi.nlm.nih.gov/Taxonomy/taxi/images/4',
                'width': None,
                'height': None,
                'asset_type': 'photo',
            }],
        }

        def fake_build_providers(args, sources, session=None):
            providers = {
                'ncbi': DummyProvider(ncbi_candidates),
            }
            return DummySession(), None, providers

        def fake_download_media(session, media_url, destination_path, cache_path=None, max_download_bytes=None, provider=None, **kwargs):
            resolved_path = destination_path[:-4] + '.jpg'
            write_valid_test_media(resolved_path)
            return {'status': 'downloaded', 'destination_path': resolved_path}

        monkeypatch.setattr('nwkit.image.build_providers', fake_build_providers)
        monkeypatch.setattr('nwkit.image.download_media', fake_download_media)

        args = make_image_args(
            infile=str(tree_path),
            out_dir=str(out_dir),
            source='ncbi',
            style='photo',
        )
        image_main(args)

        manifest_rows = read_tsv(out_dir / 'manifest.tsv')
        assert manifest_rows[0]['local_path'].endswith('.jpg')
        assert (out_dir / manifest_rows[0]['local_path']).exists()

    def test_image_main_falls_back_to_next_candidate_after_download_error(self, monkeypatch, tmp_path):
        tree_path = tmp_path / 'tree.nwk'
        tree_path.write_text('(Panthera_leo_A);')
        out_dir = tmp_path / 'out'

        inaturalist_candidates = {
            'Panthera leo': [{
                'provider': 'inaturalist',
                'provider_record_id': 'inat-fail',
                'matched_name': 'Panthera leo',
                'matched_rank': 'species',
                'license_code': 'cc-by',
                'license_url': 'https://creativecommons.org/licenses/by/4.0/',
                'attribution': 'Failing Photographer',
                'source_page_url': 'https://www.inaturalist.org/observations/fail',
                'media_url': 'https://static.inaturalist.org/photos/fail/original.jpg',
                'width': 1600,
                'height': 1200,
                'asset_type': 'photo',
            }],
        }
        wikimedia_candidates = {
            'Panthera leo': [{
                'provider': 'wikimedia',
                'provider_record_id': 'wiki-ok',
                'matched_name': 'Panthera leo',
                'matched_rank': 'species',
                'license_code': 'cc-by-sa',
                'license_url': 'https://creativecommons.org/licenses/by-sa/4.0/',
                'attribution': 'Working Photographer',
                'source_page_url': 'https://commons.wikimedia.org/wiki/File:Lion.jpg',
                'media_url': 'https://upload.wikimedia.org/lion.jpg',
                'width': 1400,
                'height': 1000,
                'asset_type': 'photo',
            }],
        }

        def fake_build_providers(args, sources, session=None):
            providers = {
                'inaturalist': DummyProvider(inaturalist_candidates),
                'wikimedia': DummyProvider(wikimedia_candidates),
            }
            return DummySession(), None, providers

        def fake_download_media(session, media_url, destination_path, cache_path=None, max_download_bytes=None, provider=None, **kwargs):
            if media_url.endswith('/fail/original.jpg'):
                raise requests.ConnectionError('transient failure')
            write_valid_test_media(destination_path)
            return 'downloaded'

        monkeypatch.setattr('nwkit.image.build_providers', fake_build_providers)
        monkeypatch.setattr('nwkit.image.download_media', fake_download_media)

        args = make_image_args(
            infile=str(tree_path),
            out_dir=str(out_dir),
            style='photo',
            source='inaturalist,wikimedia',
        )
        image_main(args)

        manifest_rows = read_tsv(out_dir / 'manifest.tsv')
        unmatched_rows = read_tsv(out_dir / 'unmatched.tsv')

        assert len(manifest_rows) == 1
        assert manifest_rows[0]['provider'] == 'wikimedia'
        assert unmatched_rows == []

    def test_materialize_value_error_falls_back_to_next_candidate(self):
        candidates = [
            {
                'provider': 'inaturalist',
                'provider_record_id': 'invalid',
                'matched_name': 'Panthera leo',
                'matched_rank': 'species',
                'license_code': 'cc-by',
                'license_url': 'https://creativecommons.org/licenses/by/4.0/',
                'attribution': 'Invalid Photographer',
                'source_page_url': 'https://example.org/invalid',
                'media_url': 'https://example.org/invalid.jpg',
                'asset_type': 'photo',
                'score': 2.0,
            },
            {
                'provider': 'wikimedia',
                'provider_record_id': 'valid',
                'matched_name': 'Panthera leo',
                'matched_rank': 'species',
                'license_code': 'cc-by',
                'license_url': 'https://creativecommons.org/licenses/by/4.0/',
                'attribution': 'Valid Photographer',
                'source_page_url': 'https://example.org/valid',
                'media_url': 'https://example.org/valid.jpg',
                'asset_type': 'photo',
                'score': 1.0,
            },
        ]

        class Materializer:
            def __init__(self):
                self.calls = list()

            def materialize(self, species_name, candidate):
                self.calls.append(candidate['provider_record_id'])
                if candidate['provider_record_id'] == 'invalid':
                    raise ValueError('decoded image is invalid')
                return {
                    'local_path': 'images/valid.jpg',
                    'download_status': 'downloaded',
                }

        materializer = Materializer()
        manifest_rows, selected_assets, unmatched_rows = (
            image_module.process_species_assets(
                species_name='Panthera leo',
                leaf_names=['Panthera_leo_A'],
                candidates=candidates,
                provider_errors=[],
                args=make_image_args(max_per_species=1),
                materializer=materializer,
            )
        )

        assert materializer.calls == ['invalid', 'valid']
        assert [row['provider'] for row in manifest_rows] == ['wikimedia']
        assert selected_assets == manifest_rows
        assert unmatched_rows == []

    def test_image_main_reports_download_error_when_all_candidates_fail(self, monkeypatch, tmp_path):
        tree_path = tmp_path / 'tree.nwk'
        tree_path.write_text('(Panthera_leo_A);')
        out_dir = tmp_path / 'out'

        inaturalist_candidates = {
            'Panthera leo': [{
                'provider': 'inaturalist',
                'provider_record_id': 'inat-fail',
                'matched_name': 'Panthera leo',
                'matched_rank': 'species',
                'license_code': 'cc-by',
                'license_url': 'https://creativecommons.org/licenses/by/4.0/',
                'attribution': 'Failing Photographer',
                'source_page_url': 'https://www.inaturalist.org/observations/fail',
                'media_url': 'https://static.inaturalist.org/photos/fail/original.jpg',
                'width': 1600,
                'height': 1200,
                'asset_type': 'photo',
            }],
        }

        def fake_build_providers(args, sources, session=None):
            providers = {
                'inaturalist': DummyProvider(inaturalist_candidates),
            }
            return DummySession(), None, providers

        def fake_download_media(session, media_url, destination_path, cache_path=None, max_download_bytes=None, provider=None, **kwargs):
            raise requests.ConnectionError('transient failure')

        monkeypatch.setattr('nwkit.image.build_providers', fake_build_providers)
        monkeypatch.setattr('nwkit.image.download_media', fake_download_media)

        args = make_image_args(
            infile=str(tree_path),
            out_dir=str(out_dir),
            style='photo',
            source='inaturalist',
        )
        image_main(args)

        manifest_rows = read_tsv(out_dir / 'manifest.tsv')
        unmatched_rows = read_tsv(out_dir / 'unmatched.tsv')

        assert manifest_rows == []
        assert unmatched_rows[0]['reason'] == 'download_error'
        assert 'inaturalist download failed for Panthera leo' in unmatched_rows[0]['details']


class TestDownloadMedia:
    def test_download_media_uses_media_accept_header(self, tmp_path):
        destination = tmp_path / 'out' / 'image.bin'
        request_kwargs = {}

        class JPEGResponse:
            headers = {'Content-Type': 'image/jpeg'}
            url = 'https://example.org/download'

            def raise_for_status(self):
                return None

            def iter_content(self, chunk_size=65536):
                yield b'\xff\xd8\xff\xe0fakejpeg'

            def close(self):
                return None

        class JPEGSession:
            def get(self, *args, **kwargs):
                request_kwargs.update(kwargs)
                return JPEGResponse()

        download_media(
            JPEGSession(),
            'https://example.org/download',
            str(destination),
        )

        assert request_kwargs['headers']['Accept'].startswith('image/*')
        assert request_kwargs['headers']['Accept'] != 'application/json'
        assert image_module.HTTP_HEADERS is image_module.API_HTTP_HEADERS
        assert image_module.API_HTTP_HEADERS['Accept'] == 'application/json'

    def test_download_media_prefers_shared_cache_when_present(self, tmp_path):
        destination = tmp_path / 'out' / 'image.svg'
        cache_path = tmp_path / 'cache' / 'image.svg'
        cache_path.parent.mkdir(parents=True)
        cache_path.write_bytes(b'<svg xmlns="http://www.w3.org/2000/svg"></svg>')

        class FailingSession:
            def get(self, *args, **kwargs):
                raise AssertionError('network should not be used when cache exists')

        result = download_media(FailingSession(), 'https://images.example.org/item.svg', str(destination), cache_path=str(cache_path))

        assert result['status'] == 'cached'
        assert destination.read_bytes().startswith(b'<svg')

    def test_cached_destination_uses_umask_and_preserves_existing_mode(
        self,
        tmp_path,
    ):
        destination = tmp_path / 'out' / 'image.svg'
        cache_path = tmp_path / 'cache' / 'image.svg'
        cache_path.parent.mkdir(parents=True)
        cache_path.write_bytes(b'<svg xmlns="http://www.w3.org/2000/svg"></svg>')

        class FailingSession:
            def get(self, *args, **kwargs):
                raise AssertionError('network should not be used when cache exists')

        previous_umask = os.umask(0o027)
        try:
            download_media(
                FailingSession(),
                'https://images.example.org/item.svg',
                str(destination),
                cache_path=str(cache_path),
            )
        finally:
            os.umask(previous_umask)
        assert stat.S_IMODE(destination.stat().st_mode) == 0o640

        destination.chmod(0o664)
        download_media(
            FailingSession(),
            'https://images.example.org/item.svg',
            str(destination),
            cache_path=str(cache_path),
            reuse_destination=False,
        )
        assert stat.S_IMODE(destination.stat().st_mode) == 0o664

    def test_download_media_removes_partial_temp_file_on_error(self, tmp_path):
        destination = tmp_path / 'out' / 'image.svg'
        destination.parent.mkdir(parents=True)

        class FailingSession:
            def get(self, *args, **kwargs):
                raise requests.ConnectionError('connection dropped')

        with pytest.raises(requests.ConnectionError, match='connection dropped'):
            download_media(FailingSession(), 'https://images.example.org/item.svg', str(destination))

        assert not destination.exists()
        assert not (tmp_path / 'out' / 'image.svg.tmp').exists()

    def test_download_media_uses_content_type_to_resolve_extension(self, tmp_path):
        destination = tmp_path / 'out' / 'image.bin'

        class JPEGResponse:
            def __init__(self):
                self.headers = {'Content-Type': 'image/jpeg'}
                self.url = 'https://example.org/download'

            def raise_for_status(self):
                return None

            def iter_content(self, chunk_size=65536):
                yield b'\xff\xd8\xff\xe0fakejpeg'

            def close(self):
                return None

        class JPEGSession:
            def get(self, *args, **kwargs):
                return JPEGResponse()

        result = download_media(JPEGSession(), 'https://example.org/download', str(destination))

        assert result['status'] == 'downloaded'
        assert result['destination_path'].endswith('.jpg')
        assert os.path.exists(result['destination_path'])
        assert not os.path.exists(destination)

    def test_direct_download_uses_umask_adjusted_mode(self, tmp_path):
        destination = tmp_path / 'out' / 'image.bin'

        class JPEGResponse:
            headers = {'Content-Type': 'image/jpeg'}
            url = 'https://example.org/download'

            def raise_for_status(self):
                return None

            def iter_content(self, chunk_size=65536):
                yield b'\xff\xd8\xff\xe0fakejpeg'

            def close(self):
                return None

        class JPEGSession:
            def get(self, *args, **kwargs):
                return JPEGResponse()

        previous_umask = os.umask(0o027)
        try:
            result = download_media(
                JPEGSession(),
                'https://example.org/download',
                str(destination),
            )
        finally:
            os.umask(previous_umask)

        assert stat.S_IMODE(os.stat(result['destination_path']).st_mode) == 0o640

    def test_download_media_reuses_existing_cache_variant_extension(self, tmp_path):
        destination = tmp_path / 'out' / 'image.bin'
        cache_path = tmp_path / 'cache' / 'image.bin'
        cache_variant = tmp_path / 'cache' / 'image.jpg'
        cache_variant.parent.mkdir(parents=True)
        cache_variant.write_bytes(b'\xff\xd8\xff\xe0jpeg-data')

        class FailingSession:
            def get(self, *args, **kwargs):
                raise AssertionError('network should not be used when cache variant exists')

        result = download_media(FailingSession(), 'https://example.org/download', str(destination), cache_path=str(cache_path))

        assert result['status'] == 'cached'
        assert result['destination_path'].endswith('.jpg')
        assert (tmp_path / 'out' / 'image.jpg').read_bytes().startswith(b'\xff\xd8\xff')

    def test_download_media_rejects_html_returned_from_image_url(self, tmp_path):
        destination = tmp_path / 'out' / 'error.jpg'

        class HTMLResponse:
            headers = {'Content-Type': 'text/html'}
            url = 'https://example.org/error.jpg'

            def raise_for_status(self):
                return None

            def iter_content(self, chunk_size=65536):
                yield b'<html><body>rate limited</body></html>'

            def close(self):
                return None

        class HTMLSession:
            def get(self, *args, **kwargs):
                return HTMLResponse()

        with pytest.raises(MediaDownloadError, match='non-image Content-Type'):
            download_media(HTMLSession(), 'https://example.org/error.jpg', str(destination))

        assert not destination.exists()
        assert list(destination.parent.glob('*.tmp')) == []

    def test_download_media_enforces_content_length_limit(self, tmp_path):
        destination = tmp_path / 'out' / 'large.jpg'

        class LargeResponse:
            headers = {'Content-Type': 'image/jpeg', 'Content-Length': '1000'}

            def raise_for_status(self):
                return None

            def iter_content(self, chunk_size=65536):
                raise AssertionError('oversized response should be rejected before streaming')

            def close(self):
                return None

        class LargeSession:
            def get(self, *args, **kwargs):
                return LargeResponse()

        with pytest.raises(MediaDownloadError, match='exceeding --max-download-bytes'):
            download_media(
                LargeSession(),
                'https://example.org/large.jpg',
                str(destination),
                max_download_bytes=100,
            )

        assert not destination.exists()
        assert list(destination.parent.glob('*.tmp')) == []


class TestImagePostprocessing:
    def test_processed_raster_uses_umask_and_preserves_existing_mode(
        self,
        tmp_path,
    ):
        Image = pytest.importorskip('PIL.Image')
        source = tmp_path / 'image.png'
        destination = tmp_path / 'image.jpg'
        Image.new('RGB', (40, 20), 'black').save(source)

        previous_umask = os.umask(0o027)
        try:
            result = postprocess_media_file(
                str(source),
                make_image_args(output_format='jpg', max_edge=16),
            )
        finally:
            os.umask(previous_umask)
        assert result == str(destination)
        assert stat.S_IMODE(destination.stat().st_mode) == 0o640

        Image.new('RGB', (40, 20), 'black').save(source)
        destination.chmod(0o664)
        postprocess_media_file(
            str(source),
            make_image_args(output_format='jpg', max_edge=12),
        )
        assert stat.S_IMODE(destination.stat().st_mode) == 0o664

    def test_postprocess_media_file_resizes_and_pads_raster(self, tmp_path):
        Image = pytest.importorskip('PIL.Image')
        source = tmp_path / 'image.png'
        Image.new('RGB', (40, 20), 'black').save(source)

        result = postprocess_media_file(
            str(source),
            make_image_args(
                output_format='png',
                max_edge=16,
                canvas='square',
                background='white',
            ),
        )

        with Image.open(result) as processed:
            assert processed.size == (16, 16)
        assert result.endswith('.png')

    def test_postprocess_media_file_trims_white_border(self, tmp_path):
        Image = pytest.importorskip('PIL.Image')
        source = tmp_path / 'trim.png'
        image = Image.new('RGB', (20, 20), 'white')
        for x in range(5, 15):
            for y in range(7, 13):
                image.putpixel((x, y), (0, 0, 0))
        image.save(source)

        result = postprocess_media_file(
            str(source),
            make_image_args(
                output_format='png',
                trim='white',
            ),
        )

        with Image.open(result) as processed:
            assert processed.size == (10, 6)

    def test_postprocess_media_file_trim_shape_square_center_crops_trimmed_result(self, tmp_path):
        Image = pytest.importorskip('PIL.Image')
        source = tmp_path / 'trim-square.png'
        image = Image.new('RGBA', (24, 24), (0, 0, 0, 0))
        for x in range(8, 14):
            for y in range(4, 20):
                image.putpixel((x, y), (0, 0, 0, 255))
        image.save(source)

        result = postprocess_media_file(
            str(source),
            make_image_args(
                output_format='png',
                trim='transparent',
                trim_shape='square',
                background='transparent',
            ),
        )

        with Image.open(result) as processed:
            assert processed.size == (6, 6)
            assert processed.mode == 'RGBA'
            assert processed.getchannel('A').getextrema() == (255, 255)

    def test_postprocess_media_file_trims_semantic_largest_rgb_component(self, tmp_path):
        Image = pytest.importorskip('PIL.Image')
        source = tmp_path / 'semantic-rgb.png'
        image = Image.new('RGB', (40, 24), 'white')
        for x in range(3, 9):
            for y in range(2, 20):
                image.putpixel((x, y), (0, 0, 0))
        for x in range(20, 35):
            for y in range(5, 17):
                image.putpixel((x, y), (0, 128, 0))
        image.save(source)

        result = postprocess_media_file(
            str(source),
            make_image_args(
                output_format='png',
                trim='semantic',
            ),
        )

        with Image.open(result) as processed:
            assert processed.size == (15, 12)
            assert processed.getpixel((0, 0)) == (0, 128, 0)

    def test_postprocess_media_file_trims_semantic_largest_alpha_component(self, tmp_path):
        Image = pytest.importorskip('PIL.Image')
        source = tmp_path / 'semantic-alpha.png'
        image = Image.new('RGBA', (30, 24), (0, 0, 0, 0))
        for x in range(2, 8):
            for y in range(3, 21):
                image.putpixel((x, y), (255, 0, 0, 255))
        for x in range(14, 26):
            for y in range(6, 16):
                image.putpixel((x, y), (0, 0, 255, 255))
        image.save(source)

        result = postprocess_media_file(
            str(source),
            make_image_args(
                output_format='png',
                trim='semantic',
                trim_shape='square',
                background='transparent',
            ),
        )

        with Image.open(result) as processed:
            assert processed.size == (10, 10)
            assert processed.mode == 'RGBA'
            assert processed.getpixel((0, 0)) == (0, 0, 255, 255)

    def test_postprocess_media_file_trim_shape_square_without_trim_center_crops_full_image(self, tmp_path):
        Image = pytest.importorskip('PIL.Image')
        source = tmp_path / 'full-square.png'
        image = Image.new('RGB', (20, 10))
        for x in range(20):
            for y in range(10):
                image.putpixel((x, y), (x, 0, 0))
        image.save(source)

        result = postprocess_media_file(
            str(source),
            make_image_args(
                output_format='png',
                trim='off',
                trim_shape='square',
                background='white',
            ),
        )

        with Image.open(result) as processed:
            assert processed.size == (10, 10)
            assert processed.getpixel((0, 0)) == (5, 0, 0)
            assert processed.getpixel((9, 0)) == (14, 0, 0)

    def test_postprocess_media_file_rasterizes_svg_when_requested(self, monkeypatch, tmp_path):
        Image = pytest.importorskip('PIL.Image')
        source = tmp_path / 'shape.svg'
        source.write_text('<svg xmlns="http://www.w3.org/2000/svg" width="20" height="10"></svg>')

        monkeypatch.setattr(
            'nwkit.image.rasterize_svg_to_image',
            lambda source_path, max_edge=None: Image.new('RGBA', (20, 10), (0, 0, 0, 255)),
        )

        result = postprocess_media_file(
            str(source),
            make_image_args(
                output_format='png',
                max_edge=12,
                canvas='square',
                background='white',
            ),
        )

        with Image.open(result) as processed:
            assert processed.size == (12, 12)
        assert result.endswith('.png')

    def test_postprocess_media_file_requires_cairosvg_for_svg_conversion(self, monkeypatch, tmp_path):
        source = tmp_path / 'shape.svg'
        source.write_text('<svg xmlns="http://www.w3.org/2000/svg" width="20" height="10"></svg>')

        monkeypatch.setattr(
            'nwkit.image.load_cairosvg_module',
            lambda: (_ for _ in ()).throw(RuntimeError('SVG image post-processing requires the optional CairoSVG dependency. Install optional image-processing dependencies with: pip install "nwkit[image]"')),
        )

        with pytest.raises(RuntimeError, match='CairoSVG dependency'):
            postprocess_media_file(
                str(source),
                make_image_args(
                    output_format='png',
                ),
            )

    def test_image_main_skips_optional_processing_deps_when_not_requested(self, monkeypatch, tmp_path):
        tree_path = tmp_path / 'tree.nwk'
        tree_path.write_text('(Apis_mellifera_A);')
        out_dir = tmp_path / 'out'

        phylopic_candidates = {
            'Apis mellifera': [{
                'provider': 'phylopic',
                'provider_record_id': 'phy-plain',
                'matched_name': 'Apis mellifera',
                'matched_rank': 'species',
                'license_code': 'cc-by',
                'license_url': 'https://creativecommons.org/licenses/by/4.0/',
                'attribution': 'PhyloPic Artist',
                'source_page_url': 'https://api.phylopic.org/images/phy-plain',
                'media_url': 'https://images.phylopic.org/apis.svg',
                'width': 1000,
                'height': 900,
                'asset_type': 'silhouette',
            }],
        }

        def fake_build_providers(args, sources, session=None):
            providers = {
                'phylopic': DummyProvider(phylopic_candidates),
            }
            return DummySession(), None, providers

        def fake_download_media(session, media_url, destination_path, cache_path=None, max_download_bytes=None, provider=None, **kwargs):
            with open(destination_path, 'wb') as handle:
                handle.write(b'<svg xmlns="http://www.w3.org/2000/svg"></svg>')
            return {'status': 'downloaded', 'destination_path': destination_path}

        monkeypatch.setattr('nwkit.image.build_providers', fake_build_providers)
        monkeypatch.setattr('nwkit.image.download_media', fake_download_media)
        monkeypatch.setattr(
            'nwkit.image.load_pillow_modules',
            lambda: (_ for _ in ()).throw(AssertionError('Pillow should not be loaded without post-processing options')),
        )

        image_main(
            make_image_args(
                infile=str(tree_path),
                out_dir=str(out_dir),
                source='phylopic',
            )
        )

        manifest_rows = read_tsv(out_dir / 'manifest.tsv')
        assert manifest_rows[0]['local_path'].endswith('.svg')

    def test_image_main_requires_pillow_when_postprocessing_is_requested(self, monkeypatch, tmp_path):
        tree_path = tmp_path / 'tree.nwk'
        tree_path.write_text('(Apis_mellifera_A);')
        out_dir = tmp_path / 'out'

        inaturalist_candidates = {
            'Apis mellifera': [{
                'provider': 'inaturalist',
                'provider_record_id': 'inat-plain',
                'matched_name': 'Apis mellifera',
                'matched_rank': 'species',
                'license_code': 'cc-by',
                'license_url': 'https://creativecommons.org/licenses/by/4.0/',
                'attribution': 'Photographer',
                'source_page_url': 'https://www.inaturalist.org/observations/1',
                'media_url': 'https://static.inaturalist.org/apis.jpg',
                'width': 1000,
                'height': 900,
                'asset_type': 'photo',
            }],
        }

        def fake_build_providers(args, sources, session=None):
            providers = {
                'inaturalist': DummyProvider(inaturalist_candidates),
            }
            return DummySession(), None, providers

        def fake_download_media(session, media_url, destination_path, cache_path=None, max_download_bytes=None, provider=None, **kwargs):
            write_valid_test_media(destination_path)
            return {'status': 'downloaded', 'destination_path': destination_path}

        monkeypatch.setattr('nwkit.image.build_providers', fake_build_providers)
        monkeypatch.setattr('nwkit.image.download_media', fake_download_media)
        monkeypatch.setattr(
            'nwkit.image.load_pillow_modules',
            lambda: (_ for _ in ()).throw(RuntimeError('Image post-processing requires the optional Pillow dependency. Install optional image-processing dependencies with: pip install "nwkit[image]"')),
        )

        with pytest.raises(RuntimeError, match='Pillow dependency'):
            image_main(
                make_image_args(
                    infile=str(tree_path),
                    out_dir=str(out_dir),
                    source='inaturalist',
                    output_format='png',
                )
            )

    def test_image_main_requires_pillow_when_trim_shape_square_is_requested(self, monkeypatch, tmp_path):
        tree_path = tmp_path / 'tree.nwk'
        tree_path.write_text('(Apis_mellifera_A);')
        out_dir = tmp_path / 'out'

        phylopic_candidates = {
            'Apis mellifera': [{
                'provider': 'phylopic',
                'provider_record_id': 'phy-plain',
                'matched_name': 'Apis mellifera',
                'matched_rank': 'species',
                'license_code': 'cc-by',
                'license_url': 'https://creativecommons.org/licenses/by/4.0/',
                'attribution': 'PhyloPic Artist',
                'source_page_url': 'https://api.phylopic.org/images/phy-plain',
                'media_url': 'https://images.phylopic.org/apis.svg',
                'width': 1000,
                'height': 900,
                'asset_type': 'silhouette',
            }],
        }

        def fake_build_providers(args, sources, session=None):
            providers = {
                'phylopic': DummyProvider(phylopic_candidates),
            }
            return DummySession(), None, providers

        def fake_download_media(session, media_url, destination_path, cache_path=None, max_download_bytes=None, provider=None, **kwargs):
            with open(destination_path, 'wb') as handle:
                handle.write(b'<svg xmlns="http://www.w3.org/2000/svg"></svg>')
            return {'status': 'downloaded', 'destination_path': destination_path}

        monkeypatch.setattr('nwkit.image.build_providers', fake_build_providers)
        monkeypatch.setattr('nwkit.image.download_media', fake_download_media)
        monkeypatch.setattr('nwkit.image.rasterize_svg_to_image', lambda source_path, max_edge=None: object())
        monkeypatch.setattr(
            'nwkit.image.load_pillow_modules',
            lambda: (_ for _ in ()).throw(RuntimeError('Image post-processing requires the optional Pillow dependency. Install optional image-processing dependencies with: pip install "nwkit[image]"')),
        )

        with pytest.raises(RuntimeError, match='Pillow dependency'):
            image_main(
                make_image_args(
                    infile=str(tree_path),
                    out_dir=str(out_dir),
                    source='phylopic',
                    trim_shape='square',
                )
            )

    def test_image_main_requires_pillow_when_trim_semantic_is_requested(self, monkeypatch, tmp_path):
        tree_path = tmp_path / 'tree.nwk'
        tree_path.write_text('(Apis_mellifera_A);')
        out_dir = tmp_path / 'out'

        inaturalist_candidates = {
            'Apis mellifera': [{
                'provider': 'inaturalist',
                'provider_record_id': 'inat-semantic',
                'matched_name': 'Apis mellifera',
                'matched_rank': 'species',
                'license_code': 'cc-by',
                'license_url': 'https://creativecommons.org/licenses/by/4.0/',
                'attribution': 'Photographer',
                'source_page_url': 'https://www.inaturalist.org/observations/1',
                'media_url': 'https://static.inaturalist.org/apis.jpg',
                'width': 1000,
                'height': 900,
                'asset_type': 'photo',
            }],
        }

        def fake_build_providers(args, sources, session=None):
            providers = {
                'inaturalist': DummyProvider(inaturalist_candidates),
            }
            return DummySession(), None, providers

        def fake_download_media(session, media_url, destination_path, cache_path=None, max_download_bytes=None, provider=None, **kwargs):
            write_valid_test_media(destination_path)
            return {'status': 'downloaded', 'destination_path': destination_path}

        monkeypatch.setattr('nwkit.image.build_providers', fake_build_providers)
        monkeypatch.setattr('nwkit.image.download_media', fake_download_media)
        monkeypatch.setattr(
            'nwkit.image.load_pillow_modules',
            lambda: (_ for _ in ()).throw(RuntimeError('Image post-processing requires the optional Pillow dependency. Install optional image-processing dependencies with: pip install "nwkit[image]"')),
        )

        with pytest.raises(RuntimeError, match='Pillow dependency'):
            image_main(
                make_image_args(
                    infile=str(tree_path),
                    out_dir=str(out_dir),
                    source='inaturalist',
                    trim='semantic',
                )
            )


class TestImageSecurityLimits:
    @pytest.mark.parametrize(
        'allowed_hosts',
        [None, ('attacker.example',)],
        ids=['arbitrary-host', 'allowlisted-host'],
    )
    def test_real_session_pins_the_single_validated_dns_result(
        self,
        monkeypatch,
        allowed_hosts,
    ):
        resolver_calls = []
        pinned_calls = []
        public_address = '93.184.216.34'

        def alternating_resolution(hostname):
            resolver_calls.append(hostname)
            if len(resolver_calls) == 1:
                return {public_address}
            return {'127.0.0.1'}

        class Response:
            status_code = 200
            headers = {}
            url = 'https://attacker.example/image.png'

            def close(self):
                return None

        def fake_pinned_request(session, method, url, address, **kwargs):
            pinned_calls.append((method, url, address))
            return Response()

        monkeypatch.setattr(
            image_module,
            'resolve_hostname_addresses',
            alternating_resolution,
        )
        monkeypatch.setattr(
            image_module,
            '_request_pinned_address',
            fake_pinned_request,
        )
        session = requests.Session()
        try:
            response = image_module.safe_external_request(
                session,
                'get',
                'https://attacker.example/image.png',
                allowed_hosts=allowed_hosts,
            )
        finally:
            session.close()

        assert response.status_code == 200
        assert resolver_calls == ['attacker.example']
        assert pinned_calls == [
            ('get', 'https://attacker.example/image.png', public_address)
        ]

    def test_pinned_adapter_connects_to_ip_with_original_tls_hostname(self):
        adapter = image_module._PinnedHTTPSAdapter(
            address='93.184.216.34',
            tls_hostname='example.com',
        )
        request = requests.Request('GET', 'https://example.com/image.png').prepare()
        try:
            pool = adapter.get_connection_with_tls_context(
                request,
                verify=True,
            )
            assert pool.host == '93.184.216.34'
            assert pool.conn_kw['server_hostname'] == 'example.com'
            assert pool.assert_hostname == 'example.com'
        finally:
            adapter.close()

    def test_pinned_request_rejects_proxy_configuration(self):
        session = requests.Session()
        session.trust_env = False
        session.proxies = {'https': 'http://proxy.example:8080'}
        try:
            with pytest.raises(MediaDownloadError, match='Proxy'):
                image_module._request_pinned_address(
                    session,
                    'get',
                    'https://example.com/image.png',
                    '93.184.216.34',
                )
        finally:
            session.close()

    def test_cross_origin_redirect_drops_session_and_request_credentials(
        self,
        monkeypatch,
    ):
        calls = list()

        class Response:
            def __init__(self, status_code, url, location=None):
                self.status_code = status_code
                self.url = url
                self.headers = {} if location is None else {'Location': location}

            def close(self):
                return None

        responses = iter([
            Response(
                302,
                'https://origin.example/start',
                'https://other.example/next',
            ),
            Response(200, 'https://other.example/next'),
        ])

        def fake_get(pinned_session, url, **kwargs):
            effective_headers = {
                str(key).lower(): value
                for key, value in pinned_session.headers.items()
            }
            effective_headers.update({
                str(key).lower(): value
                for key, value in (kwargs.get('headers') or {}).items()
            })
            calls.append({
                'url': url,
                'auth': pinned_session.auth,
                'cert': pinned_session.cert,
                'headers': effective_headers,
                'cookies': pinned_session.cookies.get_dict(),
                'params': dict(pinned_session.params),
                'request_auth': kwargs.get('auth'),
                'request_cert': kwargs.get('cert'),
                'request_cookies': kwargs.get('cookies'),
                'request_params': kwargs.get('params'),
            })
            return next(responses)

        monkeypatch.setattr(requests.Session, 'get', fake_get)
        monkeypatch.setattr(
            image_module,
            '_resolve_public_addresses',
            lambda hostname, url: ('93.184.216.34',),
        )
        session = requests.Session()
        session.trust_env = False
        session.auth = ('session-user', 'session-secret')
        session.cert = '/tmp/session-client.pem'
        session.params = {'session-api-key': 'secret'}
        session.headers.update({
            'Authorization': 'Bearer session-secret',
            'Cookie': 'session-cookie=secret',
            'Proxy-Authorization': 'Basic proxy-secret',
            'X-Api-Key': 'session-api-secret',
        })
        session.cookies.set('ambient-cookie', 'secret')
        try:
            response = image_module.safe_external_request(
                session,
                'get',
                'https://origin.example/start',
                headers={
                    'Authorization': 'Bearer request-secret',
                    'Cookie': 'request-cookie=secret',
                    'Proxy-Authorization': 'Basic request-proxy-secret',
                    'X-Api-Key': 'request-api-secret',
                    'Accept-Language': 'ja',
                },
                auth=('request-user', 'request-secret'),
                cert='/tmp/request-client.pem',
                cookies={'request-cookie': 'secret'},
                params={'request-api-key': 'secret'},
            )
            response.close()
        finally:
            session.close()

        assert calls[0]['auth'] == ('session-user', 'session-secret')
        assert calls[0]['cert'] == '/tmp/session-client.pem'
        assert calls[0]['params'] == {'session-api-key': 'secret'}
        assert calls[0]['request_auth'] == ('request-user', 'request-secret')
        assert calls[0]['request_cert'] == '/tmp/request-client.pem'
        assert calls[0]['request_cookies'] == {'request-cookie': 'secret'}
        assert calls[0]['request_params'] == {'request-api-key': 'secret'}
        assert calls[0]['cookies'] == {'ambient-cookie': 'secret'}
        assert calls[1]['auth'] is None
        assert calls[1]['cert'] is None
        assert calls[1]['params'] == {}
        assert calls[1]['request_auth'] is None
        assert calls[1]['request_cert'] is None
        assert calls[1]['request_cookies'] is None
        assert calls[1]['request_params'] is None
        assert calls[1]['cookies'] == {}
        assert calls[1]['headers']['accept-language'] == 'ja'
        assert 'authorization' not in calls[1]['headers']
        assert 'cookie' not in calls[1]['headers']
        assert 'proxy-authorization' not in calls[1]['headers']
        assert 'x-api-key' not in calls[1]['headers']

    def test_pinned_redirect_merges_response_cookies_into_caller_session(
        self,
        monkeypatch,
    ):
        cookies_seen = list()

        class Response:
            def __init__(self, status_code, url, location=None):
                self.status_code = status_code
                self.url = url
                self.headers = {} if location is None else {'Location': location}

            def close(self):
                return None

        def fake_get(pinned_session, url, **kwargs):
            cookies_seen.append(pinned_session.cookies.get_dict(
                domain='cookie.example',
                path='/',
            ))
            if len(cookies_seen) == 1:
                pinned_session.cookies.set(
                    'route',
                    'blue',
                    domain='cookie.example',
                    path='/',
                )
                return Response(302, url, '/next')
            return Response(200, url)

        monkeypatch.setattr(requests.Session, 'get', fake_get)
        monkeypatch.setattr(
            image_module,
            '_resolve_public_addresses',
            lambda hostname, url: ('93.184.216.34',),
        )
        session = requests.Session()
        session.trust_env = False
        try:
            response = image_module.safe_external_request(
                session,
                'get',
                'https://cookie.example/start',
            )
            response.close()
            assert session.cookies.get(
                'route',
                domain='cookie.example',
                path='/',
            ) == 'blue'
        finally:
            session.close()

        assert cookies_seen == [{}, {'route': 'blue'}]

    def test_custom_session_without_redirect_control_is_not_retried(self):
        calls = list()

        class LegacySession:
            def get(self, url, **kwargs):
                calls.append(dict(kwargs))
                if 'allow_redirects' in kwargs:
                    raise TypeError(
                        "get() got an unexpected keyword argument 'allow_redirects'"
                    )
                raise AssertionError('unsafe fallback request was attempted')

        with pytest.raises(MediaDownloadError, match='allow_redirects=False'):
            image_module.safe_external_request(
                LegacySession(),
                'get',
                'https://public.example/start',
            )

        assert len(calls) == 1
        assert calls[0]['allow_redirects'] is False

    def test_mixed_public_and_private_dns_answers_are_rejected(self, monkeypatch):
        monkeypatch.setattr(
            image_module,
            'resolve_hostname_addresses',
            lambda hostname: {'93.184.216.34', '127.0.0.1'},
        )
        monkeypatch.setattr(
            image_module,
            '_request_pinned_address',
            lambda *args, **kwargs: pytest.fail('request must not be sent'),
        )
        session = requests.Session()
        try:
            with pytest.raises(MediaDownloadError, match='non-public address'):
                image_module.safe_external_request(
                    session,
                    'get',
                    'https://mixed.example/image.png',
                )
        finally:
            session.close()

    def test_rejects_private_redirect_before_following_it(self):
        calls = list()

        class RedirectResponse:
            status_code = 302
            headers = {'Location': 'https://127.0.0.1/private'}
            url = 'https://example.org/start'

            def close(self):
                return None

        class RedirectSession:
            def get(self, url, **kwargs):
                calls.append(url)
                return RedirectResponse()

        with pytest.raises(MediaDownloadError, match='non-public address'):
            image_module.safe_external_request(
                RedirectSession(),
                'get',
                'https://example.org/start',
            )
        assert calls == ['https://example.org/start']

    def test_download_rejects_provider_host_outside_allowlist(self, tmp_path):
        class UnusedSession:
            def get(self, *args, **kwargs):
                raise AssertionError('request should not be sent')

        with pytest.raises(MediaDownloadError, match='not allowed'):
            download_media(
                UnusedSession(),
                'https://example.org/silhouette.svg',
                str(tmp_path / 'image.svg'),
                provider='phylopic',
            )

    def test_provider_json_body_is_bounded_while_streaming(self):
        class LargeResponse:
            headers = {}

            def iter_content(self, chunk_size):
                yield b'{"value":"'
                yield b'x' * 32
                yield b'"}'

        with pytest.raises(MediaDownloadError, match='exceeded'):
            image_module.bounded_response_json(LargeResponse(), max_bytes=16)

    def test_cached_candidate_with_insecure_url_is_discarded(self, tmp_path):
        cache_path = tmp_path / 'query.json'
        image_module.write_cached_provider_candidates(
            str(cache_path),
            [{
                'provider': 'openverse',
                'media_url': 'http://example.org/image.jpg',
            }],
            fetch_limit=10,
        )
        assert image_module.load_cached_provider_candidates(str(cache_path)) is None

    @pytest.mark.parametrize(
        'loader',
        [
            image_module.load_cached_provider_candidates,
            image_module.load_cached_bioicons_catalog,
        ],
    )
    def test_non_object_json_cache_is_treated_as_a_cache_miss(
        self,
        loader,
        tmp_path,
    ):
        cache_path = tmp_path / 'cache.json'
        cache_path.write_text('[]')

        assert loader(str(cache_path)) is None

    def test_default_svg_output_rejects_external_content(self, tmp_path):
        source = tmp_path / 'unsafe.svg'
        source.write_text(
            '<svg xmlns="http://www.w3.org/2000/svg" width="20" height="10">'
            '<image href="https://example.org/tracker.png"/>'
            '</svg>'
        )
        with pytest.raises(MediaDownloadError, match='forbidden|external'):
            postprocess_media_file(str(source), make_image_args())

    def test_default_svg_output_allows_safe_inline_and_style_tag_css(self, tmp_path):
        source = tmp_path / 'safe-css.svg'
        source.write_text(
            '<svg xmlns="http://www.w3.org/2000/svg" width="20" height="10">'
            '<defs>'
            '<filter id="shadow"><feOffset dx="1" dy="1"/></filter>'
            '<style>.safe { fill:red; stroke:#000; filter:url(#shadow); }</style>'
            '</defs>'
            '<rect class="safe" width="20" height="10" '
            'style="fill:blue;stroke:red"/>'
            '</svg>'
        )

        result = postprocess_media_file(str(source), make_image_args())

        assert result == str(source)

    @pytest.mark.parametrize(
        'unsafe_css',
        [
            '<rect style="fill:url(https://attacker.example/pattern.svg)"/>',
            '<rect style="fill:url(images/pattern.svg)"/>',
            '<style>@import "https://attacker.example/theme.css";</style>',
            '<rect style="width:expression(alert(1))"/>',
            '<rect style="fill:javascript:alert(1)"/>',
        ],
        ids=[
            'absolute-url',
            'relative-url',
            'import',
            'expression',
            'javascript',
        ],
    )
    def test_default_svg_output_rejects_unsafe_css_references(
        self,
        tmp_path,
        unsafe_css,
    ):
        source = tmp_path / 'unsafe-css.svg'
        source.write_text(
            '<svg xmlns="http://www.w3.org/2000/svg" width="20" height="10">'
            '{}'
            '</svg>'.format(unsafe_css)
        )

        with pytest.raises(MediaDownloadError, match='external|executable'):
            postprocess_media_file(str(source), make_image_args())

    def test_default_svg_output_rejects_css_escape_obfuscation(self, tmp_path):
        source = tmp_path / 'escaped-url.svg'
        source.write_text(
            '<svg xmlns="http://www.w3.org/2000/svg" width="20" height="10">'
            '<rect width="20" height="10" '
            'style="fill:u\\72l(https://attacker.example/p.svg#x)"/>'
            '</svg>'
        )

        with pytest.raises(MediaDownloadError, match='external|executable'):
            postprocess_media_file(str(source), make_image_args())

    def test_default_raster_output_rejects_truncated_png(self, tmp_path):
        source = tmp_path / 'truncated.png'
        source.write_bytes(b'\x89PNG\r\n\x1a\n\x00\x00\x00\x0d')

        with pytest.raises(MediaDownloadError, match='truncated|corrupt|unsupported'):
            postprocess_media_file(str(source), make_image_args())

    def test_default_raster_output_rejects_truncated_jpeg(self, tmp_path):
        source = tmp_path / 'truncated.jpg'
        Image, _, _ = image_module.load_pillow_modules()
        Image.new('RGB', (64, 64), 'red').save(source, format='JPEG')
        source.write_bytes(source.read_bytes()[:-2])

        with pytest.raises(MediaDownloadError, match='truncated|corrupt|unsupported'):
            postprocess_media_file(str(source), make_image_args())

    def test_svg_rejects_utf16_before_entity_parsing(self, tmp_path):
        source = tmp_path / 'utf16-entity.svg'
        source.write_bytes((
            '<?xml version="1.0" encoding="UTF-16"?>'
            '<!DOCTYPE svg [<!ENTITY local "expanded">]>'
            '<svg xmlns="http://www.w3.org/2000/svg" width="20" height="10">'
            '<text>&local;</text>'
            '</svg>'
        ).encode('utf-16'))

        with pytest.raises(MediaDownloadError, match='UTF-16|encoding'):
            postprocess_media_file(str(source), make_image_args())

    def test_default_raster_output_rejects_oversized_dimensions(self, tmp_path):
        source = tmp_path / 'oversized.png'
        Image, _, _ = image_module.load_pillow_modules()
        Image.new('RGBA', (1, 1), (0, 0, 0, 255)).save(source)
        payload = bytearray(source.read_bytes())
        payload[16:24] = struct.pack('>II', 100_000, 100_000)
        payload[29:33] = struct.pack(
            '>I',
            zlib.crc32(bytes(payload[12:29])) & 0xffffffff,
        )
        source.write_bytes(payload)

        with pytest.raises(MediaDownloadError) as error:
            postprocess_media_file(str(source), make_image_args())
        assert 'pixel' in str(error.value).lower() or 'bomb' in str(error.value).lower()

    def test_default_svg_output_strips_external_doctype(self, tmp_path):
        source = tmp_path / 'doctype.svg'
        source.write_text(
            '<?xml version="1.0"?>\n'
            '<!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.0//EN" '
            '"http://www.w3.org/TR/2001/REC-SVG-20010904/DTD/svg10.dtd">\n'
            '<svg xmlns="http://www.w3.org/2000/svg" width="20" height="10"/>'
        )
        source.chmod(0o664)

        result = postprocess_media_file(str(source), make_image_args())

        assert result == str(source)
        assert '<!DOCTYPE' not in source.read_text()
        assert stat.S_IMODE(source.stat().st_mode) == 0o664

    def test_svg_dimensions_are_read_from_sanitized_xml_root(
        self,
        monkeypatch,
        tmp_path,
    ):
        source = tmp_path / 'sanitized-dimensions.svg'
        source.write_text(
            '<?xml version="1.0"?>\n'
            '<!DOCTYPE svg SYSTEM "https://example.org/external.dtd">\n'
            '<svg xmlns="http://www.w3.org/2000/svg" width="20" height="10"/>'
        )
        monkeypatch.setattr(
            image_module.ElementTree,
            'parse',
            lambda *args, **kwargs: pytest.fail(
                'validated SVG dimensions must not reparse the original file'
            ),
        )

        result = postprocess_media_file(str(source), make_image_args())

        assert result == str(source)
        assert '<!DOCTYPE' not in source.read_text()

    def test_rasterizing_local_svg_does_not_modify_input(self, tmp_path):
        pytest.importorskip('cairosvg')
        source = tmp_path / 'local-input.svg'
        source.write_text(
            '<?xml version="1.0"?>\n'
            '<!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.0//EN" '
            '"http://www.w3.org/TR/2001/REC-SVG-20010904/DTD/svg10.dtd">\n'
            '<svg xmlns="http://www.w3.org/2000/svg" width="20" height="10">'
            '<rect width="20" height="10" fill="#000"/>'
            '</svg>'
        )
        before = source.read_bytes()

        rasterized = image_module.rasterize_svg_to_image(str(source))

        assert rasterized.size == (20, 10)
        assert source.read_bytes() == before

    def test_svg_is_scaled_before_rasterization(self, tmp_path):
        pytest.importorskip('cairosvg')
        pytest.importorskip('PIL.Image')
        source = tmp_path / 'huge.svg'
        source.write_text(
            '<svg xmlns="http://www.w3.org/2000/svg" width="100000" height="50000">'
            '<rect width="100000" height="50000" fill="#000"/>'
            '</svg>'
        )
        rasterized = image_module.rasterize_svg_to_image(str(source), max_edge=100)
        assert rasterized.size == (100, 50)

    def test_ncbi_extraction_rejects_oversized_member(self, monkeypatch, tmp_path):
        archive_path = tmp_path / 'new_taxdump.tar.gz'
        destination_path = tmp_path / 'images.dmp'
        payload = b'12345'
        with tarfile.open(archive_path, 'w:gz') as archive:
            info = tarfile.TarInfo('images.dmp')
            info.size = len(payload)
            archive.addfile(info, io.BytesIO(payload))
        monkeypatch.setattr(image_module, 'NCBI_IMAGES_TABLE_MAX_BYTES', 4)

        with pytest.raises(MediaDownloadError, match='outside the supported range'):
            image_module._extract_ncbi_images_table(
                str(archive_path),
                str(destination_path),
            )
        assert not destination_path.exists()

    def test_ncbi_cache_replaces_corrupt_images_table(self, monkeypatch, tmp_path):
        args = make_image_args(
            download_dir=str(tmp_path / 'downloads'),
            out_dir=str(tmp_path / 'out'),
        )
        cache_dir = tmp_path / 'downloads' / 'ncbi-taxonomy-images'
        cache_dir.mkdir(parents=True)
        images_path = cache_dir / 'images.dmp'
        images_path.write_text('corrupt cache')
        valid_row = (
            b'64365\t|\timage:Cyanophora paradoxa\t|\t'
            b'https://example.org/image.jpg\t|\tCC BY 4.0\t|\tAuthor\t|\t'
            b'Publisher\t|\t\t|\t2762\t|\n'
        )

        def fake_download(session, url, destination_path, **kwargs):
            with tarfile.open(destination_path, 'w:gz') as archive:
                info = tarfile.TarInfo('images.dmp')
                info.size = len(valid_row)
                archive.addfile(info, io.BytesIO(valid_row))

        monkeypatch.setattr(image_module, '_fetch_ncbi_archive_md5', lambda session: None)
        monkeypatch.setattr(image_module, '_download_to_path', fake_download)

        result = image_module.ensure_ncbi_images_table(args, DummySession())

        assert result == str(images_path)
        assert images_path.read_bytes() == valid_row

    def test_ncbi_database_validation_handles_uri_characters_and_empty_tables(
        self,
        tmp_path,
    ):
        cache_dir = tmp_path / 'cache?#%'
        cache_dir.mkdir()
        database_path = cache_dir / 'images.sqlite3'
        connection = sqlite3.connect(database_path)
        connection.execute(
            'CREATE TABLE images ('
            'taxid INTEGER, record_id TEXT, image_url TEXT)'
        )
        connection.commit()
        connection.close()

        assert not image_module._ncbi_images_database_is_valid(
            str(database_path)
        )

        connection = sqlite3.connect(database_path)
        connection.execute(
            'INSERT INTO images VALUES (?, ?, ?)',
            (1, 'record', 'https://example.org/image.jpg'),
        )
        connection.commit()
        connection.close()

        assert image_module._ncbi_images_database_is_valid(str(database_path))

    def test_ncbi_database_reloads_table_after_late_corrupt_row(
        self,
        monkeypatch,
        tmp_path,
    ):
        args = make_image_args(
            download_dir=str(tmp_path / 'downloads'),
            out_dir=str(tmp_path / 'out'),
        )
        cache_dir = tmp_path / 'downloads' / 'ncbi-taxonomy-images'
        images_path = cache_dir / 'images.dmp'
        valid_row = (
            '64365\t|\timage:Cyanophora paradoxa\t|\t'
            'https://example.org/image.jpg\t|\tCC BY 4.0\t|\tAuthor\t|\t'
            'Publisher\t|\t\t|\t2762\t|\n'
        )
        calls = {'count': 0}

        def fake_ensure_table(args, session):
            calls['count'] += 1
            cache_dir.mkdir(parents=True, exist_ok=True)
            payload = (
                valid_row + 'corrupt later row\n'
                if calls['count'] == 1
                else valid_row
            )
            images_path.write_text(payload)
            return str(images_path)

        monkeypatch.setattr(
            image_module,
            'ensure_ncbi_images_table',
            fake_ensure_table,
        )

        database_path = image_module.ensure_ncbi_images_database(
            args,
            DummySession(),
        )

        assert calls['count'] == 2
        assert image_module._ncbi_images_database_is_valid(database_path)

    def test_image_outputs_must_be_distinct(self, tmp_path):
        out_dir = tmp_path / 'out'
        with pytest.raises(ValueError, match='Output paths must be distinct'):
            image_main(make_image_args(
                out_dir=str(out_dir),
                manifest_out=str(out_dir / 'ATTRIBUTION.md'),
                source='phylopic',
            ))

import csv
import os
from argparse import Namespace

import pytest
import requests

from nwkit.image import (
    BioiconsProvider,
    EOLProvider,
    IDigBioProvider,
    NCBIProvider,
    OpenverseProvider,
    build_providers,
    candidate_score,
    collect_candidates_for_species,
    download_media,
    extract_species_mapping,
    get_provider_quality_bonus,
    get_style_priority,
    image_main,
    postprocess_media_file,
    license_allowed,
    normalize_license_code,
    parse_ncbi_images_dmp_line,
    parse_sources,
    resolve_image_cache_dir,
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
        'name_tsv': None,
        'manifest': None,
        'attribution': None,
        'fail_on_missing': False,
        'output_format': 'original',
        'max_edge': None,
        'canvas': 'none',
        'background': 'white',
        'trim': 'off',
        'trim_shape': 'bbox',
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


class TestLicenseHelpers:
    def test_normalize_license_code_from_url(self):
        assert normalize_license_code(raw_url='https://creativecommons.org/licenses/by/4.0/') == 'cc-by'
        assert normalize_license_code(raw_url='https://creativecommons.org/publicdomain/zero/1.0/') == 'public-domain'

    def test_normalize_license_code_from_code_and_attribution(self):
        assert normalize_license_code(raw_code='cc-by-nc-sa') == 'cc-by-nc-sa'
        assert normalize_license_code(raw_code='CC BY-SA 3.0') == 'cc-by-sa'
        assert normalize_license_code(raw_code='Public domain') == 'public-domain'
        assert normalize_license_code(raw_code='MIT') == 'mit'
        assert normalize_license_code(raw_code='BSD-3-Clause') == 'bsd'
        assert normalize_license_code(raw_code='by-sa') == 'cc-by-sa'
        assert normalize_license_code(raw_code='pdm') == 'public-domain'
        assert normalize_license_code(raw_code=None, attribution='(c) Someone, all rights reserved') == 'all-rights-reserved'

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

    def test_style_and_provider_quality_helpers(self):
        candidate = {'asset_type': 'silhouette', 'provider_quality': 12}
        assert get_style_priority(candidate, style='silhouette') == 2
        assert get_style_priority(candidate, style='photo') == 0
        assert get_provider_quality_bonus(candidate) == 12


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

        provider = NCBIProvider(session=DummySession(), ncbi=FakeNCBI(), args=make_image_args(out_dir=str(tmp_path / 'out')))
        candidates = provider.fetch_candidates('Cyanophora paradoxa', fallback_rank='none')

        assert len(candidates) == 1
        candidate = candidates[0]
        assert candidate['provider'] == 'ncbi'
        assert candidate['provider_record_id'] == '64365'
        assert candidate['matched_name'] == 'Cyanophora paradoxa'
        assert candidate['matched_rank'] == 'species'
        assert candidate['license_code'] == 'cc-by-sa'
        assert candidate['license_url'] == 'https://creativecommons.org/licenses/by-sa/3.0/'
        assert candidate['attribution'] == 'Wolfgang Bettighofer, Wikimedia Commons'
        assert candidate['media_url'] == 'http://www.ncbi.nlm.nih.gov/Taxonomy/taxi/images/4'

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
                    'media_url': 'https://example.org/photo.jpg',
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
                        'media_url': 'https://example.org/photo-1.jpg',
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
                        'media_url': 'https://example.org/photo-2.jpg',
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
                        'media_url': 'https://example.org/photo-3.jpg',
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
                    'media_url': 'https://example.org/photo-2.jpg',
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


class TestBioiconsProvider:
    def test_bioicons_provider_fetches_matching_svg_candidates(self, tmp_path):
        class RoutingSession:
            def get(self, url, params=None, timeout=None, headers=None):
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


class TestEOLProvider:
    def test_eol_provider_fetches_candidates_from_page_media(self):
        class RoutingSession:
            def get(self, url, params=None, timeout=None, headers=None):
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
            def post(self, url, json=None, timeout=None, headers=None):
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
            def get(self, url, params=None, timeout=None, headers=None):
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
            'details': "Expected the 'GENUS_SPECIES[_...]' convention or a matching --name_tsv entry.",
        }]

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
                'media_url': 'https://images.example.org/homo.svg',
                'width': 1200,
                'height': 800,
                'asset_type': 'silhouette',
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

        def fake_download_media(session, media_url, destination_path, cache_path=None):
            with open(destination_path, 'wb') as handle:
                handle.write(b'content')
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

        assert unmatched_rows == [{
            'leaf_name': 'BadLabel',
            'species_name': '',
            'reason': 'unparsable leaf label',
            'details': "Expected the 'GENUS_SPECIES[_...]' convention or a matching --name_tsv entry.",
        }]
        assert 'Homo sapiens' in attribution_text
        assert 'Panthera leo' in attribution_text

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
                'media_url': 'https://images.example.org/apis.svg',
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

        def fake_download_media(session, media_url, destination_path, cache_path=None):
            with open(destination_path, 'wb') as handle:
                handle.write(b'content')
            return 'downloaded'

        monkeypatch.setattr('nwkit.image.build_providers', fake_build_providers)
        monkeypatch.setattr('nwkit.image.download_media', fake_download_media)

        args = make_image_args(
            infile=str(tree_path),
            out_dir=str(out_dir),
            name_tsv=str(name_tsv),
            source='phylopic',
            fail_on_missing=True,
        )
        with pytest.raises(SystemExit, match='1'):
            image_main(args)

        manifest_rows = read_tsv(out_dir / 'manifest.tsv')
        unmatched_rows = read_tsv(out_dir / 'unmatched.tsv')

        assert manifest_rows[0]['leaf_name'] == 'Sample_1'
        assert manifest_rows[0]['species_name'] == 'Apis mellifera'
        assert unmatched_rows == [{
            'leaf_name': 'Unknown',
            'species_name': '',
            'reason': 'unparsable leaf label',
            'details': "Expected the 'GENUS_SPECIES[_...]' convention or a matching --name_tsv entry.",
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
                'media_url': 'https://images.example.org/apis-cache.svg',
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
                yield b'image-bytes'

        class CountingSession:
            def get(self, media_url, stream=True, timeout=None, headers=None):
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
                'media_url': 'https://images.example.org/apis-once.svg',
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

        def fake_download_media(session, media_url, destination_path, cache_path=None):
            with open(destination_path, 'wb') as handle:
                handle.write(b'content')
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

        def fake_download_media(session, media_url, destination_path, cache_path=None):
            call_counter['count'] += 1
            with open(destination_path, 'wb') as handle:
                handle.write(b'content')
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

        def fake_download_media(session, media_url, destination_path, cache_path=None):
            resolved_path = destination_path[:-4] + '.jpg'
            with open(resolved_path, 'wb') as handle:
                handle.write(b'jpeg-data')
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

        def fake_download_media(session, media_url, destination_path, cache_path=None):
            if media_url.endswith('/fail/original.jpg'):
                raise requests.ConnectionError('transient failure')
            with open(destination_path, 'wb') as handle:
                handle.write(b'content')
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

        def fake_download_media(session, media_url, destination_path, cache_path=None):
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
    def test_download_media_prefers_shared_cache_when_present(self, tmp_path):
        destination = tmp_path / 'out' / 'image.svg'
        cache_path = tmp_path / 'cache' / 'image.svg'
        cache_path.parent.mkdir(parents=True)
        cache_path.write_bytes(b'cached-bytes')

        class FailingSession:
            def get(self, *args, **kwargs):
                raise AssertionError('network should not be used when cache exists')

        result = download_media(FailingSession(), 'https://images.example.org/item.svg', str(destination), cache_path=str(cache_path))

        assert result['status'] == 'cached'
        assert destination.read_bytes() == b'cached-bytes'

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

    def test_download_media_reuses_existing_cache_variant_extension(self, tmp_path):
        destination = tmp_path / 'out' / 'image.bin'
        cache_path = tmp_path / 'cache' / 'image.bin'
        cache_variant = tmp_path / 'cache' / 'image.jpg'
        cache_variant.parent.mkdir(parents=True)
        cache_variant.write_bytes(b'jpeg-data')

        class FailingSession:
            def get(self, *args, **kwargs):
                raise AssertionError('network should not be used when cache variant exists')

        result = download_media(FailingSession(), 'https://example.org/download', str(destination), cache_path=str(cache_path))

        assert result['status'] == 'cached'
        assert result['destination_path'].endswith('.jpg')
        assert (tmp_path / 'out' / 'image.jpg').read_bytes() == b'jpeg-data'


class TestImagePostprocessing:
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
            lambda source_path: Image.new('RGBA', (20, 10), (0, 0, 0, 255)),
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
                'media_url': 'https://images.example.org/apis.svg',
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

        def fake_download_media(session, media_url, destination_path, cache_path=None):
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
                'media_url': 'https://images.example.org/apis.jpg',
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

        def fake_download_media(session, media_url, destination_path, cache_path=None):
            with open(destination_path, 'wb') as handle:
                handle.write(b'jpeg-data')
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
                'media_url': 'https://images.example.org/apis.svg',
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

        def fake_download_media(session, media_url, destination_path, cache_path=None):
            with open(destination_path, 'wb') as handle:
                handle.write(b'<svg xmlns="http://www.w3.org/2000/svg"></svg>')
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
                'media_url': 'https://images.example.org/apis.jpg',
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

        def fake_download_media(session, media_url, destination_path, cache_path=None):
            with open(destination_path, 'wb') as handle:
                handle.write(b'jpeg-data')
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

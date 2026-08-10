import threading

import pytest

from nwkit import image as image_module
from nwkit.image import (
    MediaDownloadError,
    NCBIProvider,
    PhylopicProvider,
    allowed_candidates_from_scored_candidates,
    build_download_session,
    build_providers,
    build_retry_config,
    candidate_score,
    classify_wikimedia_asset,
    collect_candidates_for_species,
    collect_candidates_for_species_map,
    get_aspect_fit_bonus,
    get_provider_quality_bonus,
    get_style_priority,
    license_allowed,
    normalize_license_code,
    parse_ncbi_images_dmp_line,
    parse_sources,
    resolve_download_worker_count,
    resolve_lookup_worker_count,
    resolve_ncbi_taxonomy_image_cache_dir,
    resolve_provider_fetch_limit,
    wikimedia_page_mentions_query,
)
from tests.image_test_support import (
    DummyProvider,
    DummySession,
    make_image_args,
)


class TestLicenseHelpers:
    def test_normalize_license_code_from_url(self):
        assert (
            normalize_license_code(
                raw_url="https://creativecommons.org/licenses/by/4.0/"
            )
            == "cc-by"
        )
        assert (
            normalize_license_code(
                raw_url="https://creativecommons.org/publicdomain/zero/1.0/"
            )
            == "public-domain"
        )

    def test_normalize_license_code_from_code_and_attribution(self):
        assert normalize_license_code(raw_code="cc-by-nc-sa") == "cc-by-nc-sa"
        assert normalize_license_code(raw_code="CC BY-SA 3.0") == "cc-by-sa"
        assert normalize_license_code(raw_code="Public domain") == "public-domain"
        assert normalize_license_code(raw_code="CC0 1.0") == "public-domain"
        assert normalize_license_code(raw_code="MIT") == "mit"
        assert normalize_license_code(raw_code="BSD-3-Clause") == "bsd"
        assert normalize_license_code(raw_code="by-sa") == "cc-by-sa"
        assert normalize_license_code(raw_code="pdm") == "public-domain"
        assert (
            normalize_license_code(
                raw_code=None, attribution="(c) Someone, all rights reserved"
            )
            == "all-rights-reserved"
        )

    @pytest.mark.parametrize(
        "raw_code",
        [
            "limited use only",
            "no redistribution; limited educational use",
            "submitted by user",
            "not public domain",
        ],
    )
    def test_restrictive_or_unrelated_metadata_is_not_open_license(
        self,
        raw_code,
    ):
        license_code = normalize_license_code(raw_code=raw_code)

        assert license_code == "unknown"
        assert not license_allowed(license_code, license_max="any")

    @pytest.mark.parametrize(
        ("raw_code", "expected"),
        [
            ("MIT License", "mit"),
            ("BSD-2-Clause", "bsd"),
            ("Public Domain Mark 1.0", "public-domain"),
        ],
    )
    def test_explicit_open_license_labels_remain_allowed(
        self,
        raw_code,
        expected,
    ):
        license_code = normalize_license_code(raw_code=raw_code)

        assert license_code == expected
        assert license_allowed(license_code, license_max="any")

    @pytest.mark.parametrize(
        "raw_url",
        [
            "https://evil.example/licenses/by/4.0/",
            "https://example.com/?next=/publicdomain/zero/1.0",
            "https://notcreativecommons.org/licenses/by-nc/4.0/",
            "javascript://creativecommons.org/licenses/by/4.0/",
            "https://creativecommons.org:444/licenses/by/4.0/",
        ],
    )
    def test_license_url_requires_an_official_host(self, raw_url):
        assert normalize_license_code(raw_url=raw_url) == "unknown"

    def test_license_url_accepts_official_subdomains_and_ignores_query(self):
        assert (
            normalize_license_code(
                raw_url=(
                    "https://licenses.creativecommons.org/licenses/by/4.0/"
                    "?next=/licenses/by-nc/4.0/"
                )
            )
            == "cc-by"
        )

    def test_license_allowed_respects_nd_and_ceiling(self):
        assert license_allowed("cc-by", license_max="cc-by", allow_nd=False) is True
        assert (
            license_allowed("cc-by-nc", license_max="cc-by-sa", allow_nd=False) is False
        )
        assert license_allowed("cc-by-nd", license_max="cc-by", allow_nd=False) is False
        assert license_allowed("cc-by-nd", license_max="cc-by", allow_nd=True) is True
        assert license_allowed("mit", license_max="cc-by", allow_nd=False) is True
        assert license_allowed("bsd", license_max="cc-by-sa", allow_nd=False) is True
        assert (
            license_allowed("mit", license_max="public-domain", allow_nd=False) is False
        )
        assert (
            license_allowed("all-rights-reserved", license_max="any", allow_nd=True)
            is False
        )


class TestFetchLimits:
    def test_resolve_provider_fetch_limit_scales_with_max_per_species(self):
        assert resolve_provider_fetch_limit(make_image_args(max_per_species=1)) == 10
        assert resolve_provider_fetch_limit(make_image_args(max_per_species=8)) == 12
        assert resolve_provider_fetch_limit(make_image_args(max_per_species=99)) == 30

    def test_candidate_score_prioritizes_source_order_before_license_and_quality(self):
        first_source = {
            "matched_rank": "species",
            "license_code": "cc-by",
            "asset_type": "photo",
            "width": 1200,
            "height": 900,
            "provider_quality": 0,
        }
        later_source = {
            "matched_rank": "species",
            "license_code": "public-domain",
            "asset_type": "photo",
            "width": 8000,
            "height": 8000,
            "provider_quality": 99,
        }

        assert candidate_score(
            first_source, provider_index=0, style="photo"
        ) > candidate_score(later_source, provider_index=1, style="photo")

    def test_candidate_score_uses_provider_quality_only_after_quality(self):
        lower_quality = {
            "matched_rank": "species",
            "license_code": "cc-by",
            "asset_type": "photo",
            "width": 1000,
            "height": 1000,
            "provider_quality": 90,
        }
        higher_quality = {
            "matched_rank": "species",
            "license_code": "cc-by",
            "asset_type": "photo",
            "width": 2000,
            "height": 2000,
            "provider_quality": 0,
        }

        assert candidate_score(
            higher_quality, provider_index=0, style="photo"
        ) > candidate_score(lower_quality, provider_index=0, style="photo")

    def test_candidate_score_prefers_primary_then_vector_and_drawable_aspect(self):
        base = {
            "matched_rank": "species",
            "license_code": "cc-by",
            "asset_type": "silhouette",
            "width": 1000,
            "height": 1000,
            "provider_quality": 0,
            "media_url": "https://example.org/image.png",
        }
        primary = dict(base, is_primary=True)
        vector = dict(base, is_vector=True, media_url="https://example.org/image.svg")
        raster = dict(base, is_vector=False)
        extreme_vector = dict(vector, width=10000, height=100)
        balanced_vector = dict(vector, width=1000, height=800)

        assert candidate_score(
            primary, provider_index=0, style="silhouette"
        ) > candidate_score(vector, provider_index=0, style="silhouette")
        assert candidate_score(
            vector, provider_index=0, style="silhouette"
        ) > candidate_score(raster, provider_index=0, style="silhouette")
        assert candidate_score(
            balanced_vector, provider_index=0, style="silhouette"
        ) > candidate_score(extreme_vector, provider_index=0, style="silhouette")
        assert get_aspect_fit_bonus(balanced_vector) == 100
        assert get_aspect_fit_bonus(extreme_vector) == 0

    def test_candidate_score_keeps_exact_match_ahead_of_fallback_primary(self):
        exact = {
            "matched_rank": "species",
            "license_code": "cc-by",
            "asset_type": "silhouette",
            "width": 500,
            "height": 500,
            "is_primary": False,
        }
        genus_primary = dict(exact, matched_rank="genus", is_primary=True)

        assert candidate_score(
            exact, provider_index=0, style="silhouette"
        ) > candidate_score(genus_primary, provider_index=0, style="silhouette")

    def test_disallowed_primary_falls_back_to_allowed_ranked_candidate(self):
        primary = {
            "matched_rank": "species",
            "license_code": "cc-by-nc-sa",
            "asset_type": "silhouette",
            "width": 1000,
            "height": 1000,
            "is_primary": True,
            "media_url": "https://example.org/primary.svg",
        }
        fallback = dict(
            primary,
            license_code="public-domain",
            is_primary=False,
            media_url="https://example.org/fallback.svg",
        )
        for candidate in (primary, fallback):
            candidate["score"] = candidate_score(
                candidate,
                provider_index=0,
                style="silhouette",
            )

        allowed = allowed_candidates_from_scored_candidates(
            [primary, fallback],
            args=make_image_args(license_max="public-domain"),
        )

        assert allowed == [fallback]

    def test_style_and_provider_quality_helpers(self):
        candidate = {"asset_type": "silhouette", "provider_quality": 12}
        assert get_style_priority(candidate, style="silhouette") == 2
        assert get_style_priority(candidate, style="photo") == 0
        assert get_provider_quality_bonus(candidate) == 12

    def test_candidate_score_normalizes_numeric_strings_and_nonfinite_values(
        self,
    ):
        numeric_strings = {
            "matched_rank": "species",
            "license_code": "cc-by",
            "width": "1024",
            "height": "512",
            "provider_quality": "1e3",
        }
        malformed_numbers = dict(
            numeric_strings,
            width=float("nan"),
            height=float("inf"),
            provider_quality="not-a-number",
        )

        assert isinstance(candidate_score(numeric_strings), int)
        assert isinstance(candidate_score(malformed_numbers), int)
        assert get_aspect_fit_bonus(malformed_numbers) == 0

    def test_candidate_normalization_parses_only_explicit_boolean_values(self):
        candidate = {
            "provider": "phylopic",
            "provider_record_id": "candidate-1",
            "matched_name": "Species alpha",
            "matched_rank": "species",
            "license_code": "cc-by",
            "media_url": "https://images.phylopic.org/candidate.svg",
            "is_primary": "false",
            "is_vector": "false",
        }

        normalized = image_module.normalize_provider_candidate(
            candidate,
            expected_provider="phylopic",
        )

        assert normalized["is_primary"] is False
        assert normalized["is_vector"] is False
        assert image_module.candidate_is_vector(normalized) is False

        with pytest.raises(ValueError, match="invalid 'is_vector' value"):
            image_module.normalize_provider_candidate(
                dict(candidate, is_vector="yes"),
                expected_provider="phylopic",
            )

    def test_candidate_collection_skips_malformed_records_individually(
        self,
        tmp_path,
    ):
        base = {
            "provider": "phylopic",
            "provider_record_id": "valid",
            "matched_name": "Species alpha",
            "matched_rank": "species",
            "license_code": "cc-by",
            "license_url": "https://creativecommons.org/licenses/by/4.0/",
            "attribution": "Artist",
            "source_page_url": "https://phylopic.org/images/valid",
            "media_url": "https://images.phylopic.org/valid.svg",
            "asset_type": "silhouette",
            "width": "1024",
            "height": "512",
        }
        candidates = [
            dict(base, provider_record_id="bad-rank", matched_rank={"bad": 1}),
            dict(base, provider_record_id="bad-license", license_code=["cc-by"]),
            base,
        ]
        args = make_image_args(
            out_dir=str(tmp_path / "out"),
            source="phylopic",
            style="silhouette",
        )

        collected, errors = collect_candidates_for_species(
            "Species alpha",
            args=args,
            sources=["phylopic"],
            providers={
                "phylopic": DummyProvider({"Species alpha": candidates}),
            },
        )

        assert [candidate["provider_record_id"] for candidate in collected] == ["valid"]
        assert collected[0]["width"] == 1024.0
        assert len(errors) == 2


class TestPhylopicProvider:
    @staticmethod
    def _image_item(uuid, license_url, media_url, sizes="1000x800", vector=True):
        selected_file = {
            "href": media_url,
            "sizes": sizes,
            "type": "image/svg+xml" if vector else "image/png",
        }
        links = {
            "license": {"href": license_url},
            "self": {"href": "/images/{}".format(uuid)},
            "sourceFile": selected_file,
        }
        if vector:
            links["vectorFile"] = selected_file
        return {
            "_links": links,
            "attribution": "Example Artist",
            "uuid": uuid,
        }

    def test_fetch_candidates_marks_and_fetches_linked_primary_image(self, monkeypatch):
        primary_uuid = "primary-uuid"
        fallback_item = self._image_item(
            uuid="fallback-uuid",
            license_url="https://creativecommons.org/publicdomain/zero/1.0/",
            media_url="https://images.example.org/fallback.svg",
        )
        primary_item = self._image_item(
            uuid=primary_uuid,
            license_url="https://creativecommons.org/licenses/by/4.0/",
            media_url="https://images.example.org/primary.png",
            vector=False,
        )

        class FakeNCBI:
            def get_name_translator(self, names):
                return {"Apis mellifera": [7460]}

        provider = PhylopicProvider(session=DummySession(), ncbi=FakeNCBI())
        monkeypatch.setattr(
            provider,
            "_resolve_node",
            lambda taxid: {
                "_links": {
                    "primaryImage": {
                        "href": "/images/{}?build=545".format(primary_uuid),
                    },
                },
                "build": 545,
                "uuid": "node-uuid",
            },
        )
        monkeypatch.setattr(
            provider, "_fetch_node_images", lambda node_uuid, build: [fallback_item]
        )
        linked_hrefs = list()

        def fake_fetch_linked_image(href):
            linked_hrefs.append(href)
            return primary_item

        monkeypatch.setattr(provider, "_fetch_linked_image", fake_fetch_linked_image)

        candidates = provider.fetch_candidates("Apis mellifera")
        candidates_by_uuid = {
            candidate["provider_record_id"]: candidate for candidate in candidates
        }

        assert linked_hrefs == ["/images/{}?build=545".format(primary_uuid)]
        assert candidates_by_uuid[primary_uuid]["is_primary"] is True
        assert candidates_by_uuid[primary_uuid]["is_vector"] is False
        assert candidates_by_uuid["fallback-uuid"]["is_primary"] is False
        assert candidates_by_uuid["fallback-uuid"]["is_vector"] is True
        assert candidate_score(
            candidates_by_uuid[primary_uuid], provider_index=0, style="silhouette"
        ) > candidate_score(
            candidates_by_uuid["fallback-uuid"], provider_index=0, style="silhouette"
        )


class TestSourceParsing:
    def test_parse_sources_uses_style_defaults(self):
        assert parse_sources("auto", None) == [
            "phylopic",
            "bioicons",
            "inaturalist",
            "wikimedia",
            "gbif",
            "eol",
            "idigbio",
            "openverse",
            "ncbi",
        ]
        assert parse_sources("photo", None) == [
            "inaturalist",
            "wikimedia",
            "gbif",
            "eol",
            "idigbio",
            "openverse",
            "ncbi",
        ]
        assert parse_sources("silhouette", None) == [
            "phylopic",
            "bioicons",
            "wikimedia",
        ]

    def test_parse_sources_accepts_all_implemented_sources(self):
        assert parse_sources(
            "auto",
            "phylopic,bioicons,wikimedia,gbif,inaturalist,eol,idigbio,openverse,ncbi",
        ) == [
            "phylopic",
            "bioicons",
            "wikimedia",
            "gbif",
            "inaturalist",
            "eol",
            "idigbio",
            "openverse",
            "ncbi",
        ]

    def test_parse_sources_rejects_unimplemented_sources(self):
        with pytest.raises(ValueError, match="Unsupported --source"):
            parse_sources("auto", "phylopic,example")


class TestWikimediaHelpers:
    def test_wikimedia_page_mentions_query_filters_irrelevant_results(self):
        relevant_page = {
            "title": "File:Lion (Panthera leo) old male Chobe.jpg",
            "imageinfo": [
                {
                    "extmetadata": {
                        "ObjectName": {"value": "Lion (Panthera leo) old male Chobe"},
                        "ImageDescription": {
                            "value": "Lion (<i>Panthera leo</i>), male, Chobe National Park, Botswana"
                        },
                    }
                }
            ],
        }
        irrelevant_page = {
            "title": "File:La Bohémienne endormie.jpg",
            "imageinfo": [
                {
                    "extmetadata": {
                        "ObjectName": {"value": "La Bohémienne endormie"},
                        "ImageDescription": {"value": ""},
                    }
                }
            ],
        }
        assert wikimedia_page_mentions_query(relevant_page, "Panthera leo") is True
        assert wikimedia_page_mentions_query(irrelevant_page, "Panthera leo") is False

    def test_classify_wikimedia_asset_rejects_research_figures_and_gifs(self):
        page = {
            "title": "File:Arabidopsis thaliana expression figure.gif",
            "imageinfo": [
                {
                    "url": "https://upload.wikimedia.org/expression.gif",
                    "mime": "image/gif",
                    "extmetadata": {
                        "ImageDescription": {
                            "value": "Gene-expression graph for Arabidopsis thaliana"
                        },
                    },
                }
            ],
        }

        assert classify_wikimedia_asset(page) == "illustration"

    def test_classify_wikimedia_asset_preserves_silhouettes(self):
        page = {
            "title": "File:Panthera leo silhouette.svg",
            "imageinfo": [
                {
                    "url": "https://upload.wikimedia.org/lion.svg",
                    "mime": "image/svg+xml",
                }
            ],
        }

        assert classify_wikimedia_asset(page) == "silhouette"


class TestNCBIHelpers:
    def test_parse_ncbi_images_dmp_line(self):
        line = "64373\t|\timage:Alternaria brassicae\t|\thttp://www.ncbi.nlm.nih.gov/Taxonomy/taxi/images/12\t|\tCC BY-NC (https://creativecommons.org/licenses/by-nc/4.0/)\t|\tgen_ok\t|\tiNaturalist\t|\t\t|\t29911\t|\n"
        record = parse_ncbi_images_dmp_line(line)

        assert record["record_id"] == "64373"
        assert record["title"] == "Alternaria brassicae"
        assert (
            record["image_url"] == "http://www.ncbi.nlm.nih.gov/Taxonomy/taxi/images/12"
        )
        assert record["license_code_text"] == "CC BY-NC"
        assert (
            record["license_url"] == "https://creativecommons.org/licenses/by-nc/4.0/"
        )
        assert record["attribution"] == "gen_ok"
        assert record["source_name"] == "iNaturalist"
        assert record["taxids"] == [29911]

    def test_parse_ncbi_images_dmp_line_supports_multiple_taxids(self):
        line = "64698\t|\timage:Abudefduf saxatilis\t|\thttp://www.ncbi.nlm.nih.gov/Taxonomy/taxi/images/31217\t|\tCC BY-SA 3.0 (https://creativecommons.org/licenses/by-sa/3.0/)\t|\tCralize\t|\tWikimedia Commons\t|\t\t|\t50731 1567454\t|\n"
        record = parse_ncbi_images_dmp_line(line)
        assert record["taxids"] == [50731, 1567454]

    def test_build_download_session_configures_retries_for_rate_limits(self):
        session = build_download_session()
        adapter = session.get_adapter("https://example.org/")
        retries = adapter.max_retries
        assert retries.total == 4
        assert 429 in retries.status_forcelist

    def test_build_retry_config_limits_retries_to_get_requests(self):
        retries = build_retry_config()
        allowed_methods = getattr(retries, "allowed_methods", None)
        if allowed_methods is None:
            allowed_methods = getattr(retries, "method_whitelist", None)
        assert allowed_methods == frozenset(["GET"])

    def test_resolve_lookup_worker_count_defaults_to_four_without_taxonomy(self):
        workers = resolve_lookup_worker_count(
            args=make_image_args(),
            sources=["wikimedia", "openverse"],
            species_count=20,
        )
        assert workers == 4

    def test_resolve_lookup_worker_count_defaults_to_two_with_taxonomy(self):
        workers = resolve_lookup_worker_count(
            args=make_image_args(),
            sources=["phylopic", "wikimedia"],
            species_count=20,
        )
        assert workers == 2

    def test_resolve_download_worker_count_defaults_to_four(self):
        assert resolve_download_worker_count(species_count=20) == 4

    def test_parallel_lookup_closes_taxonomy_handles_in_their_worker_threads(
        self, monkeypatch, tmp_path
    ):

        barrier = threading.Barrier(2)
        close_threads = list()

        class ThreadBoundNCBI:
            def __init__(self):
                self.created_thread = threading.get_ident()

            def close(self):
                close_threads.append((self.created_thread, threading.get_ident()))

        class BarrierProvider:
            def fetch_candidates(self, species_name, fallback_rank="none"):
                barrier.wait(timeout=5)
                return []

        def fake_build_providers(args, sources, session=None):
            return DummySession(), ThreadBoundNCBI(), {"phylopic": BarrierProvider()}

        monkeypatch.setattr("nwkit.image.build_providers", fake_build_providers)
        args = make_image_args(
            out_dir=str(tmp_path / "out"),
            download_dir=str(tmp_path / "cache"),
            source="phylopic",
        )

        results = collect_candidates_for_species_map(
            species_names=["Species alpha", "Species beta"],
            args=args,
            sources=["phylopic"],
        )

        assert set(results) == {"Species alpha", "Species beta"}
        assert len(close_threads) == 2
        assert len({created_thread for created_thread, _ in close_threads}) == 2
        assert all(
            created_thread == closed_thread
            for created_thread, closed_thread in close_threads
        )

    def test_resolve_ncbi_taxonomy_image_cache_dir(self, tmp_path):
        shared_args = make_image_args(
            download_dir=str(tmp_path / "shared"), out_dir=str(tmp_path / "out")
        )
        auto_args = make_image_args(download_dir="auto", out_dir=str(tmp_path / "out"))

        assert resolve_ncbi_taxonomy_image_cache_dir(shared_args) == str(
            tmp_path / "shared" / "ncbi-taxonomy-images"
        )
        assert resolve_ncbi_taxonomy_image_cache_dir(auto_args) == str(
            tmp_path / "out" / ".nwkit-cache" / "ncbi-taxonomy-images"
        )

    def test_ncbi_provider_fetches_candidates_from_images_table(
        self, monkeypatch, tmp_path
    ):
        images_path = tmp_path / "images.dmp"
        images_path.write_text(
            "64365\t|\timage:Cyanophora paradoxa\t|\thttp://www.ncbi.nlm.nih.gov/Taxonomy/taxi/images/4\t|\tCC BY-SA 3.0 (https://creativecommons.org/licenses/by-sa/3.0/)\t|\tWolfgang Bettighofer\t|\tWikimedia Commons\t|\t\t|\t2762\t|\n"
        )

        class FakeNCBI:
            def get_name_translator(self, names):
                mapping = {
                    "Cyanophora paradoxa": [2762],
                }
                return {name: mapping[name] for name in names if name in mapping}

        monkeypatch.setattr(
            "nwkit.image.ensure_ncbi_images_table",
            lambda args, session: str(images_path),
        )

        args = make_image_args(
            out_dir=str(tmp_path / "out"),
            download_dir=str(tmp_path / "cache"),
        )
        provider = NCBIProvider(session=DummySession(), ncbi=FakeNCBI(), args=args)
        candidates, errors = collect_candidates_for_species(
            "Cyanophora paradoxa",
            args=args,
            sources=["ncbi"],
            providers={"ncbi": provider},
        )

        assert errors == []
        assert len(candidates) == 1
        candidate = candidates[0]
        assert candidate["provider"] == "ncbi"
        assert candidate["provider_record_id"] == "64365"
        assert candidate["matched_name"] == "Cyanophora paradoxa"
        assert candidate["matched_rank"] == "species"
        assert candidate["license_code"] == "cc-by-sa"
        assert (
            candidate["license_url"]
            == "https://creativecommons.org/licenses/by-sa/3.0/"
        )
        assert candidate["attribution"] == "Wolfgang Bettighofer, Wikimedia Commons"
        assert (
            candidate["media_url"]
            == "https://www.ncbi.nlm.nih.gov/Taxonomy/taxi/images/4"
        )

    def test_ncbi_provider_does_not_upgrade_arbitrary_http_media(self):
        provider = NCBIProvider(
            session=DummySession(),
            ncbi=None,
            args=make_image_args(out_dir="unused"),
        )
        candidate = provider._candidate_from_record(
            record={
                "record_id": "1",
                "image_url": "http://example.org/Taxonomy/taxi/images/4",
                "attribution": "",
                "source_name": "",
                "license_code_text": "",
                "license_url": "",
            },
            matched_name="Example species",
            matched_rank="species",
        )

        with pytest.raises(MediaDownloadError, match="HTTPS"):
            image_module.validate_candidate_media_url(candidate)

    def test_build_providers_does_not_initialize_ncbi_eagerly(
        self, monkeypatch, tmp_path
    ):
        call_counter = {"count": 0}

        def fake_get_ete_ncbitaxa(args=None):
            call_counter["count"] += 1
            raise AssertionError(
                "NCBI taxonomy should not initialize during provider construction"
            )

        monkeypatch.setattr("nwkit.image.get_ete_ncbitaxa", fake_get_ete_ncbitaxa)

        session, ncbi, providers = build_providers(
            args=make_image_args(out_dir=str(tmp_path / "out")),
            sources=["ncbi"],
            session=DummySession(),
        )

        assert session is not None
        assert ncbi is not None
        assert "ncbi" in providers
        assert call_counter["count"] == 0

    def test_collect_candidates_skips_ncbi_when_earlier_provider_has_allowed_candidate(
        self,
    ):
        call_counter = {"ncbi": 0}

        class CountingProvider(DummyProvider):
            def fetch_candidates(self, species_name, fallback_rank="none"):
                call_counter["ncbi"] += 1
                return super().fetch_candidates(
                    species_name, fallback_rank=fallback_rank
                )

        providers = {
            "inaturalist": DummyProvider(
                {
                    "Apis mellifera": [
                        {
                            "provider": "inaturalist",
                            "provider_record_id": "inat-1",
                            "matched_name": "Apis mellifera",
                            "matched_rank": "species",
                            "license_code": "cc-by",
                            "license_url": "https://creativecommons.org/licenses/by/4.0/",
                            "attribution": "A",
                            "source_page_url": "https://example.org/obs",
                            "media_url": "https://static.inaturalist.org/photo.jpg",
                            "width": 1000,
                            "height": 900,
                            "asset_type": "photo",
                        }
                    ],
                }
            ),
            "ncbi": CountingProvider(
                {
                    "Apis mellifera": [
                        {
                            "provider": "ncbi",
                            "provider_record_id": "ncbi-1",
                            "matched_name": "Apis mellifera",
                            "matched_rank": "species",
                            "license_code": "public-domain",
                            "license_url": "https://creativecommons.org/publicdomain/zero/1.0/",
                            "attribution": "NCBI",
                            "source_page_url": "https://www.ncbi.nlm.nih.gov/Taxonomy/taxi/images/1",
                            "media_url": "https://www.ncbi.nlm.nih.gov/Taxonomy/taxi/images/1",
                            "width": None,
                            "height": None,
                            "asset_type": "photo",
                        }
                    ],
                }
            ),
        }

        candidates, provider_errors = collect_candidates_for_species(
            species_name="Apis mellifera",
            args=make_image_args(style="photo", source="inaturalist,ncbi"),
            sources=["inaturalist", "ncbi"],
            providers=providers,
        )

        assert provider_errors == []
        assert len(candidates) == 1
        assert candidates[0]["provider"] == "inaturalist"
        assert call_counter["ncbi"] == 0

    def test_collect_candidates_stops_after_exact_species_allowed_candidate(
        self, tmp_path
    ):
        call_counter = {"openverse": 0}

        class CountingProvider(DummyProvider):
            def fetch_candidates(self, species_name, fallback_rank="none"):
                call_counter["openverse"] += 1
                return super().fetch_candidates(
                    species_name, fallback_rank=fallback_rank
                )

        providers = {
            "inaturalist": DummyProvider(
                {
                    "Apis mellifera": [
                        {
                            "provider": "inaturalist",
                            "provider_record_id": "inat-early-1",
                            "matched_name": "Apis mellifera",
                            "matched_rank": "species",
                            "license_code": "cc-by",
                            "license_url": "https://creativecommons.org/licenses/by/4.0/",
                            "attribution": "A1",
                            "source_page_url": "https://example.org/obs/1",
                            "media_url": "https://static.inaturalist.org/photo-1.jpg",
                            "width": 1000,
                            "height": 900,
                            "asset_type": "photo",
                        },
                        {
                            "provider": "inaturalist",
                            "provider_record_id": "inat-early-2",
                            "matched_name": "Apis mellifera",
                            "matched_rank": "species",
                            "license_code": "cc-by",
                            "license_url": "https://creativecommons.org/licenses/by/4.0/",
                            "attribution": "A2",
                            "source_page_url": "https://example.org/obs/2",
                            "media_url": "https://static.inaturalist.org/photo-2.jpg",
                            "width": 1200,
                            "height": 950,
                            "asset_type": "photo",
                        },
                        {
                            "provider": "inaturalist",
                            "provider_record_id": "inat-early-3",
                            "matched_name": "Apis mellifera",
                            "matched_rank": "species",
                            "license_code": "cc-by",
                            "license_url": "https://creativecommons.org/licenses/by/4.0/",
                            "attribution": "A3",
                            "source_page_url": "https://example.org/obs/3",
                            "media_url": "https://static.inaturalist.org/photo-3.jpg",
                            "width": 1100,
                            "height": 920,
                            "asset_type": "photo",
                        },
                    ],
                }
            ),
            "openverse": CountingProvider(
                {
                    "Apis mellifera": [
                        {
                            "provider": "openverse",
                            "provider_record_id": "ov-1",
                            "matched_name": "Apis mellifera",
                            "matched_rank": "species",
                            "license_code": "public-domain",
                            "license_url": "https://creativecommons.org/publicdomain/zero/1.0/",
                            "attribution": "B",
                            "source_page_url": "https://example.org/ov/1",
                            "media_url": "https://static.inaturalist.org/photo-2.jpg",
                            "width": 3000,
                            "height": 2000,
                            "asset_type": "photo",
                        }
                    ],
                }
            ),
        }

        candidates, provider_errors = collect_candidates_for_species(
            species_name="Apis mellifera",
            args=make_image_args(
                out_dir=str(tmp_path / "out"),
                style="photo",
                source="inaturalist,openverse",
            ),
            sources=["inaturalist", "openverse"],
            providers=providers,
        )

        assert provider_errors == []
        assert len(candidates) == 3
        assert candidates[0]["provider"] == "inaturalist"
        assert call_counter["openverse"] == 0

    def test_collect_candidates_reuses_query_cache(self, tmp_path):
        call_counter = {"count": 0}

        class CountingProvider(DummyProvider):
            def fetch_candidates(self, species_name, fallback_rank="none"):
                call_counter["count"] += 1
                return super().fetch_candidates(
                    species_name, fallback_rank=fallback_rank
                )

        providers = {
            "eol": CountingProvider(
                {
                    "Apis mellifera": [
                        {
                            "provider": "eol",
                            "provider_record_id": "eol-1",
                            "matched_name": "Apis mellifera",
                            "matched_rank": "species",
                            "license_code": "cc-by",
                            "license_url": "https://creativecommons.org/licenses/by/4.0/",
                            "attribution": "A",
                            "source_page_url": "https://example.org/eol/1",
                            "media_url": "https://example.org/eol-photo.jpg",
                            "width": 1200,
                            "height": 800,
                            "asset_type": "photo",
                        }
                    ],
                }
            ),
        }
        args = make_image_args(
            out_dir=str(tmp_path / "out"), style="photo", source="eol"
        )

        first_candidates, first_errors = collect_candidates_for_species(
            species_name="Apis mellifera",
            args=args,
            sources=["eol"],
            providers=providers,
        )
        second_candidates, second_errors = collect_candidates_for_species(
            species_name="Apis mellifera",
            args=args,
            sources=["eol"],
            providers=providers,
        )

        assert first_errors == []
        assert second_errors == []
        assert call_counter["count"] == 1
        assert [candidate["provider_record_id"] for candidate in first_candidates] == [
            "eol-1"
        ]
        assert [candidate["provider_record_id"] for candidate in second_candidates] == [
            "eol-1"
        ]

    def test_query_cache_refetches_when_requested_capacity_increases(self, tmp_path):
        call_counter = {"count": 0}

        class CapacityProvider(DummyProvider):
            result_limit = 10

            def fetch_candidates(self, species_name, fallback_rank="none"):
                call_counter["count"] += 1
                return []

        provider = CapacityProvider({})
        providers = {"eol": provider}
        base_args = make_image_args(out_dir=str(tmp_path / "out"), source="eol")
        collect_candidates_for_species("Apis mellifera", base_args, ["eol"], providers)

        provider.result_limit = 24
        expanded_args = make_image_args(
            out_dir=str(tmp_path / "out"), source="eol", max_per_species=20
        )
        collect_candidates_for_species(
            "Apis mellifera", expanded_args, ["eol"], providers
        )

        assert call_counter["count"] == 2

    def test_query_cache_honors_expiration_and_explicit_refresh(
        self, monkeypatch, tmp_path
    ):
        call_counter = {"count": 0}
        current_time = {"value": 100.0}

        class CountingProvider(DummyProvider):
            result_limit = 10

            def fetch_candidates(self, species_name, fallback_rank="none"):
                call_counter["count"] += 1
                return []

        monkeypatch.setattr("nwkit.image.time.time", lambda: current_time["value"])
        providers = {"eol": CountingProvider({})}
        args = make_image_args(
            out_dir=str(tmp_path / "out"),
            source="eol",
            query_cache_max_age_hours=1.0,
        )
        collect_candidates_for_species("Apis mellifera", args, ["eol"], providers)
        current_time["value"] = 3801.0
        collect_candidates_for_species("Apis mellifera", args, ["eol"], providers)
        collect_candidates_for_species(
            "Apis mellifera",
            make_image_args(
                out_dir=str(tmp_path / "out"), source="eol", refresh_cache=True
            ),
            ["eol"],
            providers,
        )

        assert call_counter["count"] == 3

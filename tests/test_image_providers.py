import threading
import time
from concurrent.futures import ThreadPoolExecutor
from contextlib import contextmanager

import pytest

from nwkit import image as image_module
from nwkit.image import (
    BioiconsProvider,
    EOLProvider,
    IDigBioProvider,
    OpenverseProvider,
)
from tests.image_test_support import JSONResponse, make_image_args


class TestBioiconsProvider:
    def test_bioicons_provider_fetches_matching_svg_candidates(self, tmp_path):
        class RoutingSession:
            def get(
                self,
                url,
                params=None,
                timeout=None,
                headers=None,
                stream=None,
                allow_redirects=None,
            ):
                assert url.endswith("/git/trees/main")
                return JSONResponse(
                    {
                        "tree": [
                            {"path": "static/icons/cc-0/Animals/Ben-Murrell/Mouse.svg"},
                            {
                                "path": "static/icons/cc-by-4.0/Animals/DBCLS/Xenopus_laevis.svg"
                            },
                            {
                                "path": "static/icons/cc-by-3.0/Animals/Servier/rat-adult.svg"
                            },
                            {"path": "static/icons/categories.json"},
                        ],
                    },
                    url=url,
                )

            def close(self):
                return None

        provider = BioiconsProvider(
            session=RoutingSession(),
            args=make_image_args(out_dir=str(tmp_path / "out")),
        )
        candidates = provider.fetch_candidates("Mus musculus", fallback_rank="none")

        assert len(candidates) == 1
        candidate = candidates[0]
        assert candidate["provider"] == "bioicons"
        assert candidate["provider_record_id"] == "cc-0/Animals/Ben-Murrell/Mouse.svg"
        assert candidate["matched_name"] == "Mus musculus"
        assert candidate["matched_rank"] == "species"
        assert candidate["license_code"] == "public-domain"
        assert (
            candidate["license_url"]
            == "https://creativecommons.org/publicdomain/zero/1.0/"
        )
        assert candidate["attribution"] == "Ben Murrell"
        assert (
            candidate["media_url"]
            == "https://bioicons.com/icons/cc-0/Animals/Ben-Murrell/Mouse.svg"
        )
        assert candidate["asset_type"] == "silhouette"
        assert candidate["provider_quality"] > 0

    @pytest.mark.parametrize("parallel", [False, True], ids=["sequential", "parallel"])
    def test_bioicons_refresh_is_coalesced_across_provider_instances(
        self,
        tmp_path,
        parallel,
        monkeypatch,
    ):
        request_count = {"value": 0}
        request_count_lock = threading.Lock()
        catalog_lock = threading.Lock()

        @contextmanager
        def fast_catalog_lock(*args, **kwargs):
            with catalog_lock:
                yield

        monkeypatch.setattr(
            image_module,
            "acquire_exclusive_lock",
            fast_catalog_lock,
        )
        args = make_image_args(
            out_dir=str(tmp_path / "out"),
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
                    request_count["value"] += 1
                time.sleep(0.01)
                return JSONResponse(
                    {
                        "tree": [
                            {
                                "path": "static/icons/cc-0/Animals/Ben-Murrell/Mouse.svg",
                            }
                        ],
                    },
                    url=url,
                )

        providers = [
            BioiconsProvider(session=CountingSession(), args=args) for _ in range(4)
        ]
        if parallel:
            with ThreadPoolExecutor(max_workers=len(providers)) as executor:
                catalogs = list(
                    executor.map(
                        lambda provider: provider._load_catalog(),
                        providers,
                    )
                )
        else:
            catalogs = [provider._load_catalog() for provider in providers]

        assert request_count["value"] == 1
        assert all(catalog is catalogs[0] for catalog in catalogs)
        assert catalogs[0][0]["relative_path"].endswith("/Mouse.svg")


class TestEOLProvider:
    def test_eol_provider_fetches_candidates_from_page_media(self):
        class RoutingSession:
            def get(
                self,
                url,
                params=None,
                timeout=None,
                headers=None,
                stream=None,
                allow_redirects=None,
            ):
                if url.endswith("/search/1.0.json"):
                    return JSONResponse(
                        {
                            "results": [
                                {
                                    "id": 491995,
                                    "title": "Hapalochlaena lunulata",
                                }
                            ],
                        },
                        url=url,
                    )
                if url.endswith("/pages/1.0/491995.json"):
                    assert params["images_per_page"] == 10
                    return JSONResponse(
                        {
                            "taxonConcept": {
                                "dataObjects": [
                                    {
                                        "identifier": "EOL-media-509-demo",
                                        "dataObjectVersionID": 28895676,
                                        "dataType": "http://purl.org/dc/dcmitype/StillImage",
                                        "mediumType": "image",
                                        "dataRating": "2.5",
                                        "vettedStatus": "Trusted",
                                        "license": "http://creativecommons.org/licenses/by/3.0/",
                                        "rightsHolder": "Elias Levy",
                                        "source": "https://commons.wikimedia.org/wiki/File:Blue-Ringed_Octopus.jpg",
                                        "mediaURL": "https://upload.wikimedia.org/blue-ringed.jpg",
                                        "agents": [
                                            {
                                                "full_name": "Elias Levy",
                                                "role": "creator",
                                            }
                                        ],
                                    }
                                ],
                            },
                        },
                        url=url,
                    )
                raise AssertionError("Unexpected URL: {}".format(url))

            def close(self):
                return None

        provider = EOLProvider(
            session=RoutingSession(), args=make_image_args(max_per_species=1)
        )
        candidates = provider.fetch_candidates(
            "Hapalochlaena lunulata", fallback_rank="none"
        )

        assert len(candidates) == 1
        candidate = candidates[0]
        assert candidate["provider"] == "eol"
        assert candidate["provider_record_id"] == "28895676"
        assert candidate["matched_name"] == "Hapalochlaena lunulata"
        assert candidate["matched_rank"] == "species"
        assert candidate["license_code"] == "cc-by"
        assert candidate["license_url"] == "http://creativecommons.org/licenses/by/3.0/"
        assert candidate["attribution"] == "Elias Levy"
        assert (
            candidate["source_page_url"]
            == "https://commons.wikimedia.org/wiki/File:Blue-Ringed_Octopus.jpg"
        )
        assert candidate["media_url"] == "https://upload.wikimedia.org/blue-ringed.jpg"
        assert candidate["asset_type"] == "photo"
        assert candidate["provider_quality"] > 0


class TestIDigBioProvider:
    def test_idigbio_provider_fetches_candidates_from_media_search(self):
        class RoutingSession:
            def post(
                self,
                url,
                json=None,
                timeout=None,
                headers=None,
                stream=None,
                allow_redirects=None,
            ):
                assert url.endswith("/search/media")
                assert json["rq"] == {"scientificname": "Panthera leo"}
                assert json["limit"] == 10
                return JSONResponse(
                    {
                        "items": [
                            {
                                "uuid": "idigbio-1",
                                "data": {
                                    "ac:accessURI": "https://collections.example.org/panthera_leo.jpg",
                                    "dcterms:identifier": "https://collections.example.org/object/1",
                                    "xmpRights:UsageTerms": "https://creativecommons.org/publicdomain/zero/1.0/",
                                    "dcterms:rights": "CC0",
                                    "xmpRights:Owner": "Museum Example",
                                    "xmpRights:WebStatement": "https://collections.example.org/license",
                                    "dwc:scientificName": "Panthera leo",
                                    "exif:PixelXDimension": "2048",
                                    "exif:PixelYDimension": "1024",
                                },
                                "indexTerms": {
                                    "mediatype": "images",
                                    "dqs": 0.8,
                                },
                            }
                        ],
                    },
                    url=url,
                )

            def close(self):
                return None

        provider = IDigBioProvider(
            session=RoutingSession(), args=make_image_args(max_per_species=1)
        )
        candidates = provider.fetch_candidates("Panthera leo", fallback_rank="none")

        assert len(candidates) == 1
        candidate = candidates[0]
        assert candidate["provider"] == "idigbio"
        assert candidate["provider_record_id"] == "idigbio-1"
        assert candidate["matched_name"] == "Panthera leo"
        assert candidate["matched_rank"] == "species"
        assert candidate["license_code"] == "public-domain"
        assert (
            candidate["license_url"]
            == "https://creativecommons.org/publicdomain/zero/1.0/"
        )
        assert candidate["attribution"] == "Museum Example"
        assert (
            candidate["source_page_url"] == "https://collections.example.org/object/1"
        )
        assert (
            candidate["media_url"] == "https://collections.example.org/panthera_leo.jpg"
        )
        assert candidate["width"] == 2048
        assert candidate["height"] == 1024
        assert candidate["asset_type"] == "photo"
        assert candidate["provider_quality"] > 0


class TestOpenverseProvider:
    def test_openverse_provider_fetches_relevant_candidates(self):
        class RoutingSession:
            def get(
                self,
                url,
                params=None,
                timeout=None,
                headers=None,
                stream=None,
                allow_redirects=None,
            ):
                assert url.endswith("/images/")
                assert params["q"] == "Danio rerio"
                assert params["page_size"] == 10
                return JSONResponse(
                    {
                        "results": [
                            {
                                "id": "openverse-1",
                                "title": "Danio rerio (Peix Zebra/Zebrafish)",
                                "url": "https://cdn.example.org/danio.jpg",
                                "creator": "berarma",
                                "license": "by-sa",
                                "license_url": "https://creativecommons.org/licenses/by-sa/2.0/",
                                "foreign_landing_url": "https://www.flickr.com/photos/example/1",
                                "width": 1024,
                                "height": 683,
                                "tags": [],
                                "fields_matched": ["title"],
                            },
                            {
                                "id": "openverse-2",
                                "title": "Zebra fish in aquarium",
                                "url": "https://cdn.example.org/irrelevant.jpg",
                                "creator": "someone",
                                "license": "by",
                                "license_url": "https://creativecommons.org/licenses/by/4.0/",
                                "foreign_landing_url": "https://www.flickr.com/photos/example/2",
                                "width": 1200,
                                "height": 800,
                                "tags": [],
                                "fields_matched": ["title"],
                            },
                        ],
                    },
                    url=url,
                )

            def close(self):
                return None

        provider = OpenverseProvider(
            session=RoutingSession(), args=make_image_args(max_per_species=1)
        )
        candidates = provider.fetch_candidates("Danio rerio", fallback_rank="none")

        assert len(candidates) == 1
        candidate = candidates[0]
        assert candidate["provider"] == "openverse"
        assert candidate["provider_record_id"] == "openverse-1"
        assert candidate["matched_name"] == "Danio rerio"
        assert candidate["matched_rank"] == "species"
        assert candidate["license_code"] == "cc-by-sa"
        assert (
            candidate["license_url"]
            == "https://creativecommons.org/licenses/by-sa/2.0/"
        )
        assert candidate["attribution"] == "berarma"
        assert candidate["source_page_url"] == "https://www.flickr.com/photos/example/1"
        assert candidate["media_url"] == "https://cdn.example.org/danio.jpg"
        assert candidate["asset_type"] == "photo"
        assert candidate["provider_quality"] > 0

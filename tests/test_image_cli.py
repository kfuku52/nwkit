import io
import os

import pytest
import requests

from nwkit import image as image_module
from nwkit.image import (
    build_local_media_filename,
    extract_species_mapping,
    image_main,
    resolve_image_cache_dir,
)
from tests.image_test_support import (
    DummyProvider,
    DummySession,
    make_image_args,
    read_tsv,
    write_valid_test_media,
)


class TestImageMain:
    def test_resolve_image_cache_dir_uses_download_dir(self, tmp_path):
        args = make_image_args(download_dir=str(tmp_path / "cache"))
        assert resolve_image_cache_dir(args) == str(
            tmp_path / "cache" / "nwkit" / "image-cache"
        )

    def test_extract_species_mapping_reports_unparsable_labels(self, tmp_path):
        tree_path = tmp_path / "tree.nwk"
        tree_path.write_text("(Homo_sapiens_A,BadLabel);")

        from nwkit.util import read_tree

        tree = read_tree(
            str(tree_path), format="auto", quoted_node_names=True, quiet=True
        )
        leaf_to_species, unmatched_rows = extract_species_mapping(tree)

        assert leaf_to_species == {"Homo_sapiens_A": "Homo sapiens"}
        assert unmatched_rows == [
            {
                "leaf_name": "BadLabel",
                "species_name": "",
                "reason": "unparsable leaf label",
                "details": "Expected the configured species parser or a matching --species-name-tsv entry.",
            }
        ]

    def test_extract_species_mapping_accepts_custom_species_regex(self, tmp_path):
        tree_path = tmp_path / "tree.nwk"
        tree_path.write_text("(Homo.sapiens|A,Mus.musculus|B);")

        from nwkit.util import read_tree

        tree = read_tree(
            str(tree_path), format="auto", quoted_node_names=True, quiet=True
        )
        leaf_to_species, unmatched_rows = extract_species_mapping(
            tree,
            species_regex=r"^([A-Za-z]+)\.([A-Za-z]+)\|",
        )

        assert leaf_to_species == {
            "Homo.sapiens|A": "Homo sapiens",
            "Mus.musculus|B": "Mus musculus",
        }
        assert unmatched_rows == []

    def test_image_main_writes_manifest_attribution_and_unmatched(
        self, monkeypatch, tmp_path
    ):
        tree_path = tmp_path / "tree.nwk"
        tree_path.write_text(
            "((Homo_sapiens_A,Homo_sapiens_B),Panthera_leo_C,BadLabel);"
        )
        out_dir = tmp_path / "out"

        phylopic_candidates = {
            "Homo sapiens": [
                {
                    "provider": "phylopic",
                    "provider_record_id": "phy-1",
                    "matched_name": "Homo sapiens",
                    "matched_rank": "species",
                    "license_code": "cc-by",
                    "license_url": "https://creativecommons.org/licenses/by/4.0/",
                    "attribution": "PhyloPic Artist",
                    "source_page_url": "https://api.phylopic.org/images/phy-1",
                    "media_url": "https://images.phylopic.org/homo.svg",
                    "width": 1200,
                    "height": 800,
                    "asset_type": "silhouette",
                    "is_primary": True,
                    "is_vector": True,
                }
            ],
        }
        inaturalist_candidates = {
            "Panthera leo": [
                {
                    "provider": "inaturalist",
                    "provider_record_id": "inat-1",
                    "matched_name": "Panthera leo",
                    "matched_rank": "species",
                    "license_code": "cc-by-nc-sa",
                    "license_url": "https://creativecommons.org/licenses/by-nc-sa/4.0/",
                    "attribution": "(c) Example Photographer, some rights reserved",
                    "source_page_url": "https://www.inaturalist.org/observations/1",
                    "media_url": "https://static.inaturalist.org/photos/1/original.jpg",
                    "width": 1600,
                    "height": 1200,
                    "asset_type": "photo",
                }
            ],
        }

        def fake_build_providers(args, sources, session=None):
            providers = {
                "phylopic": DummyProvider(phylopic_candidates),
                "inaturalist": DummyProvider(inaturalist_candidates),
            }
            return DummySession(), None, providers

        def fake_download_media(
            session,
            media_url,
            destination_path,
            cache_path=None,
            max_download_bytes=None,
            provider=None,
            **kwargs,
        ):
            write_valid_test_media(destination_path)
            return "downloaded"

        monkeypatch.setattr("nwkit.image.build_providers", fake_build_providers)
        monkeypatch.setattr("nwkit.image.download_media", fake_download_media)

        args = make_image_args(
            infile=str(tree_path),
            out_dir=str(out_dir),
            source="phylopic,inaturalist",
        )
        image_main(args)

        manifest_rows = read_tsv(out_dir / "manifest.tsv")
        unmatched_rows = read_tsv(out_dir / "unmatched.tsv")
        attribution_text = (out_dir / "ATTRIBUTION.md").read_text()

        assert len(manifest_rows) == 3
        assert {row["leaf_name"] for row in manifest_rows} == {
            "Homo_sapiens_A",
            "Homo_sapiens_B",
            "Panthera_leo_C",
        }
        assert {row["species_name"] for row in manifest_rows} == {
            "Homo sapiens",
            "Panthera leo",
        }
        assert {row["status"] for row in manifest_rows} == {"downloaded"}
        assert any(row["local_path"].endswith(".svg") for row in manifest_rows)
        assert any(row["local_path"].endswith(".jpg") for row in manifest_rows)
        homo_row = next(
            row for row in manifest_rows if row["species_name"] == "Homo sapiens"
        )
        assert homo_row["is_primary"] == "yes"
        assert homo_row["is_vector"] == "yes"
        assert (
            homo_row["selection_reason"] == "exact_species_match;phylopic_primary_image"
        )

        assert unmatched_rows == [
            {
                "leaf_name": "BadLabel",
                "species_name": "",
                "reason": "unparsable leaf label",
                "details": "Expected the configured species parser or a matching --species-name-tsv entry.",
            }
        ]
        assert "Homo sapiens" in attribution_text
        assert "Panthera leo" in attribution_text

    def test_attribution_preserves_distinct_records_for_shared_media(self, tmp_path):
        from nwkit.image import write_attribution_markdown

        path = tmp_path / "ATTRIBUTION.md"
        shared_path = "images/shared.jpg"
        rows = [
            {
                "local_path": shared_path,
                "species_name": "Species alpha",
                "provider": "provider-a",
                "matched_name": "Species alpha",
                "matched_rank": "species",
                "attribution": "Alice",
                "license_code": "cc-by",
                "license_url": "https://example.org/license-a",
                "source_page_url": "https://example.org/source-a",
            },
            {
                "local_path": shared_path,
                "species_name": "Species beta",
                "provider": "provider-b",
                "matched_name": "Species beta",
                "matched_rank": "species",
                "attribution": "Bob",
                "license_code": "cc-by-sa",
                "license_url": "https://example.org/license-b",
                "source_page_url": "https://example.org/source-b",
            },
        ]

        write_attribution_markdown(str(path), rows)

        text = path.read_text()
        assert text.count("### Attribution record") == 2
        assert "Creator / attribution: Alice" in text
        assert "Creator / attribution: Bob" in text
        assert "License: cc-by\n" in text
        assert "License: cc-by-sa\n" in text

    def test_attribution_escapes_provider_controlled_markdown(self, tmp_path):
        from nwkit.image import write_attribution_markdown

        path = tmp_path / "ATTRIBUTION.md"
        rows = [
            {
                "local_path": "images/example.jpg",
                "species_name": "Species alpha\n## Forged species",
                "provider": "provider",
                "matched_name": "<b>Species alpha</b>",
                "matched_rank": "species",
                "attribution": "</p>\n## Fake license\n![track](https://evil.example/pixel)",
                "license_code": "cc-by",
                "license_url": "https://example.org/license",
                "source_page_url": "https://example.org/source",
            }
        ]

        write_attribution_markdown(str(path), rows)

        text = path.read_text()
        assert "\n## Forged species" not in text
        assert "\n## Fake license" not in text
        assert "<b>" not in text
        assert "</p>" not in text
        assert "![track](" not in text
        assert "&lt;b&gt;Species alpha&lt;/b&gt;" in text


class TestMediaFilenames:
    @pytest.mark.parametrize(
        ("url", "expected"),
        [
            ("https://example.org/image.JPG?size=large", ".jpg"),
            ("https://example.org/image.jpeg", ".jpeg"),
            ("https://example.org/image.tiff", ".tiff"),
            ("https://example.org/image.svg#preview", ".svg"),
            ("https://example.org/image.php", ".bin"),
            ("https://example.org/image.tar.gz", ".bin"),
            ("https://example.org/image.\0jpg", ".bin"),
            ("https://example.org/image\0.jpg", ".bin"),
            ("https://example.org/image." + ("x" * 4096), ".bin"),
        ],
    )
    def test_infer_extension_only_accepts_known_image_suffixes(self, url, expected):
        assert image_module.infer_extension(url) == expected

    def test_local_filename_distinguishes_urls_with_the_same_provider_record_id(self):
        candidate = {
            "provider": "gbif",
            "provider_record_id": "12345",
        }

        first = build_local_media_filename(
            "Species alpha", candidate, "https://example.org/front.jpg"
        )
        second = build_local_media_filename(
            "Species alpha", candidate, "https://example.org/back.jpg"
        )

        assert first != second
        assert first.endswith(".jpg")
        assert second.endswith(".jpg")

    def test_long_filename_components_are_bounded_and_collision_resistant(self):
        long_species = "Species_" + ("a" * 500)
        shared_prefix = "record-" + ("x" * 500)
        first_candidate = {
            "provider": "provider-" + ("p" * 500),
            "provider_record_id": shared_prefix + "-first",
        }
        second_candidate = {
            **first_candidate,
            "provider_record_id": shared_prefix + "-second",
        }

        first_local = build_local_media_filename(
            long_species,
            first_candidate,
            "https://example.org/image.jpg",
        )
        second_local = build_local_media_filename(
            long_species,
            second_candidate,
            "https://example.org/image.jpg",
        )
        cache_name = os.path.basename(
            image_module.build_media_cache_path(
                "/tmp/cache",
                "https://example.org/image.jpg",
                first_candidate["provider"],
                first_candidate["provider_record_id"],
            )
        )
        query_name = os.path.basename(
            image_module.build_query_cache_path(
                "/tmp/cache",
                first_candidate["provider"],
                long_species,
                "fallback-" + ("f" * 500),
            )
        )

        assert first_local != second_local
        assert (
            image_module.sanitize_filename_component("Apis mellifera")
            == "Apis_mellifera"
        )
        for filename in (first_local, second_local, cache_name, query_name):
            assert len(os.fsencode(filename)) < 255

    def test_image_main_uses_name_tsv_override_and_strict_mode(
        self, monkeypatch, tmp_path
    ):
        tree_path = tmp_path / "tree.nwk"
        tree_path.write_text("(Sample_1,Unknown);")
        name_tsv = tmp_path / "names.tsv"
        name_tsv.write_text("leaf_name\tspecies_name\nSample_1\tApis mellifera\n")
        out_dir = tmp_path / "out"

        phylopic_candidates = {
            "Apis mellifera": [
                {
                    "provider": "phylopic",
                    "provider_record_id": "phy-2",
                    "matched_name": "Apis mellifera",
                    "matched_rank": "species",
                    "license_code": "public-domain",
                    "license_url": "https://creativecommons.org/publicdomain/zero/1.0/",
                    "attribution": "PhyloPic Artist",
                    "source_page_url": "https://api.phylopic.org/images/phy-2",
                    "media_url": "https://images.phylopic.org/apis.svg",
                    "width": 1000,
                    "height": 900,
                    "asset_type": "silhouette",
                }
            ],
        }

        def fake_build_providers(args, sources, session=None):
            providers = {
                "phylopic": DummyProvider(phylopic_candidates),
            }
            return DummySession(), None, providers

        def fake_download_media(
            session,
            media_url,
            destination_path,
            cache_path=None,
            max_download_bytes=None,
            provider=None,
            **kwargs,
        ):
            write_valid_test_media(destination_path)
            return "downloaded"

        monkeypatch.setattr("nwkit.image.build_providers", fake_build_providers)
        monkeypatch.setattr("nwkit.image.download_media", fake_download_media)

        args = make_image_args(
            infile=str(tree_path),
            out_dir=str(out_dir),
            species_name_tsv=str(name_tsv),
            source="phylopic",
            fail_on_missing=True,
        )
        with pytest.raises(ValueError, match="could not be resolved"):
            image_main(args)

        manifest_rows = read_tsv(out_dir / "manifest.tsv")
        unmatched_rows = read_tsv(out_dir / "unmatched.tsv")

        assert manifest_rows[0]["leaf_name"] == "Sample_1"
        assert manifest_rows[0]["species_name"] == "Apis mellifera"
        assert unmatched_rows == [
            {
                "leaf_name": "Unknown",
                "species_name": "",
                "reason": "unparsable leaf label",
                "details": "Expected the configured species parser or a matching --species-name-tsv entry.",
            }
        ]

    def test_image_main_reports_filtered_by_license(self, monkeypatch, tmp_path):
        tree_path = tmp_path / "tree.nwk"
        tree_path.write_text("(Panthera_leo_A);")
        out_dir = tmp_path / "out"

        inaturalist_candidates = {
            "Panthera leo": [
                {
                    "provider": "inaturalist",
                    "provider_record_id": "inat-2",
                    "matched_name": "Panthera leo",
                    "matched_rank": "species",
                    "license_code": "all-rights-reserved",
                    "license_url": "",
                    "attribution": "(c) Example Photographer, all rights reserved",
                    "source_page_url": "https://www.inaturalist.org/observations/2",
                    "media_url": "https://static.inaturalist.org/photos/2/original.jpg",
                    "width": 1600,
                    "height": 1200,
                    "asset_type": "photo",
                }
            ],
        }

        def fake_build_providers(args, sources, session=None):
            providers = {
                "inaturalist": DummyProvider(inaturalist_candidates),
            }
            return DummySession(), None, providers

        monkeypatch.setattr("nwkit.image.build_providers", fake_build_providers)

        args = make_image_args(
            infile=str(tree_path),
            out_dir=str(out_dir),
            style="photo",
            source="inaturalist",
        )
        image_main(args)

        manifest_rows = read_tsv(out_dir / "manifest.tsv")
        unmatched_rows = read_tsv(out_dir / "unmatched.tsv")

        assert manifest_rows == []
        assert unmatched_rows == [
            {
                "leaf_name": "Panthera_leo_A",
                "species_name": "Panthera leo",
                "reason": "only disallowed-license candidates found",
                "details": "",
            }
        ]

    def test_image_main_reuses_shared_download_cache_across_runs(
        self, monkeypatch, tmp_path
    ):
        tree_path = tmp_path / "tree.nwk"
        tree_path.write_text("(Apis_mellifera_A);")
        shared_download_dir = tmp_path / "shared-cache"
        out_dir1 = tmp_path / "out1"
        out_dir2 = tmp_path / "out2"

        phylopic_candidates = {
            "Apis mellifera": [
                {
                    "provider": "phylopic",
                    "provider_record_id": "phy-cache",
                    "matched_name": "Apis mellifera",
                    "matched_rank": "species",
                    "license_code": "cc-by",
                    "license_url": "https://creativecommons.org/licenses/by/4.0/",
                    "attribution": "PhyloPic Artist",
                    "source_page_url": "https://api.phylopic.org/images/phy-cache",
                    "media_url": "https://images.phylopic.org/apis-cache.svg",
                    "width": 1000,
                    "height": 900,
                    "asset_type": "silhouette",
                }
            ],
        }
        call_counter = {"count": 0}

        class DownloadOnlyResponse:
            def raise_for_status(self):
                return None

            def iter_content(self, chunk_size=65536):
                yield b'<svg xmlns="http://www.w3.org/2000/svg"></svg>'

        class CountingSession:
            def get(
                self,
                media_url,
                stream=True,
                timeout=None,
                headers=None,
                allow_redirects=False,
            ):
                call_counter["count"] += 1
                return DownloadOnlyResponse()

            def close(self):
                return None

        def fake_build_providers(args, sources, session=None):
            providers = {
                "phylopic": DummyProvider(phylopic_candidates),
            }
            return CountingSession(), None, providers

        monkeypatch.setattr("nwkit.image.build_providers", fake_build_providers)
        monkeypatch.setattr("nwkit.image.build_download_session", CountingSession)

        args1 = make_image_args(
            infile=str(tree_path),
            out_dir=str(out_dir1),
            source="phylopic",
            download_dir=str(shared_download_dir),
        )
        image_main(args1)

        args2 = make_image_args(
            infile=str(tree_path),
            out_dir=str(out_dir2),
            source="phylopic",
            download_dir=str(shared_download_dir),
        )
        image_main(args2)

        manifest_rows_1 = read_tsv(out_dir1 / "manifest.tsv")
        manifest_rows_2 = read_tsv(out_dir2 / "manifest.tsv")

        assert call_counter["count"] == 1
        assert manifest_rows_1[0]["status"] == "downloaded"
        assert manifest_rows_2[0]["status"] == "cached"
        assert (shared_download_dir / "nwkit" / "image-cache").is_dir()

    @pytest.mark.parametrize(
        ("use_shared_cache", "expected_downloads"),
        [(False, 1), (True, 1)],
        ids=["reuse-local-raw-cache", "reuse-shared-raw-cache"],
    )
    def test_image_main_reprocesses_raw_media_when_options_change(
        self,
        monkeypatch,
        tmp_path,
        use_shared_cache,
        expected_downloads,
    ):
        Image = pytest.importorskip("PIL.Image")
        tree_path = tmp_path / "tree.nwk"
        tree_path.write_text("(Apis_mellifera_A);")
        out_dir = tmp_path / "out"
        shared_download_dir = tmp_path / "shared-cache"
        image_bytes = io.BytesIO()
        Image.new("RGB", (80, 40), "black").save(image_bytes, format="PNG")
        raw_png = image_bytes.getvalue()
        call_counter = {"count": 0}
        candidates = {
            "Apis mellifera": [
                {
                    "provider": "phylopic",
                    "provider_record_id": "phy-options",
                    "matched_name": "Apis mellifera",
                    "matched_rank": "species",
                    "license_code": "cc-by",
                    "license_url": "https://creativecommons.org/licenses/by/4.0/",
                    "attribution": "PhyloPic Artist",
                    "source_page_url": "https://api.phylopic.org/images/phy-options",
                    "media_url": "https://images.phylopic.org/apis-options.png",
                    "width": 80,
                    "height": 40,
                    "asset_type": "silhouette",
                }
            ],
        }

        class ImageResponse:
            status_code = 200
            headers = {"Content-Type": "image/png"}
            url = "https://images.phylopic.org/apis-options.png"

            def raise_for_status(self):
                return None

            def iter_content(self, chunk_size=65536):
                yield raw_png

            def close(self):
                return None

        class CountingSession:
            def get(self, url, **kwargs):
                call_counter["count"] += 1
                return ImageResponse()

            def close(self):
                return None

        def fake_build_providers(args, sources, session=None):
            return (
                DummySession(),
                None,
                {
                    "phylopic": DummyProvider(candidates),
                },
            )

        monkeypatch.setattr("nwkit.image.build_providers", fake_build_providers)
        monkeypatch.setattr("nwkit.image.build_download_session", CountingSession)
        download_dir = str(shared_download_dir) if use_shared_cache else "auto"

        image_main(
            make_image_args(
                infile=str(tree_path),
                out_dir=str(out_dir),
                source="phylopic",
                download_dir=download_dir,
                output_format="png",
                max_edge=20,
            )
        )
        first_manifest = read_tsv(out_dir / "manifest.tsv")
        first_path = out_dir / first_manifest[0]["local_path"]
        with Image.open(first_path) as first_image:
            assert first_image.size == (20, 10)

        image_main(
            make_image_args(
                infile=str(tree_path),
                out_dir=str(out_dir),
                source="phylopic",
                download_dir=download_dir,
                output_format="png",
                max_edge=40,
            )
        )
        second_manifest = read_tsv(out_dir / "manifest.tsv")
        second_path = out_dir / second_manifest[0]["local_path"]
        with Image.open(second_path) as second_image:
            assert second_image.size == (40, 20)

        assert call_counter["count"] == expected_downloads

    def test_custom_manifest_and_attribution_paths_feed_draw_with_ranked_rows(
        self,
        monkeypatch,
        tmp_path,
    ):
        tree_path = tmp_path / "tree.nwk"
        tree_path.write_text("(Apis_mellifera_A);")
        out_dir = tmp_path / "image-output"
        metadata_dir = tmp_path / "metadata"
        manifest_path = metadata_dir / "ranked-images.tsv"
        attribution_path = metadata_dir / "ATTRIBUTION.md"
        candidates = {
            "Apis mellifera": [
                {
                    "provider": "phylopic",
                    "provider_record_id": "phy-first",
                    "matched_name": "Apis mellifera",
                    "matched_rank": "species",
                    "license_code": "cc-by",
                    "license_url": "https://creativecommons.org/licenses/by/4.0/",
                    "attribution": "First Artist",
                    "source_page_url": "https://api.phylopic.org/images/phy-first",
                    "media_url": "https://images.phylopic.org/apis-first.png",
                    "width": 80,
                    "height": 40,
                    "asset_type": "silhouette",
                    "provider_quality": 20,
                },
                {
                    "provider": "phylopic",
                    "provider_record_id": "phy-second",
                    "matched_name": "Apis mellifera",
                    "matched_rank": "species",
                    "license_code": "cc-by",
                    "license_url": "https://creativecommons.org/licenses/by/4.0/",
                    "attribution": "Second Artist",
                    "source_page_url": "https://api.phylopic.org/images/phy-second",
                    "media_url": "https://images.phylopic.org/apis-second.png",
                    "width": 80,
                    "height": 40,
                    "asset_type": "silhouette",
                    "provider_quality": 10,
                },
            ],
        }

        def fake_build_providers(args, sources, session=None):
            return (
                DummySession(),
                None,
                {
                    "phylopic": DummyProvider(candidates),
                },
            )

        def fake_download_media(session, media_url, destination_path, **kwargs):
            write_valid_test_media(destination_path)
            return "downloaded"

        monkeypatch.setattr("nwkit.image.build_providers", fake_build_providers)
        monkeypatch.setattr("nwkit.image.download_media", fake_download_media)

        image_main(
            make_image_args(
                infile=str(tree_path),
                out_dir=str(out_dir),
                source="phylopic",
                max_per_species=2,
                manifest_out=str(manifest_path),
                attribution_out=str(attribution_path),
            )
        )

        manifest_rows = read_tsv(manifest_path)
        assert len(manifest_rows) == 2
        assert [row["provider_record_id"] for row in manifest_rows] == [
            "phy-first",
            "phy-second",
        ]
        for row in manifest_rows:
            assert os.path.isfile(metadata_dir / row["local_path"])
        local_file_lines = [
            line.split(": ", 1)[1]
            for line in attribution_path.read_text().splitlines()
            if line.startswith("Local file: ")
        ]
        assert len(local_file_lines) == 2
        assert all(os.path.isfile(metadata_dir / path) for path in local_file_lines)

        from nwkit.cli import main as cli_main

        draw_path = tmp_path / "ranked-images.svg"
        cli_main(
            [
                "draw",
                "-i",
                str(tree_path),
                "--species-overlap-node-plot",
                "no",
                "--tip-image-manifest",
                str(manifest_path),
                "-o",
                str(draw_path),
            ]
        )

        assert draw_path.read_text(encoding="utf-8").count("<image") == 1

    def test_image_main_does_not_rebuild_providers_for_download_stage(
        self, monkeypatch, tmp_path
    ):
        tree_path = tmp_path / "tree.nwk"
        tree_path.write_text("(Apis_mellifera_A);")
        out_dir = tmp_path / "out"
        call_counter = {"count": 0}

        phylopic_candidates = {
            "Apis mellifera": [
                {
                    "provider": "phylopic",
                    "provider_record_id": "phy-once",
                    "matched_name": "Apis mellifera",
                    "matched_rank": "species",
                    "license_code": "cc-by",
                    "license_url": "https://creativecommons.org/licenses/by/4.0/",
                    "attribution": "PhyloPic Artist",
                    "source_page_url": "https://api.phylopic.org/images/phy-once",
                    "media_url": "https://images.phylopic.org/apis-once.svg",
                    "width": 1000,
                    "height": 900,
                    "asset_type": "silhouette",
                }
            ],
        }

        def fake_build_providers(args, sources, session=None):
            call_counter["count"] += 1
            providers = {
                "phylopic": DummyProvider(phylopic_candidates),
            }
            return DummySession(), None, providers

        def fake_download_media(
            session,
            media_url,
            destination_path,
            cache_path=None,
            max_download_bytes=None,
            provider=None,
            **kwargs,
        ):
            write_valid_test_media(destination_path)
            return "downloaded"

        monkeypatch.setattr("nwkit.image.build_providers", fake_build_providers)
        monkeypatch.setattr("nwkit.image.download_media", fake_download_media)

        image_main(
            make_image_args(
                infile=str(tree_path),
                out_dir=str(out_dir),
                source="phylopic",
            )
        )

        assert call_counter["count"] == 1

    def test_image_main_reuses_same_media_across_species_within_run(
        self, monkeypatch, tmp_path
    ):
        tree_path = tmp_path / "tree.nwk"
        tree_path.write_text("(Felis_catus_A,Panthera_leo_B);")
        out_dir = tmp_path / "out"

        shared_candidates = {
            "Felis catus": [
                {
                    "provider": "wikimedia",
                    "provider_record_id": "wiki-shared-cat",
                    "matched_name": "Felis catus",
                    "matched_rank": "species",
                    "license_code": "cc-by",
                    "license_url": "https://creativecommons.org/licenses/by/4.0/",
                    "attribution": "Shared Photographer",
                    "source_page_url": "https://commons.wikimedia.org/wiki/File:SharedCatLion.jpg",
                    "media_url": "https://upload.wikimedia.org/shared-cat-lion.jpg",
                    "width": 1500,
                    "height": 1000,
                    "asset_type": "photo",
                }
            ],
            "Panthera leo": [
                {
                    "provider": "wikimedia",
                    "provider_record_id": "wiki-shared-lion",
                    "matched_name": "Panthera leo",
                    "matched_rank": "species",
                    "license_code": "cc-by",
                    "license_url": "https://creativecommons.org/licenses/by/4.0/",
                    "attribution": "Shared Photographer",
                    "source_page_url": "https://commons.wikimedia.org/wiki/File:SharedCatLion.jpg",
                    "media_url": "https://upload.wikimedia.org/shared-cat-lion.jpg",
                    "width": 1500,
                    "height": 1000,
                    "asset_type": "photo",
                }
            ],
        }
        call_counter = {"count": 0}

        def fake_build_providers(args, sources, session=None):
            providers = {
                "wikimedia": DummyProvider(shared_candidates),
            }
            return DummySession(), None, providers

        def fake_download_media(
            session,
            media_url,
            destination_path,
            cache_path=None,
            max_download_bytes=None,
            provider=None,
            **kwargs,
        ):
            call_counter["count"] += 1
            write_valid_test_media(destination_path)
            return "downloaded"

        monkeypatch.setattr("nwkit.image.build_providers", fake_build_providers)
        monkeypatch.setattr("nwkit.image.download_media", fake_download_media)

        image_main(
            make_image_args(
                infile=str(tree_path),
                out_dir=str(out_dir),
                style="photo",
                source="wikimedia",
            )
        )

        manifest_rows = read_tsv(out_dir / "manifest.tsv")
        rows_by_species = {row["species_name"]: row for row in manifest_rows}

        assert call_counter["count"] == 1
        assert rows_by_species["Felis catus"]["status"] == "downloaded"
        assert rows_by_species["Panthera leo"]["status"] == "reused"
        assert (
            rows_by_species["Felis catus"]["local_path"]
            == rows_by_species["Panthera leo"]["local_path"]
        )

    def test_image_main_records_resolved_extension_from_download(
        self, monkeypatch, tmp_path
    ):
        tree_path = tmp_path / "tree.nwk"
        tree_path.write_text("(Cyanophora_paradoxa_A);")
        out_dir = tmp_path / "out"

        ncbi_candidates = {
            "Cyanophora paradoxa": [
                {
                    "provider": "ncbi",
                    "provider_record_id": "64365",
                    "matched_name": "Cyanophora paradoxa",
                    "matched_rank": "species",
                    "license_code": "cc-by-sa",
                    "license_url": "https://creativecommons.org/licenses/by-sa/3.0/",
                    "attribution": "NCBI",
                    "source_page_url": "https://www.ncbi.nlm.nih.gov/Taxonomy/taxi/images/4",
                    "media_url": "https://www.ncbi.nlm.nih.gov/Taxonomy/taxi/images/4",
                    "width": None,
                    "height": None,
                    "asset_type": "photo",
                }
            ],
        }

        def fake_build_providers(args, sources, session=None):
            providers = {
                "ncbi": DummyProvider(ncbi_candidates),
            }
            return DummySession(), None, providers

        def fake_download_media(
            session,
            media_url,
            destination_path,
            cache_path=None,
            max_download_bytes=None,
            provider=None,
            **kwargs,
        ):
            resolved_path = destination_path[:-4] + ".jpg"
            write_valid_test_media(resolved_path)
            return {"status": "downloaded", "destination_path": resolved_path}

        monkeypatch.setattr("nwkit.image.build_providers", fake_build_providers)
        monkeypatch.setattr("nwkit.image.download_media", fake_download_media)

        args = make_image_args(
            infile=str(tree_path),
            out_dir=str(out_dir),
            source="ncbi",
            style="photo",
        )
        image_main(args)

        manifest_rows = read_tsv(out_dir / "manifest.tsv")
        assert manifest_rows[0]["local_path"].endswith(".jpg")
        assert (out_dir / manifest_rows[0]["local_path"]).exists()

    def test_image_main_falls_back_to_next_candidate_after_download_error(
        self, monkeypatch, tmp_path
    ):
        tree_path = tmp_path / "tree.nwk"
        tree_path.write_text("(Panthera_leo_A);")
        out_dir = tmp_path / "out"

        inaturalist_candidates = {
            "Panthera leo": [
                {
                    "provider": "inaturalist",
                    "provider_record_id": "inat-fail",
                    "matched_name": "Panthera leo",
                    "matched_rank": "species",
                    "license_code": "cc-by",
                    "license_url": "https://creativecommons.org/licenses/by/4.0/",
                    "attribution": "Failing Photographer",
                    "source_page_url": "https://www.inaturalist.org/observations/fail",
                    "media_url": "https://static.inaturalist.org/photos/fail/original.jpg",
                    "width": 1600,
                    "height": 1200,
                    "asset_type": "photo",
                }
            ],
        }
        wikimedia_candidates = {
            "Panthera leo": [
                {
                    "provider": "wikimedia",
                    "provider_record_id": "wiki-ok",
                    "matched_name": "Panthera leo",
                    "matched_rank": "species",
                    "license_code": "cc-by-sa",
                    "license_url": "https://creativecommons.org/licenses/by-sa/4.0/",
                    "attribution": "Working Photographer",
                    "source_page_url": "https://commons.wikimedia.org/wiki/File:Lion.jpg",
                    "media_url": "https://upload.wikimedia.org/lion.jpg",
                    "width": 1400,
                    "height": 1000,
                    "asset_type": "photo",
                }
            ],
        }

        def fake_build_providers(args, sources, session=None):
            providers = {
                "inaturalist": DummyProvider(inaturalist_candidates),
                "wikimedia": DummyProvider(wikimedia_candidates),
            }
            return DummySession(), None, providers

        def fake_download_media(
            session,
            media_url,
            destination_path,
            cache_path=None,
            max_download_bytes=None,
            provider=None,
            **kwargs,
        ):
            if media_url.endswith("/fail/original.jpg"):
                raise requests.ConnectionError("transient failure")
            write_valid_test_media(destination_path)
            return "downloaded"

        monkeypatch.setattr("nwkit.image.build_providers", fake_build_providers)
        monkeypatch.setattr("nwkit.image.download_media", fake_download_media)

        args = make_image_args(
            infile=str(tree_path),
            out_dir=str(out_dir),
            style="photo",
            source="inaturalist,wikimedia",
        )
        image_main(args)

        manifest_rows = read_tsv(out_dir / "manifest.tsv")
        unmatched_rows = read_tsv(out_dir / "unmatched.tsv")

        assert len(manifest_rows) == 1
        assert manifest_rows[0]["provider"] == "wikimedia"
        assert unmatched_rows == []

    def test_materialize_value_error_falls_back_to_next_candidate(self):
        candidates = [
            {
                "provider": "inaturalist",
                "provider_record_id": "invalid",
                "matched_name": "Panthera leo",
                "matched_rank": "species",
                "license_code": "cc-by",
                "license_url": "https://creativecommons.org/licenses/by/4.0/",
                "attribution": "Invalid Photographer",
                "source_page_url": "https://example.org/invalid",
                "media_url": "https://example.org/invalid.jpg",
                "asset_type": "photo",
                "score": 2.0,
            },
            {
                "provider": "wikimedia",
                "provider_record_id": "valid",
                "matched_name": "Panthera leo",
                "matched_rank": "species",
                "license_code": "cc-by",
                "license_url": "https://creativecommons.org/licenses/by/4.0/",
                "attribution": "Valid Photographer",
                "source_page_url": "https://example.org/valid",
                "media_url": "https://example.org/valid.jpg",
                "asset_type": "photo",
                "score": 1.0,
            },
        ]

        class Materializer:
            def __init__(self):
                self.calls = list()

            def materialize(self, species_name, candidate):
                self.calls.append(candidate["provider_record_id"])
                if candidate["provider_record_id"] == "invalid":
                    raise ValueError("decoded image is invalid")
                return {
                    "local_path": "images/valid.jpg",
                    "download_status": "downloaded",
                }

        materializer = Materializer()
        manifest_rows, selected_assets, unmatched_rows = (
            image_module.process_species_assets(
                species_name="Panthera leo",
                leaf_names=["Panthera_leo_A"],
                candidates=candidates,
                provider_errors=[],
                args=make_image_args(max_per_species=1),
                materializer=materializer,
            )
        )

        assert materializer.calls == ["invalid", "valid"]
        assert [row["provider"] for row in manifest_rows] == ["wikimedia"]
        assert selected_assets == manifest_rows
        assert unmatched_rows == []

    def test_image_main_reports_download_error_when_all_candidates_fail(
        self, monkeypatch, tmp_path
    ):
        tree_path = tmp_path / "tree.nwk"
        tree_path.write_text("(Panthera_leo_A);")
        out_dir = tmp_path / "out"

        inaturalist_candidates = {
            "Panthera leo": [
                {
                    "provider": "inaturalist",
                    "provider_record_id": "inat-fail",
                    "matched_name": "Panthera leo",
                    "matched_rank": "species",
                    "license_code": "cc-by",
                    "license_url": "https://creativecommons.org/licenses/by/4.0/",
                    "attribution": "Failing Photographer",
                    "source_page_url": "https://www.inaturalist.org/observations/fail",
                    "media_url": "https://static.inaturalist.org/photos/fail/original.jpg",
                    "width": 1600,
                    "height": 1200,
                    "asset_type": "photo",
                }
            ],
        }

        def fake_build_providers(args, sources, session=None):
            providers = {
                "inaturalist": DummyProvider(inaturalist_candidates),
            }
            return DummySession(), None, providers

        def fake_download_media(
            session,
            media_url,
            destination_path,
            cache_path=None,
            max_download_bytes=None,
            provider=None,
            **kwargs,
        ):
            raise requests.ConnectionError("transient failure")

        monkeypatch.setattr("nwkit.image.build_providers", fake_build_providers)
        monkeypatch.setattr("nwkit.image.download_media", fake_download_media)

        args = make_image_args(
            infile=str(tree_path),
            out_dir=str(out_dir),
            style="photo",
            source="inaturalist",
        )
        image_main(args)

        manifest_rows = read_tsv(out_dir / "manifest.tsv")
        unmatched_rows = read_tsv(out_dir / "unmatched.tsv")

        assert manifest_rows == []
        assert unmatched_rows[0]["reason"] == "download_error"
        assert (
            "inaturalist download failed for Panthera leo"
            in unmatched_rows[0]["details"]
        )

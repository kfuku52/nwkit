import csv
import sqlite3
from pathlib import Path
from types import SimpleNamespace

import pytest
import requests
from PIL import Image

from nwkit import image, image_workflow
from tests.image_test_support import DummySession, make_image_args, read_tsv


def workflow_args(tmp_path, leaves, **kwargs):
    tree = tmp_path / "tree.nwk"
    tree.write_text("(" + ",".join(leaves) + ");\n")
    return make_image_args(
        infile=str(tree),
        out_dir=str(tmp_path / "out"),
        download_dir=str(tmp_path / "cache"),
        **kwargs,
    )


def media_candidate(provider, query, key):
    host = "upload.wikimedia.org" if provider == "wikimedia" else "example.org"
    return {
        "provider": provider,
        "provider_record_id": key,
        "matched_name": query,
        "matched_rank": image_workflow.declared_rank(query),
        "license_code": "cc-by",
        "license_url": "https://creativecommons.org/licenses/by/4.0/",
        "attribution": "Example photographer",
        "source_page_url": "https://example.org/photo",
        "media_url": f"https://{host}/{key}.jpg",
        "asset_type": "photo",
    }


def mock_providers(monkeypatch, catalog, *, taxonomy=None):
    calls = []

    class Provider:
        def __init__(self, source):
            self.source = source

        def fetch_candidates(self, species_name, fallback_rank):
            calls.append((self.source, species_name, fallback_rank))
            return [
                media_candidate(self.source, species_name, key)
                for key in catalog.get((self.source, species_name), ())
            ]

    def build(args, sources):
        return (
            DummySession(),
            taxonomy,
            {source: Provider(source) for source in sources},
        )

    monkeypatch.setattr(image, "build_providers", build)
    monkeypatch.setenv("NWKIT_IMAGE_LOOKUP_WORKERS", "1")
    return calls


def download_stub(monkeypatch, *, broken=()):
    calls = []

    def download(*, media_url, destination_path, **kwargs):
        key = Path(media_url).stem
        calls.append(key)
        if key in broken:
            raise requests.ConnectionError("image host unreachable")
        Image.new("RGB", (8, 4), "green").save(destination_path)
        return {"status": "downloaded", "destination_path": destination_path}

    monkeypatch.setattr(image, "download_media", download)
    return calls


def verified_table(tmp_path, rows):
    path = tmp_path / "verified.tsv"
    fields = (
        "leaf_name",
        "local_path",
        "media_url",
        "source_page_url",
        "attribution",
        "license_code",
        "license_url",
        "matched_name",
        "matched_rank",
    )
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fields, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    "source_page_url": "https://example.org/paper",
                    "attribution": "Author (cropped)",
                    "license_code": "cc-by",
                    **row,
                }
            )
    return str(path)


def test_verified_overrides_are_per_leaf_and_use_the_image_pipeline(
    tmp_path, monkeypatch
):
    source = tmp_path / "figure.png"
    Image.new("RGB", (20, 10), "blue").save(source)
    path = verified_table(
        tmp_path,
        [
            {"leaf_name": "Felis_catus_sample1", "local_path": "figure.png"},
            {
                "leaf_name": "verifiedLabel",
                "media_url": "https://example.org/remote.jpg",
                "matched_name": "Oryza sativa cultivar Nipponbare",
            },
        ],
    )
    args = workflow_args(
        tmp_path,
        ["Felis_catus_sample1", "Felis_catus_sample2", "verifiedLabel"],
        source="wikimedia",
        verified_images_tsv=path,
        output_format="jpg",
        max_edge=6,
        resume=True,
        first_available=True,
    )
    calls = mock_providers(monkeypatch, {("wikimedia", "Felis catus"): ["automatic"]})
    downloads = download_stub(monkeypatch)

    image.image_main(args)

    out = Path(args.out_dir)
    rows = {row["leaf_name"]: row for row in read_tsv(out / "manifest.tsv")}
    first = rows["Felis_catus_sample1"]
    assert first["provider"] == "verified" and first["media_url"] == ""
    assert first["selection_reason"] == "exact_species_match;verified_image"
    assert first["status"] == "copied"
    with Image.open(out / first["local_path"]) as saved:
        assert saved.format == "JPEG" and saved.size == (6, 3)
    assert rows["Felis_catus_sample2"]["provider"] == "wikimedia"
    assert rows["verifiedLabel"]["matched_rank"] == "cultivar"
    assert rows["verifiedLabel"]["source_page_url"] == "https://example.org/paper"
    assert downloads == ["automatic", "remote"]
    assert calls == [("wikimedia", "Felis catus", "none")]
    assert not read_tsv(out / "unmatched.tsv")
    assert r"Author \(cropped\)" in (out / "ATTRIBUTION.md").read_text()

    image.image_main(args)
    assert {row["status"] for row in read_tsv(out / "manifest.tsv")} == {"resumed"}
    Image.new("RGB", (20, 10), "red").save(source)
    image.image_main(args)
    rows = {row["leaf_name"]: row for row in read_tsv(out / "manifest.tsv")}
    assert rows["Felis_catus_sample1"]["status"] == "copied"
    assert rows["Felis_catus_sample1"]["local_path"] != first["local_path"]
    assert rows["Felis_catus_sample2"]["status"] == "resumed"
    assert downloads == ["automatic", "remote"]


@pytest.mark.parametrize(
    ("row", "error"),
    [
        ({"media_url": "http://example.org/asset.jpg"}, "HTTPS"),
        ({"media_url": "https://127.0.0.1/a.jpg"}, "non-public"),
        (
            {"media_url": "https://example.org/a.jpg", "local_path": "also.jpg"},
            "exactly one",
        ),
        (
            {"media_url": "https://example.org/a.jpg", "license_code": "cc-by-nc"},
            "disallowed",
        ),
        (
            {
                "media_url": "https://example.org/a.jpg",
                "license_url": "https://creativecommons.org/licenses/by-nc/4.0/",
            },
            "conflicting",
        ),
    ],
)
def test_verified_metadata_validates_before_search(tmp_path, row, error):
    path = verified_table(tmp_path, [{"leaf_name": "Felis_catus", **row}])
    args = workflow_args(
        tmp_path,
        ["Felis_catus", "Panthera_leo"],
        verified_images_tsv=path,
        license_max="cc-by",
    )
    with pytest.raises((ValueError, image.MediaDownloadError), match=error):
        image.image_main(args)
    assert not (Path(args.out_dir) / "manifest.tsv").exists()


def test_first_available_falls_through_on_unusable_media(tmp_path, monkeypatch):
    args = workflow_args(
        tmp_path,
        ["Felis_catus", "Panthera_leo"],
        source="wikimedia,gbif,ncbi",
        first_available=True,
    )
    calls = mock_providers(
        monkeypatch,
        {
            ("wikimedia", "Felis catus"): ["dead"],
            ("gbif", "Felis catus"): ["safe"],
            ("wikimedia", "Panthera leo"): ["good"],
        },
    )
    downloads = download_stub(monkeypatch, broken={"dead"})
    image.image_main(args)
    assert [(source, name) for source, name, _ in calls] == [
        ("wikimedia", "Felis catus"),
        ("gbif", "Felis catus"),
        ("wikimedia", "Panthera leo"),
    ]
    assert downloads == ["dead", "safe", "good"]
    assert not read_tsv(Path(args.out_dir) / "unmatched.tsv")


def test_resume_keeps_progress_after_interruption_and_revalidates_images(
    tmp_path, monkeypatch
):
    args = workflow_args(
        tmp_path,
        ["Felis_catus", "Panthera_leo"],
        source="wikimedia",
        first_available=True,
        resume=True,
    )
    mock_providers(
        monkeypatch,
        {
            ("wikimedia", "Felis catus"): ["cat"],
            ("wikimedia", "Panthera leo"): ["lion"],
        },
    )
    downloads = download_stub(monkeypatch)
    base_download = image.download_media

    def interrupted(**kwargs):
        if kwargs["media_url"].endswith("/lion.jpg"):
            raise KeyboardInterrupt
        return base_download(**kwargs)

    monkeypatch.setattr(image, "download_media", interrupted)
    with pytest.raises(KeyboardInterrupt):
        image.image_main(args)
    out = Path(args.out_dir)
    assert not (out / "manifest.tsv").exists()
    assert len(list((out / ".nwkit-cache" / "image-resume").glob("*.json"))) == 1

    monkeypatch.setattr(image, "download_media", base_download)
    image.image_main(args)
    rows = read_tsv(out / "manifest.tsv")
    assert [row["status"] for row in rows] == ["resumed", "downloaded"]
    assert downloads == ["cat", "lion"]

    (out / rows[0]["local_path"]).write_bytes(b"corrupted")
    image.image_main(args)
    assert downloads == ["cat", "lion", "cat"]
    assert [row["status"] for row in read_tsv(out / "manifest.tsv")] == [
        "downloaded",
        "resumed",
    ]

    args.max_edge = 4
    image.image_main(args)
    assert downloads == ["cat", "lion", "cat", "cat", "lion"]
    for row in read_tsv(out / "manifest.tsv"):
        with Image.open(out / row["local_path"]) as saved:
            assert max(saved.size) == 4


def test_synonyms_use_same_unique_taxid_and_preserve_requested_name(
    tmp_path, monkeypatch
):
    db = sqlite3.connect(":memory:")
    db.executescript("""
        CREATE TABLE species (taxid INT, spname TEXT COLLATE NOCASE);
        CREATE TABLE synonym (taxid INT, spname TEXT COLLATE NOCASE);
        INSERT INTO species VALUES (1, 'Oryza glumipatula'), (2, 'Oryza ambiguous');
        INSERT INTO synonym VALUES (1, 'Oryza ambiguous'), (1, 'Oryza glumaepatula');
    """)
    assert image_workflow.unique_taxon_aliases(db, "Oryza ambiguous") == []
    taxonomy = SimpleNamespace(
        aliases=lambda name: image_workflow.unique_taxon_aliases(db, name)
    )
    args = workflow_args(
        tmp_path,
        ["Oryza_glumipatula", "Oryza_ambiguous"],
        source="wikimedia",
        synonyms=True,
        first_available=True,
    )
    calls = mock_providers(
        monkeypatch, {("wikimedia", "Oryza glumaepatula"): ["rice"]}, taxonomy=taxonomy
    )
    download_stub(monkeypatch)
    image.image_main(args)
    rows = read_tsv(Path(args.out_dir) / "manifest.tsv")
    assert len(rows) == 1
    assert rows[0]["species_name"] == "Oryza glumipatula"
    assert rows[0]["matched_name"] == "Oryza glumaepatula"
    assert rows[0]["selection_reason"].startswith("taxonomic_synonym_match;")
    assert [name for _, name, _ in calls] == [
        "Oryza ambiguous",
        "Oryza glumipatula",
        "Oryza glumaepatula",
    ]
    db.close()


def test_exact_infraspecific_queries_precede_explicitly_marked_species_images(
    tmp_path, monkeypatch
):
    args = workflow_args(
        tmp_path,
        ["Felis_silvestris_subsp._lybica", "Vigna_radiata_var._radiata"],
        source="wikimedia,gbif",
        first_available=True,
    )
    calls = mock_providers(
        monkeypatch,
        {
            ("gbif", "Felis silvestris subsp. lybica"): ["lybica"],
            ("wikimedia", "Vigna radiata"): ["mung"],
        },
    )
    download_stub(monkeypatch)
    image.image_main(args)
    rows = read_tsv(Path(args.out_dir) / "manifest.tsv")
    assert rows[0]["matched_rank"] == "subspecies"
    assert rows[0]["selection_reason"].startswith("exact_infraspecific_match;")
    assert (
        rows[1]["species_name"],
        rows[1]["matched_name"],
        rows[1]["matched_rank"],
    ) == (
        "Vigna radiata var. radiata",
        "Vigna radiata",
        "species",
    )
    assert "infraspecific_identity_unconfirmed" in rows[1]["selection_reason"]
    assert [name for _, name, _ in calls] == [
        "Felis silvestris subsp. lybica",
        "Felis silvestris subsp. lybica",
        "Vigna radiata var. radiata",
        "Vigna radiata var. radiata",
        "Vigna radiata",
    ]


def test_provider_pacer_shares_interval_and_keeps_other_providers_independent(
    monkeypatch,
):
    intervals = image_workflow.parse_provider_intervals(
        ["inaturalist=0.25", "wikimedia=0.4"], {"inaturalist", "wikimedia"}
    )
    assert intervals == {"inaturalist": 1.0, "wikimedia": 0.4}
    now = 100.0
    sleeps = []

    def sleep(delay):
        nonlocal now
        now += delay
        sleeps.append(delay)

    monkeypatch.setattr(image_workflow.time, "monotonic", lambda: now)
    monkeypatch.setattr(image_workflow.time, "sleep", sleep)
    pacer = image_workflow.RequestPacer(intervals)
    for provider in ("inaturalist", "wikimedia", "wikimedia", "inaturalist"):
        pacer.wait(provider)
    assert sleeps == pytest.approx([0.4, 0.6])


def test_wikimedia_depicts_requires_unambiguous_scientific_name_and_ranks_first(
    monkeypatch,
):
    def request(session, provider, method, url, *, params, **kwargs):
        action = params["action"]
        if action == "wbsearchentities":
            return {"search": [{"id": "Q1"}, {"id": "Q2"}]}
        if action == "wbgetentities":

            def entity(name):
                return {
                    "claims": {"P225": [{"mainsnak": {"datavalue": {"value": name}}}]}
                }

            return {
                "entities": {"Q1": entity("Felis catus"), "Q2": entity("Felis catusa")}
            }
        structured = "haswbstatement:P180=Q1" in params["gsrsearch"]
        pageid = 1 if structured else 2
        return {
            "query": {
                "pages": {
                    str(pageid): {
                        "pageid": pageid,
                        "title": "File:portrait.jpg"
                        if structured
                        else "File:Felis catus.jpg",
                        "imageinfo": [
                            {
                                "url": f"https://upload.wikimedia.org/{pageid}.jpg",
                                "extmetadata": {
                                    "LicenseShortName": {
                                        "value": "CC BY 4.0"
                                        if structured
                                        else "Public domain"
                                    }
                                },
                            }
                        ],
                    }
                }
            }
        }

    monkeypatch.setattr(image, "provider_json_request", request)
    provider = image.WikimediaProvider(
        DummySession(), args=make_image_args(wikimedia_depicts=True)
    )
    candidates = provider.fetch_candidates("Felis catus")
    assert len(candidates) == 2
    assert candidates[0]["subject_match"] and not candidates[1]["subject_match"]
    assert image.candidate_score(candidates[0]) > image.candidate_score(candidates[1])
    assert image.describe_candidate_selection(candidates[0]).endswith(
        ";wikimedia_depicts"
    )


def test_wikimedia_text_match_respects_name_boundaries():
    page = {"title": "File:Felis_catusa.jpg"}
    assert not image.wikimedia_page_mentions_query(page, "Felis catus")
    page["title"] = "File:Felis_catus_1.jpg"
    assert image.wikimedia_page_mentions_query(page, "Felis catus")


def test_validated_local_image_is_in_the_cli_audit(tmp_path, monkeypatch):
    import json

    from nwkit.cli import main

    source = tmp_path / "figure.png"
    Image.new("RGB", (20, 10), "blue").save(source)
    table = verified_table(
        tmp_path, [{"leaf_name": "Felis_catus", "local_path": "figure.png"}]
    )
    tree = tmp_path / "tree.nwk"
    tree.write_text("(Felis_catus);\n")
    audit = tmp_path / "audit.jsonl"
    out = tmp_path / "out"
    mock_providers(monkeypatch, {})

    main(
        [
            "image",
            "-i",
            str(tree),
            "--out-dir",
            str(out),
            "--audit",
            str(audit),
            "--verified-images-tsv",
            table,
            "--resume",
            "yes",
            "--first-available",
            "yes",
        ]
    )
    record = json.loads(audit.read_text())
    inputs = {entry["path"]: entry["argument"] for entry in record["inputs"]}
    assert inputs[str(source.resolve())] == "verified_images_tsv:asset"
    assert inputs[str(Path(table).resolve())] == "verified_images_tsv"


def test_partial_provider_selection_fills_remaining_slots(tmp_path, monkeypatch):
    args = workflow_args(
        tmp_path,
        ["Felis_catus", "Panthera_leo"],
        max_per_species=2,
        source="wikimedia,gbif,ncbi",
        first_available=True,
    )
    calls = mock_providers(
        monkeypatch,
        {
            ("wikimedia", "Felis catus"): ["cat1"],
            ("gbif", "Felis catus"): ["cat2", "cat3"],
            ("wikimedia", "Panthera leo"): ["lion1", "lion2"],
        },
    )
    downloads = download_stub(monkeypatch)
    image.image_main(args)
    rows = read_tsv(Path(args.out_dir) / "manifest.tsv")
    assert [row["provider"] for row in rows] == [
        "wikimedia",
        "gbif",
        "wikimedia",
        "wikimedia",
    ]
    assert len(downloads) == 4
    assert [call[0] for call in calls] == ["wikimedia", "gbif", "wikimedia"]

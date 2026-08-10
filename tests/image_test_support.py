import csv
import os
from argparse import Namespace

import requests


def make_image_args(**kwargs):
    defaults = {
        "infile": "-",
        "format": "auto",
        "quoted_node_names": True,
        "download_dir": "auto",
        "out_dir": None,
        "style": "auto",
        "source": None,
        "license_max": "cc-by-nc-sa",
        "allow_nd": False,
        "fallback_rank": "none",
        "max_per_species": 1,
        "max_download_bytes": 104857600,
        "query_cache_max_age_hours": 168.0,
        "refresh_cache": False,
        "species_name_tsv": None,
        "manifest_out": None,
        "attribution_out": None,
        "fail_on_missing": False,
        "output_format": "original",
        "max_edge": None,
        "canvas": "none",
        "background": "white",
        "trim": "off",
        "trim_shape": "bbox",
        "species_regex": r"^([^_]+_[^_]+)(?:_|$)",
    }
    defaults.update(kwargs)
    return Namespace(**defaults)


class DummySession:
    def close(self):
        return None


class JSONResponse:
    def __init__(self, payload, status_code=200, url="https://example.org"):
        self.payload = payload
        self.status_code = status_code
        self.url = url
        self.headers = {"Content-Type": "application/json"}

    def json(self):
        return self.payload

    def raise_for_status(self):
        if self.status_code >= 400:
            raise requests.HTTPError("HTTP {}".format(self.status_code))

    def close(self):
        return None


class DummyProvider:
    def __init__(self, candidates_by_species):
        self.candidates_by_species = candidates_by_species

    def fetch_candidates(self, species_name, fallback_rank="none"):
        return [
            dict(candidate)
            for candidate in self.candidates_by_species.get(species_name, [])
        ]


def read_tsv(path):
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_valid_test_media(path):
    if str(path).lower().endswith(".svg"):
        with open(path, "wb") as handle:
            handle.write(b'<svg xmlns="http://www.w3.org/2000/svg"></svg>')
        return
    from PIL import Image

    extension = os.path.splitext(str(path))[1].lower()
    image_format = {
        ".gif": "GIF",
        ".jpg": "JPEG",
        ".jpeg": "JPEG",
        ".png": "PNG",
        ".tif": "TIFF",
        ".tiff": "TIFF",
        ".webp": "WEBP",
    }.get(extension, "PNG")
    Image.new("RGB", (2, 2), "white").save(path, format=image_format)

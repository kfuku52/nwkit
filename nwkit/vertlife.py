"""Verified VertLife reference trees and species matching."""

import hashlib
import os
import stat
from dataclasses import dataclass
from pathlib import Path

import requests
from ete4.parser.newick import PARSERS

from nwkit.file_paths import validate_outputs_do_not_replace_inputs
from nwkit.output_transaction import output_lock_path, output_transaction
from nwkit.rooting_state import copy_rooting_info, require_rooted
from nwkit.species_parser import get_species_parser
from nwkit.util import (
    TREE_FORMAT_PROP,
    acquire_exclusive_lock,
    read_tree,
    resolve_download_dir,
    validate_unique_named_leaves,
)

_MAX_BYTES = 2 * 1024 * 1024


@dataclass(frozen=True)
class VertLifeDataset:
    id: str
    url: str
    sha256: str
    doi: str


_AMPHIBIANS = VertLifeDataset(
    "amphibians-2018-consensus",
    "https://storage.googleapis.com/data.vertlife.org/trees/download/"
    "amph_shl_new_Consensus_7238.tre",
    "5468a97e294f1ea8693bf9abbd85c5b81485de91bd95b835a1da42d5ec16bdf2",
    "10.5061/dryad.cc3n6j5",
)
DATASETS = {_AMPHIBIANS.id: _AMPHIBIANS}


def _dataset(args):
    dataset_id = getattr(args, "vertlife_dataset", None)
    if dataset_id not in DATASETS:
        raise ValueError(
            "Specify '--vertlife-dataset' using one of: {}.".format(", ".join(DATASETS))
        )
    return DATASETS[dataset_id]


def vertlife_input_path(args):
    dataset = _dataset(args)
    local = getattr(args, "vertlife_file", None)
    if local in ("", "-"):
        raise ValueError("'--vertlife-file' requires a file path.")
    if local:
        return Path(local)
    cache = resolve_download_dir(args) or Path.home() / ".cache" / "nwkit"
    return Path(cache) / "vertlife" / dataset.id / dataset.url.rsplit("/", 1)[-1]


def _download_lock_path(path):
    return path.with_name(f".{path.name}.download.lock")


def _cache_lock_paths(path):
    return [
        ("VertLife download lock", _download_lock_path(path)),
        ("VertLife cache output lock", output_lock_path(path)),
    ]


def vertlife_cache_lock_paths(args):
    return _cache_lock_paths(vertlife_input_path(args))


def _check_source(data, dataset):
    if len(data) > _MAX_BYTES:
        raise ValueError("VertLife source exceeds the size limit.")
    if hashlib.sha256(data).hexdigest() != dataset.sha256:
        raise ValueError(f"VertLife checksum verification failed for '{dataset.id}'.")


def _download(dataset):
    data = bytearray()
    try:
        with requests.get(dataset.url, stream=True, timeout=30) as response:
            response.raise_for_status()
            for chunk in response.iter_content(chunk_size=64 * 1024):
                data.extend(chunk)
                if len(data) > _MAX_BYTES:
                    raise ValueError("VertLife download exceeds the size limit.")
    except requests.RequestException as exc:
        raise ValueError(
            f"Failed to download VertLife dataset '{dataset.id}'."
        ) from exc
    return data


def _parse_reference(data, dataset):
    _check_source(data, dataset)
    tree = read_tree(data.decode("utf-8-sig"), "auto-strict", True, quiet=True)
    require_rooted(tree, "The VertLife reference must be rooted.")
    validate_unique_named_leaves(tree, option_name="--vertlife-dataset")
    if any(
        node.dist is None or node.dist < 0
        for node in tree.traverse()
        if not node.is_root
    ):
        raise ValueError("The VertLife reference requires non-negative branch lengths.")
    return tree


def _read_reference(path, dataset):
    def open_nonblocking(name, flags):
        return os.open(name, flags | getattr(os, "O_NONBLOCK", 0))

    with open(path, "rb", opener=open_nonblocking) as handle:
        if not stat.S_ISREG(os.fstat(handle.fileno()).st_mode):
            raise ValueError(f"VertLife source must be a regular file: '{path}'.")
        return _parse_reference(handle.read(_MAX_BYTES + 1), dataset)


def _load_cached_reference(path, dataset):
    try:
        return _read_reference(path, dataset)
    except FileNotFoundError:
        pass
    with acquire_exclusive_lock(
        _download_lock_path(path), lock_label="VertLife download", poll_seconds=0.05
    ):
        if not os.path.lexists(path):
            data = _download(dataset)
            tree = _parse_reference(data, dataset)
            with output_transaction([path], follow_symlinks=False) as staged:
                staged.write_bytes(path, data)
            return tree
    return _read_reference(path, dataset)


def load_vertlife(args):
    dataset = _dataset(args)
    path = vertlife_input_path(args)
    local = getattr(args, "vertlife_file", None) is not None
    inputs = [("VertLife source/cache", path)]
    if not local:
        inputs.extend(_cache_lock_paths(path))
    validate_outputs_do_not_replace_inputs(
        inputs,
        [
            (f"--{key}", getattr(args, key, None))
            for key in ("outfile", "report", "audit")
        ],
    )
    if local:
        tree = _read_reference(path, dataset)
    else:
        tree = _load_cached_reference(path, dataset)
        args._nwkit_vertlife_cache = str(path)
    return dataset, tree


def _name_key(name):
    return " ".join(name.replace("_", " ").split())


def match_species(reference, labels, args):
    """Return (input label, query, reference label, status) for each label."""
    by_name = {}
    for name in reference.leaf_names():
        key = _name_key(name)
        if key in by_name:
            raise ValueError(
                f"VertLife contains ambiguous normalized names for '{key}'."
            )
        by_name[key] = name
    parser = get_species_parser(args=args)
    matches = []
    for label in labels:
        override = parser.map_overrides.get(label)
        query = _name_key(label)
        status = "exact"
        if override is not None:
            query = _name_key(override.taxonomy_query or override.species_label)
            status = "mapped"
        elif query not in by_name:
            query = parser.parse(label).taxonomy_query or query
            status = "parsed"
        match = by_name.get(_name_key(query))
        matches.append(
            (label, query, match, status if match is not None else "unmatched")
        )
    return matches


def subset_reference(tree, names):
    """Consume a reference tree, retaining the selected crown and its distances."""
    retained = {}
    for node in tree.traverse(strategy="postorder"):
        if node.is_leaf:
            if node.name in names:
                retained[node] = node
            continue
        children = []
        for child in node.children:
            child.up = None
            representative = retained.pop(child, None)
            if representative is not None:
                children.append(representative)
        if len(children) == 1:
            representative = children[0]
            representative.dist += node.dist or 0.0
            retained[node] = representative
            node.children = []
        else:
            node.children = children
            if children:
                retained[node] = node
    if tree not in retained:
        raise ValueError("No requested species were found in the VertLife reference.")
    subset = retained[tree]
    subset.dist = None
    copy_rooting_info(tree, subset)
    return subset


def precise_tree_parser(format):
    """Keep numeric precision through the shared ETE serializer."""
    parser = PARSERS[int(format)]

    def write_number(value):
        return repr(float(value))

    return {
        kind: [
            {**field, "write": write_number}
            if field["pname"] in ("dist", "support")
            else field
            for field in fields
        ]
        for kind, fields in parser.items()
    }


def precise_tree_copy_format(tree, *, include_root_name=True):
    """Keep the input format when it can represent the internal node names."""
    format = tree.props.get(TREE_FORMAT_PROP, 0)
    parser = precise_tree_parser(format)
    if parser["internal"][0]["pname"] != "name" and any(
        node.name
        for node in tree.traverse()
        if not node.is_leaf and (include_root_name or not node.is_root)
    ):
        return 1, precise_tree_parser(1)
    return format, parser

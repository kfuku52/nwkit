import json as json_lib

from ete4 import Tree

import nwkit.root as root_mod
from nwkit.clade_mapping import canonical_split


def annotated_reroot_tree():
    tree = Tree(
        '(((A:2,B:3):5,C:7):11,(D:13,(E:17,F:19):23):29);',
        parser=1,
    )
    all_taxa = frozenset(tree.leaf_names())
    expected_by_split = dict()
    for node in tree.traverse():
        if node.is_root:
            node.name = 'ORIGINAL_ROOT'
            node.support = 0.99
            node.props['root_tag'] = 'root_value'
        elif node.is_leaf:
            node.support = 0.5
            node.props['tip_tag'] = 'tip_{}'.format(node.name)
        else:
            side = frozenset(node.leaf_names())
            split = canonical_split(side, all_taxa - side)
            if split not in expected_by_split:
                index = len(expected_by_split) + 1
                expected_by_split[split] = (
                    'EDGE_{}'.format(index),
                    float(10 + index),
                    'edge_value_{}'.format(index),
                )
            node.name, node.support, node.props['edge_tag'] = expected_by_split[split]
    return tree, expected_by_split


def assert_reroot_annotations(tree, expected_by_split):
    all_taxa = frozenset(tree.leaf_names())
    observed_by_split = dict()
    for node in tree.traverse():
        if node.is_root or node.is_leaf:
            continue
        side = frozenset(node.leaf_names())
        split = canonical_split(side, all_taxa - side)
        observed_by_split.setdefault(split, set()).add(
            (node.name, node.support, node.props.get('edge_tag'))
        )
    for split, expected in expected_by_split.items():
        assert observed_by_split[split] == {expected}
    assert tree.name == 'ORIGINAL_ROOT'
    assert tree.support == 0.99
    assert tree.props['root_tag'] == 'root_value'
    for leaf in tree.leaves():
        assert leaf.support == 0.5
        assert leaf.props['tip_tag'] == 'tip_{}'.format(leaf.name)


def scaled_branch_profile(tree, scale=1.0):
    all_taxa = frozenset(tree.leaf_names())
    profile = dict()
    for node in tree.traverse():
        if node.is_root:
            continue
        side = frozenset(node.leaf_names())
        split = canonical_split(side, all_taxa - side)
        profile.setdefault(split, list()).append(float(node.dist) / scale)
    return {
        split: sorted(lengths)
        for split, lengths in profile.items()
    }


def install_fake_ncbi(
    monkeypatch,
    name_to_taxid,
    lineage_by_taxid,
    taxid_to_name=None,
    rank_by_taxid=None,
):
    if taxid_to_name is None:
        taxid_to_name = {
            int(taxid): str(taxid)
            for taxid in lineage_by_taxid.keys()
        }
    if rank_by_taxid is None:
        rank_by_taxid = dict()

    class FakeNCBI:
        def __init__(self, *args, **kwargs):
            self.db = None

        def get_name_translator(self, names):
            translated = dict()
            for name in names:
                if name in name_to_taxid:
                    translated[name] = [name_to_taxid[name]]
            return translated

        def get_lineage(self, taxid):
            return list(lineage_by_taxid[int(taxid)])

        def get_taxid_translator(self, taxids):
            return {
                int(taxid): taxid_to_name[int(taxid)]
                for taxid in taxids
                if int(taxid) in taxid_to_name
            }

        def get_rank(self, taxids):
            return {
                int(taxid): rank_by_taxid.get(int(taxid), 'no rank')
                for taxid in taxids
            }

    monkeypatch.setattr(
        root_mod,
        'get_ete_ncbitaxa',
        lambda args=None: FakeNCBI(),
    )


def install_fake_timetree(
    monkeypatch,
    upload_html,
    newick_text,
    upload_status=200,
    newick_status=200,
):
    calls = list()

    class FakeResponse:
        def __init__(self, text, status_code):
            self.text = text
            self.status_code = status_code
            self.headers = {}
            self.encoding = 'utf-8'

        def iter_content(self, chunk_size):
            payload = self.text.encode(self.encoding)
            for offset in range(0, len(payload), chunk_size):
                yield payload[offset:offset + chunk_size]

        def close(self):
            pass

    class FakeSession:
        def post(self, url, files=None, data=None, timeout=None, stream=False):
            calls.append((url, stream))
            if url.endswith('/ajax/prune/load_names/'):
                return FakeResponse(upload_html, upload_status)
            if url.endswith('/ajax/newick/prunetree/download'):
                return FakeResponse(newick_text, newick_status)
            raise AssertionError('Unexpected TimeTree URL: {}'.format(url))

    monkeypatch.setattr(root_mod.requests, 'Session', FakeSession)
    return calls


def install_fake_opentree(
    monkeypatch,
    tnrs_json,
    induced_subtree_json,
    tnrs_status=200,
    induced_subtree_status=200,
):
    calls = list()

    class FakeResponse:
        def __init__(self, text, status_code, json_data=None):
            self.text = text
            self.status_code = status_code
            self._json_data = json_data
            self.headers = {}
            self.encoding = 'utf-8'

        def json(self):
            if self._json_data is None:
                raise AssertionError('JSON payload was not configured')
            return self._json_data

        def iter_content(self, chunk_size):
            payload = self.text.encode(self.encoding)
            for offset in range(0, len(payload), chunk_size):
                yield payload[offset:offset + chunk_size]

        def close(self):
            pass

    class FakeSession:
        def post(self, url, json=None, timeout=None, stream=False):
            calls.append((url, stream))
            if url.endswith('/v3/tnrs/match_names'):
                return FakeResponse(
                    json_lib.dumps(tnrs_json),
                    tnrs_status,
                    json_data=tnrs_json,
                )
            if url.endswith('/v3/tree_of_life/induced_subtree'):
                return FakeResponse(
                    json_lib.dumps(induced_subtree_json),
                    induced_subtree_status,
                    json_data=induced_subtree_json,
                )
            raise AssertionError('Unexpected OpenTree URL: {}'.format(url))

    monkeypatch.setattr(root_mod.requests, 'Session', FakeSession)
    return calls

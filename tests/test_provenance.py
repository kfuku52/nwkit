import hashlib
import io
import json
import sys
from argparse import Namespace

import pytest

from nwkit.cli import main
from nwkit.provenance import _TeeTextWriter, _argument_dict


def _sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_argument_dict_omits_private_runtime_caches():
    args = Namespace(
        command='constrain',
        species_parser='legacy',
        _nwkit_species_parser_cache=object(),
        _nwkit_species_parser_cache_key=('legacy', 'pattern', None),
        handler=object(),
    )
    assert _argument_dict(args) == {
        'command': 'constrain',
        'species_parser': 'legacy',
    }


def test_global_audit_records_arguments_hashes_and_messages(tmp_path):
    infile = tmp_path / 'input.nwk'
    outfile = tmp_path / 'output.nwk'
    audit = tmp_path / 'audit.jsonl'
    infile.write_text('((A:1,B:1):1,C:1);')
    main([
        'label',
        '-i', str(infile),
        '-o', str(outfile),
        '--target', 'intnode',
        '--audit', str(audit),
    ])
    record = json.loads(audit.read_text().strip())
    assert record['schema'] == 'nwkit-audit-v1'
    assert record['status'] == 'ok'
    assert record['command'] == 'label'
    assert record['inputs'][0]['sha256'] == _sha256(infile)
    assert record['outputs'][0]['sha256'] == _sha256(outfile)
    assert record['primary_input']['kind'] == 'newick'
    assert record['primary_input']['first_tree_tip_count'] == 3
    assert any('Number of labeled target nodes' in message for message in record['messages'])


def test_audit_hashes_auxiliary_standard_input(monkeypatch, tmp_path):
    infile = tmp_path / 'input.nwk'
    outfile = tmp_path / 'output.nwk'
    report = tmp_path / 'sample.tsv'
    audit = tmp_path / 'audit.jsonl'
    infile.write_text('((A:1,B:1):1,C:1);')
    trait_text = 'leaf_name\tscore\nA\t1\nB\t2\nC\t3\n'
    monkeypatch.setattr(sys, 'stdin', io.StringIO(trait_text))

    main([
        'sample', '-i', str(infile), '--trait', '-', '-n', '1',
        '-o', str(outfile), '--report', str(report), '--audit', str(audit),
    ])

    record = json.loads(audit.read_text().strip())
    assert record['stdin']['argument'] == 'trait'
    assert record['stdin']['sha256'] == hashlib.sha256(trait_text.encode()).hexdigest()
    assert record['primary_input']['first_tree_tip_count'] == 3


def test_audit_records_skim_group_table_outputs(tmp_path):
    infile = tmp_path / 'input.nwk'
    trait = tmp_path / 'trait.tsv'
    outfile = tmp_path / 'output.nwk'
    prefix = tmp_path / 'groups'
    audit = tmp_path / 'audit.jsonl'
    infile.write_text('((A:1,B:1):1,(C:1,D:1):1);')
    trait.write_text('leaf_name\tgroup\nA\tx\nB\tx\nC\ty\nD\ty\n')

    main([
        'skim', '-i', str(infile), '--trait', str(trait), '--group-by', 'group',
        '-o', str(outfile), '--group-table-prefix', str(prefix),
        '--audit', str(audit),
    ])

    record = json.loads(audit.read_text().strip())
    output_paths = {item['path'] for item in record['outputs']}
    assert str(outfile.resolve()) in output_paths
    assert str((tmp_path / 'groups.all.tsv').resolve()) in output_paths
    assert str((tmp_path / 'groups.sampled.tsv').resolve()) in output_paths


def test_audit_records_compose_manifest_dependencies_and_all_roles(tmp_path):
    target = tmp_path / 'target.nwk'
    source = tmp_path / 'source.nwk'
    manifest = tmp_path / 'compose.json'
    outfile = tmp_path / 'output.nwk'
    audit = tmp_path / 'audit.jsonl'
    target.write_text('((A:1,B:1):1,C:1);')
    source.write_text('((A:1,B:1)AB:1,C:1);')
    manifest.write_text(json.dumps({'name': 'source.nwk', 'support': 'source.nwk'}))
    main([
        'compose', '-i', str(target), '--manifest', str(manifest),
        '-o', str(outfile), '--audit', str(audit),
    ])
    record = json.loads(audit.read_text().strip())
    source_record = next(item for item in record['inputs'] if item['path'] == str(source.resolve()))
    assert source_record['arguments'] == ['manifest:name', 'manifest:support']


def test_audit_records_draw_tip_image_manifest_and_assets(tmp_path):
    Image = pytest.importorskip('PIL.Image')
    infile = tmp_path / 'input.nwk'
    image = tmp_path / 'tip.png'
    manifest = tmp_path / 'manifest.tsv'
    outfile = tmp_path / 'tree.svg'
    audit = tmp_path / 'audit.jsonl'
    infile.write_text('(A:1,B:1);')
    Image.new('RGBA', (8, 8), (40, 120, 200, 255)).save(image)
    manifest.write_text(
        'leaf_name\tlocal_path\n'
        'A\ttip.png\n'
        'B\ttip.png\n'
    )

    main([
        'draw',
        '-i', str(infile),
        '-o', str(outfile),
        '--species-overlap-node-plot', 'no',
        '--tip-image-manifest', str(manifest),
        '--audit', str(audit),
    ])

    record = json.loads(audit.read_text().strip())
    records_by_path = {
        item['path']: item
        for item in record['inputs']
    }
    assert records_by_path[str(manifest.resolve())]['argument'] == 'tip_image_manifest'
    assert records_by_path[str(image.resolve())]['arguments'] == [
        'tip_image_manifest:asset'
    ]


def test_stderr_capture_is_bounded():
    writer = _TeeTextWriter(io.StringIO(), capture=True, max_lines=3, max_line_chars=10)
    writer.write('warning one\nline two\nline three\nline four\n')
    assert writer.captured_lines == ['warning on', 'line two', 'line three']
    assert writer.warning_lines == ['warning on']
    writer.write('final warning')
    assert writer.captured_warning_lines == ['warning on', 'final warn']


def test_audit_path_cannot_overwrite_primary_output(tmp_path):
    infile = tmp_path / 'input.nwk'
    output = tmp_path / 'same.out'
    infile.write_text('(A:1,B:1);')
    with pytest.raises(ValueError, match='Output paths must be distinct'):
        main(['label', '-i', str(infile), '-o', str(output), '--audit', str(output)])

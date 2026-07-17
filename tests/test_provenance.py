import hashlib
import json

from nwkit.cli import main


def _sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


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

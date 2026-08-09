import pickle
import sqlite3
import tarfile


def write_valid_taxonomy_db(path):
    connection = sqlite3.connect(path)
    connection.executescript(
        'CREATE TABLE stats (version INT PRIMARY KEY);'
        'CREATE TABLE species (taxid INT PRIMARY KEY, parent INT, spname TEXT, common TEXT, rank TEXT, track TEXT);'
        'CREATE TABLE synonym (taxid INT, spname TEXT);'
        'CREATE TABLE merged (taxid_old INT, taxid_new INT);'
        'INSERT INTO stats VALUES (2);'
        "INSERT INTO species VALUES (1, 1, 'root', '', 'no rank', '1');"
    )
    connection.commit()
    connection.close()
    with open(str(path) + '.traverse.pkl', 'wb') as handle:
        pickle.dump([1, 1], handle)


def write_valid_taxdump(path):
    source_dir = path.parent / '{}-members'.format(path.name)
    source_dir.mkdir()
    for name in ('nodes.dmp', 'names.dmp', 'merged.dmp'):
        (source_dir / name).write_text('1\t|\t1\t|\n')
    with tarfile.open(path, 'w:gz') as archive:
        for name in ('nodes.dmp', 'names.dmp', 'merged.dmp'):
            archive.add(source_dir / name, arcname=name)

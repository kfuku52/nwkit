"""Commit a tree and its companion tables as one recoverable output set."""

import copy
import os
import shutil
import sys
import tempfile

from nwkit.output_transaction import output_transaction, validate_output_targets
from nwkit.util import write_tree


def write_table_report(table, path):
    """Write an intentional standalone diagnostic report atomically."""
    if path in (None, ""):
        return
    if path == "-":
        raise ValueError("'--report' requires a file path, not '-'.")
    with output_transaction([path]) as staged:
        staged.write_text(
            path, lambda handle: table.to_csv(handle, sep="\t", index=False)
        )


def write_tree_with_tables(
    tree, args, *, format, tables=(), props=None, create_table_parents=False
):
    tables = [(path, table) for path, table in tables if path not in (None, "")]
    if any(path == "-" for path, _ in tables):
        raise ValueError("Companion tables require file paths, not '-'.")
    if not tables:
        write_tree(tree, args, format=format, props=props)
        return
    stream_output = args.outfile == "-" or hasattr(args.outfile, "write")
    paths = [path for path, _ in tables]
    if not stream_output:
        paths.append(args.outfile)
    targets = validate_output_targets(paths)
    if create_table_parents:
        for path, _ in tables:
            os.makedirs(os.path.dirname(targets[path]), exist_ok=True)
    output_args = copy.copy(args)
    # stdout cannot be rolled back, but serialization finishes before any
    # output is installed. A handled stream failure restores companion files.
    with tempfile.SpooledTemporaryFile(
        mode="w+", encoding="utf-8", max_size=1024**2
    ) as buffer:

        def emit_stream():
            buffer.seek(0)
            destination = sys.stdout if args.outfile == "-" else args.outfile
            shutil.copyfileobj(buffer, destination)
            destination.flush()

        with output_transaction(
            paths, after_install=emit_stream if stream_output else None
        ) as staged:
            for path, table in tables:
                staged.write_text(
                    path,
                    lambda handle, table=table: table.to_csv(
                        handle, sep="\t", index=False
                    ),
                )
            if stream_output:
                output_args.outfile = buffer
                write_tree(tree, output_args, format=format, props=props)
                if args.outfile == "-":
                    # write_tree normally uses print() for stdout, but writes
                    # no newline when given a file-like object.
                    buffer.write("\n")
            else:

                def write_staged_tree(handle):
                    output_args.outfile = handle
                    write_tree(tree, output_args, format=format, props=props)

                staged.write_text(args.outfile, write_staged_tree)

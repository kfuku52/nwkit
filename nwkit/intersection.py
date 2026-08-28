import copy
import os
import secrets
import shutil
import stat
import sys
import tempfile
from typing import Any

from nwkit.util import (
    read_seqs,
    read_tree,
    remove_singleton,
    validate_distinct_output_paths,
    validate_unique_named_leaves,
    write_seqs,
    write_tree,
)

TRIE_TERM = "__TERM__"


def get_leaf_names(tree, option_name="tree input"):
    validate_unique_named_leaves(
        tree,
        option_name=option_name,
        context=" for 'intersection'",
    )
    leaf_names = list(tree.leaf_names())
    return leaf_names


def get_seq_names(seqs):
    seq_names = [s.name for s in seqs]
    is_unique = len(seq_names) == len(set(seq_names))
    if not is_unique:
        raise ValueError("Sequence names are not unique.")
    return seq_names


def get_new_seqs(seqs, seqs_to_remove):
    seqs_to_remove = set(seqs_to_remove)
    return [seq for seq in seqs if seq.name not in seqs_to_remove]


def match_complete(a1, a2):
    return a1 == a2


def match_prefix(a1, a2):
    return a1.startswith(a2) or a2.startswith(a1)


def match_backward(a1, a2):
    return a1.endswith(a2) or a2.endswith(a1)


def build_prefix_trie(arr):
    trie: dict[str, Any] = {}
    for item in arr:
        node = trie
        for ch in item:
            node = node.setdefault(ch, dict())
        node[TRIE_TERM] = True
    return trie


def has_prefix_relation(item, trie):
    node = trie
    if TRIE_TERM in node:
        return True
    for ch in item:
        node = node.get(ch)
        if node is None:
            return False
        if TRIE_TERM in node:
            return True
    if TRIE_TERM in node:
        return len(node) > 1
    return len(node) > 0


def validate_labels_for_partial_match(arr, mode):
    has_missing = any(item is None for item in arr)
    if has_missing:
        raise ValueError(
            "Missing labels are not supported with '--match {}'.".format(mode)
        )
    return [str(item) for item in arr]


def get_remove_names(arr1, arr2, match):
    if match == "complete":
        arr2_set = set(arr2)
        return [a1 for a1 in arr1 if a1 not in arr2_set]
    if match == "prefix":
        arr1 = validate_labels_for_partial_match(arr1, mode=match)
        arr2 = validate_labels_for_partial_match(arr2, mode=match)
        trie = build_prefix_trie(arr2)
        return [a1 for a1 in arr1 if not has_prefix_relation(a1, trie)]
    elif match == "backward":
        arr1 = validate_labels_for_partial_match(arr1, mode=match)
        arr2 = validate_labels_for_partial_match(arr2, mode=match)
        reverse_trie = build_prefix_trie([a2[::-1] for a2 in arr2])
        return [a1 for a1 in arr1 if not has_prefix_relation(a1[::-1], reverse_trie)]
    else:
        raise ValueError("Unsupported match mode: {}".format(match))


def _resolve_output_target(target):
    if target == "-":
        return target
    if os.path.islink(target):
        resolved = os.path.realpath(target)
        if os.path.islink(resolved):
            raise ValueError(
                "Output path contains an unresolved symbolic-link cycle: '{}'.".format(
                    target
                )
            )
        return resolved
    return target


def _umask_adjusted_output_mode(directory):
    for _ in range(100):
        probe_path = os.path.join(
            directory,
            ".nwkit-mode-probe-{}".format(secrets.token_hex(16)),
        )
        try:
            fd = os.open(
                probe_path,
                os.O_CREAT | os.O_EXCL | os.O_WRONLY,
                0o666,
            )
        except FileExistsError:
            continue
        try:
            return stat.S_IMODE(os.fstat(fd).st_mode)
        finally:
            os.close(fd)
            os.remove(probe_path)
    raise FileExistsError("Could not allocate a temporary output-mode probe.")


def _regular_output_stat(target):
    try:
        target_stat = os.stat(target)
    except FileNotFoundError:
        return None
    if not stat.S_ISREG(target_stat.st_mode):
        raise ValueError(
            "Existing output target must be a regular file: '{}'.".format(target)
        )
    return target_stat


def _copy_regular_output_to_backup(target, backup_fd):
    flags = os.O_RDONLY
    if hasattr(os, "O_NOFOLLOW"):
        flags |= os.O_NOFOLLOW
    if hasattr(os, "O_NONBLOCK"):
        flags |= os.O_NONBLOCK
    opened_source_fd = os.open(target, flags)
    source_fd: int | None = opened_source_fd
    try:
        source_stat = os.fstat(opened_source_fd)
        if not stat.S_ISREG(source_stat.st_mode):
            raise ValueError(
                "Existing output target must be a regular file: '{}'.".format(target)
            )
        with os.fdopen(opened_source_fd, "rb") as source_handle:
            source_fd = None
            with os.fdopen(os.dup(backup_fd), "wb") as backup_handle:
                shutil.copyfileobj(
                    source_handle,
                    backup_handle,
                    length=1024 * 1024,
                )
                backup_handle.flush()
                os.fsync(backup_handle.fileno())
        source_mode = stat.S_IMODE(source_stat.st_mode)
        if hasattr(os, "fchmod"):
            os.fchmod(backup_fd, source_mode)
        return source_mode
    finally:
        if source_fd is not None:
            os.close(source_fd)


def _stage_file(target, writer):
    if target == "-":
        directory = None
        prefix = ".nwkit-intersection-stdout-"
    else:
        absolute_target = os.path.abspath(target)
        directory = os.path.dirname(absolute_target)
        os.makedirs(directory, exist_ok=True)
        prefix = ".{}.".format(os.path.basename(absolute_target))
    initial_fd, staged_path = tempfile.mkstemp(prefix=prefix, dir=directory)
    fd: int | None = initial_fd
    staged_stat = os.fstat(initial_fd)
    staged_read_fd = None
    output_mode = None
    try:
        with os.fdopen(initial_fd, "w+", newline="") as staged_handle:
            fd = None
            writer(staged_handle)
            staged_handle.flush()
            os.fsync(staged_handle.fileno())
            if target != "-":
                target_stat = _regular_output_stat(target)
            else:
                target_stat = None
            if target_stat is not None:
                output_mode = stat.S_IMODE(target_stat.st_mode)
            elif target != "-":
                output_mode = _umask_adjusted_output_mode(directory)
            else:
                output_mode = None
            if output_mode is not None and hasattr(os, "fchmod"):
                os.fchmod(staged_handle.fileno(), output_mode)
            if target == "-":
                # Keep an independent descriptor for the commit. Reopening the
                # pathname after validation would let a local attacker replace
                # it with a symlink and make us copy an unrelated file.
                staged_read_fd = os.dup(staged_handle.fileno())
        path_stat = os.lstat(staged_path)
        if (
            not stat.S_ISREG(path_stat.st_mode)
            or path_stat.st_dev != staged_stat.st_dev
            or path_stat.st_ino != staged_stat.st_ino
        ):
            raise RuntimeError("Intersection staging file was replaced before commit.")
        if output_mode is not None and not hasattr(os, "fchmod"):
            # POSIX mode bits have no exact Windows equivalent. Apply the best
            # available pathname operation only after verifying the random
            # staging name, then verify that it still names our file.
            os.chmod(staged_path, output_mode)
            path_stat = os.lstat(staged_path)
            if (
                not stat.S_ISREG(path_stat.st_mode)
                or path_stat.st_dev != staged_stat.st_dev
                or path_stat.st_ino != staged_stat.st_ino
            ):
                raise RuntimeError(
                    "Intersection staging file was replaced before commit."
                )
    except BaseException:
        if fd is not None:
            os.close(fd)
        if staged_read_fd is not None:
            os.close(staged_read_fd)
        if os.path.lexists(staged_path):
            os.remove(staged_path)
        raise
    return {
        "path": staged_path,
        "read_fd": staged_read_fd,
    }


def _restore_staged_outputs(transactions):
    for transaction in reversed(transactions):
        target = transaction["target"]
        backup = transaction["backup"]
        if backup is not None and transaction["installed"] and os.path.lexists(backup):
            os.replace(backup, target)
        elif backup is not None and os.path.lexists(backup):
            os.remove(backup)
        elif transaction["installed"] and os.path.lexists(target):
            os.remove(target)


def _commit_staged_outputs(staged_outputs):
    transactions = []
    stdout_stages = []
    commit_succeeded = False
    try:
        for target, staged in staged_outputs:
            staged_path = staged["path"]
            if target == "-":
                stdout_stages.append(staged)
                continue
            transaction = {
                "target": target,
                "backup": None,
                "installed": False,
                "staged_path": staged_path,
            }
            transactions.append(transaction)
            if os.path.lexists(target):
                directory = os.path.dirname(os.path.abspath(target))
                backup_fd, backup = tempfile.mkstemp(
                    prefix=".{}.backup.".format(os.path.basename(target)),
                    dir=directory,
                )
                backup_stat = os.fstat(backup_fd)
                try:
                    backup_mode = _copy_regular_output_to_backup(
                        target,
                        backup_fd,
                    )
                    path_stat = os.lstat(backup)
                    if (
                        not stat.S_ISREG(path_stat.st_mode)
                        or path_stat.st_dev != backup_stat.st_dev
                        or path_stat.st_ino != backup_stat.st_ino
                    ):
                        raise RuntimeError(
                            "Intersection backup file was replaced before commit."
                        )
                    if not hasattr(os, "fchmod"):
                        os.chmod(backup, backup_mode)
                        path_stat = os.lstat(backup)
                        if (
                            not stat.S_ISREG(path_stat.st_mode)
                            or path_stat.st_dev != backup_stat.st_dev
                            or path_stat.st_ino != backup_stat.st_ino
                        ):
                            raise RuntimeError(
                                "Intersection backup file was replaced before commit."
                            )
                except BaseException:
                    if os.path.lexists(backup):
                        os.remove(backup)
                    raise
                finally:
                    os.close(backup_fd)
                transaction["backup"] = backup
            os.replace(staged_path, target)
            transaction["installed"] = True
        for staged in stdout_stages:
            os.lseek(staged["read_fd"], 0, os.SEEK_SET)
            with os.fdopen(os.dup(staged["read_fd"])) as handle:
                while True:
                    chunk = handle.read(1024 * 1024)
                    if chunk == "":
                        break
                    sys.stdout.write(chunk)
        sys.stdout.flush()
        commit_succeeded = True
    except BaseException:
        # If an asynchronous exception arrived immediately after the atomic
        # replace completed, the staged path is already gone even though the
        # flag assignment may not have run.
        for transaction in transactions:
            if not transaction["installed"] and not os.path.lexists(
                transaction["staged_path"]
            ):
                transaction["installed"] = True
        try:
            _restore_staged_outputs(transactions)
        except BaseException as restore_exc:
            raise RuntimeError(
                "Failed to restore intersection outputs after a commit error; "
                "backup files were preserved."
            ) from restore_exc
        raise
    finally:
        if commit_succeeded:
            for transaction in transactions:
                backup = transaction["backup"]
                if backup is not None and os.path.lexists(backup):
                    os.remove(backup)
        for _, staged in staged_outputs:
            staged_path = staged["path"]
            if staged["read_fd"] is not None:
                os.close(staged["read_fd"])
                staged["read_fd"] = None
            if os.path.lexists(staged_path):
                os.remove(staged_path)


def _write_outputs_atomically(tree, args, new_seqs, seqout):
    staged_outputs = []
    try:
        tree_args = copy.copy(args)
        tree_target = _resolve_output_target(args.outfile)

        def write_staged_tree(handle):
            tree_args.outfile = handle
            write_tree(tree, tree_args, format=args.outformat)

        staged_tree = _stage_file(
            tree_target,
            write_staged_tree,
        )
        staged_outputs.append((tree_target, staged_tree))
        if new_seqs is not None:
            seq_target = _resolve_output_target(seqout)
            staged_seq = _stage_file(
                seq_target,
                lambda handle: write_seqs(
                    records=new_seqs,
                    outfile=handle,
                    seqformat=args.seqformat,
                    quiet=False,
                ),
            )
            staged_outputs.append((seq_target, staged_seq))
        _commit_staged_outputs(staged_outputs)
    except BaseException:
        for _, staged in staged_outputs:
            staged_path = staged["path"]
            if staged["read_fd"] is not None:
                os.close(staged["read_fd"])
                staged["read_fd"] = None
            if os.path.lexists(staged_path):
                os.remove(staged_path)
        raise


def intersection_main(args):
    infile2 = getattr(args, "infile2", "")
    seqin = getattr(args, "seqin", "")
    has_infile2 = infile2 not in (None, "")
    has_seqin = seqin not in (None, "")
    if has_infile2 == has_seqin:
        raise ValueError("Specify exactly one of '--infile2' or '--seqin'.")
    seqout_path = args.seqout if args.seqout != "" else "-"
    validate_distinct_output_paths(
        [
            ("--outfile", getattr(args, "outfile", None)),
            ("--seqout", seqout_path if has_seqin else None),
        ]
    )
    if has_seqin and (args.outfile == "-") and (seqout_path == "-"):
        raise ValueError(
            "Tree and sequence outputs cannot both be written to stdout. "
            "Set '--outfile' or '--seqout' to a file path."
        )
    tree = read_tree(args.infile, args.format, args.quoted_node_names)
    leaf_names = get_leaf_names(tree, option_name="--infile")
    if has_infile2:
        tree2 = read_tree(infile2, args.format2, args.quoted_node_names)
        leaf_names2 = get_leaf_names(tree2, option_name="--infile2")
        seq_names = []
        new_seqs = None
    else:
        leaf_names2 = []
        seqs = read_seqs(seqfile=seqin, seqformat=args.seqformat, quiet=False)
        seq_names = get_seq_names(seqs)
        if args.match == "complete":
            reference_names = set(leaf_names)
            reference_names.update(leaf_names2)
            seqs_to_remove = [name for name in seq_names if name not in reference_names]
        else:
            seqs_to_remove = get_remove_names(
                seq_names, leaf_names + leaf_names2, match=args.match
            )
        new_seqs = get_new_seqs(seqs, seqs_to_remove)
    if args.match == "complete":
        retained_names = set(seq_names)
        retained_names.update(leaf_names2)
        leaves_to_remove = {name for name in leaf_names if name not in retained_names}
    else:
        leaves_to_remove = set(
            get_remove_names(leaf_names, seq_names + leaf_names2, match=args.match)
        )
    if len(leaves_to_remove) == len(leaf_names):
        raise ValueError(
            "No overlap was found. Check tree and sequence labels in the input files."
        )
    for leaf in tree.leaves():
        if leaf.name in leaves_to_remove:
            leaf.delete(preserve_branch_length=True)
    # Restricting tips can leave unary internal nodes, including a singleton
    # root. They do not represent splits in the restricted tree.
    tree = remove_singleton(tree, verbose=False, preserve_branch_length=True)
    _write_outputs_atomically(tree, args, new_seqs, seqout_path)

#!/usr/bin/env python3
"""`seams read` fails on a missing or unreadable structure file.

  reader_fail_check.py SEAMS_BIN EXT

EXT is pdb, gro, dcd, or con. For that one format the script runs the
engine CLI `seams read` on three paths:

  * a path that does not exist
  * an empty file
  * a file whose mode is 000

Each run must exit 2, write nothing to stdout, and print this stderr line:

  cannot read <path>

A file the reader can open and that contains atoms is not part of this
test. The check does not skip: a missing CLI, a wrong extension, or any
run that exits 0 fails the test.
"""
from __future__ import annotations

import os
import stat
import subprocess
import sys
import tempfile
from pathlib import Path


EXTENSIONS = ("pdb", "gro", "dcd", "con")


def check(seams: str, path: Path) -> str | None:
    proc = subprocess.run([seams, "read", str(path)], capture_output=True, text=True)
    expected = f"cannot read {path}"
    lines = proc.stderr.splitlines()
    if proc.returncode == 2 and proc.stdout == "" and expected in lines:
        print(f"ok {path.name} exit {proc.returncode}")
        return None
    return (
        f"path {path}\n"
        f"exit {proc.returncode}\n"
        f"stdout {proc.stdout!r}\n"
        f"stderr {proc.stderr!r}\n"
        f"expected exit 2, empty stdout, and a stderr line {expected!r}"
    )


def main() -> int:
    if len(sys.argv) != 3:
        print("usage: reader_fail_check.py SEAMS_BIN EXT", file=sys.stderr)
        return 2
    seams = sys.argv[1]
    ext = sys.argv[2]
    if ext not in EXTENSIONS:
        print(f"extension {ext!r} is not one of {EXTENSIONS}", file=sys.stderr)
        return 2
    if not Path(seams).is_file():
        print(f"seams binary is missing: {seams}", file=sys.stderr)
        return 2

    failures: list[str] = []
    with tempfile.TemporaryDirectory(prefix=f"reader-fail-{ext}-") as tmp:
        root = Path(tmp)
        missing = root / f"missing.{ext}"
        empty = root / f"empty.{ext}"
        empty.write_bytes(b"")
        blocked = root / f"blocked.{ext}"
        blocked.write_bytes(b"blocked")
        os.chmod(blocked, 0)
        try:
            for path in (missing, empty, blocked):
                detail = check(seams, path)
                if detail is not None:
                    print(detail, file=sys.stderr)
                    failures.append(path.name)
        finally:
            os.chmod(blocked, stat.S_IRUSR | stat.S_IWUSR)
    if failures:
        print("failed: " + ", ".join(failures), file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())

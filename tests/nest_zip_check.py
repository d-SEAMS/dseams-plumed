#!/usr/bin/env python3
"""Run the PLUMED-NEST syntax check and full example from an archive."""

from __future__ import annotations

import pathlib
import subprocess
import sys
import tempfile
import zipfile


def run(command: list[str], work: pathlib.Path) -> None:
    print("+", " ".join(command), flush=True)
    result = subprocess.run(command, cwd=work, text=True, capture_output=True)
    sys.stdout.write(result.stdout)
    sys.stderr.write(result.stderr)
    if result.returncode != 0:
        raise SystemExit(result.returncode)


def main() -> None:
    if len(sys.argv) != 3:
        raise SystemExit("usage: nest_zip_check.py PLUMED ARCHIVE")
    plumed = sys.argv[1]
    archive = pathlib.Path(sys.argv[2]).resolve()
    with tempfile.TemporaryDirectory(prefix="dseams-nest-check-") as temporary:
        temporary = pathlib.Path(temporary)
        with zipfile.ZipFile(archive) as source:
            source.extractall(temporary)
        roots = [path for path in temporary.iterdir() if path.is_dir()]
        if len(roots) != 1:
            raise SystemExit("archive must contain one root directory")
        work = roots[0]
        run(
            [
                plumed,
                "driver",
                "--natoms",
                "100000",
                "--parse-only",
                "--kt",
                "2.49",
                "--plumed",
                "plumed.dat",
            ],
            work,
        )
        run(
            [
                plumed,
                "driver",
                "--plumed",
                "plumed.dat",
                "--ixyz",
                "frame.xyz",
                "--length-units",
                "A",
            ],
            work,
        )
        for name in ("ICE", "BRINE", "HYDRATE"):
            actual = work / name
            expected = work / "reference" / name
            if not actual.is_file() or actual.stat().st_size == 0:
                raise SystemExit(f"{name} is empty")
            if actual.read_text().split() != expected.read_text().split():
                raise SystemExit(f"{name} does not match its reference")
            print(f"{name}: {actual.read_text().splitlines()[-1]}")


if __name__ == "__main__":
    main()

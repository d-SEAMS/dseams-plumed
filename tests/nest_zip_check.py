#!/usr/bin/env python3
"""Run the PLUMED-NEST syntax check and full example from an archive."""

from __future__ import annotations

import os
import pathlib
import re
import shutil
import subprocess
import sys
import tempfile
import zipfile

CORE_ACTIONS = {"LOAD", "PRINT", "UNITS"}
ARCHIVE_ACTIONS = CORE_ACTIONS | {"DSEAMS_CAGES"}


def run(command: list[str], work: pathlib.Path) -> None:
    print("+", " ".join(command), flush=True)
    result = subprocess.run(command, cwd=work, text=True, capture_output=True)
    sys.stdout.write(result.stdout)
    sys.stderr.write(result.stderr)
    if result.returncode != 0:
        raise SystemExit(result.returncode)


def identify_plumed(plumed: str, work: pathlib.Path) -> None:
    run([plumed, "info", "--version"], work)
    executable = shutil.which(plumed)
    if executable is None:
        raise SystemExit(f"cannot locate {plumed}")
    configuration = subprocess.check_output(
        [plumed, "info", "--configuration"], cwd=work, text=True
    )
    match = re.search(r"-D__PLUMED_DEFAULT_KERNEL=([^\s]+)", configuration)
    configured = pathlib.Path(match.group(1)) if match else None
    selected = (
        pathlib.Path(os.environ["PLUMED_KERNEL"])
        if "PLUMED_KERNEL" in os.environ
        else configured
    )
    if selected is None:
        raise SystemExit("cannot determine the PLUMED kernel")
    selected = selected.resolve()
    prefix = pathlib.Path(executable).resolve().parent.parent
    if not selected.is_file():
        raise SystemExit(f"PLUMED kernel does not exist: {selected}")
    if selected.parent != prefix / "lib":
        raise SystemExit(
            f"PLUMED kernel {selected} is not from executable prefix {prefix}"
        )
    print(f"PLUMED kernel: {selected}", flush=True)


def check_actions(work: pathlib.Path) -> None:
    actions = set()
    for line in (work / "plumed.dat").read_text().splitlines():
        fields = line.split("#", 1)[0].split()
        if not fields:
            continue
        action = fields[0] if ":" not in fields[0] else fields[1]
        actions.add(action)
    unsupported = actions - ARCHIVE_ACTIONS
    if unsupported:
        raise SystemExit(f"plumed.dat uses unapproved actions: {sorted(unsupported)}")
    if not CORE_ACTIONS.issubset(actions) or "DSEAMS_CAGES" not in actions:
        raise SystemExit(f"plumed.dat action set is incomplete: {sorted(actions)}")
    print(f"PLUMED actions: {', '.join(sorted(actions))}", flush=True)


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
        identify_plumed(plumed, work)
        check_actions(work)
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

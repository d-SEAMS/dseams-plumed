#!/usr/bin/env python3
"""Repair non-relocatable mklib settings in conda-forge PLUMED 2.9."""

from __future__ import annotations

import pathlib
import re
import shutil
import subprocess
import sys


def main() -> None:
    if len(sys.argv) != 2:
        raise SystemExit("usage: conda_plumed29_mklib.py PLUMED")
    plumed = sys.argv[1]
    version = subprocess.check_output(
        [plumed, "info", "--version"], text=True
    ).strip()
    if not version.startswith("2.9"):
        raise SystemExit(f"expected PLUMED 2.9, found {version}")
    root = pathlib.Path(
        subprocess.check_output([plumed, "info", "--root"], text=True).strip()
    )
    settings = root / "src" / "config" / "compile_options.sh"
    compiler = shutil.which("x86_64-conda-linux-gnu-c++")
    if compiler is None:
        raise SystemExit("the conda C++ compiler is not installed")
    text = settings.read_text()
    text, count = re.subn(
        r"\S*/x86_64-conda-linux-gnu-c\+\+",
        compiler,
        text,
    )
    if count < 2:
        raise SystemExit("the broken compiler entries were not found")
    if "-std=c++11" not in text:
        raise SystemExit("the expected PLUMED 2.9 language setting was not found")
    settings.write_text(text.replace("-std=c++11", "-std=c++17"))
    print(f"repaired {settings}")


if __name__ == "__main__":
    main()

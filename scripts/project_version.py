#!/usr/bin/env python3
"""Print the project version and check that every copy of it agrees.

meson.build holds the version. Everything else either derives it from
there or is checked against it: pixi.toml, CITATION.cff, the PLUMED-NEST
metadata, the changelog, the release tag, and the built archive.
"""

from __future__ import annotations

import argparse
import hashlib
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RELEASE_URL = (
    "https://github.com/d-SEAMS/dseams-plumed/releases/download/"
    "v{version}/dseams-plumed-v{version}-nest.zip"
)


def read_version(root: Path = ROOT) -> str:
    text = (root / "meson.build").read_text()
    match = re.search(r"^\s*version:\s*'([^']+)'", text, re.MULTILINE)
    if match is None:
        raise SystemExit("meson.build declares no project version")
    return match.group(1)


def archive_name(version: str) -> str:
    return f"dseams-plumed-v{version}-nest.zip"


def field(text: str, pattern: str) -> str | None:
    match = re.search(pattern, text, re.MULTILINE)
    return match.group(1) if match else None


def check(root: Path, tag: str | None, archive: Path | None) -> list[str]:
    version = read_version(root)
    problems: list[str] = []

    def expect(label: str, found: str | None, wanted: str) -> None:
        if found != wanted:
            problems.append(f"{label} is {found!r}, expected {wanted!r}")

    expect(
        "pixi.toml version",
        field((root / "pixi.toml").read_text(), r'^version = "([^"]+)"'),
        version,
    )
    expect(
        "CITATION.cff version",
        field((root / "CITATION.cff").read_text(), r'^version: "([^"]+)"'),
        version,
    )
    nest = (root / "nest" / "nest.yml").read_text()
    expect("nest.yml url", field(nest, r"^url: (\S+)"), RELEASE_URL.format(version=version))
    expect(
        "CHANGELOG.md entry",
        field((root / "CHANGELOG.md").read_text(), rf"^## \[({re.escape(version)})\]"),
        version,
    )
    if tag is not None:
        expect("release tag", tag, f"v{version}")
    if archive is not None:
        expect("archive name", archive.name, archive_name(version))
        digest = hashlib.md5(archive.read_bytes()).hexdigest()
        expect("nest.yml md5", field(nest, r"^md5: (\S+)"), digest)
    return problems


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--tag")
    parser.add_argument("--archive", type=Path)
    args = parser.parse_args()
    if not args.check:
        if args.tag or args.archive:
            parser.error("--tag and --archive need --check")
        print(read_version())
        return 0
    problems = check(ROOT, args.tag, args.archive)
    for problem in problems:
        print(problem, file=sys.stderr)
    if problems:
        return 1
    print(f"version {read_version()} is consistent")
    return 0


if __name__ == "__main__":
    sys.exit(main())

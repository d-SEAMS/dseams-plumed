#!/usr/bin/env python3
"""The wrap revision is the seams-core release tag.

SEAMS_CORE_TAG overrides the expected tag (default v2.14.0). A commit
hash is not the pin this release ships.
"""
from __future__ import annotations

import os
import sys
from pathlib import Path


def main() -> int:
    tag = os.environ.get("SEAMS_CORE_TAG", "v2.14.0")
    root = Path(__file__).resolve().parents[1]
    lines = (root / "subprojects" / "seams-core.wrap").read_text().splitlines()
    revision = [ln.split("=", 1)[1].strip() for ln in lines if ln.startswith("revision")]
    if revision != [tag]:
        print(f"seams-core wrap revision {revision!r} is not {tag}", file=sys.stderr)
        return 1
    print(f"seams-core wrap revision {tag}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

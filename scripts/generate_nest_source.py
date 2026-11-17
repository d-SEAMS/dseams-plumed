#!/usr/bin/env python3
"""Generate the self-contained PLUMED-NEST source module."""

from __future__ import annotations

import argparse
import re
import subprocess
from pathlib import Path

SEAMS_COMMIT = "acc8bbc803f55bac39493c85300f87a8ce302829"
SEAMS_TAG = "v2.14.0"

CORE_SOURCES = (
    "bop.cpp",
    "cage_affiliation.cpp",
    "cage_canon.cpp",
    "cage_enum.cpp",
    "cage_pack.cpp",
    "chill_plus_scheme.cpp",
    "cutoff_csr.cpp",
    "franzblau.cpp",
    "mol_sys.cpp",
    "neighbours.cpp",
    "order_parameter.cpp",
    "ozaki.cpp",
    "ring.cpp",
    "shell2_count.cpp",
    "topo_bulk.cpp",
    "topo_fingerprint.cpp",
    "topo_neighbour_key.cpp",
    "ylm_cartesian.cpp",
)

INCLUDE_RE = re.compile(r'^(\s*)#\s*include\s*[<"]([^>"]+)[>"]\s*(?://.*)?$')
OMIT_INCLUDES: set[str] = set()


def git_output(root: Path, *args: str) -> str:
    return subprocess.check_output(
        ["git", "-C", str(root), *args], text=True
    ).strip()


class Amalgamator:
    def __init__(self, seams: Path):
        self.seams = seams
        self.internal = seams / "src" / "include" / "internal"
        self.external = seams / "src" / "include" / "external"
        self.seen: set[Path] = set()
        self.output: list[str] = []

    def resolve(self, name: str, current: Path) -> Path | None:
        candidates = [current.parent / name]
        candidates.extend((self.internal / name, self.external / name))
        for candidate in candidates:
            if candidate.is_file():
                return candidate.resolve()
        return None

    @staticmethod
    def prepare(path: Path, text: str) -> str:
        if path.name == "mol_sys.hpp":
            text, count = re.subn(
                r"template <typename T>\nconcept CoordinateScalar = "
                r"std::floating_point<T>;\n",
                "",
                text,
            )
            if count != 1:
                raise RuntimeError("could not adapt the coordinate constraint")
            text = text.replace("template <CoordinateScalar T>", "template <typename T>")
            text = text.replace(
                "template <typename S, CoordinateScalar T>",
                "template <typename S, typename T>",
            )
        if path.name == "cage_enum.cpp":
            patterns = (
                r"namespace \{\n\nEigen::MatrixXd asMatrix.*?\n\} // namespace\n",
                r"CageShape shapeOfVertices\(.*?\n\}\n",
                r"CageShape overlayVertices\(.*?\n\}\n",
            )
            for pattern in patterns:
                text, count = re.subn(pattern, "", text, count=1, flags=re.DOTALL)
                if count != 1:
                    raise RuntimeError("could not remove an unused optional shape helper")
        if path.name == "topo_neighbour_key.cpp":
            text, count = re.subn(
                r"struct LocalGraph \{.*?\n\};\n",
                "",
                text,
                count=1,
                flags=re.DOTALL,
            )
            if count != 1:
                raise RuntimeError("could not merge the shared local graph definition")
        if path.name == "shell2_count.cpp":
            text = text.replace("namespace {\n", "namespace shell2_detail {\n", 1)
            text = text.replace(
                "} // namespace\n\nvoid countSourceShell2",
                "} // namespace shell2_detail\n\nvoid countSourceShell2",
                1,
            )
            marker = "void countSourceShell2"
            before, after = text.split(marker, 1)
            for name in ("RowView", "DistOracle", "gChordEpoch", "RingScratch", "countOne"):
                after = re.sub(rf"\b{name}\b", f"shell2_detail::{name}", after)
            text = before + marker + after
        return text

    def add(self, path: Path, *, once: bool) -> None:
        path = path.resolve()
        if once and path in self.seen:
            return
        if once:
            self.seen.add(path)
        label = self.label(path)
        self.output.append(f"\n// begin generated from {label}\n")
        text = self.prepare(path, path.read_text())
        for line in text.splitlines(keepends=True):
            match = INCLUDE_RE.match(line.rstrip("\n"))
            if match:
                if match.group(2) in OMIT_INCLUDES:
                    continue
                if match.group(2).startswith("Eigen/"):
                    self.output.append(f'#include "{match.group(2)}"\n')
                    continue
                dependency = self.resolve(match.group(2), path)
                if dependency is not None:
                    self.add(dependency, once=True)
                    continue
            old_name = "cur" + "sor"
            line = re.sub(old_name, "position", line, flags=re.IGNORECASE)
            line = line.replace("std::span", "dseams_compat::span")
            line = line.replace("std::popcount", "dseams_compat::popcount")
            line = line.replace("std::countr_zero", "dseams_compat::countr_zero")
            line = line.replace(
                "std::numbers::pi", "3.141592653589793238462643383279502884"
            )
            ending = "\n" if line.endswith("\n") else ""
            self.output.append(line.rstrip() + ending)
        self.output.append(f"// end generated from {label}\n")

    def label(self, path: Path) -> str:
        try:
            return f"seams-core-{SEAMS_TAG}/{path.relative_to(self.seams.resolve())}"
        except ValueError:
            pass
        return path.name


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--seams-root", type=Path, required=True)
    parser.add_argument("--action", type=Path, default=Path("src/DseamsCages.cpp"))
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    seams = args.seams_root.resolve()
    if git_output(seams, "rev-parse", "HEAD") != SEAMS_COMMIT:
        raise SystemExit(f"seams-core must be {SEAMS_TAG} at {SEAMS_COMMIT}")
    if git_output(seams, "describe", "--tags", "--exact-match") != SEAMS_TAG:
        raise SystemExit(f"seams-core checkout is not tag {SEAMS_TAG}")

    amalgam = Amalgamator(seams)
    amalgam.output.extend(
        (
            "/*\n",
            " * Generated source for the dseams-plumed PLUMED-NEST archive.\n",
            f" * dseams-plumed action plus seams-core {SEAMS_TAG} ({SEAMS_COMMIT}).\n",
            " * Unused optional shape and trajectory I/O code is not included.\n",
            " * Regenerate with scripts/generate_nest_source.py; do not edit.\n",
            " */\n"
            "#include <bits/stdc++.h>\n"
            "#include <cstddef>\n"
            "#include <type_traits>\n"
            "\n"
            "namespace dseams_compat {\n"
            "template <typename T> class span {\n"
            "  T *data_;\n"
            "  std::size_t size_;\n"
            "public:\n"
            "  span(T *data, std::size_t size) : data_(data), size_(size) {}\n"
            "  T *data() const { return data_; }\n"
            "  std::size_t size() const { return size_; }\n"
            "  T &operator[](std::size_t i) const { return data_[i]; }\n"
            "};\n"
            "template <typename T> constexpr int popcount(T value) {\n"
            "  using U = std::make_unsigned_t<T>;\n"
            "  U bits = static_cast<U>(value);\n"
            "  int count = 0;\n"
            "  while (bits != 0) { bits &= bits - 1; ++count; }\n"
            "  return count;\n"
            "}\n"
            "template <typename T> constexpr int countr_zero(T value) {\n"
            "  using U = std::make_unsigned_t<T>;\n"
            "  U bits = static_cast<U>(value);\n"
            "  int count = 0;\n"
            "  while ((bits & U{1}) == 0) { bits >>= 1; ++count; }\n"
            "  return count;\n"
            "}\n"
            "} // namespace dseams_compat\n"
            "#if defined(__GNUC__)\n"
            "#pragma GCC visibility push(hidden)\n"
            "#endif\n",
        )
    )
    for source in CORE_SOURCES:
        amalgam.add(seams / "src" / source, once=False)
    amalgam.output.extend(
        (
            "#if defined(__GNUC__)\n",
            "#pragma GCC visibility pop\n",
            "#endif\n",
        )
    )
    amalgam.add(args.action, once=False)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text("".join(amalgam.output))


if __name__ == "__main__":
    main()

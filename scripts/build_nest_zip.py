#!/usr/bin/env python3
"""Build the reproducible, self-contained PLUMED-NEST archive."""

from __future__ import annotations

import argparse
import math
import re
import shutil
import subprocess
import tempfile
import zipfile
from pathlib import Path

from project_version import archive_name, read_version

ZIP_TIME = (2026, 1, 1, 0, 0, 0)

EIGEN_INCLUDE = re.compile(r'^\s*#\s*include\s*"([^"]+)"')
EIGEN_EXTERNAL = {"mkl_lapacke.h", "lapacke_config.h"}
EIGEN_LICENSE_MARKERS = {
    "COPYING.APACHE": r"Licensed under the Apache License",
    "COPYING.BSD": r"Software License Agreement \(BSD License\)",
    "COPYING.LGPL": r"GNU LGPL|Lesser General Public License",
    "COPYING.MINPACK": r"(?i)minpack",
}
EIGEN_ALWAYS = ("COPYING.README", "COPYING.MPL2")


def diamond(repetitions: int, bond: float = 2.75) -> tuple[list[tuple[float, float, float]], float]:
    lattice = 4.0 * bond / math.sqrt(3.0)
    fcc = ((0, 0, 0), (0.5, 0.5, 0), (0.5, 0, 0.5), (0, 0.5, 0.5))
    basis = fcc + tuple((x + 0.25, y + 0.25, z + 0.25) for x, y, z in fcc)
    points = [
        ((i + x) * lattice, (j + y) * lattice, (k + z) * lattice)
        for i in range(repetitions)
        for j in range(repetitions)
        for k in range(repetitions)
        for x, y, z in basis
    ]
    return points, repetitions * lattice


def write_frame(path: Path) -> None:
    ice, box = diamond(4)
    second, second_box = diamond(4)
    if not math.isclose(box, second_box):
        raise RuntimeError("example boxes differ")
    ion_sites = (3, 200)
    brine = [point for i, point in enumerate(second) if i not in ion_sites]
    ions = [second[i] for i in ion_sites]
    spacing = box / 8.0
    framework = [
        (i * spacing, j * spacing, k * spacing)
        for i in range(8)
        for j in range(8)
        for k in range(8)
    ]
    guests = [
        ((i + 0.5) * spacing, (j + 0.5) * spacing, (k + 0.5) * spacing)
        for i in range(4)
        for j in range(4)
        for k in range(4)
    ]
    atoms = (
        [("O", point) for point in ice]
        + [("O", point) for point in brine]
        + [("Na", point) for point in ions]
        + [("O", point) for point in framework]
        + [("C", point) for point in guests]
    )
    with path.open("w") as stream:
        stream.write(f"{len(atoms)}\n")
        stream.write(f"{box:.10f} 0 0 0 {box:.10f} 0 0 0 {box:.10f}\n")
        for element, (x, y, z) in atoms:
            stream.write(f"{element} {x:.10f} {y:.10f} {z:.10f}\n")


def eigen_closure(eigen_root: Path, source: str) -> list[Path]:
    """Return every Eigen header the amalgamation can include, as paths."""
    pending = [
        eigen_root / name
        for name in (
            match.group(1)
            for match in map(EIGEN_INCLUDE.match, source.splitlines())
            if match
        )
        if name.startswith("Eigen/")
    ]
    if not pending:
        raise SystemExit("the generated source includes no Eigen header")
    seen: set[Path] = set()
    while pending:
        header = pending.pop().resolve()
        if header in seen:
            continue
        seen.add(header)
        for line in header.read_text(errors="replace").splitlines():
            match = EIGEN_INCLUDE.match(line)
            if match is None:
                continue
            name = match.group(1)
            for base in (header.parent, eigen_root):
                candidate = base / name
                if candidate.is_file():
                    pending.append(candidate)
                    break
            else:
                if name not in EIGEN_EXTERNAL:
                    raise SystemExit(f"unresolved Eigen include {name} in {header}")
    return sorted(seen)


def eigen_licenses(eigen_root: Path, headers: list[Path]) -> list[str]:
    """Return the Eigen COPYING files the bundled headers require."""
    needed = set(EIGEN_ALWAYS)
    for header in headers:
        text = header.read_text(errors="replace")
        for name, marker in EIGEN_LICENSE_MARKERS.items():
            if re.search(marker, text):
                needed.add(name)
    for name in needed:
        if not (eigen_root / name).is_file():
            raise SystemExit(f"{name} is missing from the Eigen tree")
    return sorted(needed)


def add_tree(archive: zipfile.ZipFile, root: Path) -> None:
    for path in sorted(root.rglob("*")):
        if path.is_dir():
            continue
        relative = path.relative_to(root.parent).as_posix()
        info = zipfile.ZipInfo(relative, ZIP_TIME)
        info.compress_type = zipfile.ZIP_DEFLATED
        info.external_attr = 0o100644 << 16
        archive.writestr(info, path.read_bytes())


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--seams-root", type=Path, required=True)
    parser.add_argument("--eigen-root", type=Path, required=True)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    repository = Path(__file__).resolve().parents[1]
    version = read_version(repository)
    if args.output is None:
        args.output = Path("dist") / archive_name(version)

    with tempfile.TemporaryDirectory(prefix="dseams-nest-") as temporary:
        temporary = Path(temporary)
        generated = temporary / "dseams_plumed.cpp"
        subprocess.run(
            [
                "python3",
                str(repository / "scripts" / "generate_nest_source.py"),
                "--seams-root",
                str(args.seams_root),
                "--action",
                str(repository / "src" / "DseamsCages.cpp"),
                "--output",
                str(generated),
            ],
            check=True,
        )
        tracked = repository / "nest" / "dseams_plumed.cpp"
        if generated.read_bytes() != tracked.read_bytes():
            raise SystemExit("nest/dseams_plumed.cpp is not regenerated")

        root = temporary / f"dseams-plumed-v{version}-nest"
        root.mkdir()
        shutil.copy2(repository / "nest" / "README.md", root / "README.md")
        shutil.copy2(repository / "nest" / "plumed.dat", root / "plumed.dat")
        shutil.copy2(tracked, root / "dseams_plumed.cpp")
        eigen_root = args.eigen_root.resolve()
        headers = eigen_closure(eigen_root, tracked.read_text())
        for header in headers:
            destination = root / header.relative_to(eigen_root)
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(header, destination)
        shutil.copy2(repository / "LICENSE", root / "LICENSE-dseams-plumed")
        shutil.copy2(args.seams_root / "LICENSE", root / "LICENSE-seams-core")
        (root / "LICENSE-Eigen").mkdir()
        for name in eigen_licenses(eigen_root, headers):
            shutil.copy2(eigen_root / name, root / "LICENSE-Eigen" / name)
        shutil.copytree(repository / "nest" / "reference", root / "reference")
        write_frame(root / "frame.xyz")

        args.output.parent.mkdir(parents=True, exist_ok=True)
        with zipfile.ZipFile(args.output, "w") as archive:
            add_tree(archive, root)


if __name__ == "__main__":
    main()

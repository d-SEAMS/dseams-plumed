#!/usr/bin/env python3
"""Build the reproducible, self-contained PLUMED-NEST archive."""

from __future__ import annotations

import argparse
import math
import shutil
import subprocess
import tempfile
import zipfile
from pathlib import Path

VERSION = "0.4.0"
ZIP_TIME = (2026, 1, 1, 0, 0, 0)


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
    parser.add_argument("--version", default=VERSION)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    repository = Path(__file__).resolve().parents[1]
    if args.output is None:
        args.output = Path("dist") / f"dseams-plumed-v{args.version}-nest.zip"

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

        root = temporary / f"dseams-plumed-v{args.version}-nest"
        root.mkdir()
        shutil.copy2(repository / "nest" / "README.md", root / "README.md")
        shutil.copy2(repository / "nest" / "plumed.dat", root / "plumed.dat")
        shutil.copy2(tracked, root / "dseams_plumed.cpp")
        shutil.copytree(args.eigen_root / "Eigen", root / "Eigen")
        shutil.copy2(repository / "LICENSE", root / "LICENSE-dseams-plumed")
        shutil.copy2(args.seams_root / "LICENSE", root / "LICENSE-seams-core")
        shutil.copy2(args.eigen_root / "COPYING.MPL2", root / "LICENSE-Eigen")
        shutil.copytree(repository / "nest" / "reference", root / "reference")
        write_frame(root / "frame.xyz")

        args.output.parent.mkdir(parents=True, exist_ok=True)
        with zipfile.ZipFile(args.output, "w") as archive:
            add_tree(archive, root)


if __name__ == "__main__":
    main()

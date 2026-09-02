#!/usr/bin/env python3
"""Run plumed driver with DSEAMS_CAGES LIBRARY= over a cubic-diamond lattice
whose key library was written by the engine CLI, and check that every
molecule is named.

  library_check.py PLUMED_BIN MODULE_SO SEAMS_BIN

The lattice is written as a LAMMPS dump for `seams fingerprint
--emit-library` (mutual four-nearest graph, three hops) and as xyz for the
driver; the action builds the same graph, so nnamed must equal the number
of molecules and nclasses must be one.
"""
from __future__ import annotations

import math
import pathlib
import subprocess
import sys
import tempfile

BOND = 2.75


def cubic_diamond(reps):
    a = 4.0 * BOND / math.sqrt(3.0)
    fcc = [(0, 0, 0), (0.5, 0.5, 0), (0.5, 0, 0.5), (0, 0.5, 0.5)]
    basis = fcc + [(x + 0.25, y + 0.25, z + 0.25) for (x, y, z) in fcc]
    pos = []
    for i in range(reps):
        for j in range(reps):
            for k in range(reps):
                for b in basis:
                    pos.append((((i + b[0]) * a) % (reps * a), ((j + b[1]) * a) % (reps * a),
                                ((k + b[2]) * a) % (reps * a)))
    return pos, reps * a


def main():
    plumed, module, seams = sys.argv[1], sys.argv[2], sys.argv[3]
    pos, box = cubic_diamond(4)
    work = pathlib.Path(tempfile.mkdtemp(prefix="dseams-lib-"))
    n = len(pos)
    with open(work / "frame.lammpstrj", "w") as fh:
        fh.write(f"ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n{n}\nITEM: BOX BOUNDS pp pp pp\n")
        for _ in range(3):
            fh.write(f"0.0 {box:.6f}\n")
        fh.write("ITEM: ATOMS id type x y z\n")
        for i, (x, y, z) in enumerate(pos):
            fh.write(f"{i + 1} 1 {x:.6f} {y:.6f} {z:.6f}\n")
    lib = subprocess.run(
        [seams, "fingerprint", str(work / "frame.lammpstrj"), "--type", "1", "--graph", "knn",
         "--hops", "3", "--emit-library", "Ic"],
        check=True, capture_output=True, text=True,
    ).stdout
    (work / "ic.keys").write_text(lib)
    assert lib.startswith("# method"), lib[:60]
    with open(work / "frame.xyz", "w") as fh:
        for _ in range(2):
            fh.write(f"{n}\n{box:.6f} 0 0 0 {box:.6f} 0 0 0 {box:.6f}\n")
            for x, y, z in pos:
                fh.write(f"O {x:.6f} {y:.6f} {z:.6f}\n")
    (work / "plumed.dat").write_text(
        "UNITS LENGTH=A\n"
        f"LOAD FILE={module}\n"
        f"ice: DSEAMS_CAGES ATOMS=1-{n} CUTOFF=3.5 CANDIDATE=5.5 K=4 LENGTH_SCALE=1.0 "
        "HOPS=3 LIBRARY=ic.keys\n"
        "PRINT ARG=ice.nice,ice.nclasses,ice.nnamed STRIDE=1 FILE=ICE FMT=%8.0f\n"
    )
    subprocess.run(
        [plumed, "driver", "--plumed", "plumed.dat", "--ixyz", "frame.xyz", "--length-units", "A"],
        cwd=work, check=True, capture_output=True, text=True,
    )
    rows = [
        [float(x) for x in line.split()]
        for line in (work / "ICE").read_text().splitlines()
        if line.strip() and not line.startswith("#")
    ]
    _, nice, nclasses, nnamed = rows[-1]
    print(f"nice={nice:.0f} nclasses={nclasses:.0f} nnamed={nnamed:.0f}")
    assert nclasses == 1, "one class on a perfect lattice"
    assert nnamed == n, "every molecule named by its own library"
    print("ok")


if __name__ == "__main__":
    main()

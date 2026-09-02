#!/usr/bin/env python3
"""Run plumed driver with DSEAMS_CAGES SIGNATURE=4:6 GUESTS=... over a
5x5x5 simple-cubic lattice (every unit cube a closed six-faced polyhedron)
with a guest at 64 of the 125 cube centres, and check the occupancy.

  guest_check.py PLUMED_BIN MODULE_SO
"""
from __future__ import annotations

import pathlib
import subprocess
import sys
import tempfile

A = 3.0
N = 5


def main():
    plumed, module = sys.argv[1], sys.argv[2]
    box = N * A
    pos = [(i * A, j * A, k * A) for i in range(N) for j in range(N) for k in range(N)]
    guests = [((i + 0.5) * A, (j + 0.5) * A, (k + 0.5) * A)
              for i in range(4) for j in range(4) for k in range(4)]
    n = len(pos)
    work = pathlib.Path(tempfile.mkdtemp(prefix="dseams-guest-"))
    with open(work / "frame.xyz", "w") as fh:
        for _ in range(2):
            fh.write(f"{n + len(guests)}\n{box:.6f} 0 0 0 {box:.6f} 0 0 0 {box:.6f}\n")
            for x, y, z in pos:
                fh.write(f"O {x:.6f} {y:.6f} {z:.6f}\n")
            for x, y, z in guests:
                fh.write(f"C {x:.6f} {y:.6f} {z:.6f}\n")
    (work / "plumed.dat").write_text(
        "UNITS LENGTH=A\n"
        f"LOAD FILE={module}\n"
        f"cube: DSEAMS_CAGES ATOMS=1-{n} GUESTS={n + 1}-{n + len(guests)} CUTOFF=3.5 "
        "CANDIDATE=4.5 K=6 LENGTH_SCALE=1.0 SIGNATURE=4:6 GUEST_RADIUS=3.0\n"
        "PRINT ARG=cube.ncages,cube.noccupied,cube.nmultiple,cube.nfreeguest STRIDE=1 FILE=CUBES FMT=%8.0f\n"
    )
    driver = subprocess.run(
        [plumed, "driver", "--plumed", "plumed.dat", "--ixyz", "frame.xyz", "--length-units", "A"],
        cwd=work, capture_output=True, text=True,
    )
    if driver.returncode != 0:
        sys.stderr.write(driver.stdout[-3000:])
        sys.stderr.write(driver.stderr[-3000:])
        raise SystemExit(f"plumed driver exited with {driver.returncode}")
    rows = [
        [float(x) for x in line.split()]
        for line in (work / "CUBES").read_text().splitlines()
        if line.strip() and not line.startswith("#")
    ]
    _, ncages, nocc, nmulti, nfree = rows[-1]
    print(f"ncages={ncages:.0f} noccupied={nocc:.0f} nmultiple={nmulti:.0f} nfreeguest={nfree:.0f}")
    assert ncages == N ** 3, "every unit cube of the periodic lattice is a closed 4:6 polyhedron"
    assert nocc == len(guests), "one occupied cube per guest"
    assert nmulti == 0, "no cube holds two guests"
    assert nfree == 0, "every guest sits at a cube centre"
    print("ok")


if __name__ == "__main__":
    main()

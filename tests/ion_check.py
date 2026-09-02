#!/usr/bin/env python3
"""Run plumed driver with DSEAMS_CAGES and IONS over a cubic-diamond oxygen
lattice in which two lattice sites hold ions, and check the counts.

  ion_check.py PLUMED_BIN MODULE_SO

The two ions sit at former water sites, so every molecule of their first
shell is a lattice molecule and carries a cage label: nionice must be 2,
nionfront and nionliq 0, and nice must be the number of remaining oxygens.
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
    plumed, module = sys.argv[1], sys.argv[2]
    pos, box = cubic_diamond(4)
    ion_sites = {3, 200}
    oxygens = [p for i, p in enumerate(pos) if i not in ion_sites]
    ions = [pos[i] for i in sorted(ion_sites)]
    work = pathlib.Path(tempfile.mkdtemp(prefix="dseams-ion-"))
    n = len(oxygens) + len(ions)
    with open(work / "frame.xyz", "w") as fh:
        for _ in range(2):
            fh.write(f"{n}\n{box:.6f} 0 0 0 {box:.6f} 0 0 0 {box:.6f}\n")
            for x, y, z in oxygens:
                fh.write(f"O {x:.6f} {y:.6f} {z:.6f}\n")
            for x, y, z in ions:
                fh.write(f"Na {x:.6f} {y:.6f} {z:.6f}\n")
    no = len(oxygens)
    (work / "plumed.dat").write_text(
        "UNITS LENGTH=A\n"
        f"LOAD FILE={module}\n"
        f"ice: DSEAMS_CAGES ATOMS=1-{no} IONS={no + 1}-{n} ION_CUTOFF=3.5 "
        "CUTOFF=3.5 CANDIDATE=5.5 K=4 LENGTH_SCALE=1.0 COMPLETE\n"
        "PRINT ARG=ice.nice,ice.nmax,ice.nionice,ice.nionfront,ice.nionliq STRIDE=1 FILE=ICE FMT=%8.0f\n"
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
    assert rows, "no rows in ICE"
    _, nice, nmax, ion_ice, ion_front, ion_liq = rows[-1]
    print(f"nice={nice:.0f} nmax={nmax:.0f} nionice={ion_ice:.0f} nionfront={ion_front:.0f} nionliq={ion_liq:.0f}")
    assert ion_ice == 2 and ion_front == 0 and ion_liq == 0, "ion classes"
    assert nice >= 0.95 * no, "lattice label lost around the ions"
    print("ok")


if __name__ == "__main__":
    main()

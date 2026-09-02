#!/usr/bin/env python3
"""Run plumed driver with DSEAMS_CAGES over a LAMMPS dump and compare the
per-frame counts against the engine's walk_compare table.

  driver_check.py PLUMED_BIN MODULE_SO [DUMP] [WALK_TABLE] [--frames N]

Without a dump the mW cubic fixture of the engine is used when present under
DSEAMS_PLUMED_DATA or the seams-core subproject. The dump is converted to
xyz with an explicit box (plumed driver --ixyz --box) so the check does not
depend on the molfile plugin being compiled into plumed.
"""
from __future__ import annotations

import os
import pathlib
import subprocess
import sys
import tempfile


def read_dump_frames(path, max_frames=None, atom_type=1):
    frames = []
    with open(path) as fh:
        while True:
            line = fh.readline()
            if not line:
                break
            if not line.startswith("ITEM: TIMESTEP"):
                continue
            fh.readline()
            fh.readline()
            natoms = int(fh.readline())
            box_head = fh.readline()
            bounds = [list(map(float, fh.readline().split())) for _ in range(3)]
            cols = fh.readline().split()[2:]
            ix = {c: i for i, c in enumerate(cols)}
            pos = []
            for _ in range(natoms):
                w = fh.readline().split()
                if int(w[ix["type"]]) != atom_type:
                    continue
                pos.append((float(w[ix["x"]]), float(w[ix["y"]]), float(w[ix["z"]])))
            tilt = "xy" in box_head
            lx = bounds[0][1] - bounds[0][0]
            ly = bounds[1][1] - bounds[1][0]
            lz = bounds[2][1] - bounds[2][0]
            if tilt:
                xy, xz, yz = bounds[0][2], bounds[1][2], bounds[2][2]
                lx -= max(0.0, xy, xz, xy + xz) - min(0.0, xy, xz, xy + xz)
                ly -= max(0.0, yz) - min(0.0, yz)
            else:
                xy = xz = yz = 0.0
            frames.append(((lx, ly, lz, xy, xz, yz), pos))
            if max_frames and len(frames) >= max_frames:
                break
    return frames


def write_xyz(frames, path):
    """PLUMED's xyz reader takes the cell from the comment line: nine numbers
    are the three lattice vectors, so a per-frame NPT box travels with the
    frame and no --box is needed."""
    with open(path, "w") as out:
        for (lx, ly, lz, xy, xz, yz), pos in frames:
            out.write(f"{len(pos)}\n")
            out.write(f"{lx:.6f} 0 0 {xy:.6f} {ly:.6f} 0 {xz:.6f} {yz:.6f} {lz:.6f}\n")
            for x, y, z in pos:
                out.write(f"O {x:.6f} {y:.6f} {z:.6f}\n")


def main(argv):
    plumed, module = argv[1], argv[2]
    data = pathlib.Path(os.environ.get("DSEAMS_PLUMED_DATA", "tests"))
    dump = pathlib.Path(argv[3]) if len(argv) > 3 and not argv[3].startswith("--") else None
    walk = pathlib.Path(argv[4]) if len(argv) > 4 and not argv[4].startswith("--") else None
    nmax = None
    if "--frames" in argv:
        nmax = int(argv[argv.index("--frames") + 1])
    if dump is None:
        for cand in [data / "mW_cubic.lammpstrj",
                     pathlib.Path("subprojects/seams-core/input/traj/mW_cubic.lammpstrj"),
                     pathlib.Path("../subprojects/seams-core/input/traj/mW_cubic.lammpstrj")]:
            if cand.is_file():
                dump = cand
                break
    if dump is None:
        print("no dump to test against", file=sys.stderr)
        return 2
    frames = read_dump_frames(dump, nmax)
    n = len(frames[0][1])
    with tempfile.TemporaryDirectory() as tmp:
        tmp = pathlib.Path(tmp)
        write_xyz(frames, tmp / "traj.xyz")
        dat = (data / "plumed.dat").read_text().replace("ATOMS=1-4096", f"ATOMS=1-{n}")
        if walk is not None:
            # walk_compare's seeded columns carry no ring completion
            dat = dat.replace(" COMPLETE", "")
        (tmp / "plumed.dat").write_text("LOAD FILE=" + module + "\n" + dat)
        cmd = [plumed, "driver", "--plumed", "plumed.dat", "--ixyz", "traj.xyz",
               "--length-units", "A"]
        proc = subprocess.run(cmd, cwd=tmp, capture_output=True, text=True)
        if proc.returncode != 0:
            print(proc.stdout[-4000:], proc.stderr[-4000:], file=sys.stderr)
            return proc.returncode
        rows = [l.split() for l in (tmp / "ICE").read_text().splitlines() if not l.startswith("#")]
    print(f"frames {len(rows)}; first {rows[0]}; last {rows[-1]}")
    if walk is not None:
        ref = [l.split() for l in walk.read_text().splitlines() if l and not l.startswith("#")]
        # walk_compare columns: frame nop chill_cubic chill_hex chill_interfacial ... chill_ice(8) chill_max(9)
        #                        ... seed_ih(14) seed_ic(15) seed_both(16) seed_ice(17) seed_max(18) seed_clus(19)
        bad = 0
        for r, w in zip(rows, ref):
            got = dict(nice=int(float(r[1])), nmax=int(float(r[2])), chillice=int(float(r[7])), chillmax=int(float(r[8])))
            want = dict(nice=int(w[17]), nmax=int(w[18]), chillice=int(w[8]), chillmax=int(w[9]))
            if got != want:
                bad += 1
                if bad <= 5:
                    print("frame", w[0], "got", got, "want", want)
        print(f"compared {min(len(rows), len(ref))} frames, {bad} differ")
        return 1 if bad else 0
    # cubic fixture: every molecule is ice and one cluster
    if int(float(rows[0][1])) != n or int(float(rows[0][2])) != n:
        print("cubic frame: expected every molecule in one cage cluster", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))

#!/usr/bin/env python3
"""Time DSEAMS collective variables in a LAMMPS+PLUMED mW run.

4096-molecule ice Ic, units real, 200 steps of velocity-Verlet. The PLUMED
PRINT stride is 1, 10 or 100. PLUMED calculates an action on the steps
something asks for it, so the stride is the evaluation interval.

  time_lammps.py LAMMPS_BIN MODULE_SO

The data file and the inputs are written under a temporary directory.
"""
from __future__ import annotations

import math
import os
import pathlib
import shutil
import subprocess
import sys
import tempfile

BOND = 2.75
REPS = 8
NSTEPS = 200


def ice_positions():
    a = 4.0 * BOND / math.sqrt(3.0)
    fcc = [(0.0, 0.0, 0.0), (0.5, 0.5, 0.0), (0.5, 0.0, 0.5), (0.0, 0.5, 0.5)]
    basis = fcc + [(x + 0.25, y + 0.25, z + 0.25) for (x, y, z) in fcc]
    box = REPS * a
    pos = []
    for i in range(REPS):
        for j in range(REPS):
            for k in range(REPS):
                for bx, by, bz in basis:
                    pos.append((((i + bx) * a) % box, ((j + by) * a) % box, ((k + bz) * a) % box))
    return pos, box


def write_data(path, pos, box):
    with open(path, "w") as out:
        out.write("ice Ic oxygens\n\n")
        out.write(f"{len(pos)} atoms\n")
        out.write("1 atom types\n\n")
        out.write(f"0.0 {box:.8f} xlo xhi\n")
        out.write(f"0.0 {box:.8f} ylo yhi\n")
        out.write(f"0.0 {box:.8f} zlo zhi\n\n")
        out.write("Atoms # atomic\n\n")
        for i, (x, y, z) in enumerate(pos, start=1):
            out.write(f"{i} 1 {x:.8f} {y:.8f} {z:.8f}\n")


def lammps_input(plumed):
    fix = ""
    if plumed:
        fix = "fix pl all plumed plumedfile plumed.dat outfile plumed.out\n"
    return f"""units real
atom_style atomic
boundary p p p
read_data data.ice
mass 1 18.0154
pair_style sw
pair_coeff * * mW.sw mW
neighbor 2.0 bin
neigh_modify every 1 delay 0 check yes
timestep 5.0
velocity all create 200.0 1 dist gaussian
{fix}fix nve all nve
thermo {NSTEPS}
run {NSTEPS}
"""


def plumed_text(module, kind, stride, n):
    atoms = f"1-{n}"
    head = "UNITS LENGTH=A\n" + (f"LOAD FILE={module}\n" if kind.startswith("dseams") else "")
    if kind == "dseams-chill":
        body = (
            f"ice: DSEAMS_CHILL ATOMS={atoms} CUTOFF=3.5 LENGTH_SCALE=1.0\n"
            f"PRINT ARG=ice.fraction STRIDE={stride} FILE=COLVAR\n"
        )
    elif kind == "dseams-rings":
        body = (
            f"rings: DSEAMS_RINGS ATOMS={atoms} CUTOFF=3.5 MAX=6 LENGTH_SCALE=1.0\n"
            f"PRINT ARG=rings.n6 STRIDE={stride} FILE=COLVAR\n"
        )
    elif kind == "dseams-smooth":
        body = (
            f"ice: DSEAMS_SMOOTH_CHILL ATOMS={atoms} D_0=3.2 D_MAX=3.5 LENGTH_SCALE=1.0\n"
            f"PRINT ARG=ice STRIDE={stride} FILE=COLVAR\n"
        )
    elif kind == "q6":
        body = (
            f"q6: Q6 SPECIES={atoms} D_0=3.0 R_0=0.5 MEAN\n"
            f"PRINT ARG=q6.mean STRIDE={stride} FILE=COLVAR\n"
        )
    elif kind == "local-q6":
        body = (
            f"q6: Q6 SPECIES={atoms} D_0=3.0 R_0=0.5\n"
            "lq6: LOCAL_Q6 SPECIES=q6 SWITCH={RATIONAL D_0=3.0 R_0=0.5} MEAN\n"
            f"PRINT ARG=lq6.mean STRIDE={stride} FILE=COLVAR\n"
        )
    elif kind == "envsim":
        a = 4.0 * BOND / math.sqrt(3.0)
        body = (
            f"es: ENVIRONMENTSIMILARITY SPECIES={atoms} SIGMA=0.5 "
            f"LATTICE_CONSTANTS={a:.6f} CRYSTAL_STRUCTURE=DIAMOND MEAN\n"
            f"PRINT ARG=es.mean STRIDE={stride} FILE=COLVAR\n"
        )
    else:
        raise SystemExit(kind)
    return head + body


def loop_time(log_text):
    for line in log_text.splitlines():
        if line.startswith("Loop time of"):
            # Loop time of 1.23 on 1 procs for 200 steps with 4096 atoms
            return float(line.split()[3]), int(line.split()[8])
    raise RuntimeError("no Loop time in lammps log\n" + log_text[-2000:])


def main():
    lmp, module = sys.argv[1], sys.argv[2]
    sw = pathlib.Path(sys.argv[3]) if len(sys.argv) > 3 else pathlib.Path("/tmp/lammps/potentials/mW.sw")
    pos, box = ice_positions()
    n = len(pos)
    kinds = ["baseline", "dseams-chill", "dseams-rings", "dseams-smooth", "q6", "local-q6", "envsim"]
    strides = [1, 10, 100]
    env = os.environ.copy()
    with tempfile.TemporaryDirectory(prefix="dseams-time-") as tmp:
        tmp = pathlib.Path(tmp)
        write_data(tmp / "data.ice", pos, box)
        shutil.copy(sw, tmp / "mW.sw")
        rows = []
        base = None
        for kind in kinds:
            use_strides = [0] if kind == "baseline" else strides
            for stride in use_strides:
                work = tmp / f"{kind}-{stride}"
                work.mkdir()
                (work / "data.ice").write_bytes((tmp / "data.ice").read_bytes())
                shutil.copy(sw, work / "mW.sw")
                (work / "in.lmp").write_text(lammps_input(kind != "baseline"))
                if kind != "baseline":
                    (work / "plumed.dat").write_text(plumed_text(module, kind, stride, n))
                proc = subprocess.run(
                    [lmp, "-in", "in.lmp", "-log", "log.lammps"],
                    cwd=work, capture_output=True, text=True, env=env,
                )
                log = (work / "log.lammps").read_text() if (work / "log.lammps").is_file() else ""
                if proc.returncode != 0 or "Loop time of" not in log:
                    sys.stderr.write(proc.stdout[-2000:])
                    sys.stderr.write(proc.stderr[-2000:])
                    sys.stderr.write(log[-2000:])
                    raise SystemExit(f"{kind} stride {stride} failed")
                seconds, steps = loop_time(log)
                per = seconds / steps
                if kind == "baseline":
                    base = per
                    overhead = 0.0
                else:
                    overhead = per - base
                rows.append((kind, stride, seconds, per, overhead))
                print(f"{kind:16} stride {stride:3}  loop {seconds:.3f}s  {per*1e3:.3f} ms/step  overhead {overhead*1e3:.3f} ms/step", flush=True)
    print(f"N={n} steps={NSTEPS} baseline {base*1e3:.3f} ms/step")


if __name__ == "__main__":
    main()

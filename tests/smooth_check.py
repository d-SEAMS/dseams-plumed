#!/usr/bin/env python3
"""Finite-difference check of DSEAMS_SMOOTH_CHILL, and the ice / liquid gap.

  smooth_check.py PLUMED_BIN MODULE_SO

The derivative frame is ideal ice Ic with one molecule displaced and a
switching window that still has slope at the hydrogen-bond length, so the
analytic force is not sitting on a flat step. plumed driver --debug-forces
compares that force with a numerical difference.
"""
from __future__ import annotations

import math
import pathlib
import subprocess
import sys
import tempfile

BOND = 2.75


def ice_ic(reps):
    a = 4.0 * BOND / math.sqrt(3.0)
    fcc = [(0.0, 0.0, 0.0), (0.5, 0.5, 0.0), (0.5, 0.0, 0.5), (0.0, 0.5, 0.5)]
    basis = fcc + [(x + 0.25, y + 0.25, z + 0.25) for (x, y, z) in fcc]
    pos = []
    box = reps * a
    for i in range(reps):
        for j in range(reps):
            for k in range(reps):
                for bx, by, bz in basis:
                    pos.append((((i + bx) * a) % box, ((j + by) * a) % box, ((k + bz) * a) % box))
    return pos, box


def write_xyz(path, pos, box, repeat=2):
    lx = ly = lz = box
    with open(path, "w") as out:
        for _ in range(repeat):
            out.write(f"{len(pos)}\n")
            out.write(f"{lx:.8f} 0 0 0 {ly:.8f} 0 0 0 {lz:.8f}\n")
            for x, y, z in pos:
                out.write(f"O {x:.8f} {y:.8f} {z:.8f}\n")


def run_driver(plumed, work, extra=None):
    cmd = [plumed, "driver", "--plumed", "plumed.dat", "--ixyz", "traj.xyz", "--length-units", "A"]
    if extra:
        cmd.extend(extra)
    proc = subprocess.run(cmd, cwd=work, capture_output=True, text=True)
    if proc.returncode != 0:
        sys.stderr.write(proc.stdout[-5000:])
        sys.stderr.write(proc.stderr[-5000:])
        raise SystemExit(f"plumed driver exited {proc.returncode}")
    return proc


def colvar_value(path):
    vals = []
    for line in pathlib.Path(path).read_text().splitlines():
        if line.startswith("#") or not line.strip():
            continue
        vals.append(float(line.split()[1]))
    return vals


def parse_debug_forces(text):
    """Worst relative |analytical - numerical| / scale in a forces.num file.

    Columns are parameter, analytical, numerical. A restraint makes the
    bias nonzero, which is what --debug-forces differentiates.
    """
    worst = 0.0
    seen = 0
    for line in text.splitlines():
        if not line or line.startswith("#") or line.startswith("#!"):
            continue
        parts = line.split()
        if len(parts) < 3:
            continue
        try:
            analytical = float(parts[1])
            numerical = float(parts[2])
        except ValueError:
            continue
        if abs(analytical) < 1e-8 and abs(numerical) < 1e-8:
            continue
        scale = max(abs(analytical), abs(numerical), 1e-8)
        worst = max(worst, abs(analytical - numerical) / scale)
        seen += 1
    return seen, worst


def main():
    plumed, module = sys.argv[1], sys.argv[2]
    pos, box = ice_ic(2)
    moved = list(pos)
    moved[0] = ((moved[0][0] + 0.35) % box, (moved[0][1] - 0.22) % box, (moved[0][2] + 0.18) % box)
    liquidish = [((x + 1.15 * math.sin(3 * i)) % box,
                  (y + 0.95 * math.cos(2 * i)) % box,
                  (z + 1.05 * math.sin(i + 1.7)) % box)
                 for i, (x, y, z) in enumerate(pos)]
    with tempfile.TemporaryDirectory(prefix="dseams-smooth-") as tmp:
        tmp = pathlib.Path(tmp)
        n = len(pos)
        # Value gap: default switch, ideal ice against a strongly disordered copy.
        write_xyz(tmp / "ideal.xyz", pos, box)
        write_xyz(tmp / "liquidish.xyz", liquidish, box)
        dat = (
            "UNITS LENGTH=A\n"
            f"LOAD FILE={module}\n"
            f"ice: DSEAMS_SMOOTH_CHILL ATOMS=1-{n} LENGTH_SCALE=1.0\n"
            "PRINT ARG=ice STRIDE=1 FILE=COLVAR FMT=%16.8f\n"
        )
        (tmp / "plumed.dat").write_text(dat)
        (tmp / "traj.xyz").write_bytes((tmp / "ideal.xyz").read_bytes())
        run_driver(plumed, tmp)
        ideal = colvar_value(tmp / "COLVAR")[-1]
        (tmp / "traj.xyz").write_bytes((tmp / "liquidish.xyz").read_bytes())
        (tmp / "COLVAR").unlink()
        run_driver(plumed, tmp)
        disordered = colvar_value(tmp / "COLVAR")[-1]
        print(f"smooth ideal {ideal:.4f} disordered {disordered:.4f}")
        if ideal < 0.85:
            raise SystemExit("smooth CV is not near 1 on ideal ice Ic")
        if disordered > ideal - 0.4:
            raise SystemExit("smooth CV does not drop on a disordered frame")

        # Derivatives, with the hydrogen bond inside the sloping part of the switch.
        deriv = tmp / "deriv"
        deriv.mkdir()
        write_xyz(deriv / "traj.xyz", moved, box, repeat=1)
        (deriv / "plumed.dat").write_text(
            "UNITS LENGTH=A\n"
            f"LOAD FILE={module}\n"
            f"ice: DSEAMS_SMOOTH_CHILL ATOMS=1-{n} D_0=2.4 D_MAX=3.3 WIDTH=0.15 "
            "BETA=0.4 LENGTH_SCALE=1.0\n"
            "res: RESTRAINT ARG=ice AT=0.2 KAPPA=20.0\n"
            "PRINT ARG=ice,res.bias STRIDE=1 FILE=COLVAR FMT=%16.8f\n"
        )
        proc = run_driver(plumed, deriv, ["--debug-forces", "forces.num"])
        body = (deriv / "forces.num").read_text()
        if "bias for configuration appears to be zero" in (proc.stdout + proc.stderr):
            raise SystemExit("restraint did not apply a force")
        seen, worst = parse_debug_forces(body)
        print(f"debug-forces comparisons {seen} worst relative {worst:.3e}")
        if seen == 0:
            sys.stderr.write(body[-4000:])
            raise SystemExit("debug-forces produced no comparison")
        if worst > 1e-3:
            sys.stderr.write(body[-4000:])
            raise SystemExit(f"analytic derivative mismatch {worst}")

        # A bias accepts the smooth CV. The discrete fraction has no derivative
        # and is refused.
        bias = tmp / "bias"
        bias.mkdir()
        write_xyz(bias / "traj.xyz", pos, box)
        (bias / "plumed.dat").write_text(
            "UNITS LENGTH=A\n"
            f"LOAD FILE={module}\n"
            f"ice: DSEAMS_SMOOTH_CHILL ATOMS=1-{n} LENGTH_SCALE=1.0\n"
            "metad: METAD ARG=ice SIGMA=0.05 HEIGHT=0.1 PACE=1 "
            "GRID_MIN=0 GRID_MAX=1 GRID_BIN=50 FILE=HILLS\n"
            "PRINT ARG=ice,metad.bias STRIDE=1 FILE=COLVAR FMT=%16.8f\n"
        )
        run_driver(plumed, bias)
        biased = colvar_value(bias / "COLVAR")
        print(f"metad on smooth CV ran, value {biased[-1]:.4f}")

        refused = tmp / "refused"
        refused.mkdir()
        write_xyz(refused / "traj.xyz", pos, box)
        (refused / "plumed.dat").write_text(
            "UNITS LENGTH=A\n"
            f"LOAD FILE={module}\n"
            f"ice: DSEAMS_CHILL ATOMS=1-{n} LENGTH_SCALE=1.0\n"
            "metad: METAD ARG=ice.fraction SIGMA=0.05 HEIGHT=0.1 PACE=1 "
            "GRID_MIN=0 GRID_MAX=1 GRID_BIN=50 FILE=HILLS\n"
            "PRINT ARG=ice.fraction STRIDE=1 FILE=COLVAR\n"
        )
        proc = subprocess.run(
            [plumed, "driver", "--plumed", "plumed.dat", "--ixyz", "traj.xyz", "--length-units", "A"],
            cwd=refused, capture_output=True, text=True,
        )
        if proc.returncode == 0:
            sys.stderr.write(proc.stdout[-2000:])
            raise SystemExit("METAD accepted the non-differentiable fraction")
        print("metad refused the discrete fraction")
    print("ok")


if __name__ == "__main__":
    main()

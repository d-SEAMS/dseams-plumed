#!/usr/bin/env python3
"""Generate ice Ic, a liquid and an interface, and check that plumed driver
matches offline d-SEAMS on both CHILL+ rules and on the primitive ring census.

  chill_match.py PLUMED_BIN MODULE_SO OFFLINE_BIN

Frames are written in a temporary directory and are not kept.
"""
from __future__ import annotations

import math
import pathlib
import random
import subprocess
import sys
import tempfile

BOND = 2.75
CUTOFF = 3.5


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


def liquid(n, seed, min_dist=2.35):
    density = 0.0334
    length = (n / density) ** (1.0 / 3.0)
    rng = random.Random(seed)
    pts = []
    misses = 0
    while len(pts) < n:
        trial = (rng.random() * length, rng.random() * length, rng.random() * length)
        if all(_mic_dist(trial, other, length) >= min_dist for other in pts):
            pts.append(trial)
            misses = 0
        else:
            misses += 1
            if misses > 200000:
                raise RuntimeError("liquid packing stalled")
    return pts, length


def _mic_dist(a, b, length):
    d2 = 0.0
    for u, v in zip(a, b):
        d = u - v
        d -= length * round(d / length)
        d2 += d * d
    return math.sqrt(d2)


def interface(reps, seed):
    """Ice Ic slab, a gap inside the cutoff, a hard-sphere liquid, then a
    gap wider than the cutoff so the periodic image does not touch."""
    ice, ice_box = ice_ic(reps)
    # Keep the lateral box of the crystal and give the liquid its own height.
    n_liq = len(ice)
    lx = ly = ice_box
    density = 0.0334
    lz_liq = n_liq / (density * lx * ly)
    rng = random.Random(seed)
    liquid_pts = []
    misses = 0
    while len(liquid_pts) < n_liq:
        trial = (rng.random() * lx, rng.random() * ly, rng.random() * lz_liq)
        ok = True
        for other in liquid_pts:
            dx = trial[0] - other[0]
            dy = trial[1] - other[1]
            dz = trial[2] - other[2]
            dx -= lx * round(dx / lx)
            dy -= ly * round(dy / ly)
            dz -= lz_liq * round(dz / lz_liq)
            if dx * dx + dy * dy + dz * dz < 2.35 ** 2:
                ok = False
                break
        if ok:
            liquid_pts.append(trial)
            misses = 0
        else:
            misses += 1
            if misses > 400000:
                raise RuntimeError("interface liquid stalled")
    gap_contact = 2.0
    gap_isolate = CUTOFF + 0.6
    pos = list(ice)
    z0 = ice_box + gap_contact
    for x, y, z in liquid_pts:
        pos.append((x, y, z + z0))
    box_z = z0 + lz_liq + gap_isolate
    return pos, (lx, ly, box_z)


def write_xyz(path, frames):
    with open(path, "w") as out:
        for box, pos in frames:
            if isinstance(box, tuple):
                lx, ly, lz = box
            else:
                lx = ly = lz = box
            out.write(f"{len(pos)}\n")
            out.write(f"{lx:.8f} 0 0 0 {ly:.8f} 0 0 0 {lz:.8f}\n")
            for x, y, z in pos:
                out.write(f"O {x:.8f} {y:.8f} {z:.8f}\n")


def plumed_input(module, n):
    atoms = f"1-{n}"
    return (
        "UNITS LENGTH=A\n"
        f"LOAD FILE={module}\n"
        f"ice: DSEAMS_CHILL ATOMS={atoms} CUTOFF=3.5 LENGTH_SCALE=1.0 RULE=chill+\n"
        "PRINT ARG=ice.fraction,ice.ncubic,ice.nhex,ice.ninterfacial,ice.nclathrate,"
        "ice.ninterclath,ice.nwater,ice.nmax STRIDE=1 FILE=CHILL FMT=%16.8f\n"
        f"ov: DSEAMS_CHILL ATOMS={atoms} CUTOFF=3.5 LENGTH_SCALE=1.0 RULE=ovito-chill+\n"
        "PRINT ARG=ov.fraction,ov.ncubic,ov.nhex,ov.ninterfacial,ov.nclathrate,"
        "ov.ninterclath,ov.nwater,ov.nmax STRIDE=1 FILE=OVITO FMT=%16.8f\n"
        f"rings: DSEAMS_RINGS ATOMS={atoms} CUTOFF=3.5 MAX=6 LENGTH_SCALE=1.0\n"
        "PRINT ARG=rings.n3,rings.n4,rings.n5,rings.n6,rings.n7,rings.n8 "
        "STRIDE=1 FILE=RINGS FMT=%8.0f\n"
    )


def data_rows(path):
    rows = []
    for line in pathlib.Path(path).read_text().splitlines():
        if line.startswith("#") or not line.strip():
            continue
        rows.append([float(x) for x in line.split()])
    return rows


def close(got, want, tol, label):
    if abs(got - want) > tol:
        raise SystemExit(f"{label}: driver {got} offline {want}")


def main():
    plumed, module, offline = sys.argv[1], sys.argv[2], sys.argv[3]
    ice, ice_box = ice_ic(3)
    liq, liq_box = liquid(192, seed=11)
    face, face_box = interface(3, seed=19)
    systems = {
        "ice": [(ice_box, ice)],
        "liquid": [(liq_box, liq)],
        "interface": [(face_box, face)],
    }
    with tempfile.TemporaryDirectory(prefix="dseams-chill-") as tmp:
        tmp = pathlib.Path(tmp)
        for name, frames in systems.items():
            xyz = tmp / f"{name}.xyz"
            write_xyz(xyz, frames)
            ref = subprocess.run([offline, str(xyz), "--max-ring", "6"], check=True,
                                 capture_output=True, text=True)
            offline_rows = [list(map(float, line.split())) for line in ref.stdout.splitlines() if line.strip()]
            work = tmp / name
            work.mkdir()
            n = len(frames[0][1])
            (work / "plumed.dat").write_text(plumed_input(module, n))
            (work / "traj.xyz").write_bytes(xyz.read_bytes())
            driver = subprocess.run(
                [plumed, "driver", "--plumed", "plumed.dat", "--ixyz", "traj.xyz", "--length-units", "A"],
                cwd=work, capture_output=True, text=True,
            )
            if driver.returncode != 0:
                sys.stderr.write(driver.stdout[-4000:])
                sys.stderr.write(driver.stderr[-4000:])
                raise SystemExit(f"plumed driver failed on {name}")
            chill = data_rows(work / "CHILL")
            ovito = data_rows(work / "OVITO")
            rings = data_rows(work / "RINGS")
            if not (len(chill) == len(ovito) == len(rings) == len(offline_rows)):
                raise SystemExit(f"{name}: frame count mismatch")
            for row, ovi, ring, off in zip(chill, ovito, rings, offline_rows):
                # offline: step nop pub_f pub_c pub_h pub_i pub_cl pub_ic pub_w pub_max
                #          ovi_f ... ovi_max n3 n4 n5 n6 n7 n8
                close(row[1], off[2], 1e-6, f"{name} published fraction")
                for k, label in enumerate(("ncubic", "nhex", "ninterfacial", "nclathrate",
                                           "ninterclath", "nwater", "nmax")):
                    close(row[2 + k], off[3 + k], 1e-6, f"{name} published {label}")
                close(ovi[1], off[10], 1e-6, f"{name} ovito fraction")
                for k, label in enumerate(("ncubic", "nhex", "ninterfacial", "nclathrate",
                                           "ninterclath", "nwater", "nmax")):
                    close(ovi[2 + k], off[11 + k], 1e-6, f"{name} ovito {label}")
                for k, label in enumerate(("n3", "n4", "n5", "n6", "n7", "n8")):
                    close(ring[1 + k], off[18 + k], 1e-6, f"{name} {label}")
            print(
                f"{name}: N={n} published fraction {chill[-1][1]:.4f} "
                f"ovito fraction {ovito[-1][1]:.4f} n6 {rings[-1][4]:.0f}"
            )
        ylm = subprocess.run([offline, str(tmp / "ice.xyz"), "--check-ylm", "--check-correl"],
                             check=False, capture_output=True, text=True)
        sys.stdout.write(ylm.stdout)
        sys.stderr.write(ylm.stderr)
        if ylm.returncode != 0:
            raise SystemExit("q3 correlation check failed")
        ice_fraction = data_rows(tmp / "ice" / "CHILL")[-1][1]
        liq_fraction = data_rows(tmp / "liquid" / "CHILL")[-1][1]
        face_fraction = data_rows(tmp / "interface" / "CHILL")[-1][1]
        if ice_fraction != 1.0:
            raise SystemExit(f"ideal ice Ic fraction is {ice_fraction}, want 1")
        if liq_fraction > 0.05:
            raise SystemExit(f"liquid fraction {liq_fraction} is not a liquid")
        if not (liq_fraction < face_fraction < ice_fraction):
            raise SystemExit(
                f"interface fraction {face_fraction} is not between liquid {liq_fraction} and ice"
            )
        ice_rings = data_rows(tmp / "ice" / "RINGS")[-1]
        if ice_rings[1] != 0 or ice_rings[2] != 0 or ice_rings[3] != 0 or ice_rings[4] <= 0:
            raise SystemExit(f"ideal ice Ic ring census {ice_rings[1:]} is not pure six-rings")
    print("ok")


if __name__ == "__main__":
    main()

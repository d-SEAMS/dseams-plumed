"""Fifty steps of TIP4P-Ew with the smooth ice fraction.

Requires OpenMM and the openmm-plumed plugin (PlumedForce). Point
PLUMED_KERNEL at libplumedKernel.so and put that library, the module's
dependencies, and libOpenMMPlumed on LD_LIBRARY_PATH. No trajectory is
written.
"""

import ctypes
import os
import subprocess
from pathlib import Path

from openmm import LangevinMiddleIntegrator, Platform, Vec3, unit
from openmm.app import CutoffPeriodic, ForceField, Modeller, Simulation, Topology
from openmmplumed import PlumedForce

HERE = Path(__file__).resolve().parent
MODULE = Path(os.environ.get("DSEAMS_MODULE", HERE.parents[1] / "build" / "libdseams_plumed.so"))


def attach(system, force):
    """Add a PlumedForce. Fall back to the C++ call when SWIG will not cast it."""
    try:
        system.addForce(force)
        return
    except TypeError:
        pass
    helper = HERE / "addforce.so"
    if not helper.is_file():
        include = os.environ.get("OPENMM_INCLUDE")
        lib = os.environ.get("OPENMM_LIB")
        if not include or not lib:
            raise SystemExit("system.addForce rejected PlumedForce; set OPENMM_INCLUDE and OPENMM_LIB")
        cxx = os.environ.get("CXX", "g++")
        subprocess.check_call(
            [
                cxx,
                "-shared",
                "-fPIC",
                "-O2",
                f"-I{include}",
                f"-L{lib}",
                "-lOpenMM",
                str(HERE / "addforce.cpp"),
                "-o",
                str(helper),
            ]
        )
    blob = ctypes.CDLL(str(helper))
    blob.dseams_openmm_add_force.argtypes = [ctypes.c_void_p, ctypes.c_void_p]
    force.thisown = 0
    blob.dseams_openmm_add_force(ctypes.c_void_p(int(system.this)), ctypes.c_void_p(int(force.this)))


def main():
    if not MODULE.is_file():
        raise SystemExit(f"module not found: {MODULE}")
    forcefield = ForceField("tip4pew.xml")
    modeller = Modeller(Topology(), [])
    modeller.addSolvent(
        forcefield,
        boxSize=Vec3(1.9, 1.9, 1.9) * unit.nanometer,
        model="tip4pew",
    )
    system = forcefield.createSystem(
        modeller.topology,
        nonbondedMethod=CutoffPeriodic,
        nonbondedCutoff=0.9 * unit.nanometer,
        rigidWater=True,
        constraints=None,
    )
    n_mol = modeller.topology.getNumResidues()
    colvar = HERE / "COLVAR"
    script = (
        "UNITS LENGTH=nm\n"
        f"LOAD FILE={MODULE}\n"
        f"ice: DSEAMS_SMOOTH_CHILL ATOMS=1-{4 * n_mol}:4 D_0=3.2 D_MAX=3.5 LENGTH_SCALE=10.0\n"
        f"PRINT ARG=ice STRIDE=10 FILE={colvar}\n"
    )
    attach(system, PlumedForce(script))
    integrator = LangevinMiddleIntegrator(
        300 * unit.kelvin, 1.0 / unit.picosecond, 0.002 * unit.picoseconds
    )
    platform = Platform.getPlatformByName("CPU")
    simulation = Simulation(modeller.topology, system, integrator, platform)
    simulation.context.setPositions(modeller.positions)
    simulation.context.setVelocitiesToTemperature(300 * unit.kelvin, 1)
    print(f"molecules {n_mol}")
    simulation.step(50)
    if not colvar.is_file():
        raise SystemExit("PLUMED did not write COLVAR")
    lines = [ln for ln in colvar.read_text().splitlines() if ln.strip() and not ln.startswith("#")]
    print("first", lines[0])
    print("last", lines[-1])
    colvar.unlink()


if __name__ == "__main__":
    main()

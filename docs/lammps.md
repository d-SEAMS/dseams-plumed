# Inside LAMMPS

Use `fix plumed` to evaluate `DSEAMS_CAGES` on the fly during a LAMMPS run.
Because `DSEAMS_CAGES` registers integer counts with no derivatives, it adds
no bias forces and does not perturb the integrator, thermostat, or barostat.

## Fix ordering

Place `fix plumed` **before** `fix npt` or `fix press/berendsen`. Those
barostat fixes compute the pressure internally, and LAMMPS refuses a
reversed order.

A minimal mW nucleation input logs the cage score every 500 steps and stops
when the largest ice cluster enters a committor basin:

```lammps
fix pl  all plumed plumedfile plumed.dat outfile plumed.log
fix npt all npt temp 215 215 1000 iso 1 1 10000
```

## Matching `plumed.dat`

LAMMPS `units real` and `units metal` both pass lengths in Ångströms, which
matches the engine's internal unit. Set `UNITS LENGTH=A` and
`LENGTH_SCALE=1.0` in `plumed.dat`:

```text
UNITS LENGTH=A
LOAD FILE=./build/libdseams_plumed.so

ice: DSEAMS_CAGES ...
  ATOMS=1-4096
  CUTOFF=3.5
  CANDIDATE=5.5
  K=4
  LENGTH_SCALE=1.0
  COMPLETE
...

PRINT ...
  ARG=ice.nice,ice.nmax,ice.nic,ice.nih,ice.chillice,ice.chillmax
  STRIDE=500
  FILE=ICE
...

COMMITTOR ...
  ARG=ice.nmax
  STRIDE=500
  BASIN_LL1=0 BASIN_UL1=20
  BASIN_LL2=800 BASIN_UL2=100000
...
```

When running a simulation engine that passes coordinates in nanometres (such
as GROMACS), omit `UNITS LENGTH=A` and `LENGTH_SCALE=1.0` so the default
`LENGTH_SCALE=10.0` converts nanometres to Ångströms.

## Multi-site water models

`ATOMS` must list only the oxygen atoms (one atom per water molecule). For
models with explicit hydrogens or dummy sites, select the oxygen indices
with PLUMED's stride syntax. For example, in a 3000-molecule TIP4P/2005
system where each molecule occupies four consecutive IDs starting with
oxygen (`O, H, H, M`), followed by 60 ions:

```text
UNITS LENGTH=A
LOAD FILE=./build/libdseams_plumed.so

brine: DSEAMS_CAGES ...
  ATOMS=1-12000:4
  IONS=12001-12060
  ION_CUTOFF=3.5
  LENGTH_SCALE=1.0
  COMPLETE
...

PRINT ...
  ARG=brine.nice,brine.nmax,brine.nionice,brine.nionfront,brine.nionliq
  STRIDE=500
  FILE=BRINE
...
```

## Threading

In PLUMED 2.10 and newer, control PLUMED's thread team with
`PLUMED_NUM_THREADS` before launching LAMMPS.

## Bundled multi-action example

The PLUMED-NEST egg
([plumID:26.018](https://www.plumed-nest.org/eggs/26/018/), `v0.4.1`;
`nest/` in the repository) includes a complete `plumed.dat` with one
`DSEAMS_CAGES` action each for ice nucleation in water (`ICE`), a growing
front in NaCl brine (`BRINE`), and filled cages in a periodic cage network
(`HYDRATE`).

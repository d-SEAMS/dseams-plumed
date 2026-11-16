# LAMMPS: bias the smooth ice fraction

64 mW molecules in a periodic cell, a distorted ice Ic lattice. `fix plumed`
applies a restraint that pulls `DSEAMS_SMOOTH_CHILL` toward 0.4. Nothing
written here is a trajectory; `COLVAR` and `log.lammps` are produced by the run.

```
lmp -in in.lammps
```

`plumed.dat` loads `../../build/libdseams_plumed.so`. Change that path if the
module lives somewhere else. The mW pair style reads `mW.sw` (Molinero and
Moore, J. Phys. Chem. B 113, 4008, 2009).

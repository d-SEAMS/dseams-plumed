# GROMACS: bias the smooth ice fraction

32 TIP4P/Ice molecules. `gmx mdrun -plumed plumed.dat` restrains
`DSEAMS_SMOOTH_CHILL` toward 0.3. The oxygen of each molecule is the first of
four sites, so the group is `1-128:4`. Coordinates in `conf.gro` are the
starting box, not a trajectory.

```
gmx grompp -f md.mdp -c conf.gro -p topol.top -o run.tpr
gmx mdrun -deffnm run -plumed plumed.dat -nsteps 20
```

GROMACS has to be built with PLUMED (`-DGMX_USE_PLUMED=ON` on 2025). Point
`LOAD` in `plumed.dat` at `libdseams_plumed.so`.

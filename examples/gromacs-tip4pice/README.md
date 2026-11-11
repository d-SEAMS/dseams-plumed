# GROMACS: grow ice on a TIP4P/Ice slab

`conf.gro` is a 432-molecule ice Ih / liquid slab. GenIce built the
crystal (`genice2 1h --rep 3 3 3 --water tip4p --seed 1`). The lower
half was held in place and the upper half was heated, so the start is
about 40 percent ice. `D_0` and `D_MAX` are in Angstrom;
`LENGTH_SCALE=10` converts GROMACS nanometres.

```
gmx grompp -f md.mdp -c conf.gro -p topol.top -o run.tpr
gmx mdrun -deffnm run -plumed plumed.dat
```

Point `LOAD` at `libdseams_plumed.so`. GROMACS has to be built with
PLUMED. The run is 100 ps at 210 K and 1 bar, with a restraint of the
smooth fraction toward 0.90 (`KAPPA` is kJ/mol). No trajectory is
stored. One such run:

| time / ps | smooth fraction | restraint / kJ mol^-1 | published fraction |
| ---: | ---: | ---: | ---: |
| 0 | 0.427 | 2240 | 0.403 |
| 1 | 0.767 | 178 | 0.766 |
| 10 | 0.859 | 17.1 | 0.919 |
| 20 | 0.850 | 24.8 | 0.891 |
| 50 | 0.886 | 2.07 | 0.970 |
| 100 | 0.880 | 3.83 | 0.975 |

The published CHILL+ fraction is monitored. It has no derivative. The
restraint acts on the smooth fraction.

`plumed-metad.dat` and `md-metad.mdp` are 2 ns of well-tempered
metadynamics on the same slab at 210 K (height 1.5 kJ/mol, sigma 0.03,
bias factor 15, a hill every 1 ps). `fes.dat` is the free energy from
`plumed sum_hills`, shifted so the minimum is zero.

| smooth fraction | free energy / kJ mol^-1 |
| ---: | ---: |
| 0.35 | 107 |
| 0.46 | 97 |
| 0.58 | 61 |
| 0.69 | 34 |
| 0.81 | 22 |
| 0.92 | 0 |
| 0.98 | 60 |

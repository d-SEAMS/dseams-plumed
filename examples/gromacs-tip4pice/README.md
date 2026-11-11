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

`plumed-metad.dat` and `md-metad.mdp` are well-tempered metadynamics
on the same slab (height 1.5 kJ/mol, sigma 0.03, bias factor 15, a
hill every 1 ps). `fes.dat` is the 2 ns free energy from
`plumed sum_hills` at 210 K, shifted so the minimum is zero.

The minimum moves as the hills accumulate. Each column below is a
`sum_hills` profile shifted so its own minimum is zero. The last two
columns are a block estimate: the 2 ns trajectory split into four
500 ps blocks, reweighted by exp(V/kT), with the standard error of
the mean across blocks that visited the bin.

| smooth fraction | 0.5 ns | 1 ns | 2 ns | block mean | stderr |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 0.40 | 72 | 81 | 106 |  |  |
| 0.70 | 0 | 9 | 33 | 0.7 |  |
| 0.80 | 51 | 4 | 24 | 20 | 11 |
| 0.84 | 69 | 0 | 13 | 7 | 3 |
| 0.92 |  | 51 | 0 | 0.6 |  |
| 0.98 |  | 81 | 54 |  |  |

The minimum sits at 0.69 after 0.5 ns, 0.83 after 1 ns and 0.92 after
2 ns. The 1.0-1.5 ns block and the 1.5-2.0 ns block both put it near
0.91. `fes-conv.dat` and `docs/smooth-fes-conv.png` are the full curves.

A second walk of the same slab at 330 K, above the TIP4P/Ice melting
point, is the liquid control. After 1.5 ns the collective variable
sits at 0.001 and the `sum_hills` minimum is at fraction 0.005
(`fes-330K.dat`). The ice well is gone.

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

## Eight walkers, 28 ns at 210 K

The single-walker minimum was still moving at 2 ns. The follow-up is
eight walkers sharing hills, bias factor 200, sigma 0.06, height
1.5 kJ/mol, one hill every 0.5 ps, 3.5 ns each (28 ns in total). The
starts are spread from fraction 0.17 to 0.97. A fully liquid box was
also run; the smooth fraction has no gradient there, so that walker
stays at zero and is left out of the profile below.

`fes-mw.dat` and `docs/fes-mw.png` are the free energy with the
minimum set to zero. The grey band is the range of the profiles at
3.0, 3.2, 3.35 and 3.5 ns.

| smooth fraction | F at 3.5 ns / kJ mol^-1 | range over 3.0–3.5 ns |
| ---: | ---: | ---: |
| 0.20 | 277 | 4.2 |
| 0.40 | 188 | 3.5 |
| 0.60 | 114 | 5.3 |
| 0.80 | 47 | 3.6 |
| 0.90 | 14 | 1.7 |
| 0.930 | 0 | 0 |

The only minimum is at fraction 0.930. The free energy falls
monotonically from low fraction to that minimum, so the uphill cost
of leaving the ice basin is the ice-liquid difference and there is
no separate barrier on this coordinate. From fraction 0.40 to the
minimum the difference is 188 kJ/mol for the 432-molecule slab
(0.44 kJ/mol per molecule). The minimum does not move across the last
0.5 ns. The free energies above it still change by 2 to 5 kJ/mol
over that window, so the 2 kJ/mol target is met at the minimum and
at fraction 0.90, and the residual on the slope is the range in the
table.

A second ensemble with the same bias factor and sigma 0.03, four
walkers, 1.2 ns, puts the minimum at 0.944. At the same 1.2 ns the
sigma 0.06 run has its minimum at 0.929. The narrower kernel leaves
the wing about 50 to 80 kJ/mol higher, because those hills do not
overlap the gaps between walkers as quickly.

## Smooth fraction against Q6, equal wall time

Both start from the same slab (published fraction 0.403) with one
thread and a restraint toward the ice value of that collective
variable. The initial force, kappa times the distance to the target,
is 4000 kJ/mol per unit of the coordinate. The wall clock is 603 s,
which is 0.167 CPU-hour.

| collective variable | simulated time | published fraction at the end | time to published 0.8 |
| --- | ---: | ---: | ---: |
| smooth fraction | 1.437 ns | 0.891 | 3 ps |
| Q6 MEAN | 5 ps | 0.623 | not reached |

The smooth run delivers 207 ns/day and Q6 delivers 0.844 ns/day.
Per CPU-hour the smooth fraction crosses the freezing threshold and
Q6 does not. A static ENVIRONMENTSIMILARITY to cubic diamond is 0.127
on the ice end of this slab and 0.248 on the liquid, so that
coordinate points the wrong way for ice Ih and was not used as a bias.

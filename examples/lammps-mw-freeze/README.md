# mW: drive freezing with the smooth ice fraction

`data.seed` is 512 mW particles. It starts as ice Ic. A sphere of radius
9 A about the centre (99 particles) is held fixed while the rest is
heated to 480 K for 15 ps, which leaves a crystalline core in liquid
water. The smooth fraction on that frame is 0.134 and the published
CHILL+ fraction is 0.074.

`in.smooth` is 20 ps of NPT at 210 K and 1 atm. A moving restraint takes
the smooth fraction from 0.184 toward 0.95. `KAPPA` is kcal/mol.
`fix plumed` is before `fix npt`. One run:

| time / ps | smooth fraction | restraint / kcal mol^-1 | published fraction |
| ---: | ---: | ---: | ---: |
| 0 | 0.134 | 7.4 | 0.074 |
| 2 | 0.144 | 41 | 0.078 |
| 6 | 0.208 | 127 | 0.133 |
| 10 | 0.287 | 236 | 0.148 |
| 14 | 0.379 | 349 | 0.252 |
| 20 | 0.484 | 652 | 0.320 |

The published count is only printed. The restraint acts on the smooth
fraction.

Wall time per step on one thread, same 2000 steps, from the LAMMPS
loop timer. The smooth row has no extra printed count.

| run | ms/step |
| --- | ---: |
| no collective variable | 0.779 |
| smooth fraction, moving restraint | 1.49 |
| ENVIRONMENTSIMILARITY diamond, moving restraint | 81.1 |
| Q6 MEAN, moving restraint | 323 |

`in.plain`, `in.q6` and `in.envsim` are those three runs. No
trajectories are stored.

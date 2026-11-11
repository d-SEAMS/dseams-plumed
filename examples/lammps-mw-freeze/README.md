# mW: drive freezing with the smooth ice fraction

`data.seed` is 512 mW particles. It starts as ice Ic. A sphere of radius
9 A about the centre (99 particles) is held fixed while the rest is
heated to 480 K for 15 ps, which leaves a crystalline core in liquid
water. The smooth fraction on that frame is 0.075 and the published
CHILL+ fraction is 0.074.

`in.smooth` is NPT at 230 K and 1 atm with a 10 fs step. A moving
restraint takes the smooth fraction from 0.20 to 0.98 over 80 ps
(`KAPPA` is kcal/mol) and then holds it. `fix plumed` is before
`fix npt`. Three velocity seeds, same coordinates:

| velocity seed | time to published fraction 0.9 | final published fraction |
| ---: | ---: | ---: |
| 2 | 64 ps | 1.00 |
| 3 | 62 ps | 0.95 |
| 4 | 54 ps | 0.98 |

The published count is only printed. The restraint acts on the smooth
fraction. No trajectories are stored.

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

## Liquid slope

The published count is flat in a fully liquid shell, so a restraint
on that count alone does not nucleate. The smooth score adds two
terms that are already one on a frozen molecule and still sloping in
the liquid. `TAIL_C=-0.65` and `TAIL_W=0.20` put a second staggered
weight between a liquid correlation and the published -0.8 edge, and
the weight is still climbing at that edge, so four saturated soft
bonds do not by themselves reach 1 (`TAIL_N=4`). `COORD_W=0.45` is a
Gaussian in the smooth neighbour count, equal to 1 at the published
shell of four.

On the ideal ice Ic frame the score is 0.99544 against a published
fraction of 1, residual -0.0046. On the thermal TIP4P/Ice crystal the
score is 0.937 against a published fraction of 0.991, residual
-0.054. The root mean square of those two frozen residuals is 0.038.
The hot liquid scores 0.096, with published fraction 0. A restraint
from that liquid expands the box toward the ice density and raises
the published fraction to about 0.04; it does not cross into the
crystal basin on a 100 ps timescale.

## Eight walkers, 24 ns at 210 K, liquid included

Eight walkers share hills. Bias factor 200, sigma 0.06, height
1.5 kJ/mol, one hill every 1 ps, 3 ns each (24 ns in total). One
walker starts from the hot liquid and moves through 0.02–0.24. The
others start from restrained intermediates and from the crystal. The
water model is TIP4P/Ice.

`fes-mw.dat` and `docs/fes-mw.png` are the free energy. Each column
is zeroed at that column's minimum for fraction above 0.85. The band
on the figure is the range of the 2.4, 2.7 and 3.0 ns profiles. No
walker entered 0.34–0.41, and a restraint held at the low edge of
that interval did not cross it in 40 ps, so the height of a barrier
there is not measured. The figure leaves the interval blank. The
rows of `fes-mw.dat` inside it are only the overlap of the two
kernels.

| smooth fraction | F at 3.0 ns / kJ mol^-1 | range over 2.4–3.0 ns | per molecule / kJ mol^-1 |
| ---: | ---: | ---: | ---: |
| 0.10 | 129 | 12 | 0.30 |
| 0.20 | 57 | 21 | 0.13 |
| 0.28 | 22 | 11 | 0.05 |
| 0.45 | 277 | 5 | 0.64 |
| 0.70 | 117 | 1 | 0.27 |
| 0.90 | 12 | 1 | 0.028 |
| 0.935 | 0 | 0 | 0 |

The ice minimum stays at 0.935 from 2.4 ns onward. The low side has a
well near 0.28. At 3.0 ns that well is 22 kJ/mol above the ice
minimum. F(0.26) still moved by 16 kJ/mol over the preceding 0.6 ns,
so the depth of the low well is not settled. F(0.45), on the crystal
side of the empty interval, moved by 5 kJ/mol in that window.

Per molecule of the 432-molecule box, the hot-liquid value (fraction
0.10) sits 0.30 kJ/mol above the ice minimum, with a block range of
0.03 kJ/mol. TIP4P/Ice (Abascal, Sanz, García Fernández and Vega,
J. Chem. Phys. 122, 234511, 2005) melts at 272.2 K with an enthalpy of
fusion of 1.29 kcal/mol (5.40 kJ/mol). Neglecting the heat-capacity
difference, the bulk driving force at 210 K is about
5.40*(1-210/272.2) = 1.23 kJ/mol per molecule. The 0.30 kJ/mol is
smaller because fraction 0.10 to about 0.30 is mostly the soft tail:
a restrained frame at 0.33 still has a published fraction near 0.04.
From fraction 0.45, where a restrained frame has published fraction
0.38, to the minimum is 277 kJ/mol for the box, 0.64 kJ/mol per
molecule. Spread over the molecules that join the crystal between
those two frames (published fraction 0.38 to 0.99) that is about
1.05 kJ/mol per molecule that freezes, the same order as the 1.23
kJ/mol estimate. The slab at fraction 0.45 is an interface, not a
bulk liquid, and the heat-capacity term is not in the estimate.

## Equal wall time, best setup of each coordinate

Each run starts from the same slab (published fraction 0.403) and is
stopped after 300 s of wall time. The restraint targets the thermal
crystal value of that coordinate, with the initial force
kappa*|target - slab| equal to 4000 kJ/mol. One job runs at a time.

The smooth target is 0.937 (kappa 9313). Q6 uses
`SWITCH={CUBIC D_0=0.32 D_MAX=0.35}`, `NL_CUTOFF=0.50`, `NL_STRIDE=10`
and `MEAN`. The target is 0.584 (kappa 250000). `LOWMEM` is accepted
and ignored: this PLUMED warns that the flag is deprecated and no
longer changes the calculation. ENVIRONMENTSIMILARITY uses
`SIGMA=0.03` and four oxygen shells taken from this crystal
(`iceh/env1.pdb` through `env4.pdb`), which is hexagonal ice, not
cubic diamond. The slab value is 0.340 and the crystal value is
0.399 (kappa 67797). On the hot liquid the same coordinate is 0.096,
so it points toward this ice.

| coordinate | threads | ns/day | simulated | published at the end | time to published 0.8 | ns per CPU-hour |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| smooth | 1 | 194 | 0.674 ns | 0.896 | 2 ps | 8.09 |
| Q6 | 1 | 7.20 | 25 ps | 0.722 | not reached | 0.30 |
| environment similarity | 1 | 2.01 | 7 ps | 0.769 | not reached | 0.084 |
| smooth | 8 | 535 | 1.86 ns | 0.928 | 12 ps | 2.79 |
| Q6 | 8 | 49.5 | 172 ps | 0.889 | 29 ps | 0.26 |
| environment similarity | 8 | 15.0 | 52 ps | 0.720 | not reached | 0.078 |

A CPU-hour is wall time times the thread count. At one thread the
smooth fraction is 27 times the Q6 throughput and is the only run
that reaches a published fraction of 0.8. At eight threads Q6 does
reach 0.8, at 29 ps, and the throughput ratio falls to 11. The smooth
fraction still produces more simulated time per CPU-hour at either
width. One thread is the better smooth rate per CPU-hour (8.1 ns
against 2.8). Environment similarity to this hexagonal ice moves the
published fraction in the right direction and is the most expensive
of the three.

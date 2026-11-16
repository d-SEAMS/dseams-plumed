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

## Liquid slope, capped

The published count is flat in a fully liquid shell. A second
staggered weight (`TAIL_C=-0.65`, `TAIL_W=0.20`, `TAIL_N=4`) still
has a slope there. Without a ceiling that weight let a liquid with
no crystal outscore a seeded nucleus: the same restrained frame
scored 0.34 at a published fraction of 0.04, while every seeded
frame with a published fraction above 0.20 scored at most 0.29. The
free-energy gap at 0.34–0.41 was that jump. `TAIL_AMP=0.65` caps the
soft channel. `BETA=0.06` and `COORD_W=0.70` keep the crystal anchor.

On ideal ice Ic the score is 0.99851 against a published fraction of
1, residual -0.0015. On the thermal crystal it is 0.932 against
0.991, residual -0.059. The root mean square of those two frozen
residuals is 0.042. The hot liquid scores 0.083. The restrained
frame that used to score 0.34 now scores 0.233, below the seeded
frames with a published fraction above 0.20.

A clash-free spherical seed, 0 to 100% of the 432 molecules kept
from the thermal crystal, was scored outside this tree. The smooth
score and the published fraction on that series have correlation
0.989. A seed that keeps 85% of the crystal scores 0.398 at a
published fraction of 0.306, inside the old gap, and the series is
ordered through it.

## Eight MPI walkers, shared OPES bias, 13 ns

PLUMED 2.10.2 and GROMACS 2025.2 were built against the same Open
MPI. `OPES_METAD` uses `WALKERS_MPI`, `SIGMA=0.06`, `BARRIER=300`,
`PACE=500`, `TEMP=210`. `METAD` in the same binary accepts
`WALKERS_MPI` as well. Eight ranks, one thread each, started at
smooth scores 0.083, 0.161, 0.232, 0.309, 0.398, 0.531, 0.676 and
0.932. The run is 13 ns. The eight COLVAR files (time, score, bias,
rct) stay outside git. `scripts/fes_blocks.py` rebuilds the block
table from them.

The reweighted free energy on the last 20% (10.4–13.0 ns) is split
into five equal blocks. A bin needs at least 10 samples in a block.
Each column is zeroed at its own minimum.
Twelve bins are present in every block. Eight of those twelve move
by less than 1 kJ/mol. The four that do not are fraction 0.250
(span 1.79 kJ/mol), 0.370 (1.74), 0.890 (1.68) and 0.410 (1.10).
The median span on the twelve bins is 0.53 kJ/mol. Bins that a block
does not visit are left blank; the walkers do not put 10 samples in
every 0.02-wide bin in every 520 ps block, so the grid from 0.10 to
0.98 is not fully covered.

The ice well (lowest free energy in 0.90–0.98) minus the low well
(lowest in 0.12–0.25) in the five blocks is -2.82, -2.94, -2.81,
-3.19 and -3.23 kJ/mol. The range is 0.42 kJ/mol, under 1 kJ/mol.
The mean is -3.00 kJ/mol for the box, -0.0069 kJ/mol per molecule.
The ice well is the lower of the two.

TIP4P/Ice (Abascal, Sanz, García Fernández and Vega, J. Chem. Phys.
122, 234511, 2005) melts at 272.2 K with an enthalpy of fusion of
1.29 kcal/mol (5.40 kJ/mol). Neglecting the heat-capacity difference,
the bulk driving force at 210 K is about 1.23 kJ/mol per molecule.
The collective-variable offset is about 180 times smaller. Three
things sit between those numbers. The box has 432 molecules, so a
single global fraction cannot separate bulk ice, bulk liquid and the
interface; the free energy at a given score is the minimum over all
structures with that score, and an interface-sized fluctuation costs
far less than 432 times the bulk value. The low well is still a
nearly unlabeled liquid (the seeded series puts a published fraction
near 0.02 at smooth 0.19), not a bulk phase boundary. And the soft
channel, even capped, still gives the liquid a slope, so part of the
motion in the score is local bond rearrangement rather than a change
of phase.

## Equal wall time, capped score against Q6

Both start from the same slab. The smooth target is 0.932 with kappa
8909, so the initial force is 4000 kJ/mol. Q6 uses
`SWITCH={CUBIC D_0=0.32 D_MAX=0.35}`, `NL_CUTOFF=0.50` and
`NL_STRIDE=10`. Its target is 0.584 with kappa 250000, the same
initial force. The neighbour list is active: the log says the cutoff
is 0.50 nm and the list is rebuilt every 10 steps. One job at a time,
300 s of wall time, `mpirun -np 1`.

| coordinate | threads | ns/day | simulated |
| --- | ---: | ---: | ---: |
| smooth, capped | 1 | 183 | 0.640 ns |
| Q6 with a neighbour list | 1 | 7.44 | 26 ps |
| smooth, capped | 8 | 542 | 1.90 ns |
| Q6 with a neighbour list | 8 | 43.2 | 151 ps |

At one thread the capped score is 25 times the Q6 throughput. At
eight threads the ratio is 13. One thread gives 7.6 ns per CPU-hour
for the smooth score and 0.31 for Q6. Eight threads give 2.8 and
0.23 ns per CPU-hour.

## Equal wall time on the uncapped tail

The table below is the previous collective variable, with
`TAIL_AMP=1`. It is not a timing of the capped score.
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

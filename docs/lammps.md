# Inside LAMMPS

Place `fix plumed` before `fix npt` or `fix press/berendsen`. Those
fixes compute the pressure internally. LAMMPS refuses a reversed
order. A minimal mW run logs the score every 500 steps. It stops on
a committor basin:

```text
fix pl  all plumed plumedfile plumed.dat outfile plumed.log
fix npt all npt temp 215 215 1000 iso 1 1 10000
```

Use a `plumed.dat` as on the [home page](index.md) and `LENGTH_SCALE=1.0`
after `UNITS LENGTH=A`. LAMMPS `units real` matches the A unit of the
engine.

The PLUMED-NEST egg under `nest/` carries a `plumed.dat` with one action
each for ice nucleation in water, a growing front in an NaCl brine, and
the filled cages of a methane hydrate.

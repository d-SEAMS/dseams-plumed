# OpenMM: smooth ice fraction on TIP4P-Ew

`run.py` fills a 1.9 nm box with TIP4P-Ew, loads the module through
`PlumedForce`, and runs 50 steps at 300 K. OpenMM's built-in four-site
model is TIP4P-Ew. `D_0` and `D_MAX` are Angstrom; `LENGTH_SCALE=10`
converts the nanometres OpenMM passes to PLUMED. Oxygens are every
fourth atom (`ATOMS=1-N:4`).

`plumed.dat` is the same input for 64 molecules. The script writes the
atom list for however many waters the box holds.

```
export PLUMED_KERNEL=/path/to/libplumedKernel.so
export LD_LIBRARY_PATH=/path/to/plumed/lib:/path/to/libOpenMMPlumed:$LD_LIBRARY_PATH
python3 run.py
```

The plugin is [openmm-plumed](https://github.com/openmm/openmm-plumed)
(`from openmmplumed import PlumedForce`). No trajectory is stored.
`run.py` deletes the COLVAR it prints. If the plugin's Python wrapper
will not pass `PlumedForce` through `System.addForce`, the script
builds `addforce.cpp` and adds the force from C++. That path needs
`OPENMM_INCLUDE` and `OPENMM_LIB`.

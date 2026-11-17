# Ions

Ions stay off the hydrogen-bond network. `IONS` builds the graph on
the oxygen atoms and reads each ion against that assignment. Water
inside `ION_CUTOFF` of the ion forms the first shell. An ion sits in
ice if every shell molecule carries a cage label. It sits in liquid if
no shell molecule does. Mixed shells put the ion at the front. The
three ion components report brine rejection at a growing front in an
electrolyte, for example TIP4P/2005 with the Madrid 2019 ions.

```text
LOAD FILE=./libdseams_plumed.so
brine: DSEAMS_CAGES ATOMS=1-3000 IONS=3001-3060 ION_CUTOFF=3.5 COMPLETE
PRINT ARG=brine.nice,brine.nmax,brine.nionice,brine.nionfront,brine.nionliq STRIDE=500 FILE=BRINE
```

| Component | Ion state |
| --- | --- |
| `nionice` | Every first-shell molecule carries a cage label. |
| `nionfront` | The first shell is mixed. |
| `nionliq` | No first-shell molecule carries a cage label. |

Without `IONS` the three components stay 0.

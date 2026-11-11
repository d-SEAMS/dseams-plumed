/*
** DSEAMS_SMOOTH_CHILL: a differentiable ice fraction for biasing.
**
** SPDX-License-Identifier: MIT
*/
#include "colvar/Colvar.h"
#include "core/ActionRegister.h"

#include "smooth_chill.hpp"

#include <cmath>
#include <vector>

namespace PLMD {
namespace colvar {

//+PLUMEDOC COLVAR DSEAMS_SMOOTH_CHILL
/*
Smooth CHILL+-like ice fraction, with derivatives.

DSEAMS_CHILL is the published Nguyen and Molinero classifier. Its labels
jump when a correlation crosses -0.8 or a neighbour enters the cutoff, so
the derivative is zero almost everywhere and undefined on the jump. This
action is the quantity to bias.

Each pair inside D_MAX contributes with the cubic switching function that
is 1 at D_0 and 0 at D_MAX (the PLUMED CUBIC switch). q_3 is the
switching-function average of Y_3m, in the same spherical-harmonic
convention as d-SEAMS. A bond is staggered to the degree
1 / (1 + exp((c_ij + 0.8) / WIDTH)). The molecule score is a sigmoid of
its staggered-bond count, centred at MID (default 2.2) so that cubic ice
(about four staggered bonds) and hexagonal ice (about three) both score
near one and a liquid scores near zero. The collective variable is the
mean molecule score.

On an ideal ice Ic lattice, with D_0 above the hydrogen-bond length and
D_MAX below the second shell, the shell is the four nearest neighbours and
the correlations match the published rule. The value is still not the
published fraction: the steps are soft.

Derivatives are analytic and include the angular dependence of q_3 and the
radial switching function. Use the value in METAD, RESTRAINT or
MOVINGRESTRAINT. STRIDE is not accepted: a bias needs the force every step.

\par Examples

\plumedfile
LOAD FILE=libdseams_plumed.so
UNITS LENGTH=A
ice: DSEAMS_SMOOTH_CHILL ATOMS=1-4096 D_0=3.2 D_MAX=3.5 LENGTH_SCALE=1.0
metad: METAD ARG=ice SIGMA=0.05 HEIGHT=0.2 PACE=500 GRID_MIN=0 GRID_MAX=1 GRID_BIN=200 FILE=HILLS
PRINT ARG=ice,metad.bias STRIDE=100 FILE=COLVAR
\endplumedfile

*/
//+ENDPLUMEDOC

class DseamsSmooth : public Colvar {
  dseams::SmoothParams par_;
  double lengthScale_ = 10.0;
  int nAtom_ = 0;

public:
  static void registerKeywords(Keywords &keys);
  explicit DseamsSmooth(const ActionOptions &);
  void calculate() override;
};

PLUMED_REGISTER_ACTION(DseamsSmooth, "DSEAMS_SMOOTH_CHILL")

void DseamsSmooth::registerKeywords(Keywords &keys) {
  Colvar::registerKeywords(keys);
  keys.add("atoms", "ATOMS", "the oxygen atoms, one per water molecule");
  keys.add("compulsory", "D_0", "3.2",
           "inner distance of the cubic neighbour switch, in Angstrom");
  keys.add("compulsory", "D_MAX", "3.5",
           "outer distance of the cubic neighbour switch, in Angstrom");
  keys.add("compulsory", "STAGGERED", "-0.8",
           "centre of the staggered sigmoid, the published c_ij edge");
  keys.add("compulsory", "WIDTH", "0.05", "width of the staggered sigmoid");
  keys.add("compulsory", "MID", "2.2",
           "staggered-count midpoint of the per-molecule sigmoid");
  keys.add("compulsory", "BETA", "0.30", "width of the per-molecule sigmoid");
  keys.add("compulsory", "LENGTH_SCALE", "10.0",
           "factor from the PLUMED length unit to Angstrom (10 for nm, 1 after UNITS LENGTH=A)");
}

DseamsSmooth::DseamsSmooth(const ActionOptions &ao) : PLUMED_COLVAR_INIT(ao) {
  std::vector<AtomNumber> atoms;
  parseAtomList("ATOMS", atoms);
  if (atoms.empty()) {
    error("ATOMS must name at least one oxygen");
  }
  parse("D_0", par_.d0);
  parse("D_MAX", par_.dmax);
  parse("STAGGERED", par_.staggered);
  parse("WIDTH", par_.width);
  parse("MID", par_.mid);
  parse("BETA", par_.beta);
  parse("LENGTH_SCALE", lengthScale_);
  checkRead();
  if (!(par_.dmax > par_.d0) || !(par_.d0 > 0.0)) {
    error("require 0 < D_0 < D_MAX");
  }
  if (!(par_.width > 0.0) || !(par_.beta > 0.0) || !(lengthScale_ > 0.0)) {
    error("WIDTH, BETA and LENGTH_SCALE must be positive");
  }
  nAtom_ = static_cast<int>(atoms.size());
  addValueWithDerivatives();
  setNotPeriodic();
  requestAtoms(atoms);
  log.printf("  %d atoms, cubic switch %.3f to %.3f A, staggered sigmoid at %.3f width %.3f\n",
             nAtom_, par_.d0, par_.dmax, par_.staggered, par_.width);
  log.printf("  molecule sigmoid midpoint %.3f width %.3f\n", par_.mid, par_.beta);
  log.printf("  smooth surrogate of CHILL+ with derivatives; not the published label. "
             "Monitor the label with DSEAMS_CHILL\n");
}

void DseamsSmooth::calculate() {
  // Pairs come from PLUMED's minimum image, so the force and the virial
  // see the same box tensor the driver differentiates.
  std::vector<dseams::PairBond> pairs;
  const double dmax2 = par_.dmax * par_.dmax;
  for (int i = 0; i < nAtom_; i++) {
    for (int j = i + 1; j < nAtom_; j++) {
      Vector d = pbcDistance(getPosition(i), getPosition(j)) * lengthScale_;
      if (d.modulo2() < dmax2) {
        pairs.push_back({i, j, d[0], d[1], d[2]});
      }
    }
  }
  const bool deriv = !doNotCalculateDerivatives();
  const auto value = dseams::smoothChill(nAtom_, pairs, par_, deriv);
  setValue(value.fraction);
  if (!deriv) {
    return;
  }
  for (int i = 0; i < nAtom_; i++) {
    const auto &g = value.grad[static_cast<std::size_t>(i)];
    setAtomsDerivatives(i, lengthScale_ * Vector(g[0], g[1], g[2]));
  }
  Tensor virial;
  for (int a = 0; a < 3; a++) {
    for (int b = 0; b < 3; b++) {
      virial(a, b) = value.virial[a][b];
    }
  }
  setBoxDerivatives(virial);
}

} // namespace colvar
} // namespace PLMD

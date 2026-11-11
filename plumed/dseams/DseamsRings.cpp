/*
** DSEAMS_RINGS: primitive ring census on the cutoff neighbour graph.
**
** SPDX-License-Identifier: MIT
*/
#include "colvar/Colvar.h"
#include "core/ActionRegister.h"

#include "dseams_frame.hpp"

#include <franzblau.hpp>
#include <neighbours.hpp>

#include <string>
#include <vector>

namespace PLMD {
namespace colvar {

//+PLUMEDOC COLVAR DSEAMS_RINGS
/*
Primitive ring counts of d-SEAMS on a group of oxygen atoms.

The neighbour graph is the cutoff shell (the same cutoff CHILL+ uses, 3.5 A
by default). primitive::ringNetwork enumerates cycles up to MAX and drops
every cycle that fails the Franzblau shortest-path test. n3 through n8 are
the number of primitive rings of that size. On ideal ice Ic or Ih the cutoff
graph is the hydrogen-bond network and the census is six-rings only.

This is not the six-ring count of DSEAMS_CAGES. That count is taken on the
union of the four-nearest graphs. The two graphs coincide on an undistorted
tetrahedral lattice and differ in the liquid.

Ring membership is discrete. No derivative is defined, and none is returned.
The counts are for PRINT, COMMITTOR and analysis, not for a bias.

\par Examples

\plumedfile
LOAD FILE=libdseams_plumed.so
UNITS LENGTH=A
rings: DSEAMS_RINGS ATOMS=1-4096 CUTOFF=3.5 MAX=6 LENGTH_SCALE=1.0
PRINT ARG=rings.n5,rings.n6,rings.n7 STRIDE=100 FILE=RINGS
\endplumedfile

*/
//+ENDPLUMEDOC

class DseamsRings : public Colvar {
  double cutoff_ = 3.5;
  double lengthScale_ = 10.0;
  int maxRing_ = 6;
  int nAtom_ = 0;

public:
  static void registerKeywords(Keywords &keys);
  explicit DseamsRings(const ActionOptions &);
  void calculate() override;
};

PLUMED_REGISTER_ACTION(DseamsRings, "DSEAMS_RINGS")

void DseamsRings::registerKeywords(Keywords &keys) {
  Colvar::registerKeywords(keys);
  keys.add("atoms", "ATOMS", "the oxygen atoms, one per water molecule");
  keys.add("compulsory", "CUTOFF", "3.5",
           "neighbour cutoff in Angstrom of the graph the rings are counted on");
  keys.add("compulsory", "MAX", "6", "largest primitive ring to enumerate (3 to 8)");
  keys.add("compulsory", "LENGTH_SCALE", "10.0",
           "factor from the PLUMED length unit to Angstrom (10 for nm, 1 after UNITS LENGTH=A)");
  keys.addOutputComponent("n3", "default", "primitive 3-rings");
  keys.addOutputComponent("n4", "default", "primitive 4-rings");
  keys.addOutputComponent("n5", "default", "primitive 5-rings");
  keys.addOutputComponent("n6", "default", "primitive 6-rings");
  keys.addOutputComponent("n7", "default", "primitive 7-rings; 0 when MAX < 7");
  keys.addOutputComponent("n8", "default", "primitive 8-rings; 0 when MAX < 8");
}

DseamsRings::DseamsRings(const ActionOptions &ao) : PLUMED_COLVAR_INIT(ao) {
  std::vector<AtomNumber> atoms;
  parseAtomList("ATOMS", atoms);
  if (atoms.empty()) {
    error("ATOMS must name at least one oxygen");
  }
  parse("CUTOFF", cutoff_);
  parse("MAX", maxRing_);
  parse("LENGTH_SCALE", lengthScale_);
  checkRead();
  if (!(cutoff_ > 0.0) || !(lengthScale_ > 0.0)) {
    error("CUTOFF and LENGTH_SCALE must be positive");
  }
  if (maxRing_ < 3 || maxRing_ > 8) {
    error("MAX must lie between 3 and 8");
  }
  nAtom_ = static_cast<int>(atoms.size());
  for (const char *name : {"n3", "n4", "n5", "n6", "n7", "n8"}) {
    addComponent(name);
    componentIsNotPeriodic(name);
  }
  requestAtoms(atoms);
  log.printf("  %d atoms, cutoff %.3f A, primitive rings through size %d\n", nAtom_, cutoff_,
             maxRing_);
  log.printf("  the graph is the cutoff shell, not the four-nearest union of DSEAMS_CAGES\n");
  log.printf("  ring counts have no derivatives; use them for PRINT, COMMITTOR and analysis\n");
}

void DseamsRings::calculate() {
  if (!doNotCalculateDerivatives()) {
    error("DSEAMS_RINGS counts are discrete and have no derivatives. "
          "Use them for PRINT and COMMITTOR, not for a bias.");
  }
  std::vector<Vector> pos(static_cast<std::size_t>(nAtom_));
  for (int i = 0; i < nAtom_; i++) {
    pos[static_cast<std::size_t>(i)] = getPosition(i);
  }
  auto cloud = dseams::cloudFromAtoms(nAtom_, pos, getBox(), lengthScale_, getStep());
  auto rows = nneigh::neighListO(cutoff_, cloud, 1);
  auto idx = nneigh::neighbourListByIndex(cloud, rows);
  const auto rings = primitive::ringNetwork(idx, maxRing_);
  int hist[9] = {};
  for (const auto &ring : rings) {
    const int s = static_cast<int>(ring.size());
    if (s >= 3 && s <= 8) {
      hist[s] += 1;
    }
  }
  getPntrToComponent("n3")->set(hist[3]);
  getPntrToComponent("n4")->set(hist[4]);
  getPntrToComponent("n5")->set(hist[5]);
  getPntrToComponent("n6")->set(hist[6]);
  getPntrToComponent("n7")->set(hist[7]);
  getPntrToComponent("n8")->set(hist[8]);
}

} // namespace colvar
} // namespace PLMD

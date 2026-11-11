/*
** DSEAMS_CHILL: published CHILL+ ice fraction, and the OVITO alternative.
**
** SPDX-License-Identifier: MIT
*/
#include "colvar/Colvar.h"
#include "core/ActionRegister.h"

#include "dseams_frame.hpp"

#include <bop.hpp>
#include <chill_plus_scheme.hpp>
#include <mol_sys.hpp>
#include <neighbours.hpp>

#include <algorithm>
#include <string>
#include <vector>

namespace PLMD {
namespace colvar {

//+PLUMEDOC COLVAR DSEAMS_CHILL
/*
Nguyen and Molinero CHILL+ ice fraction on a group of oxygen atoms.

The default rule is the published CHILL+ classifier (Nguyen and Molinero,
J. Phys. Chem. B 119, 9369, 2015). A neighbour inside the cutoff is staggered
when the q_3 correlation satisfies c_ij <= -0.8 and eclipsed when
-0.35 <= c_ij <= 0.25, boundaries included. q_3 is averaged over the four
nearest neighbours of that shell. Cubic ice is four staggered bonds,
hexagonal ice is three staggered and one eclipsed, and both rows keep the
published neighbour test for interfacial ice. fraction is
(cubic + hexagonal) / N, the bulk-ice fraction used for nucleation.

RULE=ovito-chill+ selects the OVITO 3.17.2 Chill+ modifier instead: q_3 is
averaged over the whole cutoff shell, the c_ij edges are open, and the
interfacial row has no neighbour test. Cubic and hexagonal ice on an ideal
lattice are unchanged. The interfacial count is not.

The labels are a step function of the coordinates. No derivative is defined,
and none is returned. Use the components for PRINT, COMMITTOR and analysis.
A bias or metadynamics needs DSEAMS_SMOOTH_CHILL, which is a smooth
surrogate and is not this classifier.

\par Examples

\plumedfile
LOAD FILE=libdseams_plumed.so
UNITS LENGTH=A
ice: DSEAMS_CHILL ATOMS=1-4096 CUTOFF=3.5 LENGTH_SCALE=1.0
PRINT ARG=ice.fraction,ice.nice,ice.ncubic,ice.nhex,ice.ninterfacial STRIDE=100 FILE=ICE
ovito: DSEAMS_CHILL ATOMS=1-4096 RULE=ovito-chill+ LENGTH_SCALE=1.0
PRINT ARG=ovito.fraction,ovito.ninterfacial STRIDE=100 FILE=OVITO
\endplumedfile

*/
//+ENDPLUMEDOC

class DseamsChill : public Colvar {
  double cutoff_ = 3.5;
  double lengthScale_ = 10.0;
  chill::ChillPlusRule rule_ = chill::ChillPlusRule::Published;
  int nAtom_ = 0;
  std::vector<std::string> names_;

public:
  static void registerKeywords(Keywords &keys);
  explicit DseamsChill(const ActionOptions &);
  void calculate() override;
};

PLUMED_REGISTER_ACTION(DseamsChill, "DSEAMS_CHILL")

void DseamsChill::registerKeywords(Keywords &keys) {
  Colvar::registerKeywords(keys);
  keys.add("atoms", "ATOMS", "the oxygen atoms, one per water molecule");
  keys.add("compulsory", "CUTOFF", "3.5", "neighbour cutoff in Angstrom");
  keys.add("compulsory", "LENGTH_SCALE", "10.0",
           "factor from the PLUMED length unit to Angstrom (10 for nm, 1 after UNITS LENGTH=A)");
  keys.add("compulsory", "RULE", "chill+",
           "chill+ is the published Nguyen and Molinero 2015 rule; "
           "ovito-chill+ is the OVITO 3.17.2 Chill+ modifier. A file whose "
           "text contains rule: or chill-plus: is read the same way");
  keys.addOutputComponent("fraction", "default",
                          "bulk ice fraction, (cubic + hexagonal) / N");
  keys.addOutputComponent("nice", "default", "cubic plus hexagonal molecules");
  keys.addOutputComponent("ncubic", "default", "cubic ice molecules");
  keys.addOutputComponent("nhex", "default", "hexagonal ice molecules");
  keys.addOutputComponent("ninterfacial", "default", "interfacial ice molecules");
  keys.addOutputComponent("nclathrate", "default", "clathrate molecules");
  keys.addOutputComponent("ninterclath", "default", "interfacial clathrate molecules");
  keys.addOutputComponent("nwater", "default", "molecules classified as liquid water");
  keys.addOutputComponent("nmax", "default", "largest connected cluster of bulk ice");
}

namespace {

int findRoot(std::vector<int> &parent, int a) {
  while (parent[static_cast<std::size_t>(a)] != a) {
    parent[static_cast<std::size_t>(a)] =
        parent[static_cast<std::size_t>(parent[static_cast<std::size_t>(a)])];
    a = parent[static_cast<std::size_t>(a)];
  }
  return a;
}

int largestCluster(const std::vector<char> &flag, const std::vector<std::vector<int>> &idx) {
  const int n = static_cast<int>(flag.size());
  std::vector<int> parent(static_cast<std::size_t>(n));
  std::vector<int> sz(static_cast<std::size_t>(n), 1);
  for (int i = 0; i < n; i++) {
    parent[static_cast<std::size_t>(i)] = i;
  }
  for (int i = 0; i < n; i++) {
    if (!flag[static_cast<std::size_t>(i)] || static_cast<int>(idx.size()) <= i) {
      continue;
    }
    for (std::size_t m = 1; m < idx[static_cast<std::size_t>(i)].size(); m++) {
      const int j = idx[static_cast<std::size_t>(i)][m];
      if (j > i && j < n && flag[static_cast<std::size_t>(j)]) {
        int a = findRoot(parent, i);
        int b = findRoot(parent, j);
        if (a != b) {
          if (sz[static_cast<std::size_t>(a)] < sz[static_cast<std::size_t>(b)]) {
            std::swap(a, b);
          }
          parent[static_cast<std::size_t>(b)] = a;
          sz[static_cast<std::size_t>(a)] += sz[static_cast<std::size_t>(b)];
        }
      }
    }
  }
  int nMax = 0;
  for (int i = 0; i < n; i++) {
    if (flag[static_cast<std::size_t>(i)] && findRoot(parent, i) == i) {
      nMax = std::max(nMax, sz[static_cast<std::size_t>(i)]);
    }
  }
  return nMax;
}

} // namespace

DseamsChill::DseamsChill(const ActionOptions &ao) : PLUMED_COLVAR_INIT(ao) {
  std::vector<AtomNumber> atoms;
  parseAtomList("ATOMS", atoms);
  if (atoms.empty()) {
    error("ATOMS must name at least one oxygen");
  }
  parse("CUTOFF", cutoff_);
  parse("LENGTH_SCALE", lengthScale_);
  std::string ruleSpec;
  parse("RULE", ruleSpec);
  checkRead();
  try {
    rule_ = chill::chillPlusRuleFromFile(ruleSpec);
  } catch (const std::invalid_argument &ex) {
    error(ex.what());
  }
  if (!(cutoff_ > 0.0) || !(lengthScale_ > 0.0)) {
    error("CUTOFF and LENGTH_SCALE must be positive");
  }
  nAtom_ = static_cast<int>(atoms.size());
  names_ = {"fraction", "nice", "ncubic", "nhex", "ninterfacial",
            "nclathrate", "ninterclath", "nwater", "nmax"};
  for (const auto &name : names_) {
    addComponent(name);
    componentIsNotPeriodic(name);
  }
  requestAtoms(atoms);
  log.printf("  %d atoms, cutoff %.3f A, rule %s\n", nAtom_, cutoff_,
             chill::chillPlusRuleName(rule_));
  log.printf("  published rule: four nearest inside the cutoff, closed c_ij edges, "
             "interfacial neighbour test\n");
  log.printf("  ovito-chill+: whole cutoff shell, open c_ij edges, no neighbour test\n");
  log.printf("  fraction = (cubic + hexagonal) / N. These labels have no derivatives "
             "and cannot be biased; use DSEAMS_SMOOTH_CHILL for metadynamics\n");
}

void DseamsChill::calculate() {
  if (!doNotCalculateDerivatives()) {
    error("DSEAMS_CHILL is a step function and has no derivatives. "
          "Use it for PRINT and COMMITTOR. Bias DSEAMS_SMOOTH_CHILL.");
  }
  std::vector<Vector> pos(static_cast<std::size_t>(nAtom_));
  for (int i = 0; i < nAtom_; i++) {
    pos[static_cast<std::size_t>(i)] = getPosition(i);
  }
  auto cloud = dseams::cloudFromAtoms(nAtom_, pos, getBox(), lengthScale_, getStep());
  auto rows = nneigh::neighListO(cutoff_, cloud, 1);
  auto idx = nneigh::neighbourListByIndex(cloud, rows);
  chill::applyChillPlus(cloud, rows, rule_, false);

  int cubic = 0;
  int hex = 0;
  int interfacial = 0;
  int clathrate = 0;
  int interClath = 0;
  int water = 0;
  std::vector<char> bulk(static_cast<std::size_t>(nAtom_), 0);
  for (int i = 0; i < nAtom_; i++) {
    switch (cloud.pts[static_cast<std::size_t>(i)].iceType) {
    case molSys::atom_state_type::cubic:
    case molSys::atom_state_type::reCubic:
      ++cubic;
      bulk[static_cast<std::size_t>(i)] = 1;
      break;
    case molSys::atom_state_type::hexagonal:
    case molSys::atom_state_type::reHex:
      ++hex;
      bulk[static_cast<std::size_t>(i)] = 1;
      break;
    case molSys::atom_state_type::interfacial:
      ++interfacial;
      break;
    case molSys::atom_state_type::clathrate:
      ++clathrate;
      break;
    case molSys::atom_state_type::interClathrate:
      ++interClath;
      break;
    default:
      ++water;
      break;
    }
  }
  const int nice = cubic + hex;
  const double fraction = static_cast<double>(nice) / static_cast<double>(nAtom_);
  getPntrToComponent("fraction")->set(fraction);
  getPntrToComponent("nice")->set(nice);
  getPntrToComponent("ncubic")->set(cubic);
  getPntrToComponent("nhex")->set(hex);
  getPntrToComponent("ninterfacial")->set(interfacial);
  getPntrToComponent("nclathrate")->set(clathrate);
  getPntrToComponent("ninterclath")->set(interClath);
  getPntrToComponent("nwater")->set(water);
  getPntrToComponent("nmax")->set(largestCluster(bulk, idx));
}

} // namespace colvar
} // namespace PLMD

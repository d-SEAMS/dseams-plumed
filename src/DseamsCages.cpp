/*
** dseams-plumed: the d-SEAMS topological ice score as a PLUMED action.
**
** SPDX-License-Identifier: MIT
*/
#include "colvar/Colvar.h"
#include "core/ActionRegister.h"
#include "tools/Pbc.h"

#include <bop.hpp>
#include <cage_affiliation.hpp>
#include <franzblau.hpp>
#include <mol_sys.hpp>
#include <neighbours.hpp>

#include <algorithm>
#include <string>
#include <vector>

namespace PLMD {
namespace colvar {

//+PLUMEDOC COLVAR DSEAMS_CAGES
/*
Per-frame topological ice score of d-SEAMS on a group of oxygen atoms.

Builds the mutual and union four-nearest-neighbour graphs, the primitive
six-membered rings and the seeded cage assignment (hexagonal and
double-diamond cages), optionally completed over edge-sharing rings, and
reports molecule counts and the largest connected ice cluster. CHILL+ on
the cutoff graph is reported next to it. The counts are integers with no
derivatives: they serve PRINT, COMMITTOR basins and analysis, not biasing.

Ions are not part of the hydrogen-bond network. With IONS the graph is
built on the oxygens alone and each ion is read against it: the water
within ION_CUTOFF of the ion is its first shell, and the ion counts as in
ice when every shell molecule carries a cage label, in liquid when none
does, and at the front otherwise. That is the brine-rejection observable
of a growing front in an electrolyte (TIP4P/2005 with the Madrid 2019
ions, for instance).

\par Examples

\plumedfile
LOAD FILE=libdseams_plumed.so
ice: DSEAMS_CAGES ATOMS=1-4096 CUTOFF=3.5 COMPLETE
PRINT ARG=ice.nice,ice.nmax,ice.nic,ice.nih,ice.chillice,ice.chillmax STRIDE=100 FILE=ICE
COMMITTOR ARG=ice.nmax STRIDE=100 BASIN_LL1=0 BASIN_UL1=20 BASIN_LL2=800 BASIN_UL2=100000
brine: DSEAMS_CAGES ATOMS=1-3000 IONS=3001-3060 ION_CUTOFF=3.5 COMPLETE
PRINT ARG=brine.nice,brine.nmax,brine.nionice,brine.nionfront,brine.nionliq STRIDE=500 FILE=BRINE
\endplumedfile

*/
//+ENDPLUMEDOC

class DseamsCages : public Colvar {
  double cutoff_ = 3.5;     // Angstrom, cutoff graph for CHILL+ and the fallback
  double candidate_ = 5.5;  // Angstrom, candidate shell for the k-nearest walk
  int k_ = 4;
  double lengthScale_ = 10.0;  // internal length unit -> Angstrom (nm default)
  bool complete_ = false;
  int nOxygen_ = 0;
  int nIon_ = 0;
  double ionCutoff_ = 3.5;  // Angstrom, first water shell of an ion
  std::vector<std::string> names_;

public:
  static void registerKeywords(Keywords &keys);
  explicit DseamsCages(const ActionOptions &);
  void calculate() override;
};

PLUMED_REGISTER_ACTION(DseamsCages, "DSEAMS_CAGES")

void DseamsCages::registerKeywords(Keywords &keys) {
  Colvar::registerKeywords(keys);
  keys.add("atoms", "ATOMS", "the oxygen atoms of the water molecules");
  keys.add("compulsory", "CUTOFF", "3.5",
           "neighbour cutoff in Angstrom for the CHILL+ graph");
  keys.add("compulsory", "CANDIDATE", "5.5",
           "candidate shell in Angstrom for the four-nearest walk");
  keys.add("compulsory", "K", "4", "neighbours per molecule in the bonded graph");
  keys.add("compulsory", "LENGTH_SCALE", "10.0",
           "factor from the PLUMED length unit to Angstrom (10 for nm)");
  keys.addFlag("COMPLETE", false,
               "fill the last vertex of six-rings whose other vertices carry a label");
  keys.add("atoms", "IONS", "ions read against the water assignment; not part of the graph");
  keys.add("compulsory", "ION_CUTOFF", "3.5",
           "radius in Angstrom of an ion's first water shell");
  keys.addOutputComponent("nice", "default", "molecules in a hexagonal or double-diamond cage");
  keys.addOutputComponent("nmax", "default", "largest connected cluster of cage molecules");
  keys.addOutputComponent("nclus", "default", "number of connected cage clusters");
  keys.addOutputComponent("nic", "default", "molecules in a double-diamond cage only");
  keys.addOutputComponent("nih", "default", "molecules in a hexagonal cage only");
  keys.addOutputComponent("nmixed", "default", "molecules in both cage types");
  keys.addOutputComponent("chillice", "default", "CHILL+ cubic plus hexagonal molecules");
  keys.addOutputComponent("chillmax", "default", "largest CHILL+ bulk-ice cluster");
  keys.addOutputComponent("chillinterfacial", "default", "CHILL+ interfacial molecules");
  keys.addOutputComponent("sixrings", "default", "primitive six-membered rings on the union graph");
  keys.addOutputComponent("nionice", "default", "ions whose first water shell is all cage molecules");
  keys.addOutputComponent("nionfront", "default", "ions with a mixed first shell");
  keys.addOutputComponent("nionliq", "default", "ions with no cage molecule in the first shell");
}

DseamsCages::DseamsCages(const ActionOptions &ao) : PLUMED_COLVAR_INIT(ao) {
  std::vector<AtomNumber> atoms;
  parseAtomList("ATOMS", atoms);
  if (atoms.empty()) {
    error("ATOMS must name at least one oxygen");
  }
  parse("CUTOFF", cutoff_);
  parse("CANDIDATE", candidate_);
  parse("K", k_);
  parse("LENGTH_SCALE", lengthScale_);
  parseFlag("COMPLETE", complete_);
  std::vector<AtomNumber> ions;
  parseAtomList("IONS", ions);
  parse("ION_CUTOFF", ionCutoff_);
  checkRead();
  nOxygen_ = static_cast<int>(atoms.size());
  nIon_ = static_cast<int>(ions.size());
  // PLUMED reserves the underscore in component names
  names_ = {"nice", "nmax", "nclus", "nic", "nih", "nmixed",
            "chillice", "chillmax", "chillinterfacial", "sixrings",
            "nionice", "nionfront", "nionliq"};
  for (const auto &n : names_) {
    addComponent(n);
    componentIsNotPeriodic(n);
  }
  atoms.insert(atoms.end(), ions.begin(), ions.end());
  requestAtoms(atoms);
  log.printf("  %d oxygens, cutoff %.3f A, k=%d within %.3f A, completion %s\n",
             nOxygen_, cutoff_, k_, candidate_, complete_ ? "on" : "off");
  if (nIon_ > 0) {
    log.printf("  %d ions read against the assignment, first shell %.3f A\n", nIon_,
               ionCutoff_);
  }
  log.printf("  counts carry no derivatives; use them for PRINT, COMMITTOR and analysis\n");
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

// Largest connected component and component count of the flagged atoms
// over the index graph (rows lead with the atom itself)
void clusterFlags(const std::vector<char> &flag,
                  const std::vector<std::vector<int>> &idx, int &nMax,
                  int &nClus) {
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
  nMax = 0;
  nClus = 0;
  for (int i = 0; i < n; i++) {
    if (flag[static_cast<std::size_t>(i)] && findRoot(parent, i) == i) {
      ++nClus;
      nMax = std::max(nMax, sz[static_cast<std::size_t>(i)]);
    }
  }
}

} // namespace

void DseamsCages::calculate() {
  const int nop = nOxygen_;
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  cloud.nop = nop;
  cloud.currentFrame = getStep();
  cloud.pts.resize(static_cast<std::size_t>(nop));
  cloud.idIndexMap.reserve(static_cast<std::size_t>(nop));
  for (int i = 0; i < nop; i++) {
    const Vector p = getPosition(i);
    auto &pt = cloud.pts[static_cast<std::size_t>(i)];
    pt.type = 1;
    pt.atomID = i + 1;
    pt.molID = i + 1;
    pt.x = p[0] * lengthScale_;
    pt.y = p[1] * lengthScale_;
    pt.z = p[2] * lengthScale_;
    cloud.idIndexMap[i + 1] = i;
  }
  // PLUMED keeps a along x and b in the xy plane, the LAMMPS convention:
  // spans lx, ly, lz then tilts xy, xz, yz. The engine's minimum image
  // reads that six-vector.
  const Tensor box = getBox();
  const double lx = box(0, 0) * lengthScale_;
  const double ly = box(1, 1) * lengthScale_;
  const double lz = box(2, 2) * lengthScale_;
  const double xy = box(1, 0) * lengthScale_;
  const double xz = box(2, 0) * lengthScale_;
  const double yz = box(2, 1) * lengthScale_;
  const bool tilted = xy != 0.0 || xz != 0.0 || yz != 0.0;
  cloud.box = tilted ? std::vector<double>{lx, ly, lz, xy, xz, yz}
                     : std::vector<double>{lx, ly, lz};
  cloud.boxLow = {0.0, 0.0, 0.0};

  // CHILL+ on the cutoff graph
  auto cutRows = nneigh::neighListO(cutoff_, cloud, 1);
  auto idxC = nneigh::neighbourListByIndex(cloud, cutRows);
  chill::getCorrelPlus(cloud, cutRows, false);
  chill::getIceTypePlusNoPrint(cloud, cutRows, false);
  std::vector<char> chillIce(static_cast<std::size_t>(nop), 0);
  int chillN = 0;
  int chillInterfacial = 0;
  for (int i = 0; i < nop; i++) {
    switch (cloud.pts[static_cast<std::size_t>(i)].iceType) {
    case molSys::atom_state_type::cubic:
    case molSys::atom_state_type::reCubic:
    case molSys::atom_state_type::hexagonal:
    case molSys::atom_state_type::reHex:
      chillIce[static_cast<std::size_t>(i)] = 1;
      ++chillN;
      break;
    case molSys::atom_state_type::interfacial:
      ++chillInterfacial;
      break;
    default:
      break;
    }
  }
  int chillMax = 0;
  int chillClus = 0;
  clusterFlags(chillIce, idxC, chillMax, chillClus);

  // Seeded cage assignment on the four-nearest graphs
  auto graphs = nneigh::kNearestNeighbourPair(cloud, k_, candidate_, 1);
  auto idxS = nneigh::neighbourListByIndex(cloud, graphs.first);
  auto idxU = nneigh::neighbourListByIndex(cloud, graphs.second);
  std::vector<std::vector<int>> sixS;
  std::vector<std::vector<int>> sixU;
  for (auto &r : primitive::ringNetwork(idxS, 6)) {
    if (r.size() == 6) {
      sixS.push_back(std::move(r));
    }
  }
  for (auto &r : primitive::ringNetwork(idxU, 6)) {
    if (r.size() == 6) {
      sixU.push_back(std::move(r));
    }
  }
  const auto seeded =
      ring::seededCageAffiliation(sixS, idxS, sixU, idxU, complete_);
  int ih = 0;
  int ic = 0;
  int mixed = 0;
  std::vector<char> ice(static_cast<std::size_t>(nop), 0);
  for (int i = 0; i < nop; i++) {
    const bool h = seeded.hc[static_cast<std::size_t>(i)];
    const bool d = seeded.ddc[static_cast<std::size_t>(i)];
    if (h && d) {
      ++mixed;
    } else if (h) {
      ++ih;
    } else if (d) {
      ++ic;
    }
    ice[static_cast<std::size_t>(i)] = (h || d) ? 1 : 0;
  }
  int nMax = 0;
  int nClus = 0;
  clusterFlags(ice, idxU, nMax, nClus);

  // Ions against the assignment: the water within ionCutoff_ is the first
  // shell; all labelled is in ice, none labelled is liquid, else front
  int ionIce = 0;
  int ionFront = 0;
  int ionLiq = 0;
  const double shell2 = ionCutoff_ * ionCutoff_ / (lengthScale_ * lengthScale_);
  for (int a = 0; a < nIon_; a++) {
    const Vector ionPos = getPosition(nop + a);
    int inShell = 0;
    int labelled = 0;
    for (int i = 0; i < nop; i++) {
      const Vector d = pbcDistance(getPosition(i), ionPos);
      if (d.modulo2() < shell2) {
        ++inShell;
        labelled += ice[static_cast<std::size_t>(i)] ? 1 : 0;
      }
    }
    if (inShell > 0 && labelled == inShell) {
      ++ionIce;
    } else if (labelled > 0) {
      ++ionFront;
    } else {
      ++ionLiq;
    }
  }

  getPntrToComponent("nice")->set(ih + ic + mixed);
  getPntrToComponent("nmax")->set(nMax);
  getPntrToComponent("nclus")->set(nClus);
  getPntrToComponent("nic")->set(ic);
  getPntrToComponent("nih")->set(ih);
  getPntrToComponent("nmixed")->set(mixed);
  getPntrToComponent("chillice")->set(chillN);
  getPntrToComponent("chillmax")->set(chillMax);
  getPntrToComponent("chillinterfacial")->set(chillInterfacial);
  getPntrToComponent("sixrings")->set(static_cast<double>(sixU.size()));
  getPntrToComponent("nionice")->set(ionIce);
  getPntrToComponent("nionfront")->set(ionFront);
  getPntrToComponent("nionliq")->set(ionLiq);
}

} // namespace colvar
} // namespace PLMD

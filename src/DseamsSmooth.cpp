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

The neighbour list is a Verlet list: it is built out to D_MAX plus a skin
and kept until an atom has moved by half that skin or the box has changed.
Pair vectors are refreshed from the minimum image every step, so the force
and the virial see the current box. Y_3m and its Cartesian gradient are
closed form, and the scratch used by the sum is kept between steps.

On an ideal ice Ic lattice, with D_0 above the hydrogen-bond length and
D_MAX below the second shell, the shell is the four nearest neighbours and
the correlations match the published rule. The value is still not the
published fraction: the steps are soft.

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
  // Half a skin of motion cannot bring the second ice shell (about 4.5 A)
  // inside D_MAX before the list is rebuilt.
  double skin_ = 0.9;
  int nAtom_ = 0;
  bool haveList_ = false;
  double listBox_[9] = {};
  std::vector<double> anchor_; // positions in Angstrom at the last rebuild
  std::vector<int> listI_;
  std::vector<int> listJ_;
  std::vector<int> cellHead_;
  std::vector<int> cellNext_;
  std::vector<dseams::PairBond> pairs_;
  dseams::SmoothChill engine_;

public:
  static void registerKeywords(Keywords &keys);
  explicit DseamsSmooth(const ActionOptions &);
  void calculate() override;

private:
  bool listStale(const Tensor &box) const;
  void rebuildList(const Tensor &box);
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
  anchor_.assign(static_cast<std::size_t>(nAtom_) * 3, 0.0);
  log.printf("  %d atoms, cubic switch %.3f to %.3f A, staggered sigmoid at %.3f width %.3f\n",
             nAtom_, par_.d0, par_.dmax, par_.staggered, par_.width);
  log.printf("  molecule sigmoid midpoint %.3f width %.3f\n", par_.mid, par_.beta);
  log.printf("  Verlet neighbour list with skin %.2f A, rebuilt when an atom moves "
             "by half of that; Y_3m gradient is analytic\n",
             skin_);
  log.printf("  smooth surrogate of CHILL+ with derivatives; not the published label. "
             "Monitor the label with DSEAMS_CHILL\n");
}

bool DseamsSmooth::listStale(const Tensor &box) const {
  if (!haveList_) {
    return true;
  }
  const double scale = lengthScale_;
  const double cur[9] = {box(0, 0) * scale, box(0, 1) * scale, box(0, 2) * scale,
                         box(1, 0) * scale, box(1, 1) * scale, box(1, 2) * scale,
                         box(2, 0) * scale, box(2, 1) * scale, box(2, 2) * scale};
  for (int k = 0; k < 9; k++) {
    if (std::abs(cur[k] - listBox_[k]) > 1e-3) {
      return true;
    }
  }
  const double lim2 = (0.5 * skin_) * (0.5 * skin_);
  for (int i = 0; i < nAtom_; i++) {
    const Vector p = getPosition(i) * scale;
    const double dx = p[0] - anchor_[static_cast<std::size_t>(3 * i)];
    const double dy = p[1] - anchor_[static_cast<std::size_t>(3 * i + 1)];
    const double dz = p[2] - anchor_[static_cast<std::size_t>(3 * i + 2)];
    if (dx * dx + dy * dy + dz * dz > lim2) {
      return true;
    }
  }
  return false;
}

void DseamsSmooth::rebuildList(const Tensor &box) {
  const double scale = lengthScale_;
  const double cut = par_.dmax + skin_;
  const double cut2 = cut * cut;
  listBox_[0] = box(0, 0) * scale;
  listBox_[1] = box(0, 1) * scale;
  listBox_[2] = box(0, 2) * scale;
  listBox_[3] = box(1, 0) * scale;
  listBox_[4] = box(1, 1) * scale;
  listBox_[5] = box(1, 2) * scale;
  listBox_[6] = box(2, 0) * scale;
  listBox_[7] = box(2, 1) * scale;
  listBox_[8] = box(2, 2) * scale;
  const double lx = listBox_[0];
  const double ly = listBox_[4];
  const double lz = listBox_[8];
  for (int i = 0; i < nAtom_; i++) {
    const Vector p = getPosition(i) * scale;
    anchor_[static_cast<std::size_t>(3 * i)] = p[0];
    anchor_[static_cast<std::size_t>(3 * i + 1)] = p[1];
    anchor_[static_cast<std::size_t>(3 * i + 2)] = p[2];
  }
  listI_.clear();
  listJ_.clear();
  const bool ortho = std::abs(listBox_[1]) < 1e-8 && std::abs(listBox_[2]) < 1e-8 &&
                     std::abs(listBox_[3]) < 1e-8 && std::abs(listBox_[5]) < 1e-8 &&
                     std::abs(listBox_[6]) < 1e-8 && std::abs(listBox_[7]) < 1e-8;
  const bool cells = ortho && lx > 3.0 * cut && ly > 3.0 * cut && lz > 3.0 * cut;
  auto accept = [&](int i, int j) {
    const Vector d = pbcDistance(getPosition(i), getPosition(j)) * scale;
    if (d.modulo2() < cut2) {
      listI_.push_back(i);
      listJ_.push_back(j);
    }
  };
  if (!cells) {
    for (int i = 0; i < nAtom_; i++) {
      for (int j = i + 1; j < nAtom_; j++) {
        accept(i, j);
      }
    }
  } else {
    const int nx = std::max(1, static_cast<int>(lx / cut));
    const int ny = std::max(1, static_cast<int>(ly / cut));
    const int nz = std::max(1, static_cast<int>(lz / cut));
    const int ncell = nx * ny * nz;
    cellHead_.assign(static_cast<std::size_t>(ncell), -1);
    cellNext_.assign(static_cast<std::size_t>(nAtom_), -1);
    auto cellOf = [&](int i) {
      double x = anchor_[static_cast<std::size_t>(3 * i)];
      double y = anchor_[static_cast<std::size_t>(3 * i + 1)];
      double z = anchor_[static_cast<std::size_t>(3 * i + 2)];
      x -= lx * std::floor(x / lx);
      y -= ly * std::floor(y / ly);
      z -= lz * std::floor(z / lz);
      int ix = std::min(nx - 1, static_cast<int>(x / lx * nx));
      int iy = std::min(ny - 1, static_cast<int>(y / ly * ny));
      int iz = std::min(nz - 1, static_cast<int>(z / lz * nz));
      if (ix < 0) {
        ix = 0;
      }
      if (iy < 0) {
        iy = 0;
      }
      if (iz < 0) {
        iz = 0;
      }
      return (iz * ny + iy) * nx + ix;
    };
    for (int i = 0; i < nAtom_; i++) {
      const int c = cellOf(i);
      cellNext_[static_cast<std::size_t>(i)] = cellHead_[static_cast<std::size_t>(c)];
      cellHead_[static_cast<std::size_t>(c)] = i;
    }
    auto at = [&](int ix, int iy, int iz) {
      ix = (ix % nx + nx) % nx;
      iy = (iy % ny + ny) % ny;
      iz = (iz % nz + nz) % nz;
      return (iz * ny + iy) * nx + ix;
    };
    const int off[14][3] = {{0, 0, 0},  {1, 0, 0},  {-1, 1, 0}, {0, 1, 0},  {1, 1, 0},
                            {-1, -1, 1}, {0, -1, 1}, {1, -1, 1}, {-1, 0, 1}, {0, 0, 1},
                            {1, 0, 1},   {-1, 1, 1}, {0, 1, 1},  {1, 1, 1}};
    for (int iz = 0; iz < nz; iz++) {
      for (int iy = 0; iy < ny; iy++) {
        for (int ix = 0; ix < nx; ix++) {
          const int home = at(ix, iy, iz);
          for (int o = 0; o < 14; o++) {
            const int other = at(ix + off[o][0], iy + off[o][1], iz + off[o][2]);
            for (int a = cellHead_[static_cast<std::size_t>(home)]; a >= 0;
                 a = cellNext_[static_cast<std::size_t>(a)]) {
              for (int b = cellHead_[static_cast<std::size_t>(other)]; b >= 0;
                   b = cellNext_[static_cast<std::size_t>(b)]) {
                if (o == 0 && b <= a) {
                  continue;
                }
                if (o != 0 && a == b) {
                  continue;
                }
                accept(a, b);
              }
            }
          }
        }
      }
    }
  }
  haveList_ = true;
}

void DseamsSmooth::calculate() {
  const Tensor box = getBox();
  if (listStale(box)) {
    rebuildList(box);
  }
  pairs_.clear();
  if (pairs_.capacity() < listI_.size()) {
    pairs_.reserve(listI_.size());
  }
  const double dmax2 = par_.dmax * par_.dmax;
  const double scale = lengthScale_;
  // Minimum image from PLUMED, including a tilted box. A hand-rolled
  // orthorhombic image disagrees with that derivative at the 1e-4 level.
  for (std::size_t n = 0; n < listI_.size(); n++) {
    const int i = listI_[n];
    const int j = listJ_[n];
    const Vector d = pbcDistance(getPosition(i), getPosition(j)) * scale;
    if (d.modulo2() < dmax2) {
      pairs_.push_back({i, j, d[0], d[1], d[2]});
    }
  }
  const bool deriv = !doNotCalculateDerivatives();
  engine_.evaluate(nAtom_, pairs_, par_, deriv);
  setValue(engine_.fraction);
  if (!deriv) {
    return;
  }
  for (int i = 0; i < nAtom_; i++) {
    const auto &g = engine_.grad[static_cast<std::size_t>(i)];
    setAtomsDerivatives(i, scale * Vector(g[0], g[1], g[2]));
  }
  Tensor virial;
  for (int a = 0; a < 3; a++) {
    for (int b = 0; b < 3; b++) {
      virial(a, b) = engine_.virial[a][b];
    }
  }
  setBoxDerivatives(virial);
}

} // namespace colvar
} // namespace PLMD

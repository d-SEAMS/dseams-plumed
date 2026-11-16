/*
** DSEAMS_SMOOTH_CHILL: a differentiable ice fraction for biasing.
**
** SPDX-License-Identifier: MIT
*/
#include "colvar/Colvar.h"
#include "core/ActionRegister.h"
#include "tools/OpenMP.h"

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
1 / (1 + exp((c_ij + 0.8) / WIDTH)). An eclipsed bond is the published
window from -0.35 to 0.25, with edges of width ECL_WIDTH. A molecule
scores by a soft OR of the published counts: four staggered neighbours
(cubic), or three staggered and one eclipsed (hexagonal). BETA is the
width of those count sigmoids. A second, softer staggered weight
(TAIL_C, TAIL_W) turns on between a liquid correlation and the
published -0.8 edge, and the molecule score includes
1 - exp(-n / TAIL_N) for that count. The centre is close enough to
-0.8 that the weight is still climbing at the published threshold, so
a restraint cannot be satisfied by bonds that are only mildly
staggered. TAIL_AMP caps that channel so the liquid score stays below
a small crystal nucleus. A Gaussian gate in the smooth neighbour count is 1 at the
published shell of four (COORD_W is its width) and is below one when
a liquid carries an extra neighbour. On a frozen molecule the
published channel is already one and the shell is already four, so
neither term moves the value. In a liquid both terms still have a
slope. The collective variable is the mean molecule score.

The neighbour list is a Verlet list: it is built out to D_MAX plus a skin
and kept until an atom has moved by half that skin or the box has changed.
Pair vectors are refreshed from the minimum image every step, so the force
and the virial see the current box. Y_3m and its Cartesian gradient are
closed form, and the scratch used by the sum is kept between steps.
The loops run on PLUMED's NUM_THREADS. Y_3m is summed from a
molecule-major copy of the neighbour list, still in pair-list order.
The molecule scores and the virial are each reduced as eight partials
folded from the first chunk to the last, so one thread and many threads
write the same bits.

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
  std::vector<int> candI_;
  std::vector<int> candJ_;
  std::vector<double> candX_;
  std::vector<double> candY_;
  std::vector<double> candZ_;
  std::vector<char> candKeep_;
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
  keys.add("compulsory", "WIDTH", "0.06", "width of the staggered sigmoid in c_ij");
  keys.add("compulsory", "ECL_WIDTH", "0.06",
           "width of each edge of the published eclipsed window");
  keys.add("compulsory", "BETA", "0.06",
           "width of the staggered-count sigmoids (cubic at 4, hexagonal at 3)");
  keys.add("compulsory", "BETA_E", "0.10", "width of the eclipsed-count sigmoid");
  keys.add("compulsory", "TAIL_C", "-0.65",
           "centre of the soft staggered weight, between a liquid correlation and -0.8");
  keys.add("compulsory", "TAIL_W", "0.20", "width of the soft staggered weight in c_ij");
  keys.add("compulsory", "TAIL_N", "4.0",
           "soft-count scale. The tail is 1 - exp(-n / TAIL_N)");
  keys.add("compulsory", "TAIL_AMP", "0.65",
           "ceiling of the soft channel, so a liquid cannot outscore a small nucleus");
  keys.add("compulsory", "COORD_W", "0.70",
           "width of the Gaussian that is 1 when the smooth neighbour count is 4");
  keys.add("compulsory", "LENGTH_SCALE", "10.0",
           "factor from the PLUMED length unit to Angstrom (10 for nm, 1 after UNITS LENGTH=A)");
  keys.setValueDescription("mean molecule score, between 0 and 1");
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
  parse("ECL_WIDTH", par_.widthE);
  parse("BETA", par_.beta);
  parse("BETA_E", par_.betaE);
  parse("TAIL_C", par_.tailCenter);
  parse("TAIL_W", par_.tailWidth);
  parse("TAIL_N", par_.tailN);
  parse("TAIL_AMP", par_.tailAmp);
  parse("COORD_W", par_.coordWidth);
  parse("LENGTH_SCALE", lengthScale_);
  checkRead();
  if (!(par_.dmax > par_.d0) || !(par_.d0 > 0.0)) {
    error("require 0 < D_0 < D_MAX");
  }
  if (!(par_.width > 0.0) || !(par_.widthE > 0.0) || !(par_.beta > 0.0) || !(par_.betaE > 0.0) ||
      !(par_.tailWidth > 0.0) || !(par_.tailN > 0.0) || !(par_.tailAmp > 0.0) || !(par_.coordWidth > 0.0) ||
      !(lengthScale_ > 0.0)) {
    error("WIDTH, ECL_WIDTH, BETA, BETA_E, TAIL_W, TAIL_N, TAIL_AMP, COORD_W and LENGTH_SCALE must be positive");
  }
  nAtom_ = static_cast<int>(atoms.size());
  addValueWithDerivatives();
  setNotPeriodic();
  requestAtoms(atoms);
  anchor_.assign(static_cast<std::size_t>(nAtom_) * 3, 0.0);
  log.printf("  %d atoms, cubic switch %.3f to %.3f A, staggered sigmoid at %.3f width %.3f\n",
             nAtom_, par_.d0, par_.dmax, par_.staggered, par_.width);
  log.printf("  published counts, softened: cubic at 4 staggered, hexagonal at 3 staggered "
             "plus 1 eclipsed; count widths %.3f and %.3f\n",
             par_.beta, par_.betaE);
  log.printf("  soft staggered tail centred at %.3f width %.3f, count scale %.3f, "
             "ceiling %.3f; neighbour-count gate width %.3f around 4\n",
             par_.tailCenter, par_.tailWidth, par_.tailN, par_.tailAmp, par_.coordWidth);
  log.printf("  Verlet neighbour list with skin %.2f A, rebuilt when an atom moves "
             "by half of that; Y_3m gradient is analytic\n",
             skin_);
  log.printf("  OpenMP width is PLUMED NUM_THREADS (currently %u). One thread and many "
             "threads evaluate the same value and the same forces\n",
             OpenMP::getNumThreads());
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
  candI_.clear();
  candJ_.clear();
  const bool ortho = std::abs(listBox_[1]) < 1e-8 && std::abs(listBox_[2]) < 1e-8 &&
                     std::abs(listBox_[3]) < 1e-8 && std::abs(listBox_[5]) < 1e-8 &&
                     std::abs(listBox_[6]) < 1e-8 && std::abs(listBox_[7]) < 1e-8;
  const bool cells = ortho && lx > 3.0 * cut && ly > 3.0 * cut && lz > 3.0 * cut;
  auto consider = [&](int i, int j) {
    candI_.push_back(i);
    candJ_.push_back(j);
  };
  if (!cells) {
    for (int i = 0; i < nAtom_; i++) {
      for (int j = i + 1; j < nAtom_; j++) {
        consider(i, j);
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
                consider(a, b);
              }
            }
          }
        }
      }
    }
  }
  const int nc = static_cast<int>(candI_.size());
  candX_.resize(static_cast<std::size_t>(nc));
  candY_.resize(static_cast<std::size_t>(nc));
  candZ_.resize(static_cast<std::size_t>(nc));
  candKeep_.resize(static_cast<std::size_t>(nc));
  const int nt = std::max(1, static_cast<int>(OpenMP::getNumThreads()));
#pragma omp parallel for num_threads(nt) schedule(static)
  for (int n = 0; n < nc; n++) {
    const Vector d =
        pbcDistance(getPosition(candI_[static_cast<std::size_t>(n)]),
                    getPosition(candJ_[static_cast<std::size_t>(n)])) *
        scale;
    candX_[static_cast<std::size_t>(n)] = d[0];
    candY_[static_cast<std::size_t>(n)] = d[1];
    candZ_[static_cast<std::size_t>(n)] = d[2];
    candKeep_[static_cast<std::size_t>(n)] = d.modulo2() < cut2 ? 1 : 0;
  }
  for (int n = 0; n < nc; n++) {
    if (candKeep_[static_cast<std::size_t>(n)]) {
      listI_.push_back(candI_[static_cast<std::size_t>(n)]);
      listJ_.push_back(candJ_[static_cast<std::size_t>(n)]);
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
  const int nlist = static_cast<int>(listI_.size());
  candX_.resize(static_cast<std::size_t>(nlist));
  candY_.resize(static_cast<std::size_t>(nlist));
  candZ_.resize(static_cast<std::size_t>(nlist));
  candKeep_.resize(static_cast<std::size_t>(nlist));
  const int nt = std::max(1, static_cast<int>(OpenMP::getNumThreads()));
  // Minimum image from PLUMED, including a tilted box. A hand-rolled
  // orthorhombic image disagrees with that derivative at the 1e-4 level.
  // Distances are written by pair index and then gathered in list order.
#pragma omp parallel for num_threads(nt) schedule(static)
  for (int n = 0; n < nlist; n++) {
    const Vector d = pbcDistance(getPosition(listI_[static_cast<std::size_t>(n)]),
                                 getPosition(listJ_[static_cast<std::size_t>(n)])) *
                     scale;
    candX_[static_cast<std::size_t>(n)] = d[0];
    candY_[static_cast<std::size_t>(n)] = d[1];
    candZ_[static_cast<std::size_t>(n)] = d[2];
    candKeep_[static_cast<std::size_t>(n)] = d.modulo2() < dmax2 ? 1 : 0;
  }
  // Pack survivors in list order. Chunk offsets make the result independent
  // of how many threads run the copy.
  std::vector<int> chunkBase(static_cast<std::size_t>(nt + 1));
  std::vector<int> chunkCount(static_cast<std::size_t>(nt));
#pragma omp parallel for num_threads(nt) schedule(static)
  for (int t = 0; t < nt; t++) {
    const int begin = nlist * t / nt;
    const int end = nlist * (t + 1) / nt;
    int count = 0;
    for (int n = begin; n < end; n++) {
      count += candKeep_[static_cast<std::size_t>(n)];
    }
    chunkCount[static_cast<std::size_t>(t)] = count;
  }
  int packed = 0;
  for (int t = 0; t < nt; t++) {
    chunkBase[static_cast<std::size_t>(t)] = packed;
    packed += chunkCount[static_cast<std::size_t>(t)];
  }
  pairs_.resize(static_cast<std::size_t>(packed));
#pragma omp parallel for num_threads(nt) schedule(static)
  for (int t = 0; t < nt; t++) {
    const int begin = nlist * t / nt;
    const int end = nlist * (t + 1) / nt;
    int dest = chunkBase[static_cast<std::size_t>(t)];
    for (int n = begin; n < end; n++) {
      if (!candKeep_[static_cast<std::size_t>(n)]) {
        continue;
      }
      pairs_[static_cast<std::size_t>(dest++)] = {
          listI_[static_cast<std::size_t>(n)], listJ_[static_cast<std::size_t>(n)],
          candX_[static_cast<std::size_t>(n)], candY_[static_cast<std::size_t>(n)],
          candZ_[static_cast<std::size_t>(n)]};
    }
  }
  const bool deriv = !doNotCalculateDerivatives();
  engine_.evaluate(nAtom_, pairs_, par_, deriv, nt);
  setValue(engine_.fraction);
  if (!deriv) {
    return;
  }
#pragma omp parallel for num_threads(nt) schedule(static)
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

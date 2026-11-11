/*
** Shared frame packing for the d-SEAMS PLUMED actions.
**
** SPDX-License-Identifier: MIT
*/
#ifndef DSEAMS_FRAME_HPP_
#define DSEAMS_FRAME_HPP_

#include "tools/Tensor.h"

#include <mol_sys.hpp>

#include <vector>

namespace dseams {

// PLUMED stores lattice vectors as the rows of the box tensor. The engine
// reads the LAMMPS tilt form: spans lx, ly, lz then xy, xz, yz.
inline molSys::PointCloud<molSys::Point<double>, double>
cloudFromAtoms(int nop, const std::vector<PLMD::Vector> &pos,
               const PLMD::Tensor &box, double lengthScale, int step) {
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  cloud.nop = nop;
  cloud.currentFrame = step;
  cloud.pts.resize(static_cast<std::size_t>(nop));
  cloud.idIndexMap.reserve(static_cast<std::size_t>(nop));
  for (int i = 0; i < nop; i++) {
    auto &pt = cloud.pts[static_cast<std::size_t>(i)];
    pt.type = 1;
    pt.atomID = i + 1;
    pt.molID = i + 1;
    pt.x = pos[static_cast<std::size_t>(i)][0] * lengthScale;
    pt.y = pos[static_cast<std::size_t>(i)][1] * lengthScale;
    pt.z = pos[static_cast<std::size_t>(i)][2] * lengthScale;
    cloud.idIndexMap[i + 1] = i;
  }
  const double lx = box(0, 0) * lengthScale;
  const double ly = box(1, 1) * lengthScale;
  const double lz = box(2, 2) * lengthScale;
  const double xy = box(1, 0) * lengthScale;
  const double xz = box(2, 0) * lengthScale;
  const double yz = box(2, 1) * lengthScale;
  const bool tilted = xy != 0.0 || xz != 0.0 || yz != 0.0;
  cloud.box = tilted ? std::vector<double>{lx, ly, lz, xy, xz, yz}
                     : std::vector<double>{lx, ly, lz};
  cloud.boxLow = {0.0, 0.0, 0.0};
  return cloud;
}

} // namespace dseams

#endif

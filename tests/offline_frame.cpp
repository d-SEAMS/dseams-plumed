/*
** Offline d-SEAMS counts for the same frames the PLUMED driver reads.
**
** SPDX-License-Identifier: MIT
*/
#include "smooth_chill.hpp"

#include <bop.hpp>
#include <chill_plus_scheme.hpp>
#include <franzblau.hpp>
#include <mol_sys.hpp>
#include <neighbours.hpp>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <random>
#include <sstream>
#include <string>
#include <vector>

namespace {

struct Frame {
  double lx = 0, ly = 0, lz = 0, xy = 0, xz = 0, yz = 0;
  std::vector<std::array<double, 3>> pos;
};

std::vector<Frame> readXyz(const std::string &path) {
  std::ifstream in(path);
  if (!in) {
    throw std::runtime_error("cannot read " + path);
  }
  std::vector<Frame> frames;
  std::string line;
  while (std::getline(in, line)) {
    if (line.empty()) {
      continue;
    }
    Frame frame;
    const int n = std::stoi(line);
    std::getline(in, line);
    std::istringstream box(line);
    double a[9] = {};
    for (int i = 0; i < 9; i++) {
      box >> a[i];
    }
    frame.lx = a[0];
    frame.xy = a[3];
    frame.ly = a[4];
    frame.xz = a[6];
    frame.yz = a[7];
    frame.lz = a[8];
    frame.pos.resize(static_cast<std::size_t>(n));
    for (int i = 0; i < n; i++) {
      std::getline(in, line);
      std::istringstream row(line);
      std::string elem;
      row >> elem >> frame.pos[static_cast<std::size_t>(i)][0] >>
          frame.pos[static_cast<std::size_t>(i)][1] >> frame.pos[static_cast<std::size_t>(i)][2];
    }
    frames.push_back(std::move(frame));
  }
  return frames;
}

molSys::PointCloud<molSys::Point<double>, double> toCloud(const Frame &frame, int step) {
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  cloud.nop = static_cast<int>(frame.pos.size());
  cloud.currentFrame = step;
  cloud.pts.resize(frame.pos.size());
  for (int i = 0; i < cloud.nop; i++) {
    auto &pt = cloud.pts[static_cast<std::size_t>(i)];
    pt.type = 1;
    pt.atomID = i + 1;
    pt.molID = i + 1;
    pt.x = frame.pos[static_cast<std::size_t>(i)][0];
    pt.y = frame.pos[static_cast<std::size_t>(i)][1];
    pt.z = frame.pos[static_cast<std::size_t>(i)][2];
    cloud.idIndexMap[i + 1] = i;
  }
  const bool tilted = frame.xy != 0.0 || frame.xz != 0.0 || frame.yz != 0.0;
  cloud.box = tilted ? std::vector<double>{frame.lx, frame.ly, frame.lz, frame.xy, frame.xz, frame.yz}
                     : std::vector<double>{frame.lx, frame.ly, frame.lz};
  cloud.boxLow = {0.0, 0.0, 0.0};
  return cloud;
}

struct Counts {
  int cubic = 0, hex = 0, interfacial = 0, clathrate = 0, interClath = 0, water = 0, nmax = 0;
  double fraction = 0;
};

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

Counts tally(molSys::PointCloud<molSys::Point<double>, double> &cloud,
             const std::vector<std::vector<int>> &idx) {
  Counts c;
  std::vector<char> bulk(static_cast<std::size_t>(cloud.nop), 0);
  for (int i = 0; i < cloud.nop; i++) {
    switch (cloud.pts[static_cast<std::size_t>(i)].iceType) {
    case molSys::atom_state_type::cubic:
    case molSys::atom_state_type::reCubic:
      ++c.cubic;
      bulk[static_cast<std::size_t>(i)] = 1;
      break;
    case molSys::atom_state_type::hexagonal:
    case molSys::atom_state_type::reHex:
      ++c.hex;
      bulk[static_cast<std::size_t>(i)] = 1;
      break;
    case molSys::atom_state_type::interfacial:
      ++c.interfacial;
      break;
    case molSys::atom_state_type::clathrate:
      ++c.clathrate;
      break;
    case molSys::atom_state_type::interClathrate:
      ++c.interClath;
      break;
    default:
      ++c.water;
      break;
    }
  }
  c.nmax = largestCluster(bulk, idx);
  c.fraction = static_cast<double>(c.cubic + c.hex) / static_cast<double>(cloud.nop);
  return c;
}

int checkYlm() {
  std::mt19937 rng(7);
  std::uniform_real_distribution<double> dist(-2.0, 2.0);
  double worst = 0.0;
  int tested = 0;
  for (int n = 0; n < 48; n++) {
    const double x = dist(rng);
    const double y = dist(rng);
    const double z = dist(rng);
    if (x * x + y * y + z * z < 1e-4) {
      continue;
    }
    const auto got = dseams::smooth_detail::ylm3(x, y, z);
    const auto ref = sph::spheriHarmo(3, sph::radialCoord({x, y, z}));
    for (int m = 0; m < 7; m++) {
      worst = std::max(worst, std::abs(got.re[m] - ref[static_cast<std::size_t>(m)].real()));
      worst = std::max(worst, std::abs(got.im[m] - ref[static_cast<std::size_t>(m)].imag()));
    }
    ++tested;
  }
  std::cout << "ylm directions " << tested << " max abs error " << worst << "\n";
  return worst < 1e-8 ? 0 : 1;
}

int checkCorrel(const Frame &frame) {
  auto cloud = toCloud(frame, 0);
  auto rows = nneigh::neighListO(3.5, cloud, 1);
  chill::applyChillPlus(cloud, rows, chill::ChillPlusRule::Published, false);
  const int n = cloud.nop;
  std::vector<std::vector<int>> neigh(static_cast<std::size_t>(n));
  for (int i = 0; i < n; i++) {
    std::vector<std::pair<double, int>> cand;
    for (int j = 0; j < n; j++) {
      if (j == i) {
        continue;
      }
      const double dx = dseams::minimumImage(frame.pos[static_cast<std::size_t>(j)][0] -
                                                  frame.pos[static_cast<std::size_t>(i)][0],
                                              frame.lx);
      const double dy = dseams::minimumImage(frame.pos[static_cast<std::size_t>(j)][1] -
                                                  frame.pos[static_cast<std::size_t>(i)][1],
                                              frame.ly);
      const double dz = dseams::minimumImage(frame.pos[static_cast<std::size_t>(j)][2] -
                                                  frame.pos[static_cast<std::size_t>(i)][2],
                                              frame.lz);
      const double d2 = dx * dx + dy * dy + dz * dz;
      if (d2 < 3.5 * 3.5) {
        cand.emplace_back(d2, j);
      }
    }
    std::sort(cand.begin(), cand.end());
    const std::size_t keep = std::min<std::size_t>(4, cand.size());
    for (std::size_t k = 0; k < keep; k++) {
      neigh[static_cast<std::size_t>(i)].push_back(cand[k].second);
    }
  }
  const auto ours = dseams::hardShellCorrelations(frame.pos, neigh, frame.lx, frame.ly, frame.lz);
  double worst = 0.0;
  int bonds = 0;
  for (int i = 0; i < n; i++) {
    std::vector<double> engine;
    for (const auto &bond : cloud.pts[static_cast<std::size_t>(i)].c_ij) {
      engine.push_back(bond.c_value);
    }
    auto mine = ours[static_cast<std::size_t>(i)];
    if (engine.size() != mine.size()) {
      std::cerr << "atom " << i << " bond count engine " << engine.size() << " smooth " << mine.size()
                << "\n";
      return 1;
    }
    std::sort(engine.begin(), engine.end());
    std::sort(mine.begin(), mine.end());
    for (std::size_t k = 0; k < engine.size(); k++) {
      worst = std::max(worst, std::abs(engine[k] - mine[k]));
      ++bonds;
    }
  }
  std::cout << "correl bonds " << bonds << " max abs error " << worst << "\n";
  return worst < 1e-6 ? 0 : 1;
}

void printCounts(const std::vector<Frame> &frames, int maxRing) {
  int step = 0;
  for (const auto &frame : frames) {
    auto cloud = toCloud(frame, step);
    auto rows = nneigh::neighListO(3.5, cloud, 1);
    auto idx = nneigh::neighbourListByIndex(cloud, rows);
    auto published = cloud;
    auto ovito = cloud;
    chill::applyChillPlus(published, rows, chill::ChillPlusRule::Published, false);
    chill::applyChillPlus(ovito, rows, chill::ChillPlusRule::Ovito, false);
    const Counts pub = tally(published, idx);
    const Counts ovi = tally(ovito, idx);
    const auto rings = primitive::ringNetwork(idx, maxRing);
    int hist[9] = {};
    for (const auto &ring : rings) {
      const int s = static_cast<int>(ring.size());
      if (s >= 3 && s <= 8) {
        hist[s] += 1;
      }
    }
    std::cout << step << " " << cloud.nop << " " << pub.fraction << " " << pub.cubic << " " << pub.hex
              << " " << pub.interfacial << " " << pub.clathrate << " " << pub.interClath << " "
              << pub.water << " " << pub.nmax << " " << ovi.fraction << " " << ovi.cubic << " "
              << ovi.hex << " " << ovi.interfacial << " " << ovi.clathrate << " " << ovi.interClath
              << " " << ovi.water << " " << ovi.nmax << " " << hist[3] << " " << hist[4] << " "
              << hist[5] << " " << hist[6] << " " << hist[7] << " " << hist[8] << "\n";
    ++step;
  }
}

} // namespace

int main(int argc, char **argv) {
  try {
    if (argc < 2) {
      std::cerr << "usage: dseams_offline traj.xyz [--max-ring N] [--check-ylm] [--check-correl]\n";
      return 2;
    }
    int maxRing = 6;
    bool ylm = false;
    bool correl = false;
    std::string path;
    for (int i = 1; i < argc; i++) {
      const std::string arg = argv[i];
      if (arg == "--check-ylm") {
        ylm = true;
      } else if (arg == "--check-correl") {
        correl = true;
      } else if (arg == "--max-ring" && i + 1 < argc) {
        maxRing = std::stoi(argv[++i]);
      } else {
        path = arg;
      }
    }
    int status = 0;
    if (ylm) {
      status |= checkYlm();
    }
    if (path.empty()) {
      return (ylm && status == 0) ? 0 : (status == 0 ? 2 : status);
    }
    const auto frames = readXyz(path);
    if (frames.empty()) {
      std::cerr << "no frames in " << path << "\n";
      return 2;
    }
    if (correl) {
      status |= checkCorrel(frames.front());
    }
    if (!ylm && !correl) {
      printCounts(frames, maxRing);
    }
    return status;
  } catch (const std::exception &ex) {
    std::cerr << ex.what() << "\n";
    return 1;
  }
}

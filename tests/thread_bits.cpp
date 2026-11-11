/*
** One OpenMP thread and eight threads must write the same bits.
**
** SPDX-License-Identifier: MIT
*/
#include "smooth_chill.hpp"

#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>

namespace {

double iceA() { return 4.0 * 2.75 / std::sqrt(3.0); }

void ice(int reps, std::vector<std::array<double, 3>> &pos, double &box) {
  const double a = iceA();
  const double fcc[4][3] = {{0, 0, 0}, {0.5, 0.5, 0}, {0.5, 0, 0.5}, {0, 0.5, 0.5}};
  box = reps * a;
  pos.clear();
  for (int i = 0; i < reps; i++) {
    for (int j = 0; j < reps; j++) {
      for (int k = 0; k < reps; k++) {
        for (int b = 0; b < 8; b++) {
          const double bx = fcc[b % 4][0] + (b >= 4 ? 0.25 : 0.0);
          const double by = fcc[b % 4][1] + (b >= 4 ? 0.25 : 0.0);
          const double bz = fcc[b % 4][2] + (b >= 4 ? 0.25 : 0.0);
          double x = std::fmod((i + bx) * a, box);
          double y = std::fmod((j + by) * a, box);
          double z = std::fmod((k + bz) * a, box);
          if (x < 0) {
            x += box;
          }
          if (y < 0) {
            y += box;
          }
          if (z < 0) {
            z += box;
          }
          pos.push_back({x, y, z});
        }
      }
    }
  }
  // A few molecules off their lattice sites, so the switch has slope.
  if (pos.size() > 3) {
    pos[0][0] += 0.12;
    pos[1][1] -= 0.07;
    pos[2][2] += 0.04;
  }
}

double mic(double d, double length) {
  return d - length * std::round(d / length);
}

std::vector<dseams::PairBond> pairsOf(const std::vector<std::array<double, 3>> &pos, double box,
                                     double cut) {
  std::vector<dseams::PairBond> pairs;
  const double cut2 = cut * cut;
  const int n = static_cast<int>(pos.size());
  for (int i = 0; i < n; i++) {
    for (int j = i + 1; j < n; j++) {
      const double dx = mic(pos[static_cast<std::size_t>(j)][0] - pos[static_cast<std::size_t>(i)][0], box);
      const double dy = mic(pos[static_cast<std::size_t>(j)][1] - pos[static_cast<std::size_t>(i)][1], box);
      const double dz = mic(pos[static_cast<std::size_t>(j)][2] - pos[static_cast<std::size_t>(i)][2], box);
      if (dx * dx + dy * dy + dz * dz < cut2) {
        pairs.push_back({i, j, dx, dy, dz});
      }
    }
  }
  return pairs;
}

bool sameBits(double a, double b) {
  std::uint64_t ia = 0;
  std::uint64_t ib = 0;
  std::memcpy(&ia, &a, sizeof ia);
  std::memcpy(&ib, &b, sizeof ib);
  return ia == ib;
}

} // namespace

int main() {
  dseams::SmoothParams par;
  par.d0 = 2.6;
  par.dmax = 3.4;
  par.width = 0.12;
  par.beta = 0.35;
  std::vector<std::array<double, 3>> pos;
  double box = 0.0;
  ice(2, pos, box);
  const auto pairs = pairsOf(pos, box, par.dmax + 0.9);
  dseams::SmoothChill one;
  dseams::SmoothChill many;
  one.evaluate(static_cast<int>(pos.size()), pairs, par, true, 1);
  many.evaluate(static_cast<int>(pos.size()), pairs, par, true, 8);
  int bad = 0;
  if (!sameBits(one.fraction, many.fraction)) {
    std::fprintf(stderr, "fraction bits differ\n");
    bad++;
  }
  if (one.grad.size() != many.grad.size()) {
    std::fprintf(stderr, "grad size differs\n");
    return 1;
  }
  for (std::size_t i = 0; i < one.grad.size(); i++) {
    for (int k = 0; k < 3; k++) {
      if (!sameBits(one.grad[i][static_cast<std::size_t>(k)], many.grad[i][static_cast<std::size_t>(k)])) {
        std::fprintf(stderr, "grad %zu %d differs\n", i, k);
        bad++;
        if (bad > 8) {
          break;
        }
      }
    }
  }
  for (int a = 0; a < 3; a++) {
    for (int b = 0; b < 3; b++) {
      if (!sameBits(one.virial[a][b], many.virial[a][b])) {
        std::fprintf(stderr, "virial %d %d differs\n", a, b);
        bad++;
      }
    }
  }
  // A second width, still the same bits.
  dseams::SmoothChill four;
  four.evaluate(static_cast<int>(pos.size()), pairs, par, true, 4);
  if (!sameBits(one.fraction, four.fraction)) {
    std::fprintf(stderr, "4-thread fraction bits differ\n");
    bad++;
  }
  if (bad) {
    return 1;
  }
  std::printf("thread bits match, fraction %a, atoms %zu, pairs %zu\n", one.fraction, pos.size(),
              pairs.size());
  return 0;
}

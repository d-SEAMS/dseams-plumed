/*
** The smooth fraction tracks the published CHILL+ label from liquid to ice,
** and it does not decrease as that label increases along a slab cut.
**
** SPDX-License-Identifier: MIT
*/
#include "smooth_chill.hpp"

#include <chill_plus_scheme.hpp>
#include <mol_sys.hpp>
#include <neighbours.hpp>

#include <cmath>
#include <cstdio>
#include <random>
#include <string>
#include <vector>

namespace {

struct Cloud {
  std::vector<std::array<double, 3>> pos;
  double box = 0;
};

Cloud iceIc(int reps) {
  const double bond = 2.75;
  const double a = 4.0 * bond / std::sqrt(3.0);
  const double fcc[4][3] = {{0, 0, 0}, {0.5, 0.5, 0}, {0.5, 0, 0.5}, {0, 0.5, 0.5}};
  Cloud out;
  out.box = reps * a;
  for (int i = 0; i < reps; i++) {
    for (int j = 0; j < reps; j++) {
      for (int k = 0; k < reps; k++) {
        for (int b = 0; b < 8; b++) {
          const double bx = fcc[b % 4][0] + (b >= 4 ? 0.25 : 0.0);
          const double by = fcc[b % 4][1] + (b >= 4 ? 0.25 : 0.0);
          const double bz = fcc[b % 4][2] + (b >= 4 ? 0.25 : 0.0);
          double x = std::fmod((i + bx) * a, out.box);
          double y = std::fmod((j + by) * a, out.box);
          double z = std::fmod((k + bz) * a, out.box);
          if (x < 0) x += out.box;
          if (y < 0) y += out.box;
          if (z < 0) z += out.box;
          out.pos.push_back({x, y, z});
        }
      }
    }
  }
  return out;
}

double mic(double d, double box) { return d - box * std::round(d / box); }

std::vector<dseams::PairBond> pairsOf(const Cloud &c, double cut) {
  std::vector<dseams::PairBond> pairs;
  const double cut2 = cut * cut;
  const int n = static_cast<int>(c.pos.size());
  for (int i = 0; i < n; i++) {
    for (int j = i + 1; j < n; j++) {
      const double dx = mic(c.pos[static_cast<std::size_t>(j)][0] - c.pos[static_cast<std::size_t>(i)][0], c.box);
      const double dy = mic(c.pos[static_cast<std::size_t>(j)][1] - c.pos[static_cast<std::size_t>(i)][1], c.box);
      const double dz = mic(c.pos[static_cast<std::size_t>(j)][2] - c.pos[static_cast<std::size_t>(i)][2], c.box);
      if (dx * dx + dy * dy + dz * dz < cut2) {
        pairs.push_back({i, j, dx, dy, dz});
      }
    }
  }
  return pairs;
}

double smoothOf(const Cloud &c, const dseams::SmoothParams &par) {
  dseams::SmoothChill eng;
  const auto pairs = pairsOf(c, par.dmax + 0.9);
  eng.evaluate(static_cast<int>(c.pos.size()), pairs, par, false, 1);
  return eng.fraction;
}

double publishedOf(const Cloud &c) {
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  cloud.nop = static_cast<int>(c.pos.size());
  cloud.pts.resize(c.pos.size());
  cloud.box = {c.box, c.box, c.box};
  cloud.boxLow = {0.0, 0.0, 0.0};
  for (int i = 0; i < cloud.nop; i++) {
    auto &pt = cloud.pts[static_cast<std::size_t>(i)];
    pt.type = 1;
    pt.atomID = i + 1;
    pt.molID = i + 1;
    pt.x = c.pos[static_cast<std::size_t>(i)][0];
    pt.y = c.pos[static_cast<std::size_t>(i)][1];
    pt.z = c.pos[static_cast<std::size_t>(i)][2];
    cloud.idIndexMap[i + 1] = i;
  }
  auto rows = nneigh::neighListO(3.5, cloud, 1);
  chill::applyChillPlus(cloud, rows, chill::ChillPlusRule::Published, false);
  int ice = 0;
  for (const auto &pt : cloud.pts) {
    if (pt.iceType == molSys::atom_state_type::cubic || pt.iceType == molSys::atom_state_type::hexagonal) {
      ++ice;
    }
  }
  return static_cast<double>(ice) / static_cast<double>(cloud.nop);
}

Cloud slice(const Cloud &crystal, double crystalFrac, std::mt19937 &rng) {
  Cloud out = crystal;
  std::uniform_real_distribution<double> shift(-1.6, 1.6);
  const double zcut = crystalFrac * crystal.box;
  for (auto &p : out.pos) {
    if (p[2] <= zcut) {
      continue;
    }
    p[0] = std::fmod(p[0] + shift(rng) + crystal.box, crystal.box);
    p[1] = std::fmod(p[1] + shift(rng) + crystal.box, crystal.box);
    p[2] = std::fmod(p[2] + shift(rng) + crystal.box, crystal.box);
  }
  return out;
}

} // namespace

int main(int argc, char **argv) {
  const bool fit = argc > 1 && std::string(argv[1]) == "--fit";
  const auto base = iceIc(4);
  std::vector<double> levels;
  for (int s = 0; s <= 10; s++) {
    levels.push_back(s / 10.0);
  }
  dseams::SmoothParams par;
  auto score = [&](const dseams::SmoothParams &p, double &rmse, double &worst, bool &mono) {
    std::vector<double> pub, sm;
    std::mt19937 rng(11);
    for (double level : levels) {
      const auto frame = slice(base, level, rng);
      pub.push_back(publishedOf(frame));
      sm.push_back(smoothOf(frame, p));
    }
    rmse = 0;
    worst = 0;
    for (std::size_t i = 0; i < pub.size(); i++) {
      const double e = sm[i] - pub[i];
      rmse += e * e;
      worst = std::max(worst, std::abs(e));
    }
    rmse = std::sqrt(rmse / static_cast<double>(pub.size()));
    mono = true;
    for (std::size_t i = 1; i < sm.size(); i++) {
      if (sm[i] + 1e-3 < sm[i - 1]) {
        mono = false;
      }
    }
    if (fit) {
      std::printf("beta %.3f w %.3f we %.3f be %.3f rmse %.4f worst %.4f mono %d ice %.4f liq %.4f pub_ice %.4f pub_liq %.4f\n",
                  p.beta, p.width, p.widthE, p.betaE, rmse, worst, mono ? 1 : 0, sm.back(), sm.front(), pub.back(),
                  pub.front());
    }
    return std::pair<std::vector<double>, std::vector<double>>{pub, sm};
  };

  if (fit) {
    for (double beta : {0.08, 0.10, 0.12, 0.16}) {
      for (double width : {0.04, 0.06, 0.10}) {
        for (double widthE : {0.06, 0.10}) {
          for (double betaE : {0.10, 0.16}) {
            dseams::SmoothParams p;
            p.beta = beta;
            p.width = width;
            p.widthE = widthE;
            p.betaE = betaE;
            double rmse = 0, worst = 0;
            bool mono = false;
            score(p, rmse, worst, mono);
          }
        }
      }
    }
    return 0;
  }

  double rmse = 0, worst = 0;
  bool mono = false;
  const auto both = score(par, rmse, worst, mono);
  const auto &pub = both.first;
  const auto &sm = both.second;
  std::printf("level published smooth residual\n");
  for (std::size_t i = 0; i < pub.size(); i++) {
    std::printf("%.1f %.5f %.5f %.5f\n", levels[i], pub[i], sm[i], sm[i] - pub[i]);
  }
  std::printf("rmse %.5f worst %.5f mono %d\n", rmse, worst, mono ? 1 : 0);
  int bad = 0;
  if (!(sm.back() > 0.97)) {
    std::fprintf(stderr, "crystal smooth %.5f is below 0.97\n", sm.back());
    bad++;
  }
  double frozen = 0.0;
  int nfrozen = 0;
  for (std::size_t i = 0; i < pub.size(); i++) {
    if (pub[i] > 0.9) {
      const double e = sm[i] - pub[i];
      frozen += e * e;
      nfrozen++;
    }
  }
  if (nfrozen > 0) {
    frozen = std::sqrt(frozen / static_cast<double>(nfrozen));
  }
  std::printf("frozen_rmse %.5f over %d frames with published > 0.9\n", frozen, nfrozen);
  dseams::SmoothChill eng;
  std::mt19937 liquidRng(11);
  const auto liquid = slice(base, 0.0, liquidRng);
  const auto pairs = pairsOf(liquid, par.dmax + 0.9);
  eng.evaluate(static_cast<int>(liquid.pos.size()), pairs, par, true, 1);
  double gsum = 0.0;
  for (const auto &g : eng.grad) {
    gsum += std::sqrt(g[0] * g[0] + g[1] * g[1] + g[2] * g[2]);
  }
  const double gmean = gsum / static_cast<double>(eng.grad.size());
  std::printf("liquid_grad %.6e value %.5f\n", gmean, eng.fraction);
  if (!(frozen < 0.03)) {
    std::fprintf(stderr, "frozen rmse %.5f\n", frozen);
    bad++;
  }
  if (!(gmean > 1e-4)) {
    std::fprintf(stderr, "liquid gradient %.3e is flat\n", gmean);
    bad++;
  }
  if (!mono) {
    std::fprintf(stderr, "smooth fraction decreased along the slab cut\n");
    bad++;
  }
  (void)rmse;
  return bad == 0 ? 0 : 1;
}

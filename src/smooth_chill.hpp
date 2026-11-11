/*
** Smooth CHILL+-like ice fraction.
**
** The published Nguyen and Molinero classifier is a step function of the
** bond correlation c_ij, so it has no force. This surrogate keeps the same
** q_3 correlation, replaces the hard cutoff and the hard staggered threshold
** by smooth steps, and returns a fraction in (0, 1) with Cartesian
** derivatives.
**
** q_3 uses the seams-core convention (azimuth atan2(x, y)). The value and
** the Y_3m gradient are closed form. Scratch buffers live on SmoothChill
** and are reused across calls.
**
** SPDX-License-Identifier: MIT
*/
#ifndef DSEAMS_SMOOTH_CHILL_HPP_
#define DSEAMS_SMOOTH_CHILL_HPP_

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <vector>

namespace dseams {

struct SmoothParams {
  double d0 = 3.2;         // Angstrom, cubic switch is 1 inside this
  double dmax = 3.5;       // Angstrom, cubic switch is 0 outside this
  double staggered = -0.8; // Nguyen and Molinero staggered edge
  double width = 0.05;     // width of the staggered sigmoid
  double mid = 2.2;        // staggered-count midpoint of a molecule
  double beta = 0.30;      // width of that molecule sigmoid
};

struct PairBond {
  int i = 0;
  int j = 0;
  double dx = 0.0; // minimum image of r_j - r_i, Angstrom
  double dy = 0.0;
  double dz = 0.0;
};

namespace smooth_detail {

struct Ylm3 {
  double re[7] = {};
  double im[7] = {};
  double dre[7][3] = {};
  double dim[7][3] = {};
};

inline void ylmNorms(double n[4]) {
  constexpr double pi = 3.14159265358979323846;
  n[0] = 0.25 * std::sqrt(7.0 / pi);
  n[1] = 0.125 * std::sqrt(21.0 / pi);
  n[2] = 0.25 * std::sqrt(105.0 / (2.0 * pi));
  n[3] = 0.125 * std::sqrt(35.0 / pi);
}

// Value of Y_3m, m = -3..+3. Matches seams-core sph::spheriHarmo for l = 3.
inline void ylm3Value(double x, double y, double z, double re[7], double im[7]) {
  double nrm[4];
  ylmNorms(nrm);
  const double r = std::sqrt(x * x + y * y + z * z);
  const double rho2 = x * x + y * y;
  const double rho = std::sqrt(rho2 + 1e-30);
  const double ct = z / r;
  const double st = rho / r;
  const bool pole = rho2 < 1e-24;
  const double cp = pole ? 1.0 : y / rho;
  const double sp = pole ? 0.0 : x / rho;
  const double ct2 = ct * ct;
  const double ct3 = ct2 * ct;
  const double st2 = st * st;
  const double a0 = nrm[0] * (5.0 * ct3 - 3.0 * ct);
  const double a1 = nrm[1] * st * (5.0 * ct2 - 1.0);
  const double a2 = nrm[2] * st2 * ct;
  const double a3 = nrm[3] * st2 * st;
  const double pr2 = cp * cp - sp * sp;
  const double pi2 = 2.0 * cp * sp;
  const double pr3 = pr2 * cp - pi2 * sp;
  const double pi3 = pr2 * sp + pi2 * cp;
  const double amp[4] = {a0, a1, a2, a3};
  const double pr[4] = {1.0, cp, pr2, pr3};
  const double pi[4] = {0.0, sp, pi2, pi3};
  for (int m = 0; m <= 3; m++) {
    const double negRe = amp[m] * pr[m];
    const double negIm = -amp[m] * pi[m];
    const double sign = (m % 2 == 0) ? 1.0 : -1.0;
    const int in = 3 - m;
    const int ip = 3 + m;
    re[in] = negRe;
    im[in] = negIm;
    re[ip] = sign * negRe;
    im[ip] = sign * amp[m] * pi[m];
  }
}

// Value and dY/d(x, y, z). The chain rule follows ylm3Value, including the
// 1e-30 floor under the cylindrical radius and a zero azimuthal derivative
// on the z axis.
inline void ylm3ValueGrad(double x, double y, double z, Ylm3 &out) {
  double nrm[4];
  ylmNorms(nrm);
  const double r2 = x * x + y * y + z * z;
  const double r = std::sqrt(r2);
  const double rho2 = x * x + y * y;
  const double rho = std::sqrt(rho2 + 1e-30);
  const double ct = z / r;
  const double st = rho / r;
  const bool pole = rho2 < 1e-24;
  const double cp = pole ? 1.0 : y / rho;
  const double sp = pole ? 0.0 : x / rho;
  const double ct2 = ct * ct;
  const double st2 = st * st;
  const double a[4] = {nrm[0] * (5.0 * ct2 * ct - 3.0 * ct), nrm[1] * st * (5.0 * ct2 - 1.0),
                       nrm[2] * st2 * ct, nrm[3] * st2 * st};
  const double dact[4] = {nrm[0] * (15.0 * ct2 - 3.0), nrm[1] * st * (10.0 * ct), nrm[2] * st2, 0.0};
  const double dast[4] = {0.0, nrm[1] * (5.0 * ct2 - 1.0), nrm[2] * (2.0 * st * ct),
                          nrm[3] * (3.0 * st2)};
  const double pr2 = cp * cp - sp * sp;
  const double pi2 = 2.0 * cp * sp;
  const double pr[4] = {1.0, cp, pr2, pr2 * cp - pi2 * sp};
  const double pi[4] = {0.0, sp, pi2, pr2 * sp + pi2 * cp};
  const double r3 = r2 * r;
  const double rho3 = rho * rho * rho;
  const double dct[3] = {-z * x / r3, -z * y / r3, rho2 / r3};
  const double dst[3] = {x * (r / rho - rho / r) / r2, y * (r / rho - rho / r) / r2,
                         -rho * z / r3};
  double dcp[3] = {0.0, 0.0, 0.0};
  double dsp[3] = {0.0, 0.0, 0.0};
  if (!pole) {
    dcp[0] = -y * x / rho3;
    dcp[1] = (rho * rho - y * y) / rho3;
    dsp[0] = (rho * rho - x * x) / rho3;
    dsp[1] = -x * y / rho3;
  }
  for (int m = 0; m <= 3; m++) {
    const int in = 3 - m;
    const int ip = 3 + m;
    const double sign = (m % 2 == 0) ? 1.0 : -1.0;
    out.re[in] = a[m] * pr[m];
    out.im[in] = -a[m] * pi[m];
    out.re[ip] = sign * out.re[in];
    out.im[ip] = sign * a[m] * pi[m];
  }
  for (int k = 0; k < 3; k++) {
    const double dpr2 = 2.0 * cp * dcp[k] - 2.0 * sp * dsp[k];
    const double dpi2 = 2.0 * (dcp[k] * sp + cp * dsp[k]);
    const double dpr3 = dpr2 * cp + pr2 * dcp[k] - dpi2 * sp - pi2 * dsp[k];
    const double dpi3 = dpr2 * sp + pr2 * dsp[k] + dpi2 * cp + pi2 * dcp[k];
    const double dpr[4] = {0.0, dcp[k], dpr2, dpr3};
    const double dpi[4] = {0.0, dsp[k], dpi2, dpi3};
    for (int m = 0; m <= 3; m++) {
      const double da = dact[m] * dct[k] + dast[m] * dst[k];
      const double dNegRe = da * pr[m] + a[m] * dpr[m];
      const double dNegIm = -(da * pi[m] + a[m] * dpi[m]);
      const double sign = (m % 2 == 0) ? 1.0 : -1.0;
      const int in = 3 - m;
      const int ip = 3 + m;
      out.dre[in][k] = dNegRe;
      out.dim[in][k] = dNegIm;
      out.dre[ip][k] = sign * dNegRe;
      out.dim[ip][k] = sign * (da * pi[m] + a[m] * dpi[m]);
    }
  }
}

inline double cubicSwitch(double r, double d0, double dmax, double &dsdr) {
  if (r <= d0) {
    dsdr = 0.0;
    return 1.0;
  }
  if (r >= dmax) {
    dsdr = 0.0;
    return 0.0;
  }
  const double span = dmax - d0;
  const double y = (r - d0) / span;
  const double s = (y - 1.0) * (y - 1.0) * (1.0 + 2.0 * y);
  dsdr = (6.0 * y * (y - 1.0)) / span;
  return s;
}

inline double sigmoid(double x) {
  if (x > 40.0) {
    return 1.0;
  }
  if (x < -40.0) {
    return 0.0;
  }
  return 1.0 / (1.0 + std::exp(-x));
}

struct Qacc {
  double re[7] = {};
  double im[7] = {};
  void clear() {
    for (int m = 0; m < 7; m++) {
      re[m] = 0.0;
      im[m] = 0.0;
    }
  }
};

struct Edge {
  int i = 0;
  int j = 0;
  double dx = 0.0;
  double dy = 0.0;
  double dz = 0.0;
  double sigma = 0.0;
  double dsdr = 0.0;
  double c = 0.0;
  double s = 0.0;
  double dsdc = 0.0;
  Ylm3 y;
};

} // namespace smooth_detail

// Reused scratch for one group of atoms. evaluate() does not allocate
// after the first call of a given size, apart from growth of the edge list.
class SmoothChill {
public:
  double fraction = 0.0;
  std::vector<std::array<double, 3>> grad;
  double virial[3][3] = {};

  void evaluate(int n, const std::vector<PairBond> &pairs, const SmoothParams &par, bool wantDeriv) {
    using namespace smooth_detail;
    ensure(n);
    fraction = 0.0;
    if (n <= 0) {
      return;
    }
    edges_.clear();
    if (edges_.capacity() < pairs.size() * 2) {
      edges_.reserve(pairs.size() * 2);
    }
    for (const auto &p : pairs) {
      push(p.i, p.j, p.dx, p.dy, p.dz, par, wantDeriv);
      push(p.j, p.i, -p.dx, -p.dy, -p.dz, par, wantDeriv);
    }
    std::fill(coord_.begin(), coord_.end(), 0.0);
    for (auto &u : U_) {
      u.clear();
    }
    for (const auto &e : edges_) {
      coord_[static_cast<std::size_t>(e.i)] += e.sigma;
      auto &u = U_[static_cast<std::size_t>(e.i)];
      for (int m = 0; m < 7; m++) {
        u.re[m] += e.sigma * e.y.re[m];
        u.im[m] += e.sigma * e.y.im[m];
      }
    }
    std::fill(norm2_.begin(), norm2_.end(), 0.0);
    for (int i = 0; i < n; i++) {
      const double z = coord_[static_cast<std::size_t>(i)];
      auto &q = Qlm_[static_cast<std::size_t>(i)];
      q.clear();
      if (z <= 0.0) {
        continue;
      }
      const auto &u = U_[static_cast<std::size_t>(i)];
      double n2 = 0.0;
      for (int m = 0; m < 7; m++) {
        q.re[m] = u.re[m] / z;
        q.im[m] = u.im[m] / z;
        n2 += q.re[m] * q.re[m] + q.im[m] * q.im[m];
      }
      norm2_[static_cast<std::size_t>(i)] = n2;
    }
    std::fill(ns_.begin(), ns_.end(), 0.0);
    for (auto &e : edges_) {
      const double ni2 = norm2_[static_cast<std::size_t>(e.i)];
      const double nj2 = norm2_[static_cast<std::size_t>(e.j)];
      if (ni2 < 1e-24 || nj2 < 1e-24) {
        e.c = 0.0;
        e.s = 0.0;
        e.dsdc = 0.0;
        continue;
      }
      const auto &qi = Qlm_[static_cast<std::size_t>(e.i)];
      const auto &qj = Qlm_[static_cast<std::size_t>(e.j)];
      double dot = 0.0;
      for (int m = 0; m < 7; m++) {
        dot += qi.re[m] * qj.re[m] + qi.im[m] * qj.im[m];
      }
      e.c = dot / std::sqrt(ni2 * nj2);
      const double x = (e.c - par.staggered) / par.width;
      e.s = sigmoid(-x);
      e.dsdc = wantDeriv ? -e.s * (1.0 - e.s) / par.width : 0.0;
      ns_[static_cast<std::size_t>(e.i)] += e.sigma * e.s;
    }
    double sum = 0.0;
    for (int i = 0; i < n; i++) {
      const double x = (ns_[static_cast<std::size_t>(i)] - par.mid) / par.beta;
      const double f = sigmoid(x);
      dfdn_[static_cast<std::size_t>(i)] = wantDeriv ? f * (1.0 - f) / par.beta : 0.0;
      sum += f;
    }
    fraction = sum / static_cast<double>(n);
    if (!wantDeriv) {
      return;
    }
    for (auto &g : grad) {
      g = {0.0, 0.0, 0.0};
    }
    for (int a = 0; a < 3; a++) {
      for (int b = 0; b < 3; b++) {
        virial[a][b] = 0.0;
      }
    }
    for (auto &g : gQ_) {
      g.clear();
    }
    const double invN = 1.0 / static_cast<double>(n);
    for (const auto &e : edges_) {
      const double gc = invN * dfdn_[static_cast<std::size_t>(e.i)] * e.sigma * e.dsdc;
      if (gc == 0.0) {
        continue;
      }
      const double ni2 = norm2_[static_cast<std::size_t>(e.i)];
      const double nj2 = norm2_[static_cast<std::size_t>(e.j)];
      if (ni2 < 1e-24 || nj2 < 1e-24) {
        continue;
      }
      const double scale = 1.0 / std::sqrt(ni2 * nj2);
      const auto &qi = Qlm_[static_cast<std::size_t>(e.i)];
      const auto &qj = Qlm_[static_cast<std::size_t>(e.j)];
      auto &gi = gQ_[static_cast<std::size_t>(e.i)];
      auto &gj = gQ_[static_cast<std::size_t>(e.j)];
      for (int m = 0; m < 7; m++) {
        gi.re[m] += gc * (qj.re[m] * scale - e.c * qi.re[m] / ni2);
        gi.im[m] += gc * (qj.im[m] * scale - e.c * qi.im[m] / ni2);
        gj.re[m] += gc * (qi.re[m] * scale - e.c * qj.re[m] / nj2);
        gj.im[m] += gc * (qi.im[m] * scale - e.c * qj.im[m] / nj2);
      }
    }
    std::fill(gZ_.begin(), gZ_.end(), 0.0);
    for (int i = 0; i < n; i++) {
      const double z = coord_[static_cast<std::size_t>(i)];
      auto &gu = gU_[static_cast<std::size_t>(i)];
      gu.clear();
      if (z <= 0.0) {
        continue;
      }
      const auto &gq = gQ_[static_cast<std::size_t>(i)];
      const auto &q = Qlm_[static_cast<std::size_t>(i)];
      double dot = 0.0;
      for (int m = 0; m < 7; m++) {
        gu.re[m] = gq.re[m] / z;
        gu.im[m] = gq.im[m] / z;
        dot += gq.re[m] * q.re[m] + gq.im[m] * q.im[m];
      }
      gZ_[static_cast<std::size_t>(i)] = -dot / z;
    }
    for (const auto &e : edges_) {
      const double gns = invN * dfdn_[static_cast<std::size_t>(e.i)];
      double gsig = gns * e.s + gZ_[static_cast<std::size_t>(e.i)];
      const auto &gu = gU_[static_cast<std::size_t>(e.i)];
      for (int m = 0; m < 7; m++) {
        gsig += gu.re[m] * e.y.re[m] + gu.im[m] * e.y.im[m];
      }
      double g[3] = {0.0, 0.0, 0.0};
      const double r = std::sqrt(e.dx * e.dx + e.dy * e.dy + e.dz * e.dz);
      if (r > 1e-8) {
        const double gr = gsig * e.dsdr;
        g[0] += gr * e.dx / r;
        g[1] += gr * e.dy / r;
        g[2] += gr * e.dz / r;
      }
      for (int m = 0; m < 7; m++) {
        const double gyRe = gu.re[m] * e.sigma;
        const double gyIm = gu.im[m] * e.sigma;
        for (int k = 0; k < 3; k++) {
          g[k] += gyRe * e.y.dre[m][k] + gyIm * e.y.dim[m][k];
        }
      }
      auto &gi = grad[static_cast<std::size_t>(e.i)];
      auto &gj = grad[static_cast<std::size_t>(e.j)];
      gi[0] -= g[0];
      gi[1] -= g[1];
      gi[2] -= g[2];
      gj[0] += g[0];
      gj[1] += g[1];
      gj[2] += g[2];
      const double delta[3] = {e.dx, e.dy, e.dz};
      for (int a = 0; a < 3; a++) {
        for (int b = 0; b < 3; b++) {
          virial[a][b] -= delta[a] * g[b];
        }
      }
    }
  }

private:
  void ensure(int n) {
    if (n == n_) {
      return;
    }
    n_ = n;
    const auto nn = static_cast<std::size_t>(std::max(n, 0));
    coord_.assign(nn, 0.0);
    norm2_.assign(nn, 0.0);
    ns_.assign(nn, 0.0);
    dfdn_.assign(nn, 0.0);
    gZ_.assign(nn, 0.0);
    U_.assign(nn, {});
    Qlm_.assign(nn, {});
    gQ_.assign(nn, {});
    gU_.assign(nn, {});
    grad.assign(nn, {0.0, 0.0, 0.0});
  }

  void push(int i, int j, double dx, double dy, double dz, const SmoothParams &par, bool wantDeriv) {
    const double r = std::sqrt(dx * dx + dy * dy + dz * dz);
    if (!(r > 1e-8) || r >= par.dmax) {
      return;
    }
    smooth_detail::Edge e;
    e.i = i;
    e.j = j;
    e.dx = dx;
    e.dy = dy;
    e.dz = dz;
    e.sigma = smooth_detail::cubicSwitch(r, par.d0, par.dmax, e.dsdr);
    if (e.sigma <= 0.0 && e.dsdr == 0.0) {
      return;
    }
    if (wantDeriv) {
      smooth_detail::ylm3ValueGrad(dx, dy, dz, e.y);
    } else {
      smooth_detail::ylm3Value(dx, dy, dz, e.y.re, e.y.im);
    }
    edges_.push_back(e);
  }

  int n_ = -1;
  std::vector<smooth_detail::Edge> edges_;
  std::vector<double> coord_;
  std::vector<double> norm2_;
  std::vector<double> ns_;
  std::vector<double> dfdn_;
  std::vector<double> gZ_;
  std::vector<smooth_detail::Qacc> U_;
  std::vector<smooth_detail::Qacc> Qlm_;
  std::vector<smooth_detail::Qacc> gQ_;
  std::vector<smooth_detail::Qacc> gU_;
};

inline double minimumImage(double d, double length) {
  if (!(length > 0.0)) {
    return d;
  }
  return d - length * std::round(d / length);
}

// Correlations against a hard neighbour shell. Weights are 1, so the q_3 sum
// matches the published four-nearest average when `neigh` is those four.
// The box is orthorhombic. Correlations are grouped by the central atom, in
// `neigh` order.
inline std::vector<std::vector<double>>
hardShellCorrelations(const std::vector<std::array<double, 3>> &pos,
                      const std::vector<std::vector<int>> &neigh, double lx, double ly,
                      double lz) {
  const int n = static_cast<int>(pos.size());
  std::vector<std::array<double, 7>> re(static_cast<std::size_t>(n));
  std::vector<std::array<double, 7>> im(static_cast<std::size_t>(n));
  for (int i = 0; i < n; i++) {
    re[static_cast<std::size_t>(i)].fill(0.0);
    im[static_cast<std::size_t>(i)].fill(0.0);
  }
  double yre[7];
  double yim[7];
  for (int i = 0; i < n; i++) {
    const auto &o = pos[static_cast<std::size_t>(i)];
    for (int j : neigh[static_cast<std::size_t>(i)]) {
      const auto &p = pos[static_cast<std::size_t>(j)];
      smooth_detail::ylm3Value(minimumImage(p[0] - o[0], lx), minimumImage(p[1] - o[1], ly),
                               minimumImage(p[2] - o[2], lz), yre, yim);
      for (int m = 0; m < 7; m++) {
        re[static_cast<std::size_t>(i)][static_cast<std::size_t>(m)] += yre[m];
        im[static_cast<std::size_t>(i)][static_cast<std::size_t>(m)] += yim[m];
      }
    }
  }
  std::vector<std::vector<double>> cij(static_cast<std::size_t>(n));
  for (int i = 0; i < n; i++) {
    double ni = 0.0;
    for (int m = 0; m < 7; m++) {
      ni += re[static_cast<std::size_t>(i)][static_cast<std::size_t>(m)] *
                re[static_cast<std::size_t>(i)][static_cast<std::size_t>(m)] +
            im[static_cast<std::size_t>(i)][static_cast<std::size_t>(m)] *
                im[static_cast<std::size_t>(i)][static_cast<std::size_t>(m)];
    }
    for (int j : neigh[static_cast<std::size_t>(i)]) {
      double nj = 0.0;
      double dot = 0.0;
      for (int m = 0; m < 7; m++) {
        dot += re[static_cast<std::size_t>(i)][static_cast<std::size_t>(m)] *
                   re[static_cast<std::size_t>(j)][static_cast<std::size_t>(m)] +
               im[static_cast<std::size_t>(i)][static_cast<std::size_t>(m)] *
                   im[static_cast<std::size_t>(j)][static_cast<std::size_t>(m)];
        nj += re[static_cast<std::size_t>(j)][static_cast<std::size_t>(m)] *
                  re[static_cast<std::size_t>(j)][static_cast<std::size_t>(m)] +
              im[static_cast<std::size_t>(j)][static_cast<std::size_t>(m)] *
                  im[static_cast<std::size_t>(j)][static_cast<std::size_t>(m)];
      }
      cij[static_cast<std::size_t>(i)].push_back(dot / std::sqrt(ni * nj));
    }
  }
  return cij;
}

} // namespace dseams

#endif

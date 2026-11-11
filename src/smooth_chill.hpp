/*
** Smooth CHILL+-like ice fraction.
**
** The published Nguyen and Molinero classifier is a step function of the
** bond correlation c_ij and of the neighbour count, so it has no force.
** This surrogate keeps the same q_3 correlation, replaces the hard cutoff
** and the hard staggered threshold by smooth steps, and returns a fraction
** in (0, 1) together with its Cartesian derivatives.
**
** It is not the published label. Cubic ice (four staggered bonds) and
** hexagonal ice (three) both sit above the molecule midpoint; liquid water
** does not. q_3 uses the same Condon-Shortley convention as seams-core
** (azimuth atan2(x, y)), so on a shell that is exactly the four nearest
** neighbours the correlations match chill::applyChillPlus.
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

struct SmoothValue {
  double fraction = 0.0;
  // derivative of the fraction w.r.t. each Cartesian coordinate, Angstrom
  std::vector<std::array<double, 3>> grad;
  // box derivative matching PLUMED setBoxDerivatives: -sum (delta ⊗ dV/ddelta)
  double virial[3][3] = {};
};

namespace smooth_detail {

struct Dual {
  double v = 0.0;
  double d[3] = {0.0, 0.0, 0.0};

  Dual() = default;
  explicit Dual(double x) : v(x) {}
  static Dual var(double x, int axis) {
    Dual a(x);
    a.d[axis] = 1.0;
    return a;
  }
};

inline Dual operator+(Dual a, Dual b) {
  Dual c;
  c.v = a.v + b.v;
  for (int k = 0; k < 3; k++) {
    c.d[k] = a.d[k] + b.d[k];
  }
  return c;
}
inline Dual operator-(Dual a, Dual b) {
  Dual c;
  c.v = a.v - b.v;
  for (int k = 0; k < 3; k++) {
    c.d[k] = a.d[k] - b.d[k];
  }
  return c;
}
inline Dual operator-(Dual a) {
  Dual c;
  c.v = -a.v;
  for (int k = 0; k < 3; k++) {
    c.d[k] = -a.d[k];
  }
  return c;
}
inline Dual operator*(Dual a, Dual b) {
  Dual c;
  c.v = a.v * b.v;
  for (int k = 0; k < 3; k++) {
    c.d[k] = a.d[k] * b.v + a.v * b.d[k];
  }
  return c;
}
inline Dual operator*(double a, Dual b) {
  Dual c;
  c.v = a * b.v;
  for (int k = 0; k < 3; k++) {
    c.d[k] = a * b.d[k];
  }
  return c;
}
inline Dual operator/(Dual a, Dual b) {
  Dual c;
  c.v = a.v / b.v;
  const double inv = 1.0 / (b.v * b.v);
  for (int k = 0; k < 3; k++) {
    c.d[k] = (a.d[k] * b.v - a.v * b.d[k]) * inv;
  }
  return c;
}
inline Dual sqrtDual(Dual a) {
  Dual c;
  c.v = std::sqrt(a.v);
  const double inv = 0.5 / c.v;
  for (int k = 0; k < 3; k++) {
    c.d[k] = a.d[k] * inv;
  }
  return c;
}

struct CDual {
  Dual re;
  Dual im;
};

inline CDual operator*(CDual a, CDual b) {
  return {a.re * b.re - a.im * b.im, a.re * b.im + a.im * b.re};
}
inline CDual operator*(Dual a, CDual b) { return {a * b.re, a * b.im}; }
inline CDual conj(CDual a) { return {a.re, -a.im}; }

struct Ylm3 {
  double re[7] = {};
  double im[7] = {};
  double dre[7][3] = {};
  double dim[7][3] = {};
};

// Y_3m, m = -3..+3, in the seams-core hard-coded convention.
inline Ylm3 ylm3(double x, double y, double z) {
  constexpr double pi = 3.14159265358979323846;
  const double norm[4] = {0.25 * std::sqrt(7.0 / pi), 0.125 * std::sqrt(21.0 / pi),
                          0.25 * std::sqrt(105.0 / (2.0 * pi)),
                          0.125 * std::sqrt(35.0 / pi)};
  Dual dx = Dual::var(x, 0);
  Dual dy = Dual::var(y, 1);
  Dual dz = Dual::var(z, 2);
  const Dual r = sqrtDual(dx * dx + dy * dy + dz * dz);
  const Dual rho = sqrtDual(dx * dx + dy * dy + Dual(1e-30));
  const Dual cosT = dz / r;
  const Dual sinT = rho / r;
  // phi = atan2(x, y): cos phi = y/rho, sin phi = x/rho. A bond along z
  // has no azimuth; the m != 0 harmonics vanish with sin theta.
  const bool pole = (x * x + y * y) < 1e-24;
  const Dual cosP = pole ? Dual(1.0) : dy / rho;
  const Dual sinP = pole ? Dual(0.0) : dx / rho;
  const CDual phase1{cosP, sinP};
  const CDual phase2 = phase1 * phase1;
  const CDual phase3 = phase2 * phase1;
  const CDual phase[4] = {CDual{Dual(1.0), Dual(0.0)}, phase1, phase2, phase3};

  const Dual c1 = cosT;
  const Dual c2 = cosT * cosT;
  const Dual c3 = c2 * cosT;
  const Dual s1 = sinT;
  const Dual s2 = sinT * sinT;
  const Dual s3 = s2 * sinT;
  const Dual amp[4] = {norm[0] * (5.0 * c3 - 3.0 * c1), norm[1] * (s1 * (5.0 * c2 - Dual(1.0))),
                       norm[2] * (s2 * c1), norm[3] * s3};

  Ylm3 out;
  for (int m = 0; m <= 3; m++) {
    const CDual neg = amp[m] * conj(phase[m]);
    const Dual sign(m % 2 == 0 ? 1.0 : -1.0);
    const CDual pos = (sign * amp[m]) * phase[m];
    const int in = 3 - m;
    const int ip = 3 + m;
    out.re[in] = neg.re.v;
    out.im[in] = neg.im.v;
    out.re[ip] = pos.re.v;
    out.im[ip] = pos.im.v;
    for (int k = 0; k < 3; k++) {
      out.dre[in][k] = neg.re.d[k];
      out.dim[in][k] = neg.im.d[k];
      out.dre[ip][k] = pos.re.d[k];
      out.dim[ip][k] = pos.im.d[k];
    }
  }
  return out;
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
  bool live = false;
};

} // namespace smooth_detail

inline SmoothValue smoothChill(int n, const std::vector<PairBond> &pairs,
                               const SmoothParams &par, bool wantDeriv) {
  using namespace smooth_detail;
  SmoothValue out;
  out.grad.assign(static_cast<std::size_t>(n), {0.0, 0.0, 0.0});
  if (n <= 0) {
    return out;
  }

  std::vector<Edge> edges;
  edges.reserve(pairs.size() * 2);
  auto push = [&](int i, int j, double dx, double dy, double dz) {
    const double r = std::sqrt(dx * dx + dy * dy + dz * dz);
    if (!(r > 1e-8) || r >= par.dmax) {
      return;
    }
    Edge e;
    e.i = i;
    e.j = j;
    e.dx = dx;
    e.dy = dy;
    e.dz = dz;
    e.sigma = cubicSwitch(r, par.d0, par.dmax, e.dsdr);
    if (e.sigma <= 0.0 && e.dsdr == 0.0) {
      return;
    }
    e.y = ylm3(dx, dy, dz);
    e.live = true;
    edges.push_back(e);
  };
  for (const auto &p : pairs) {
    push(p.i, p.j, p.dx, p.dy, p.dz);
    push(p.j, p.i, -p.dx, -p.dy, -p.dz);
  }

  struct Q {
    double re[7] = {};
    double im[7] = {};
  };
  std::vector<double> coord(static_cast<std::size_t>(n), 0.0);
  std::vector<Q> U(static_cast<std::size_t>(n));
  for (const auto &e : edges) {
    coord[static_cast<std::size_t>(e.i)] += e.sigma;
    auto &u = U[static_cast<std::size_t>(e.i)];
    for (int m = 0; m < 7; m++) {
      u.re[m] += e.sigma * e.y.re[m];
      u.im[m] += e.sigma * e.y.im[m];
    }
  }
  std::vector<Q> Qlm(static_cast<std::size_t>(n));
  std::vector<double> norm2(static_cast<std::size_t>(n), 0.0);
  for (int i = 0; i < n; i++) {
    const double z = coord[static_cast<std::size_t>(i)];
    if (z <= 0.0) {
      continue;
    }
    auto &q = Qlm[static_cast<std::size_t>(i)];
    const auto &u = U[static_cast<std::size_t>(i)];
    double n2 = 0.0;
    for (int m = 0; m < 7; m++) {
      q.re[m] = u.re[m] / z;
      q.im[m] = u.im[m] / z;
      n2 += q.re[m] * q.re[m] + q.im[m] * q.im[m];
    }
    norm2[static_cast<std::size_t>(i)] = n2;
  }

  std::vector<double> ns(static_cast<std::size_t>(n), 0.0);
  for (auto &e : edges) {
    const double ni2 = norm2[static_cast<std::size_t>(e.i)];
    const double nj2 = norm2[static_cast<std::size_t>(e.j)];
    if (ni2 < 1e-24 || nj2 < 1e-24) {
      e.c = 0.0;
      e.s = 0.0;
      e.dsdc = 0.0;
      continue;
    }
    const auto &qi = Qlm[static_cast<std::size_t>(e.i)];
    const auto &qj = Qlm[static_cast<std::size_t>(e.j)];
    double dot = 0.0;
    for (int m = 0; m < 7; m++) {
      dot += qi.re[m] * qj.re[m] + qi.im[m] * qj.im[m];
    }
    e.c = dot / std::sqrt(ni2 * nj2);
    const double x = (e.c - par.staggered) / par.width;
    e.s = sigmoid(-x);
    e.dsdc = -e.s * (1.0 - e.s) / par.width;
    ns[static_cast<std::size_t>(e.i)] += e.sigma * e.s;
  }

  double sum = 0.0;
  std::vector<double> molecule(static_cast<std::size_t>(n), 0.0);
  std::vector<double> dfdn(static_cast<std::size_t>(n), 0.0);
  for (int i = 0; i < n; i++) {
    const double x = (ns[static_cast<std::size_t>(i)] - par.mid) / par.beta;
    const double f = sigmoid(x);
    molecule[static_cast<std::size_t>(i)] = f;
    dfdn[static_cast<std::size_t>(i)] = f * (1.0 - f) / par.beta;
    sum += f;
  }
  out.fraction = sum / static_cast<double>(n);
  if (!wantDeriv) {
    return out;
  }

  const double invN = 1.0 / static_cast<double>(n);
  std::vector<Q> gQ(static_cast<std::size_t>(n));
  for (const auto &e : edges) {
    const double gns = invN * dfdn[static_cast<std::size_t>(e.i)];
    const double gc = gns * e.sigma * e.dsdc;
    if (gc == 0.0) {
      continue;
    }
    const double ni2 = norm2[static_cast<std::size_t>(e.i)];
    const double nj2 = norm2[static_cast<std::size_t>(e.j)];
    if (ni2 < 1e-24 || nj2 < 1e-24) {
      continue;
    }
    const double ni = std::sqrt(ni2);
    const double nj = std::sqrt(nj2);
    const auto &qi = Qlm[static_cast<std::size_t>(e.i)];
    const auto &qj = Qlm[static_cast<std::size_t>(e.j)];
    auto &gi = gQ[static_cast<std::size_t>(e.i)];
    auto &gj = gQ[static_cast<std::size_t>(e.j)];
    const double scale = 1.0 / (ni * nj);
    for (int m = 0; m < 7; m++) {
      gi.re[m] += gc * (qj.re[m] * scale - e.c * qi.re[m] / ni2);
      gi.im[m] += gc * (qj.im[m] * scale - e.c * qi.im[m] / ni2);
      gj.re[m] += gc * (qi.re[m] * scale - e.c * qj.re[m] / nj2);
      gj.im[m] += gc * (qi.im[m] * scale - e.c * qj.im[m] / nj2);
    }
  }

  std::vector<Q> gU(static_cast<std::size_t>(n));
  std::vector<double> gZ(static_cast<std::size_t>(n), 0.0);
  for (int i = 0; i < n; i++) {
    const double z = coord[static_cast<std::size_t>(i)];
    if (z <= 0.0) {
      continue;
    }
    const auto &gq = gQ[static_cast<std::size_t>(i)];
    const auto &q = Qlm[static_cast<std::size_t>(i)];
    auto &gu = gU[static_cast<std::size_t>(i)];
    double dot = 0.0;
    for (int m = 0; m < 7; m++) {
      gu.re[m] = gq.re[m] / z;
      gu.im[m] = gq.im[m] / z;
      dot += gq.re[m] * q.re[m] + gq.im[m] * q.im[m];
    }
    gZ[static_cast<std::size_t>(i)] = -dot / z;
  }

  for (const auto &e : edges) {
    const double gns = invN * dfdn[static_cast<std::size_t>(e.i)];
    double gsig = gns * e.s + gZ[static_cast<std::size_t>(e.i)];
    const auto &gu = gU[static_cast<std::size_t>(e.i)];
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
    auto &gi = out.grad[static_cast<std::size_t>(e.i)];
    auto &gj = out.grad[static_cast<std::size_t>(e.j)];
    gi[0] -= g[0];
    gi[1] -= g[1];
    gi[2] -= g[2];
    gj[0] += g[0];
    gj[1] += g[1];
    gj[2] += g[2];
    const double delta[3] = {e.dx, e.dy, e.dz};
    for (int a = 0; a < 3; a++) {
      for (int b = 0; b < 3; b++) {
        out.virial[a][b] -= delta[a] * g[b];
      }
    }
  }
  return out;
}

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
  std::vector<smooth_detail::Ylm3> acc(static_cast<std::size_t>(n));
  for (int i = 0; i < n; i++) {
    const auto &o = pos[static_cast<std::size_t>(i)];
    for (int j : neigh[static_cast<std::size_t>(i)]) {
      const auto &p = pos[static_cast<std::size_t>(j)];
      const auto y = smooth_detail::ylm3(minimumImage(p[0] - o[0], lx),
                                         minimumImage(p[1] - o[1], ly),
                                         minimumImage(p[2] - o[2], lz));
      auto &a = acc[static_cast<std::size_t>(i)];
      for (int m = 0; m < 7; m++) {
        a.re[m] += y.re[m];
        a.im[m] += y.im[m];
      }
    }
  }
  std::vector<std::vector<double>> cij(static_cast<std::size_t>(n));
  for (int i = 0; i < n; i++) {
    double ni = 0.0;
    for (int m = 0; m < 7; m++) {
      ni += acc[static_cast<std::size_t>(i)].re[m] * acc[static_cast<std::size_t>(i)].re[m] +
            acc[static_cast<std::size_t>(i)].im[m] * acc[static_cast<std::size_t>(i)].im[m];
    }
    for (int j : neigh[static_cast<std::size_t>(i)]) {
      double nj = 0.0;
      double dot = 0.0;
      for (int m = 0; m < 7; m++) {
        dot += acc[static_cast<std::size_t>(i)].re[m] * acc[static_cast<std::size_t>(j)].re[m] +
               acc[static_cast<std::size_t>(i)].im[m] * acc[static_cast<std::size_t>(j)].im[m];
        nj += acc[static_cast<std::size_t>(j)].re[m] * acc[static_cast<std::size_t>(j)].re[m] +
              acc[static_cast<std::size_t>(j)].im[m] * acc[static_cast<std::size_t>(j)].im[m];
      }
      cij[static_cast<std::size_t>(i)].push_back(dot / std::sqrt(ni * nj));
    }
  }
  return cij;
}

} // namespace dseams

#endif

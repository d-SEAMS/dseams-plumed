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

inline const double *ylmNorms() {
  constexpr double pi = 3.14159265358979323846;
  static const double n[4] = {0.25 * std::sqrt(7.0 / pi), 0.125 * std::sqrt(21.0 / pi),
                              0.25 * std::sqrt(105.0 / (2.0 * pi)),
                              0.125 * std::sqrt(35.0 / pi)};
  return n;
}

// Value of Y_3m, m = -3..+3. Matches seams-core sph::spheriHarmo for l = 3.
inline void ylm3Value(double x, double y, double z, double re[7], double im[7]) {
  const double *nrm = ylmNorms();
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
  const double *nrm = ylmNorms();
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

} // namespace smooth_detail

// Reused scratch for one group of atoms. evaluate() keeps its buffers
// between calls. nthreads is PLUMED's NUM_THREADS. Every sum that a
// molecule owns walks that molecule's bonds in pair-list order, and the
// fraction and the virial are left folds, so one thread and many threads
// write the same bits.
class SmoothChill {
public:
  double fraction = 0.0;
  std::vector<std::array<double, 3>> grad;
  double virial[3][3] = {};

  void evaluate(int n, const std::vector<PairBond> &pairs, const SmoothParams &par, bool wantDeriv,
                int nthreads = 1) {
    using namespace smooth_detail;
    ensure(n);
    fraction = 0.0;
    if (n <= 0) {
      return;
    }
    const int nt = std::max(nthreads, 1);
    const int np = static_cast<int>(pairs.size());
    slot_.resize(static_cast<std::size_t>(np));
    keep_.resize(static_cast<std::size_t>(np));
    const SmoothParams parm = par;

#pragma omp parallel num_threads(nt)
    {
#pragma omp for schedule(static)
      for (int p = 0; p < np; p++) {
        fillSlot(pairs[static_cast<std::size_t>(p)], slot_[static_cast<std::size_t>(p)],
                 keep_[static_cast<std::size_t>(p)], parm, wantDeriv);
      }
#pragma omp single
      { compact(n); }
#pragma omp for schedule(static)
      for (int a = 0; a < n; a++) {
        coord_[static_cast<std::size_t>(a)] = 0.0;
        U_[static_cast<std::size_t>(a)].clear();
        for (int k = off_[static_cast<std::size_t>(a)]; k < off_[static_cast<std::size_t>(a + 1)]; k++) {
          const Bond &b = bondAt(adj_[static_cast<std::size_t>(k)]);
          const double sign = (a == b.i) ? 1.0 : -1.0;
          coord_[static_cast<std::size_t>(a)] += b.sigma;
          auto &u = U_[static_cast<std::size_t>(a)];
          for (int m = 0; m < 7; m++) {
            const double yr = sign * b.re[m];
            const double yi = sign * b.im[m];
            u.re[m] += b.sigma * yr;
            u.im[m] += b.sigma * yi;
          }
        }
      }
#pragma omp for schedule(static)
      for (int a = 0; a < n; a++) {
        const double z = coord_[static_cast<std::size_t>(a)];
        auto &q = Qlm_[static_cast<std::size_t>(a)];
        q.clear();
        norm2_[static_cast<std::size_t>(a)] = 0.0;
        if (z <= 0.0) {
          continue;
        }
        const auto &u = U_[static_cast<std::size_t>(a)];
        double n2 = 0.0;
        for (int m = 0; m < 7; m++) {
          q.re[m] = u.re[m] / z;
          q.im[m] = u.im[m] / z;
          n2 += q.re[m] * q.re[m] + q.im[m] * q.im[m];
        }
        norm2_[static_cast<std::size_t>(a)] = n2;
      }
#pragma omp for schedule(static)
      for (int p = 0; p < nb_; p++) {
        bondCorrelation(bondAt(p), parm, wantDeriv);
      }
#pragma omp for schedule(static)
      for (int a = 0; a < n; a++) {
        double ns = 0.0;
        for (int k = off_[static_cast<std::size_t>(a)]; k < off_[static_cast<std::size_t>(a + 1)]; k++) {
          const Bond &b = bondAt(adj_[static_cast<std::size_t>(k)]);
          if (b.live) {
            ns += b.sigma * b.s;
          }
        }
        ns_[static_cast<std::size_t>(a)] = ns;
      }
    }

    double sum = 0.0;
    for (int i = 0; i < n; i++) {
      const double x = (ns_[static_cast<std::size_t>(i)] - parm.mid) / parm.beta;
      const double f = sigmoid(x);
      dfdn_[static_cast<std::size_t>(i)] = wantDeriv ? f * (1.0 - f) / parm.beta : 0.0;
      sum += f;
    }
    fraction = sum / static_cast<double>(n);
    if (!wantDeriv) {
      return;
    }
    const double invN = 1.0 / static_cast<double>(n);
    gFwd_.resize(static_cast<std::size_t>(nb_));
    gRev_.resize(static_cast<std::size_t>(nb_));

#pragma omp parallel num_threads(nt)
    {
#pragma omp for schedule(static)
      for (int a = 0; a < n; a++) {
        auto &gq = gQ_[static_cast<std::size_t>(a)];
        gq.clear();
        const double n2a = norm2_[static_cast<std::size_t>(a)];
        if (!(n2a >= 1e-24)) {
          continue;
        }
        const auto &qa = Qlm_[static_cast<std::size_t>(a)];
        for (int k = off_[static_cast<std::size_t>(a)]; k < off_[static_cast<std::size_t>(a + 1)]; k++) {
          const Bond &b = bondAt(adj_[static_cast<std::size_t>(k)]);
          if (!b.live) {
            continue;
          }
          const int other = (a == b.i) ? b.j : b.i;
          const double n2b = norm2_[static_cast<std::size_t>(other)];
          if (n2b < 1e-24) {
            continue;
          }
          const double scale = 1.0 / std::sqrt(n2a * n2b);
          const auto &qo = Qlm_[static_cast<std::size_t>(other)];
          const double gc1 = invN * dfdn_[static_cast<std::size_t>(b.i)] * b.sigma * b.dsdc;
          const double gc2 = invN * dfdn_[static_cast<std::size_t>(b.j)] * b.sigma * b.dsdc;
          addGQ(gq, gc1, qa, qo, n2a, scale, b.c);
          addGQ(gq, gc2, qa, qo, n2a, scale, b.c);
        }
      }
#pragma omp for schedule(static)
      for (int a = 0; a < n; a++) {
        const double z = coord_[static_cast<std::size_t>(a)];
        auto &gu = gU_[static_cast<std::size_t>(a)];
        gu.clear();
        gZ_[static_cast<std::size_t>(a)] = 0.0;
        if (z <= 0.0) {
          continue;
        }
        const auto &gq = gQ_[static_cast<std::size_t>(a)];
        const auto &q = Qlm_[static_cast<std::size_t>(a)];
        double dot = 0.0;
        for (int m = 0; m < 7; m++) {
          gu.re[m] = gq.re[m] / z;
          gu.im[m] = gq.im[m] / z;
          dot += gq.re[m] * q.re[m] + gq.im[m] * q.im[m];
        }
        gZ_[static_cast<std::size_t>(a)] = -dot / z;
      }
#pragma omp for schedule(static)
      for (int p = 0; p < nb_; p++) {
        forceOnBond(bondAt(p), invN, gFwd_[static_cast<std::size_t>(p)],
                    gRev_[static_cast<std::size_t>(p)]);
      }
#pragma omp for schedule(static)
      for (int a = 0; a < n; a++) {
        auto &ga = grad[static_cast<std::size_t>(a)];
        ga = {0.0, 0.0, 0.0};
        for (int k = off_[static_cast<std::size_t>(a)]; k < off_[static_cast<std::size_t>(a + 1)]; k++) {
          const int p = adj_[static_cast<std::size_t>(k)];
          const Bond &b = bondAt(p);
          const auto &gf = gFwd_[static_cast<std::size_t>(p)];
          const auto &gr = gRev_[static_cast<std::size_t>(p)];
          if (a == b.i) {
            ga[0] -= gf[0];
            ga[1] -= gf[1];
            ga[2] -= gf[2];
            ga[0] += gr[0];
            ga[1] += gr[1];
            ga[2] += gr[2];
          } else {
            ga[0] += gf[0];
            ga[1] += gf[1];
            ga[2] += gf[2];
            ga[0] -= gr[0];
            ga[1] -= gr[1];
            ga[2] -= gr[2];
          }
        }
      }
    }

    for (int a = 0; a < 3; a++) {
      for (int b = 0; b < 3; b++) {
        virial[a][b] = 0.0;
      }
    }
    for (int p = 0; p < nb_; p++) {
      const Bond &b = bondAt(p);
      const auto &gf = gFwd_[static_cast<std::size_t>(p)];
      const auto &gr = gRev_[static_cast<std::size_t>(p)];
      const double d[3] = {b.dx, b.dy, b.dz};
      for (int a = 0; a < 3; a++) {
        for (int c = 0; c < 3; c++) {
          virial[a][c] -= d[a] * gf[c];
        }
      }
      for (int a = 0; a < 3; a++) {
        for (int c = 0; c < 3; c++) {
          virial[a][c] -= (-d[a]) * gr[c];
        }
      }
    }
  }

private:
  struct Bond {
    int i = 0;
    int j = 0;
    char live = 0;
    double dx = 0.0;
    double dy = 0.0;
    double dz = 0.0;
    double sigma = 0.0;
    double dsdr = 0.0;
    double c = 0.0;
    double s = 0.0;
    double dsdc = 0.0;
    double re[7] = {};
    double im[7] = {};
    double dre[7][3] = {};
    double dim[7][3] = {};
  };

  void fillSlot(const PairBond &p, Bond &b, char &keep, const SmoothParams &par, bool wantDeriv) const {
    using namespace smooth_detail;
    keep = 0;
    const double dx = p.dx;
    const double dy = p.dy;
    const double dz = p.dz;
    const double r = std::sqrt(dx * dx + dy * dy + dz * dz);
    if (!(r > 1e-8) || r >= par.dmax) {
      return;
    }
    double dsdr = 0.0;
    const double sigma = cubicSwitch(r, par.d0, par.dmax, dsdr);
    if (sigma <= 0.0 && dsdr == 0.0) {
      return;
    }
    keep = 1;
    b.i = p.i;
    b.j = p.j;
    b.live = 0;
    b.dx = dx;
    b.dy = dy;
    b.dz = dz;
    b.sigma = sigma;
    b.dsdr = dsdr;
    b.c = 0.0;
    b.s = 0.0;
    b.dsdc = 0.0;
    if (wantDeriv) {
      Ylm3 y;
      ylm3ValueGrad(dx, dy, dz, y);
      for (int m = 0; m < 7; m++) {
        b.re[m] = y.re[m];
        b.im[m] = y.im[m];
        for (int k = 0; k < 3; k++) {
          b.dre[m][k] = y.dre[m][k];
          b.dim[m][k] = y.dim[m][k];
        }
      }
    } else {
      ylm3Value(dx, dy, dz, b.re, b.im);
    }
  }

  Bond &bondAt(int p) { return slot_[static_cast<std::size_t>(kept_[static_cast<std::size_t>(p)])]; }
  const Bond &bondAt(int p) const {
    return slot_[static_cast<std::size_t>(kept_[static_cast<std::size_t>(p)])];
  }

  void compact(int n) {
    kept_.clear();
    if (kept_.capacity() < slot_.size()) {
      kept_.reserve(slot_.size());
    }
    for (std::size_t p = 0; p < slot_.size(); p++) {
      if (keep_[p]) {
        kept_.push_back(static_cast<int>(p));
      }
    }
    nb_ = static_cast<int>(kept_.size());
    counts_.assign(static_cast<std::size_t>(n), 0);
    for (int p = 0; p < nb_; p++) {
      const Bond &b = bondAt(p);
      counts_[static_cast<std::size_t>(b.i)]++;
      counts_[static_cast<std::size_t>(b.j)]++;
    }
    off_.resize(static_cast<std::size_t>(n + 1));
    off_[0] = 0;
    for (int a = 0; a < n; a++) {
      off_[static_cast<std::size_t>(a + 1)] =
          off_[static_cast<std::size_t>(a)] + counts_[static_cast<std::size_t>(a)];
    }
    adj_.resize(static_cast<std::size_t>(off_[static_cast<std::size_t>(n)]));
    std::vector<int> mark(off_.begin(), off_.end() - 1);
    for (int p = 0; p < nb_; p++) {
      const Bond &b = bondAt(p);
      adj_[static_cast<std::size_t>(mark[static_cast<std::size_t>(b.i)]++)] = p;
      adj_[static_cast<std::size_t>(mark[static_cast<std::size_t>(b.j)]++)] = p;
    }
  }

  void bondCorrelation(Bond &b, const SmoothParams &par, bool wantDeriv) const {
    const double ni2 = norm2_[static_cast<std::size_t>(b.i)];
    const double nj2 = norm2_[static_cast<std::size_t>(b.j)];
    if (ni2 < 1e-24 || nj2 < 1e-24) {
      b.live = 0;
      b.c = 0.0;
      b.s = 0.0;
      b.dsdc = 0.0;
      return;
    }
    const auto &qi = Qlm_[static_cast<std::size_t>(b.i)];
    const auto &qj = Qlm_[static_cast<std::size_t>(b.j)];
    double dot = 0.0;
    for (int m = 0; m < 7; m++) {
      dot += qi.re[m] * qj.re[m] + qi.im[m] * qj.im[m];
    }
    b.c = dot / std::sqrt(ni2 * nj2);
    const double x = (b.c - par.staggered) / par.width;
    b.s = smooth_detail::sigmoid(-x);
    b.dsdc = wantDeriv ? -b.s * (1.0 - b.s) / par.width : 0.0;
    b.live = 1;
  }

  static void addGQ(smooth_detail::Qacc &g, double gc, const smooth_detail::Qacc &qSelf,
                    const smooth_detail::Qacc &qOther, double n2self, double scale, double c) {
    if (gc == 0.0) {
      return;
    }
    for (int m = 0; m < 7; m++) {
      g.re[m] += gc * (qOther.re[m] * scale - c * qSelf.re[m] / n2self);
      g.im[m] += gc * (qOther.im[m] * scale - c * qSelf.im[m] / n2self);
    }
  }

  void forceOnBond(const Bond &b, double invN, std::array<double, 3> &gf, std::array<double, 3> &gr) const {
    directedForce(b, b.i, b.re, b.im, 1.0, invN, gf);
    double nre[7];
    double nim[7];
    for (int m = 0; m < 7; m++) {
      nre[m] = -b.re[m];
      nim[m] = -b.im[m];
    }
    directedForce(b, b.j, nre, nim, -1.0, invN, gr);
  }

  // sign is +1 along (dx, dy, dz) and -1 along the opposite bond. The
  // Cartesian gradient of Y_3m is the same in both directions because l = 3
  // is odd: dY(-r)/d(-r) = dY(r)/dr.
  void directedForce(const Bond &b, int central, const double re[7], const double im[7], double sign,
                     double invN, std::array<double, 3> &g) const {
    const double gns = invN * dfdn_[static_cast<std::size_t>(central)];
    double gsig = gns * b.s + gZ_[static_cast<std::size_t>(central)];
    const auto &gu = gU_[static_cast<std::size_t>(central)];
    for (int m = 0; m < 7; m++) {
      gsig += gu.re[m] * re[m] + gu.im[m] * im[m];
    }
    g = {0.0, 0.0, 0.0};
    const double r = std::sqrt(b.dx * b.dx + b.dy * b.dy + b.dz * b.dz);
    if (r > 1e-8) {
      const double gr = gsig * b.dsdr;
      g[0] += gr * (sign * b.dx) / r;
      g[1] += gr * (sign * b.dy) / r;
      g[2] += gr * (sign * b.dz) / r;
    }
    for (int m = 0; m < 7; m++) {
      const double gyRe = gu.re[m] * b.sigma;
      const double gyIm = gu.im[m] * b.sigma;
      for (int k = 0; k < 3; k++) {
        g[k] += gyRe * b.dre[m][k] + gyIm * b.dim[m][k];
      }
    }
  }

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

  int n_ = -1;
  int nb_ = 0;
  std::vector<Bond> slot_;
  std::vector<char> keep_;
  std::vector<int> kept_;
  std::vector<int> counts_;
  std::vector<int> off_;
  std::vector<int> adj_;
  std::vector<double> coord_;
  std::vector<double> norm2_;
  std::vector<double> ns_;
  std::vector<double> dfdn_;
  std::vector<double> gZ_;
  std::vector<smooth_detail::Qacc> U_;
  std::vector<smooth_detail::Qacc> Qlm_;
  std::vector<smooth_detail::Qacc> gQ_;
  std::vector<smooth_detail::Qacc> gU_;
  std::vector<std::array<double, 3>> gFwd_;
  std::vector<std::array<double, 3>> gRev_;
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

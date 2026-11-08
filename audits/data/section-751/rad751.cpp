// Audit 751: semiclassical (Crisp-Jaynes type) radiation of the relative-motion matter-wave field of hydrogen (3D, modes 1s, 2s, 2px, 2py, 2pz;
// atomic units, mu = 1).  i dc/dt = H0 c - F.(r c), F = (2/3) K alpha^3 d^3<r>/dt^3, d^3<r>/dt^3 = c^+ M c, M_nm = (i w_nm)^3 r_nm.
// argv: K s0 tmax_in_tau [dt=0.02] [every=200]; output: t P2p P1s P2s dx dy
#include <cmath>
#include <complex>
#include <cstdio>
#include <cstdlib>
typedef std::complex<double> C;
int main(int argc, char** argv) {
  const double K = std::atof(argv[1]), s0 = std::atof(argv[2]), tmaxTau = std::atof(argv[3]);
  const double dt = argc > 4 ? std::atof(argv[4]) : 0.02; const long every = argc > 5 ? std::atol(argv[5]) : 200;
  const double alpha = 7.2973525643e-3, a3 = alpha * alpha * alpha;
  const double E[5] = {-0.5, -0.125, -0.125, -0.125, -0.125};          // 1s 2s 2px 2py 2pz
  const double R1 = 128 * std::sqrt(6.0) / 243 / std::sqrt(3.0), R2 = -3 * std::sqrt(3.0) / std::sqrt(3.0);
  double r[3][5][5] = {};                                                // x, y, z
  for (int k = 0; k < 3; ++k) { r[k][0][2 + k] = r[k][2 + k][0] = R1; r[k][1][2 + k] = r[k][2 + k][1] = R2; }
  C M[3][5][5];
  for (int k = 0; k < 3; ++k) for (int n = 0; n < 5; ++n) for (int m = 0; m < 5; ++m) { const C iw(0, E[n] - E[m]); M[k][n][m] = iw * iw * iw * r[k][n][m]; }
  const double w = 0.375, d2 = R1 * R1; const double GF = 4.0 / 3 * a3 * w * w * w * d2, tmax = tmaxTau / (K * GF);
  C c[5] = {std::sqrt(s0), 0, -std::sqrt((1 - s0) / 2), C(0, -std::sqrt((1 - s0) / 2)), 0};   // 2p m=+1 = -(px + i py)/sqrt2
  auto rhs = [&](const C* y, C* f) {
    double F[3]; for (int k = 0; k < 3; ++k) { C s = 0; for (int n = 0; n < 5; ++n) for (int m = 0; m < 5; ++m) s += std::conj(y[n]) * M[k][n][m] * y[m]; F[k] = 2.0 / 3 * K * a3 * s.real(); }
    for (int n = 0; n < 5; ++n) { C h = E[n] * y[n]; for (int k = 0; k < 3; ++k) for (int m = 0; m < 5; ++m) h -= F[k] * r[k][n][m] * y[m]; f[n] = C(0, -1) * h; } };
  std::printf("# 751 K %.3e s0 %.6e GammaFermi(mu=1) %.10e tmax %.6e dt %.4f\n# t P2p P1s P2s dx dy\n", K, s0, GF, tmax, dt);
  C k1[5], k2[5], k3[5], k4[5], y[5]; long step = 0;
  for (double t = 0; t <= tmax; t += dt, ++step) {
    if (step % every == 0) { double dx = 0, dy = 0; for (int n = 0; n < 5; ++n) for (int m = 0; m < 5; ++m) { dx += (std::conj(c[n]) * r[0][n][m] * c[m]).real(); dy += (std::conj(c[n]) * r[1][n][m] * c[m]).real(); }
      std::printf("%.6f %.12e %.12e %.6e %.9e %.9e\n", t, std::norm(c[2]) + std::norm(c[3]) + std::norm(c[4]), std::norm(c[0]), std::norm(c[1]), dx, dy); }
    rhs(c, k1); for (int n = 0; n < 5; ++n) y[n] = c[n] + 0.5 * dt * k1[n]; rhs(y, k2); for (int n = 0; n < 5; ++n) y[n] = c[n] + 0.5 * dt * k2[n];
    rhs(y, k3); for (int n = 0; n < 5; ++n) y[n] = c[n] + dt * k3[n]; rhs(y, k4);
    for (int n = 0; n < 5; ++n) c[n] += dt / 6 * (k1[n] + 2.0 * k2[n] + 2.0 * k3[n] + k4[n]); }
  return 0; }

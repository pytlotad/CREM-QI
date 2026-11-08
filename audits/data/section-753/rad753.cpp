// Audit 753: rad751 (semiclassical radiation of the relative-motion field, modes 1s 2s 2px 2py 2pz, a.u., mu = 1) plus a classical random
// zero-point field in the dipole approximation, H' = +E_zp(t).r.  Each component: sum of N cosines, omega uniform in [w21 - B/2, w21 + B/2],
// random phase, A_k = sqrt(2 S(omega_k) dOmega), S(omega) = K (2/(3 pi)) alpha^3 omega^3 (one-sided, per component).  Each run deterministic,
// generator seeded with the run number.
// argv: K start(0 = 2p m=+1, 1 = 1s) zpf(0/1) s0 runs tmax_in_tau [dt=0.05] [N=1500] [band=0.1] [samples_per_tau=20]
// output: "# TE run t_e(tau; -1 = censored)" and rows "run t(tau) P2p P1s P2s"
#include <cmath>
#include <complex>
#include <cstdio>
#include <cstdlib>
#include <random>
#include <string>
#include <vector>
typedef std::complex<double> C;
int main(int argc, char** argv) {
  const double K = std::atof(argv[1]); const int start = std::atoi(argv[2]), zpf = std::atoi(argv[3]); const double s0 = std::atof(argv[4]);
  const int runs = std::atoi(argv[5]); const double tmaxTau = std::atof(argv[6]);
  const double dt = argc > 7 ? std::atof(argv[7]) : 0.05; const int N = argc > 8 ? std::atoi(argv[8]) : 1500;
  const double band = argc > 9 ? std::atof(argv[9]) : 0.1; const int spt = argc > 10 ? std::atoi(argv[10]) : 20;
  const double alpha = 7.2973525643e-3, a3 = alpha * alpha * alpha, PI = 3.14159265358979323846;
  const double E[5] = {-0.5, -0.125, -0.125, -0.125, -0.125};
  const double R1 = 128 * std::sqrt(6.0) / 243 / std::sqrt(3.0), R2 = -3 * std::sqrt(3.0) / std::sqrt(3.0);
  double r[3][5][5] = {};
  for (int k = 0; k < 3; ++k) { r[k][0][2 + k] = r[k][2 + k][0] = R1; r[k][1][2 + k] = r[k][2 + k][1] = R2; }
  C M[3][5][5];
  for (int k = 0; k < 3; ++k) for (int n = 0; n < 5; ++n) for (int m = 0; m < 5; ++m) { const C iw(0, E[n] - E[m]); M[k][n][m] = iw * iw * iw * r[k][n][m]; }
  const double w = 0.375, GF = 4.0 / 3 * a3 * w * w * w * R1 * R1, GK = K * GF, tmax = tmaxTau / GK;
  const long nsteps = long(tmax / dt + 0.5); const double sampleDt = 1.0 / (GK * spt);
  std::printf("# 753 K %.3e start %d zpf %d s0 %.3e runs %d tmax %.6e tau_K %.6e dt %.4f N %d band %.4f\n# run t/tau P2p P1s P2s\n",
              K, start, zpf, s0, runs, tmax, 1 / GK, dt, N, band);
  std::vector<std::string> outs(runs);
#pragma omp parallel for schedule(dynamic)
  for (int run = 0; run < runs; ++run) {
    std::mt19937_64 g(run); std::uniform_real_distribution<double> U(0, 1);
    std::vector<double> zr(3 * N), zi(3 * N), cr(3 * N), sr(3 * N);
    const double dOm = band / N;
    for (int j = 0; j < 3 * N; ++j) {
      const double om = w - band / 2 + band * U(g), ph = 2 * PI * U(g);
      const double A = zpf ? std::sqrt(2 * K * 2.0 / (3 * PI) * a3 * om * om * om * dOm) : 0.0;
      zr[j] = A * std::cos(ph); zi[j] = A * std::sin(ph); cr[j] = std::cos(om * dt / 2); sr[j] = std::sin(om * dt / 2); }
    auto field = [&](double* Ez) { for (int k = 0; k < 3; ++k) { double s = 0; const double* p = &zr[k * N]; for (int j = 0; j < N; ++j) s += p[j]; Ez[k] = s; } };
    auto rotate = [&]() { for (int j = 0; j < 3 * N; ++j) { const double a = zr[j] * cr[j] - zi[j] * sr[j], b = zr[j] * sr[j] + zi[j] * cr[j]; zr[j] = a; zi[j] = b; } };
    C c[5];
    if (start == 0) { c[0] = std::sqrt(s0); c[1] = 0; c[2] = -std::sqrt((1 - s0) / 2); c[3] = C(0, -std::sqrt((1 - s0) / 2)); c[4] = 0; }
    else { c[0] = 1; c[1] = c[2] = c[3] = c[4] = 0; }
    auto rhs = [&](const C* y, const double* Ez, C* f) {
      double F[3]; for (int k = 0; k < 3; ++k) { C s = 0; for (int n = 0; n < 5; ++n) for (int m = 0; m < 5; ++m) s += std::conj(y[n]) * M[k][n][m] * y[m]; F[k] = 2.0 / 3 * K * a3 * s.real() - Ez[k]; }
      for (int n = 0; n < 5; ++n) { C h = E[n] * y[n]; for (int k = 0; k < 3; ++k) for (int m = 0; m < 5; ++m) h -= F[k] * r[k][n][m] * y[m]; f[n] = C(0, -1) * h; } };
    std::string o; char buf[160]; double te = -1, nextS = 0; double E0[3], Eh[3], E1[3]; field(E0);
    C k1[5], k2[5], k3[5], k4[5], y[5];
    for (long st = 0; st <= nsteps; ++st) {
      const double t = st * dt; const double P2p = std::norm(c[2]) + std::norm(c[3]) + std::norm(c[4]);
      if (te < 0 && start == 0 && P2p < 0.5) te = t * GK;
      if (t >= nextS - 1e-9) { std::snprintf(buf, sizeof buf, "%d %.4f %.6e %.6e %.6e\n", run, t * GK, P2p, std::norm(c[0]), std::norm(c[1])); o += buf; nextS += sampleDt; }
      if (st == nsteps) break;
      rotate(); field(Eh); rotate(); field(E1);
      rhs(c, E0, k1); for (int n = 0; n < 5; ++n) y[n] = c[n] + 0.5 * dt * k1[n]; rhs(y, Eh, k2); for (int n = 0; n < 5; ++n) y[n] = c[n] + 0.5 * dt * k2[n];
      rhs(y, Eh, k3); for (int n = 0; n < 5; ++n) y[n] = c[n] + dt * k3[n]; rhs(y, E1, k4);
      for (int n = 0; n < 5; ++n) c[n] += dt / 6 * (k1[n] + 2.0 * k2[n] + 2.0 * k3[n] + k4[n]);
      for (int k = 0; k < 3; ++k) E0[k] = E1[k]; }
    std::snprintf(buf, sizeof buf, "# TE %d %.6f\n", run, te); outs[run] = std::string(buf) + o; }
  for (auto& s : outs) std::fputs(s.c_str(), stdout);
  return 0; }

// Audit 754 (run 2): Bohmian field trajectories, two-level atom + N quantized modes (RWA, single excitation; omega_0 = 1).  Same model as
// bohm754.py (cross-check), adaptive step h = min(dt_max, f |c_e|^2) common to all runs (prereg754, uzupelnienie 2).  Amplitudes: RK4 in the
// interaction picture with two half steps per step (values at t, t + h/2, t + h); configurations: RK4 per run (OpenMP over runs).
// Guidance: dq_k/dt = [Im(S* u_k) - b_k Im(c_e* S)] / (|c_e|^2 + |S|^2 + p_g0),  S = sum u_k q_k,  u_k = c_k sqrt(2 w_k),  b_k = 2 g / sqrt(2 w_k).
// argv: N M start(e|g) dt_max f_ad [m0 = 0: runs m0 .. m0+M-1, generator seed = run number; single-threaded processes in parallel] ; output (stdout): "# ..." header, rows "t/tau |c_e|^2 w_0 ... w_{M-1}" every 0.05 tau, "# TE i t_e", "# CR i n"
#include <cmath>
#include <complex>
#include <cstdio>
#include <algorithm>
#include <cstdlib>
#include <random>
#include <vector>
typedef std::complex<double> C;
int main(int argc, char** argv) {
  const int N = std::atoi(argv[1]), M = std::atoi(argv[2]); const bool E = argv[3][0] == 'e'; const double dtmax = std::atof(argv[4]), fad = std::atof(argv[5]); const int m0 = argc > 6 ? std::atoi(argv[6]) : 0;
  const double Gm = 0.01, W = 40 * Gm, w0 = 1.0, tau = 1 / Gm, T = 5 * tau, PI = 3.14159265358979323846;
  std::vector<double> om(N), sq(N), bk(N); const double dom = W / N, g = std::sqrt(Gm * dom / (2 * PI)), pg0 = E ? 0.0 : 1.0;
  for (int k = 0; k < N; ++k) { om[k] = w0 - W / 2 + W * (k + 0.5) / N; sq[k] = std::sqrt(2 * om[k]); bk[k] = 2 * g / sq[k]; }
  C ae = E ? 1.0 : 0.0; std::vector<C> ak(N, 0.0);
  auto famp = [&](double t, C a, const std::vector<C>& b, C& da, std::vector<C>& db) {
    C s = 0; for (int k = 0; k < N; ++k) { const C ph = std::polar(1.0, -(om[k] - w0) * t); s += b[k] * ph; db[k] = C(0, -1) * g * a * std::conj(ph); } da = C(0, -1) * g * s; };
  auto stepAmp = [&](double t, double h) {   // one RK4 step of the amplitudes
    std::vector<C> d1(N), d2(N), d3(N), d4(N), y(N); C a1, a2, a3, a4;
    famp(t, ae, ak, a1, d1); for (int k = 0; k < N; ++k) y[k] = ak[k] + h / 2 * d1[k]; famp(t + h / 2, ae + h / 2 * a1, y, a2, d2);
    for (int k = 0; k < N; ++k) y[k] = ak[k] + h / 2 * d2[k]; famp(t + h / 2, ae + h / 2 * a2, y, a3, d3);
    for (int k = 0; k < N; ++k) y[k] = ak[k] + h * d3[k]; famp(t + h, ae + h * a3, y, a4, d4);
    ae += h / 6 * (a1 + 2.0 * a2 + 2.0 * a3 + a4); for (int k = 0; k < N; ++k) ak[k] += h / 6 * (d1[k] + 2.0 * d2[k] + 2.0 * d3[k] + d4[k]); };
  struct Snap { std::vector<double> ur, ui; double cr, ci; };   // real/imag parts (vectorizable; no complex-multiply library calls)
  auto snap = [&](double t) { Snap s; s.ur.resize(N); s.ui.resize(N);
    for (int k = 0; k < N; ++k) { const C u = ak[k] * std::polar(1.0, -om[k] * t) * sq[k]; s.ur[k] = u.real(); s.ui[k] = u.imag(); }
    const C ce = ae * std::polar(1.0, -w0 * t); s.cr = ce.real(); s.ci = ce.imag(); return s; };
  std::vector<std::vector<double>> q(M, std::vector<double>(N));
  for (int m = 0; m < M; ++m) { std::mt19937_64 gen(m0 + m); std::normal_distribution<double> nd(0, 1); for (int k = 0; k < N; ++k) q[m][k] = nd(gen) / std::sqrt(2 * om[k]); }
  auto vel = [&](const Snap& s, const double* x, double* v) {
    const double* ur = s.ur.data(); const double* ui = s.ui.data(); double Sr = 0, Si = 0;
    for (int k = 0; k < N; ++k) { Sr += ur[k] * x[k]; Si += ui[k] * x[k]; }
    const double S2 = Sr * Sr + Si * Si, den = s.cr * s.cr + s.ci * s.ci + S2 + pg0, a = s.cr * Si - s.ci * Sr, id = 1 / den;   // a = Im(c_e* S)
    if (v) { const double* b = bk.data(); for (int k = 0; k < N; ++k) v[k] = (Sr * ui[k] - Si * ur[k] - b[k] * a) * id; }      // Im(S* u_k)
    return (S2 + pg0) * id; };
  std::vector<double> wprev(M), te(M, -1.0); std::vector<int> cr(M, 0);
  Snap s0 = snap(0); for (int m = 0; m < M; ++m) wprev[m] = vel(s0, q[m].data(), nullptr);
  std::printf("# 754 m0 %d N %d M %d start %c dt_max %.4f f %.4f Gamma %.4f band %.4f\n", m0, N, M, E ? 'e' : 'g', dtmax, fad, Gm, W);
  auto row = [&](double tt, double pe, const std::vector<double>& w) { std::printf("%.4f %.8e", tt, pe); for (double x : w) std::printf(" %.6f", x); std::printf("\n"); };
  row(0, std::norm(ae), wprev);
  double t = 0, nextG = 0.05; long nst = 0; std::vector<double> wnew(M); double vmax = 0;
  while (t < T - 1e-12) {
    const double h = std::min({dtmax, E ? fad * std::norm(ae) : dtmax, T - t});
    Snap sa = snap(t); stepAmp(t, h / 2); Snap sb = snap(t + h / 2); stepAmp(t + h / 2, h / 2); Snap sc = snap(t + h);
    double vm = 0;
#pragma omp parallel reduction(max:vm)
    { thread_local std::vector<double> k1, k2, k3, k4, y; if ((int)k1.size() != N) { k1.resize(N); k2.resize(N); k3.resize(N); k4.resize(N); y.resize(N); }
#pragma omp for schedule(static)
    for (int m = 0; m < M; ++m) { double* x = q[m].data();
      vel(sa, x, k1.data()); for (int k = 0; k < N; ++k) { y[k] = x[k] + h / 2 * k1[k]; vm = std::max(vm, std::fabs(k1[k])); }
      vel(sb, y.data(), k2.data()); for (int k = 0; k < N; ++k) y[k] = x[k] + h / 2 * k2[k];
      vel(sb, y.data(), k3.data()); for (int k = 0; k < N; ++k) y[k] = x[k] + h * k3[k];
      vel(sc, y.data(), k4.data()); for (int k = 0; k < N; ++k) x[k] += h / 6 * (k1[k] + 2 * k2[k] + 2 * k3[k] + k4[k]);
      wnew[m] = vel(sc, x, nullptr); } }
    vmax = std::max(vmax, vm); t += h; ++nst;
    for (int m = 0; m < M; ++m) { const bool pa = wprev[m] >= 0.5, pb = wnew[m] >= 0.5; if (pa != pb) ++cr[m];
      if (pb && !pa && te[m] < 0) te[m] = (t - h + h * (0.5 - wprev[m]) / (wnew[m] - wprev[m])) / tau; }
    while (nextG <= 5 + 1e-9 && nextG * tau <= t + 1e-9) { const double x = (nextG * tau - (t - h)) / h; std::vector<double> wi(M);
      for (int m = 0; m < M; ++m) wi[m] = wprev[m] + x * (wnew[m] - wprev[m]); row(nextG, std::norm(ae), wi); nextG += 0.05; }
    wprev = wnew; }
  for (int m = 0; m < M; ++m) std::printf("# TE %d %.6f\n# CR %d %d\n", m0 + m, te[m], m0 + m, cr[m]);
  std::printf("# steps %ld vmax %.4e ce2_end %.8e\n", nst, vmax, std::norm(ae));
  return 0; }

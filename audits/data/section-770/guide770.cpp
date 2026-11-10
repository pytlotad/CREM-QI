// Audit 770: kinematic condition (first-order guidance).  The particle has no actions of its own: its torus (J_g, L_g) is read from the
// field occupations, J_g = sum w_n (n_r + 1/2), L_g = sum w_n m, w_n = |c_n|^2 / sum |c|^2, and it moves on that torus in angle variables
// (p = grad S_+-(q; J_g, L_g); branch = sign sin u).  Field of the relative motion projected on the exact bound modes of 2D -1/r + beta/r^2
// (hbar = mu = 1, modes as ../section-757/lock757.cpp): dc_n/dt = -(i E_n + Gamma) c_n - i A psi_n*(x_p) exp(i phi_s) (exact linear step).
// Source phase: phi_s = int (v^2/2 - V) dt + CAUST at every crossing of u through a multiple of pi (caustic counter, 769).
// argv: periods f P0 pair(A|B) ; env CAUST (0) A (1) BETA (0.2) GT (0.1) NRMAX (12) MMAX (6) STEPS (2000)
// Initial field: c_a = sqrt((1-f) P0), c_b = sqrt(f P0); pair A: a = (0,1), b = (1,1); pair B: a = (0,1), b = (0,2).
// Output per T_r0 (STEPS steps): k t J_g L_g P w_max nr_max m_max
#include <cmath>
#include <complex>
#include <cstdio>
#include <cstdlib>
#include <vector>
typedef std::complex<double> C;
static double envd(const char* k, double d){ const char* v = std::getenv(k); return v ? std::atof(v) : d; }
struct Mode { int nr, m; double nu, N, E, al, lognorm; };
int main(int argc, char** argv){
  setvbuf(stdout, nullptr, _IOLBF, 0);
  const int periods = std::atoi(argv[1]); const double f = std::atof(argv[2]), P0 = std::atof(argv[3]); const char pair = argv[4][0];
  const double caust = envd("CAUST", 0.0), A = envd("A", 1.0), beta = envd("BETA", 0.2), GT = envd("GT", 0.1);
  const int NRMAX = (int)envd("NRMAX", 12), MMAX = (int)envd("MMAX", 6), STEPS = (int)envd("STEPS", 2000);
  const double N0 = 0.5 + std::sqrt(1.0 + 2 * beta), Tr0 = 2 * M_PI * N0 * N0 * N0, Gam = GT / Tr0, h = Tr0 / STEPS;
  std::vector<Mode> md;
  for (int m = -MMAX; m <= MMAX; ++m) { const double nu = std::sqrt(m * m + 2 * beta);
    for (int nr = 0; nr <= NRMAX; ++nr) { Mode q; q.nr = nr; q.m = m; q.nu = nu; q.N = nr + nu + 0.5; q.E = -0.5 / (q.N * q.N); q.al = 2 * nu;
      q.lognorm = 0.5 * (-2 * std::log(q.N / 2) - (std::lgamma(nr + q.al + 1) - std::lgamma(nr + 1.0)) - std::log(2 * q.N)); md.push_back(q); } }
  const int K = md.size(); auto idx = [&](int nr, int m) { return (m + MMAX) * (NRMAX + 1) + nr; };
  std::vector<C> c(K, 0.0), prop(K), stepS(K), psiV(K);
  for (int i = 0; i < K; ++i) { const C z(Gam, md[i].E); prop[i] = std::exp(-z * h); stepS[i] = (1.0 - prop[i]) / z; }
  const int ia = idx(0, 1), ib = (pair == 'A') ? idx(1, 1) : idx(0, 2);
  c[ia] = std::sqrt((1 - f) * P0); c[ib] += std::sqrt(f * P0);
  auto modes = [&](double r, double ph) {
    for (int i = 0; i < K; ++i) {
      const Mode& q = md[i]; const double rho = 2 * r / q.N;
      double l0 = 1, l1 = 1 + q.al - rho, L = (q.nr == 0) ? 1 : l1;
      for (int k = 1; k < q.nr; ++k) { const double l2 = ((2 * k + 1 + q.al - rho) * l1 - (k + q.al) * l0) / (k + 1); l0 = l1; l1 = l2; L = l2; }
      const double pre = std::exp(q.lognorm + q.nu * std::log(rho) - rho / 2) / std::sqrt(2 * M_PI);
      psiV[i] = pre * L * std::exp(C(0, q.m * ph)); } };
  auto actions = [&](double& J, double& L, double& P, double& wmax, int& imax) {
    P = 0; J = 0; L = 0; wmax = 0; imax = 0;
    for (int i = 0; i < K; ++i) { const double n = std::norm(c[i]); P += n; J += n * (md[i].nr + 0.5); L += n * md[i].m; if (n > wmax) { wmax = n; imax = i; } }
    J /= P; L /= P; wmax /= P; };
  double w = 0, phi = 0, phis = 0, t = 0; long cau = 0;
  std::printf("# 770 beta %.3f GT %.3f Tr0 %.6f caust %.6f A %.3f f %.3f P0 %.1f pair %c modes %d steps %d\n", beta, GT, Tr0, caust, A, f, P0, pair, K, STEPS);
  std::printf("# k t J_g L_g P w_max nr_max m_max\n");
  for (int k = 1; k <= periods; ++k) {
    double J = 0, Lg = 0, P = 0, wmax = 0; int imax = 0;
    for (int s = 0; s < STEPS; ++s) {
      actions(J, Lg, P, wmax, imax);
      const double Lam = std::sqrt(Lg * Lg + 2 * beta), N = J + Lam, e = std::sqrt(std::max(0.0, 1 - Lam * Lam / (N * N))), om = 1 / (N * N * N);
      // Kepler: u from w (continuous)
      const double kk = std::floor(w / (2 * M_PI)), wp = w - 2 * M_PI * kk; double u = wp + e * std::sin(wp);
      for (int it = 0; it < 50; ++it) { const double du = (u - e * std::sin(u) - wp) / (1 - e * std::cos(u)); u -= du; if (std::fabs(du) < 1e-15) break; }
      u += 2 * M_PI * kk;
      const long cn = (long)std::floor(u / M_PI); if (cn != cau) { phis += caust * (double)(cn - cau); cau = cn; }
      const double r = N * N * (1 - e * std::cos(u)), dudt = om / (1 - e * std::cos(u)), vr = N * N * e * std::sin(u) * dudt;
      const double lag = 0.5 * (vr * vr + Lg * Lg / (r * r)) + 1 / r - beta / (r * r);
      modes(r, phi);
      const C src = C(0, -A) * std::exp(C(0, phis));
      for (int i = 0; i < K; ++i) c[i] = prop[i] * c[i] + stepS[i] * src * std::conj(psiV[i]);
      w += h * om; phi += h * Lg / (r * r); phis += h * lag; t += h; }
    actions(J, Lg, P, wmax, imax);
    std::printf("%d %.4f %.9f %.9f %.6e %.6f %d %d\n", k, t, J, Lg, P, wmax, md[imax].nr, md[imax].m);
  }
  return 0; }

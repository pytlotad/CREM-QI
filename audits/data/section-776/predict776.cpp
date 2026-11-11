// Audit 776: drive spectrum on the FROZEN start torus (evaluation tool; model ../section-772/guide772.cpp unchanged).
// Start torus J0 = sum w (n_r + 1/2), L0 = sum w m from INIT weights; orbit, angle variables, source phase with caustic counter, modes, step -- as in guide772.
// S_n = sum_k h psi_n*(x_k) exp(i phi_s,k) exp(i E_n t_k) over `periods` periods (Gamma = 0); D_n = |S_n|^2 / T.
// argv: periods ; env INIT (required) CAUST (0) BETA (0.2) NRMAX (12) MMAX (6) STEPS (2000).  Output: header + one line per mode "nr m D", sorted by D descending.
#include <cmath>
#include <complex>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <algorithm>
typedef std::complex<double> C;
static double envd(const char* k, double d){ const char* v = std::getenv(k); return v ? std::atof(v) : d; }
struct Mode { int nr, m; double nu, N, E, al, lognorm; };
int main(int argc, char** argv){
  const int periods = std::atoi(argv[1]);
  const double caust = envd("CAUST", 0.0), beta = envd("BETA", 0.2);
  const int NRMAX = (int)envd("NRMAX", 12), MMAX = (int)envd("MMAX", 6), STEPS = (int)envd("STEPS", 2000);
  const double N0 = 0.5 + std::sqrt(1.0 + 2 * beta), Tr0 = 2 * M_PI * N0 * N0 * N0, h = Tr0 / STEPS;
  std::vector<Mode> md;
  for (int m = -MMAX; m <= MMAX; ++m) { const double nu = std::sqrt(m * m + 2 * beta);
    for (int nr = 0; nr <= NRMAX; ++nr) { Mode q; q.nr = nr; q.m = m; q.nu = nu; q.N = nr + nu + 0.5; q.E = -0.5 / (q.N * q.N); q.al = 2 * nu;
      q.lognorm = 0.5 * (-2 * std::log(q.N / 2) - (std::lgamma(nr + q.al + 1) - std::lgamma(nr + 1.0)) - std::log(2 * q.N)); md.push_back(q); } }
  const int K = md.size();
  std::vector<C> S(K, 0.0), psiV(K);
  double J = 0, Lg = 0, W = 0;
  { const char* p = std::getenv("INIT"); int nr, m, nch; double w, ph;
    while (p && std::sscanf(p, "%d:%d:%lf:%lf%n", &nr, &m, &w, &ph, &nch) == 4) { J += w * (nr + 0.5); Lg += w * m; W += w; p += nch; if (*p == ';') ++p; else break; } }
  J /= W; Lg /= W;
  auto modes = [&](double r, double ph) {
    for (int i = 0; i < K; ++i) {
      const Mode& q = md[i]; const double rho = 2 * r / q.N;
      double l0 = 1, l1 = 1 + q.al - rho, L = (q.nr == 0) ? 1 : l1;
      for (int k = 1; k < q.nr; ++k) { const double l2 = ((2 * k + 1 + q.al - rho) * l1 - (k + q.al) * l0) / (k + 1); l0 = l1; l1 = l2; L = l2; }
      const double pre = std::exp(q.lognorm + q.nu * std::log(rho) - rho / 2) / std::sqrt(2 * M_PI);
      psiV[i] = pre * L * std::exp(C(0, q.m * ph)); } };
  const double Lam = std::sqrt(Lg * Lg + 2 * beta), N = J + Lam, e = std::sqrt(std::max(0.0, 1 - Lam * Lam / (N * N))), om = 1 / (N * N * N);
  double w = 0, phi = 0, phis = 0, t = 0; long cau = 0;
  for (int k = 1; k <= periods; ++k)
    for (int s = 0; s < STEPS; ++s) {
      const double kk = std::floor(w / (2 * M_PI)), wp = w - 2 * M_PI * kk; double u = wp + e * std::sin(wp);
      for (int it = 0; it < 50; ++it) { const double du = (u - e * std::sin(u) - wp) / (1 - e * std::cos(u)); u -= du; if (std::fabs(du) < 1e-15) break; }
      u += 2 * M_PI * kk;
      const long cn = (long)std::floor(u / M_PI); if (cn != cau) { phis += caust * (double)(cn - cau); cau = cn; }
      const double r = N * N * (1 - e * std::cos(u)), dudt = om / (1 - e * std::cos(u)), vr = N * N * e * std::sin(u) * dudt;
      const double lag = 0.5 * (vr * vr + Lg * Lg / (r * r)) + 1 / r - beta / (r * r);
      modes(r, phi);
      const C src = std::exp(C(0, phis));
      for (int i = 0; i < K; ++i) S[i] += h * std::conj(psiV[i]) * src * std::exp(C(0, md[i].E * t));
      w += h * om; phi += h * Lg / (r * r); phis += h * lag; t += h; }
  std::vector<int> o(K); for (int i = 0; i < K; ++i) o[i] = i;
  std::sort(o.begin(), o.end(), [&](int a, int b) { return std::norm(S[a]) > std::norm(S[b]); });
  std::printf("# 776 beta %.3f caust %.6f J0 %.6f L0 %.6f periods %d T %.6f\n", beta, caust, J, Lg, periods, t);
  for (int i : o) std::printf("%d %d %.6e\n", md[i].nr, md[i].m, std::norm(S[i]) / t);
  return 0; }

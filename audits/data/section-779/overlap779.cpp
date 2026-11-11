// Audit 779: one-period drive integral I_n on a FROZEN torus (evaluation tool; model ../section-772/guide772.cpp unchanged).
// Modes, orbit (Kepler angle variables), caustic counter, Lagrangian source phase -- as in ../section-776/predict776.cpp; here the step h = T / STEPS is commensurate with
// the radial period T = 2 pi N^3 and exactly ONE period is integrated (left Riemann sum, start at periapsis u = 0 with the periapsis jump already counted, as in predict776 from period 2 on).
// stdin: lines "key beta caust J L".  stdout per line: "## key beta caust J L N T int_lag int_phidot jumps_inside" then per mode "nr m ReI ImI |I|^2 Phi_num",
// Phi_num = int_lag + 2 caust + E_n T - m int_phidot (wrapped to (-pi, pi]).   env NRMAX (12) MMAX (6) STEPS (2000).
#include <cmath>
#include <complex>
#include <cstdio>
#include <cstdlib>
#include <vector>
typedef std::complex<double> C;
static double envd(const char* k, double d){ const char* v = std::getenv(k); return v ? std::atof(v) : d; }
struct Mode { int nr, m; double nu, N, E, al, lognorm; };
static double wrap(double x){ x = std::fmod(x + M_PI, 2 * M_PI); if (x <= 0) x += 2 * M_PI; return x - M_PI; }
int main(){
  const int NRMAX = (int)envd("NRMAX", 12), MMAX = (int)envd("MMAX", 6), STEPS = (int)envd("STEPS", 2000);
  char key[256]; double beta, caust, J, Lg;
  while (std::scanf("%255s %lf %lf %lf %lf", key, &beta, &caust, &J, &Lg) == 5) {
    std::vector<Mode> md;
    for (int m = -MMAX; m <= MMAX; ++m) { const double nu = std::sqrt(m * m + 2 * beta);
      for (int nr = 0; nr <= NRMAX; ++nr) { Mode q; q.nr = nr; q.m = m; q.nu = nu; q.N = nr + nu + 0.5; q.E = -0.5 / (q.N * q.N); q.al = 2 * nu;
        q.lognorm = 0.5 * (-2 * std::log(q.N / 2) - (std::lgamma(nr + q.al + 1) - std::lgamma(nr + 1.0)) - std::log(2 * q.N)); md.push_back(q); } }
    const int K = md.size();
    std::vector<C> I(K, 0.0), psiV(K);
    auto modes = [&](double r, double ph) {
      for (int i = 0; i < K; ++i) {
        const Mode& q = md[i]; const double rho = 2 * r / q.N;
        double l0 = 1, l1 = 1 + q.al - rho, L = (q.nr == 0) ? 1 : l1;
        for (int k = 1; k < q.nr; ++k) { const double l2 = ((2 * k + 1 + q.al - rho) * l1 - (k + q.al) * l0) / (k + 1); l0 = l1; l1 = l2; L = l2; }
        const double pre = std::exp(q.lognorm + q.nu * std::log(rho) - rho / 2) / std::sqrt(2 * M_PI);
        psiV[i] = pre * L * std::exp(C(0, q.m * ph)); } };
    const double Lam = std::sqrt(Lg * Lg + 2 * beta), N = J + Lam, e = std::sqrt(std::max(0.0, 1 - Lam * Lam / (N * N))), om = 1 / (N * N * N);
    const double T = 2 * M_PI * N * N * N, h = T / STEPS;
    double w = 0, phi = 0, ilag = 0, iphi = 0, t = 0, jumps = 0; long cau = 0;
    for (int s = 0; s < STEPS; ++s) {
      const double wp = w; double u = wp + e * std::sin(wp);
      for (int it = 0; it < 50; ++it) { const double du = (u - e * std::sin(u) - wp) / (1 - e * std::cos(u)); u -= du; if (std::fabs(du) < 1e-15) break; }
      const long cn = (long)std::floor(u / M_PI); if (cn != cau) { jumps += caust * (double)(cn - cau); cau = cn; }
      const double r = N * N * (1 - e * std::cos(u)), dudt = om / (1 - e * std::cos(u)), vr = N * N * e * std::sin(u) * dudt;
      const double lag = 0.5 * (vr * vr + Lg * Lg / (r * r)) + 1 / r - beta / (r * r);
      modes(r, phi);
      const C src = std::exp(C(0, ilag + jumps));
      for (int i = 0; i < K; ++i) I[i] += h * std::conj(psiV[i]) * src * std::exp(C(0, md[i].E * t));
      w += h * om; phi += h * Lg / (r * r); ilag += h * lag; iphi += h * Lg / (r * r); t += h; }
    std::printf("## %s %.4f %.6f %.6f %.6f %.8f %.6f %.10f %.10f %.6f\n", key, beta, caust, J, Lg, N, T, ilag, iphi, jumps);
    for (int i = 0; i < K; ++i)
      std::printf("%d %d %.8e %.8e %.8e %.8f\n", md[i].nr, md[i].m, I[i].real(), I[i].imag(), std::norm(I[i]), wrap(ilag + 2 * caust + md[i].E * T - md[i].m * iphi));
  }
  return 0; }

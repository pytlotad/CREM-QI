// Audit 769 (copy of ../section-757/lock757.cpp + env CAUST: phi_s += CAUST at every sign change of the radial velocity).
// Audit 757 (copy of ../section-744/lock744.cpp + work columns): W_chi = sum Re[e^{-i(phi_s+chi)} grad psi_mem . v] h per radial period,
// chi = 0 and pi/2 -- the work of the engine-type coupling V = -Re[e^{-i(phi_s+chi)} psi_mem(x_p)] on a test particle (eps = 0: no back-reaction).
// Audit 744: source takes over the phase of the (retarded) memory field.  Field of the relative motion projected on the exact bound
// modes of 2D -1/r + beta/r^2 (hbar = mu = 1): psi = sum c_n psi_n, dc_n/dt = -(i E_n + Gamma) c_n - i A psi_n*(x_p) exp(i phi_s).
// Memory (retarded) field: emitted up to t - tau, i.e. c_n^mem(t) = exp(-(i E_n + Gamma) tau) c_n(t - tau) (exact, ring buffer).
// Source phase: dphi_s/dt = (v^2/2 - V) + kappa sin(arg psi_mem(x_p) + sgn pi/2 - phi_s)   (sgn +1 harmony, -1 anti-harmony).
// Particle: Verlet (SUB substeps per field step), force eps Re[exp(-i(phi_s + chi)) grad psi_mem(x_p)].
// argv: periods kappaT tauFrac sgn eps chi ; env JR0 (0.5) L0 (1.0) BETA (0.2) GT (0.1) GABS (absolute Gamma) A (1) NRMAX (12) MMAX (6) STEPS (2000) SUB (4)
// Output per radial period (periapsis to periapsis): k t J_r L slip <P> <|psi_mem(x_p)|>   slip = Dphi_s - Dphi_gauge (wrapped).
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
  const int periods = std::atoi(argv[1]); const double kT = std::atof(argv[2]), tauF = std::atof(argv[3]), sgn = std::atof(argv[4]);
  const double eps = std::atof(argv[5]), chi = std::atof(argv[6]);
  const double Jr0 = envd("JR0", 0.5), L0 = envd("L0", 1.0), beta = envd("BETA", 0.2), GT = envd("GT", 0.1), A = envd("A", 1.0); const double caust = envd("CAUST", 0.0);
  const int NRMAX = (int)envd("NRMAX", 12), MMAX = (int)envd("MMAX", 6), STEPS = (int)envd("STEPS", 2000), SUB = (int)envd("SUB", 4);
  const double Lam0 = std::sqrt(L0 * L0 + 2 * beta), N0 = Jr0 + Lam0, Tr0 = 2 * M_PI * N0 * N0 * N0, Gam = std::getenv("GABS") ? envd("GABS", 0) : GT / Tr0, h = Tr0 / STEPS;
  const double kappa = kT / Tr0; const int ntau = (int)std::lround(tauF * STEPS); const double tau = ntau * h;
  std::vector<Mode> md;
  for (int m = -MMAX; m <= MMAX; ++m) { const double nu = std::sqrt(m * m + 2 * beta);
    for (int nr = 0; nr <= NRMAX; ++nr) { Mode q; q.nr = nr; q.m = m; q.nu = nu; q.N = nr + nu + 0.5; q.E = -0.5 / (q.N * q.N); q.al = 2 * nu;
      q.lognorm = 0.5 * (-2 * std::log(q.N / 2) - (std::lgamma(nr + q.al + 1) - std::lgamma(nr + 1.0)) - std::log(2 * q.N)); md.push_back(q); } }
  const int K = md.size();
  std::vector<C> c(K, 0.0), prop(K), propTau(K), stepS(K);
  for (int i = 0; i < K; ++i) { const C z(Gam, md[i].E); prop[i] = std::exp(-z * h); stepS[i] = (1.0 - prop[i]) / z; propTau[i] = std::exp(-z * tau); }
  std::vector<std::vector<C>> ring(ntau + 1, std::vector<C>(K, 0.0));   // ring[j % (ntau+1)] = c(t_j)
  // mode values and gradients at (x, y)
  std::vector<C> psiV(K), gx(K), gy(K);
  auto modes = [&](double x, double y, bool grad) {
    const double r = std::sqrt(x * x + y * y), ph = std::atan2(y, x), cs = std::cos(ph), sn = std::sin(ph);
    for (int i = 0; i < K; ++i) {
      const Mode& q = md[i]; const double rho = 2 * r / q.N;
      // generalized Laguerre L_n^al and L_{n-1}^{al+1}
      double l0 = 1, l1 = 1 + q.al - rho, L = (q.nr == 0) ? 1 : l1;
      for (int k = 1; k < q.nr; ++k) { const double l2 = ((2 * k + 1 + q.al - rho) * l1 - (k + q.al) * l0) / (k + 1); l0 = l1; l1 = l2; L = l2; }
      double dL = 0; if (q.nr >= 1) { const double a1 = q.al + 1; double m0 = 1, m1 = 1 + a1 - rho; dL = (q.nr - 1 == 0) ? 1 : m1;
        for (int k = 1; k < q.nr - 1; ++k) { const double m2 = ((2 * k + 1 + a1 - rho) * m1 - (k + a1) * m0) / (k + 1); m0 = m1; m1 = m2; dL = m2; } dL = -dL; }
      const double pre = std::exp(q.lognorm + q.nu * std::log(rho) - rho / 2) / std::sqrt(2 * M_PI);
      const double R = pre * L; const C e = std::exp(C(0, q.m * ph)); psiV[i] = R * e;
      if (grad) { const double dR = pre * (L * (q.nu / r - 1 / q.N) + dL * 2 / q.N); const C gr = dR * e, gp = C(0, q.m / r) * R * e;
        gx[i] = gr * cs - gp * sn; gy[i] = gr * sn + gp * cs; } } };
  auto Vp = [&](double x, double y) { const double r = std::sqrt(x * x + y * y); return -1 / r + beta / (r * r); };
  auto ac = [&](double x, double y, double& ax, double& ay) { const double r2 = x * x + y * y, r = std::sqrt(r2), g = 1 / (r2 * r) - 2 * beta / (r2 * r2); ax = -x * g; ay = -y * g; };
  const double a0 = N0 * N0, e0 = std::sqrt(std::max(0.0, 1 - Lam0 * Lam0 / (N0 * N0))), rp = a0 * (1 - e0);
  double x = rp, y = 0, vx = 0, vy = L0 / rp, phis = 0, phig = 0, t = 0, vrPrev = 0, phis0 = 0, phig0 = 0; long nstep = 0;
  std::printf("# 744 beta %.3f GT %.3f Jr0 %.4f L0 %.4f Tr0 %.6f kappaT %.3f tau/Tr0 %.4f sgn %+.0f eps %.3e chi %.6f A %.2f modes %d steps %d\n",
              beta, GT, Jr0, L0, Tr0, kT, tauF, sgn, eps, chi, A, K, STEPS);
  std::printf("# k t J_r L slip <P> <|psi_mem|> W0 W90\n");
  int k = 0; double Ps = 0, Ms = 0, W0 = 0, W90 = 0; long ns = 0;
  while (k < periods) {
    modes(x, y, true);
    // memory field at the particle
    const std::vector<C>& cpast = ring[(nstep + 1) % (ntau + 1)];   C pm = 0, gmx = 0, gmy = 0;   // cpast = c(t_{j - ntau}) C pm = 0, gmx = 0, gmy = 0;
    if (nstep >= ntau) for (int i = 0; i < K; ++i) { const C cm = propTau[i] * cpast[i]; pm += cm * psiV[i]; gmx += cm * gx[i]; gmy += cm * gy[i]; }
    // field update (source from current point, exact linear step)
    const C src = C(0, -A) * std::exp(C(0, phis)); double P = 0;
    for (int i = 0; i < K; ++i) { c[i] = prop[i] * c[i] + stepS[i] * src * std::conj(psiV[i]); P += std::norm(c[i]); }
    // phase
    const double am = std::abs(pm); const double lag = 0.5 * (vx * vx + vy * vy) - Vp(x, y);
    double rate = lag; if (kappa != 0 && am > 1e-8) rate += kappa * std::sin(std::arg(pm) + sgn * M_PI / 2 - phis);
    { const C gv = gmx * vx + gmy * vy, e0 = std::exp(C(0, -phis)); W0 += std::real(e0 * gv) * h; W90 += std::real(C(0, -1) * e0 * gv) * h; }
    double fx = 0, fy = 0; if (eps != 0) { const C ph1 = std::exp(C(0, -(phis + chi))); fx = eps * std::real(ph1 * gmx); fy = eps * std::real(ph1 * gmy); }
    // particle
    const double hs = h / SUB; double ax, ay; ac(x, y, ax, ay); double g0 = 0;
    for (int s = 0; s < SUB; ++s) { const double l1 = 0.5 * (vx * vx + vy * vy) - Vp(x, y);
      vx += 0.5 * hs * (ax + fx); vy += 0.5 * hs * (ay + fy); x += hs * vx; y += hs * vy; ac(x, y, ax, ay); vx += 0.5 * hs * (ax + fx); vy += 0.5 * hs * (ay + fy);
      const double l2 = 0.5 * (vx * vx + vy * vy) - Vp(x, y); g0 += 0.5 * (l1 + l2) * hs; }
    phis += h * (rate - lag) + g0; phig += g0;   // gauge part integrated with the particle, relaxation part explicit
    // ring buffer: store c(t+h)
    ring[(nstep + 1) % (ntau + 1)] = c;
    t += h; ++nstep; Ps += P; Ms += am; ++ns;
    const double vr = x * vx + y * vy;
    if (caust != 0 && ((vrPrev < 0 && vr >= 0) || (vrPrev > 0 && vr <= 0))) phis += caust;   // audit 769: caustic counter
    if (vrPrev < 0 && vr >= 0) { ++k;
      const double E = 0.5 * (vx * vx + vy * vy) + Vp(x, y); if (!(E < 0)) { std::printf("# unbound\n"); break; }
      const double L = x * vy - y * vx, Nn = 1 / std::sqrt(-2 * E), Jr = Nn - std::sqrt(L * L + 2 * beta);
      double sl = std::remainder((phis - phis0) - (phig - phig0), 2 * M_PI);
      std::printf("%d %.4f %.9f %.9f %+.6f %.6e %.6e %+.9e %+.9e\n", k, t, Jr, L, sl, Ps / ns, Ms / ns, W0, W90);
      phis0 = phis; phig0 = phig; Ps = Ms = W0 = W90 = 0; ns = 0; }
    vrPrev = vr; }
  return 0; }

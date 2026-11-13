// Audit 783: copy of ../section-782/guide782.cpp + env GW (radiative widths: uniform damping replaced by Gamma_n = sum_f Gamma_F(n -> f), amplitude damping Gamma_n / 2;
// ground 0) and GREF (reference unit Gamma_ref = GREF / Tr0 for kappa, default GT).  Without GW (or GW=0) and GREF the output is byte-identical to guide782.
// Audit 782: copy of ../section-781/guide781.cpp + env DIP (Abraham-Lorentz reaction on the FIELD dipole, semiclassical Crisp-Jaynes as in ../section-751; without DIP or DIP=0 skipped).
// DIP: d = <psi|r|psi>/P, d'''_+ = sum c_n'* c_n (R+)_n'n (i(E_n'-E_n))^3 / P (free evolution), F = (2/3) kappa d''', H_int = -(F_- R+ + F_+ R-)/2 (Hermitian; norm conserved),
// c += h G c + h^2/2 G G c, G = (i/2)(F_- R+ + F_+ R-).  (R+)_n'n = I(n',n) delta(m', m+1), I = int R_n' R_n r^2 dr (numerical at start).  kappa = DIP Gamma / ((4/3) w^3 I^2/2) for (0,1)->(0,0).
// With DIP: column 9 = energy radiated in the period (drop of <E> by H_int), column 10 = weight |m| >= 5; env WDUMP=1: every 10 periods "#W k" + weights of 13 lowest modes.
// Audit 781: copy of ../section-772/guide772.cpp + env RAD (classical radiation reaction taken from the field; without RAD or RAD=0 the block is skipped -> output byte-identical).
// RAD: Larmor power P = k a^2 (a = 1/r^2 - 2 beta/r^3) and torque dL/dt = -k L_g (1/r^3 - 2 beta/r^4) of the actual orbit are removed from the field by the smallest tilt of the weights,
// c_n *= exp(-(l1 (E_n - <E>) + l2 (m - <m>))/2), norm restored; l1, l2 from linear response; |exponent| <= 1e-3 per step (else scaled down, shortfall counted as deficit).
// k = RAD * 0.178 * (E_(0,1) - E_(0,0)) * Gamma / <a^2> on EBK torus (0,1)  (0.178 = Larmor n = 2 circular / hbar omega A(2p), hydrogen).  With RAD: extra column 'deficit' per period.
// Audit 772 (copy of ../section-770/guide770.cpp + env INIT).  Audit 770: kinematic condition (first-order guidance).  The particle has no actions of its own: its torus (J_g, L_g) is read from the
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
#include <algorithm>
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
  const char* radS = std::getenv("RAD"); const double RAD = radS ? std::atof(radS) : 0.0; const bool rad = RAD != 0.0; double kR = 0;
  if (rad) { const double Lm = std::sqrt(1 + 2 * beta), Nn = 0.5 + Lm, ee = std::sqrt(1 - Lm * Lm / (Nn * Nn)); double s2 = 0; const int M = 20000;
    for (int i = 0; i < M; ++i) { const double wp = 2 * M_PI * (i + 0.5) / M; double u = wp; for (int it = 0; it < 60; ++it) u -= (u - ee * std::sin(u) - wp) / (1 - ee * std::cos(u));
      const double r = Nn * Nn * (1 - ee * std::cos(u)), a = 1 / (r * r) - 2 * beta / (r * r * r); s2 += a * a / M; }
    const double N00 = 0.5 + std::sqrt(2 * beta); kR = RAD * 0.178 * (0.5 / (N00 * N00) - 0.5 / (Nn * Nn)) * Gam / s2; }
  double defT = 0, defA = 0;
  const char* dipS = std::getenv("DIP"); const double DIP = dipS ? std::atof(dipS) : 0.0; const bool dip = DIP != 0.0; const bool wdump = envd("WDUMP", 0) != 0;
  const char* gwS = std::getenv("GW"); const double GWs = gwS ? std::atof(gwS) : 0.0; const bool gw = GWs != 0.0, ext = dip || gw; const double GamR = envd("GREF", GT) / Tr0; double kapW = 0;
  const int NR1 = NRMAX + 1; std::vector<double> Ip;   // Ip[((mp + MMAX - 1) * NR1 + nrp) * NR1 + nr] = I((nrp, mp), (nr, mp - 1)), mp = -MMAX+1..MMAX
  double kap = 0, normdev = 0;
  if (dip || gw) {
    const int Km = (2 * MMAX + 1) * NR1; Ip.assign((size_t)(2 * MMAX) * NR1 * NR1, 0.0); std::vector<double> nrm(Km, 0.0), Rv(Km);
    const int G = 200000; const double smax = std::sqrt(3000.0), ds = smax / G;
    for (int g = 0; g < G; ++g) { const double sg = (g + 0.5) * ds, r = sg * sg, wgt = 2 * sg * ds * r;   // r dr
      for (int i = 0; i < Km; ++i) { const Mode& q = md[i]; const double rho = 2 * r / q.N;
        double l0 = 1, l1 = 1 + q.al - rho, L = (q.nr == 0) ? 1 : l1;
        for (int k = 1; k < q.nr; ++k) { const double l2 = ((2 * k + 1 + q.al - rho) * l1 - (k + q.al) * l0) / (k + 1); l0 = l1; l1 = l2; L = l2; }
        Rv[i] = std::exp(q.lognorm + q.nu * std::log(rho) - rho / 2) * L; nrm[i] += wgt * Rv[i] * Rv[i]; }
      for (int mp = -MMAX + 1; mp <= MMAX; ++mp) for (int a = 0; a < NR1; ++a) { const double ra = Rv[(mp + MMAX) * NR1 + a] * wgt * r;
        double* row = &Ip[((size_t)(mp + MMAX - 1) * NR1 + a) * NR1]; const double* rb = &Rv[(mp - 1 + MMAX) * NR1]; for (int b = 0; b < NR1; ++b) row[b] += ra * rb[b]; } }
    for (int i = 0; i < Km; ++i) normdev = std::max(normdev, std::fabs(nrm[i] - 1));
    const double I01 = Ip[((size_t)(1 + MMAX - 1) * NR1 + 0) * NR1 + 0], w01 = md[(1 + MMAX) * NR1].E - md[MMAX * NR1].E;
    const double kden = (4.0 / 3) * w01 * w01 * w01 * I01 * I01 / 2; kap = DIP * GamR / kden; kapW = GWs * GamR / kden; }
  auto Iof = [&](int nrp, int mp, int nr, int m) -> double {   // <(nrp,mp)| r e^{+-i phi} |(nr,m)> radial integral, |mp - m| = 1
    if (mp == m + 1) return Ip[((size_t)(mp + MMAX - 1) * NR1 + nrp) * NR1 + nr];
    if (m == mp + 1) return Ip[((size_t)(m + MMAX - 1) * NR1 + nr) * NR1 + nrp];
    return 0.0; };
  double Erad = 0;
  const int K = md.size(); auto idx = [&](int nr, int m) { return (m + MMAX) * (NRMAX + 1) + nr; };
  struct Pr { int ip, in; double I, om3; }; std::vector<Pr> prs;   // pairs m(ip) = m(in) + 1
  if (dip || gw) for (int mp = -MMAX + 1; mp <= MMAX; ++mp) for (int a = 0; a < NR1; ++a) for (int b = 0; b < NR1; ++b) {
    const int ip = idx(a, mp), in = idx(b, mp - 1); const double om = md[ip].E - md[in].E; prs.push_back({ip, in, Ip[((size_t)(mp + MMAX - 1) * NR1 + a) * NR1 + b], om * om * om}); }
  std::vector<C> g1(K), g2(K);
  std::vector<C> c(K, 0.0), prop(K), stepS(K), psiV(K);
  std::vector<double> Gn(K, 0.0);   // audit 783: radiative widths (weight decay rates)
  if (gw) for (const Pr& p : prs) { const double om = std::fabs(md[p.ip].E - md[p.in].E), g = (4.0 / 3) * kapW * om * om * om * p.I * p.I / 2; Gn[md[p.ip].E > md[p.in].E ? p.ip : p.in] += g; }
  for (int i = 0; i < K; ++i) { const C z(gw ? Gn[i] / 2 : Gam, md[i].E); prop[i] = std::exp(-z * h); stepS[i] = (1.0 - prop[i]) / z; }
  const int ia = idx(0, 1), ib = (pair == 'A') ? idx(1, 1) : idx(0, 2);
  c[ia] = std::sqrt((1 - f) * P0); c[ib] += std::sqrt(f * P0);
  if (const char* ini = std::getenv("INIT")) {   // audit 772: "nr:m:w:phase;..." overrides the pair start
    std::fill(c.begin(), c.end(), C(0.0)); const char* p = ini; int nr, m, nch; double w, ph;
    while (std::sscanf(p, "%d:%d:%lf:%lf%n", &nr, &m, &w, &ph, &nch) == 4) { c[idx(nr, m)] = std::sqrt(w * P0) * std::exp(C(0, ph)); p += nch; if (*p == ';') ++p; else break; } }
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
  if (ext) {
    std::printf("# 783 GW %.3f kappaW %.6e GREF %.6f\n", GWs, kapW, GamR * Tr0);
    std::printf("# 782 DIP %.3f kappa %.6e radial-norm max dev %.2e ; extra columns: E_rad in period, w(|m|>=5)\n", DIP, kap, normdev);
    const int LN[13][2] = {{0,0},{0,1},{0,-1},{1,0},{0,2},{0,-2},{1,1},{1,-1},{2,0},{1,2},{1,-2},{2,1},{2,-1}};
    std::printf("# 782 Fermi rates / Gamma (i -> f, E_f < E_i, |dm| = 1) among 13 lowest modes:\n");
    for (int a = 0; a < 13; ++a) for (int b = 0; b < 13; ++b) { const int ia2 = idx(LN[a][0], LN[a][1]), ib2 = idx(LN[b][0], LN[b][1]);
      if (std::abs(LN[a][1] - LN[b][1]) != 1 || md[ib2].E >= md[ia2].E) continue; const double om2 = md[ia2].E - md[ib2].E, I2 = Iof(LN[b][0], LN[b][1], LN[a][0], LN[a][1]);
      std::printf("#F %d %d -> %d %d  omega %.6f  I %.6f  GammaF/Gamma %.6e\n", LN[a][0], LN[a][1], LN[b][0], LN[b][1], om2, I2, (4.0 / 3) * (dip ? kap : kapW) * om2 * om2 * om2 * I2 * I2 / 2 / GamR); }
    if (gw) for (int a = 0; a < 13; ++a) std::printf("#G %d %d  Gamma_n/Gamma_ref %.6e\n", LN[a][0], LN[a][1], Gn[idx(LN[a][0], LN[a][1])] / GamR); }
  if (rad) std::printf("# 781 RAD %.3f k %.6e ; extra column: deficit = 1 - (energy taken from field) / (Larmor energy) in the period\n", RAD, kR);
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
      if (rad) {
        double Pt = 0, Em = 0, Mm = 0; for (int i = 0; i < K; ++i) { const double n = std::norm(c[i]); Pt += n; Em += n * md[i].E; Mm += n * md[i].m; } Em /= Pt; Mm /= Pt;
        double VE = 0, CEm = 0, VM = 0; for (int i = 0; i < K; ++i) { const double n = std::norm(c[i]) / Pt, dE = md[i].E - Em, dm = md[i].m - Mm; VE += n * dE * dE; CEm += n * dE * dm; VM += n * dm * dm; }
        const double ar = 1 / (r * r) - 2 * beta / (r * r * r), tE = kR * ar * ar * h, tM = -kR * Lg * (1 / (r * r * r) - 2 * beta / (r * r * r * r)) * h;
        VE += 1e-24; VM += 1e-24; const double det = VE * VM - CEm * CEm;
        double l1 = (tE * VM - CEm * (-tM)) / det, l2 = (VE * (-tM) - CEm * tE) / det, xm = 0;
        for (int i = 0; i < K; ++i) xm = std::max(xm, std::fabs(l1 * (md[i].E - Em) + l2 * (md[i].m - Mm)));
        if (xm > 1e-3) { l1 *= 1e-3 / xm; l2 *= 1e-3 / xm; }
        double P2 = 0, E2 = 0; for (int i = 0; i < K; ++i) { c[i] *= std::exp(-0.5 * (l1 * (md[i].E - Em) + l2 * (md[i].m - Mm))); P2 += std::norm(c[i]); }
        const double sc = std::sqrt(Pt / P2); for (int i = 0; i < K; ++i) { c[i] *= sc; E2 += std::norm(c[i]) * md[i].E; }
        defT += tE; defA += Em - E2 / Pt; }
      if (dip) {
        double Pt = 0, E0 = 0; for (int i = 0; i < K; ++i) { const double n = std::norm(c[i]); Pt += n; E0 += n * md[i].E; }
        C d3(0.0);   // d'''_+ = sum conj(c_n') c_n I(n',n) (-i w^3), w = E_n' - E_n, m' = m + 1
        for (const Pr& p : prs) d3 += std::conj(c[p.ip]) * c[p.in] * (p.I * p.om3);
        d3 *= C(0, -1.0 / Pt); const C Fp = (2.0 / 3) * kap * d3, Fm = std::conj(Fp);
        auto applyG = [&](const std::vector<C>& x, std::vector<C>& y) {   // y = G x = (i/2)(F_- R+ x + F_+ R- x)
          std::fill(y.begin(), y.end(), C(0.0));
          for (const Pr& p : prs) { y[p.ip] += (Fm * p.I) * x[p.in]; y[p.in] += (Fp * p.I) * x[p.ip]; }
          for (auto& v : y) v *= C(0, 0.5); };
        applyG(c, g1); applyG(g1, g2);
        double E1 = 0; for (int i = 0; i < K; ++i) { c[i] += h * g1[i] + 0.5 * h * h * g2[i]; E1 += std::norm(c[i]) * md[i].E; }
        Erad += (E0 - E1) / Pt; }
      w += h * om; phi += h * Lg / (r * r); phis += h * lag; t += h; }
    actions(J, Lg, P, wmax, imax);
    if (ext) { double P5 = 0; for (int i = 0; i < K; ++i) if (std::abs(md[i].m) >= 5) P5 += std::norm(c[i]);
      std::printf("%d %.4f %.9f %.9f %.6e %.6f %d %d %.6e %.6e\n", k, t, J, Lg, P, wmax, md[imax].nr, md[imax].m, Erad, P5 / P); Erad = 0;
      if (wdump && k % 10 == 0) { const int LN[13][2] = {{0,0},{0,1},{0,-1},{1,0},{0,2},{0,-2},{1,1},{1,-1},{2,0},{1,2},{1,-2},{2,1},{2,-1}};
        std::printf("#W %d", k); for (int a = 0; a < 13; ++a) std::printf(" %.5e", std::norm(c[idx(LN[a][0], LN[a][1])]) / P); std::printf("\n"); } }
    else if (!rad) std::printf("%d %.4f %.9f %.9f %.6e %.6f %d %d\n", k, t, J, Lg, P, wmax, md[imax].nr, md[imax].m);
    else { std::printf("%d %.4f %.9f %.9f %.6e %.6f %d %d %.4e\n", k, t, J, Lg, P, wmax, md[imax].nr, md[imax].m, defT > 0 ? 1 - defA / defT : 0.0); defT = defA = 0; }
  }
  return 0; }

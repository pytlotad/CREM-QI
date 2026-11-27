// Audit 871: local-field prototype (no mode basis, no potential in the field).  Two charges (1: +1, 2: -1), centre-of-mass frame, hbar = e^2 = 1,
// masses in m_e.  Local field psi on a grid in the configuration space of r = x2 - x1 (no propagation, own frequency 0, damping GAM); the source
// at r(t) writes |v_rel| h/(sqrt(2 pi) W) G_W(x - r) e^{i phi_s}; a write reaches the field after TAUEX (+ |x - r_w|/C per cell with WDELAY).
// Source phase: CLOCK 3 (default) d phi/dt = sum p_i.v_i; CLOCK 1 = sum T_i - V (only MODE 0); minus sum m_i Omega_z,i (Thomas-BMT) when MS1/MS2 != 0.
// Guidance: F = KAPPA grad Re[psi(r) e^{-i phi_s}] on the relative motion (particle 2 +F, particle 1 -F).
// MODE 0: non-relativistic, instantaneous Coulomb.  MODE 1: relativistic momenta, Lienard-Wiechert fields of the partner at the retarded time.
// Output: one line per cycle (pericentre if E0 >= 0.1, else angular revolution of r).
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <deque>
#include <complex>
#include <algorithm>
typedef std::complex<double> Cx;
static double envd(const char* k, double d) { const char* v = std::getenv(k); return v ? std::atof(v) : d; }
struct V { double x, y; };
static inline V operator+(V a, V b) { return {a.x + b.x, a.y + b.y}; }
static inline V operator-(V a, V b) { return {a.x - b.x, a.y - b.y}; }
static inline V operator*(double s, V a) { return {s * a.x, s * a.y}; }
static inline double dot(V a, V b) { return a.x * b.x + a.y * b.y; }
static inline double crs(V a, V b) { return a.x * b.y - a.y * b.x; }
static inline double nrm(V a) { return std::sqrt(dot(a, a)); }
struct St { V x1, x2, p1, p2; double phi; };
struct Hist { double t; V x, v, a; };
int main() {
  setvbuf(stdout, nullptr, _IOLBF, 0);
  const int MODE = (int)envd("MODE", 0), CLOCK = (int)envd("CLOCK", 3), WDELAY = (int)envd("WDELAY", 0);
  const double M1 = envd("M1", 1836.15267343), M2 = envd("M2", 1.0), C = envd("C", 10.0), G1 = envd("G1", 5.5856946893), G2 = envd("G2", 2.00231930436);
  const double MS1 = envd("MS1", 0.0), MS2 = envd("MS2", 0.0), N0 = envd("N0", 1.3), E0 = envd("E0", 0.0), TEND = envd("TEND", 3e4), H = envd("H", 0.005);
  const double W = envd("W", 0.25), DX = envd("DX", 0.05), GAM = envd("GAM", 1e-3), TAUF = envd("TAUF", 0.3);
  const double MT = M1 + M2, MU = M1 * M2 / MT, KAPPA = envd("KAPPA", -0.2 * GAM * GAM / MU);
  const double A0 = N0 * N0 / MU, T0 = 2 * M_PI * N0 * N0 * N0 / MU, TAUEX = TAUF * T0, BOX = std::ceil(A0 * (1 + E0) + 3.0);
  const bool spin = (MS1 != 0.0 || MS2 != 0.0), peri = E0 >= 0.1;
  if (MODE == 0 && spin) { std::fprintf(stderr, "spin needs MODE 1\n"); return 1; }
  const int NG = (int)std::lround(2 * BOX / DX) + 1;
  std::vector<Cx> F((size_t)NG * NG, Cx(0.0));
  // kernel offsets sorted by distance
  const int KR = (int)std::ceil(3 * W / DX) + 1; std::vector<std::pair<int, int>> off;
  for (int i = -KR; i <= KR; ++i) for (int j = -KR; j <= KR; ++j) if ((i * i + j * j) * DX * DX <= (3 * W + DX) * (3 * W + DX)) off.push_back({i, j});
  std::sort(off.begin(), off.end(), [](auto a, auto b) { return a.first * a.first + a.second * a.second < b.first * b.first + b.second * b.second; });
  // initial state: relative orbit at pericentre
  const double rp = A0 * (1 - E0), vp = std::sqrt((2 / rp - 1 / A0) / MU);
  V r0 = {rp, 0}, v0 = {0, vp};
  St s; s.x1 = (-M2 / MT) * r0; s.x2 = (M1 / MT) * r0; V v1 = (-M2 / MT) * v0, v2 = (M1 / MT) * v0;
  auto gam = [&](V v) { return MODE ? 1 / std::sqrt(1 - dot(v, v) / (C * C)) : 1.0; };
  s.p1 = (M1 * gam(v1)) * v1; s.p2 = (M2 * gam(v2)) * v2; s.phi = 0;
  auto vel = [&](V p, double m) { if (!MODE) return (1 / m) * p; const double g = std::sqrt(1 + dot(p, p) / (m * m * C * C)); return (1 / (g * m)) * p; };
  // history for retarded fields
  const size_t HR = 1 << 17; std::vector<Hist> h1(HR), h2(HR); size_t hn = 0;
  auto hget = [&](std::vector<Hist>& hh, double tr, V& x, V& v, V& a) {
    const size_t n = std::min(hn, HR); const Hist& last = hh[(hn - 1) % HR];
    if (tr >= last.t) { const double d = tr - last.t; x = last.x + d * last.v + (0.5 * d * d) * last.a; v = last.v + d * last.a; a = last.a; return; }
    size_t lo = hn - n, hi = hn - 1;   // absolute indices
    if (tr < hh[lo % HR].t) { std::fprintf(stderr, "history too short\n"); std::exit(2); }
    while (hi - lo > 1) { const size_t m = (lo + hi) / 2; if (hh[m % HR].t <= tr) lo = m; else hi = m; }
    const Hist &A = hh[lo % HR], &B = hh[hi % HR]; const double dt = B.t - A.t, u = (tr - A.t) / dt;
    const double h00 = 2 * u * u * u - 3 * u * u + 1, h10 = u * u * u - 2 * u * u + u, h01 = -2 * u * u * u + 3 * u * u, h11 = u * u * u - u * u;
    x = h00 * A.x + (h10 * dt) * A.v + h01 * B.x + (h11 * dt) * B.v; v = (1 - u) * A.v + u * B.v; a = (1 - u) * A.a + u * B.a; };
  // LW field of particle j (charge q) at point xo, time to
  auto lw = [&](std::vector<Hist>& hh, double q, V xo, double to, V& E, double& Bz) {
    V xj, vj, aj; double tr = to - nrm(xo - hh[(hn - 1) % HR].x) / C;
    for (int it = 0; it < 60; ++it) { hget(hh, tr, xj, vj, aj); const double tn = to - nrm(xo - xj) / C; if (std::fabs(tn - tr) < 1e-14) { tr = tn; break; } tr = tn; }
    hget(hh, tr, xj, vj, aj);
    const V Rv = xo - xj; const double R = nrm(Rv); const V n = (1 / R) * Rv, b = (1 / C) * vj, bd = (1 / C) * aj;
    const double k = 1 - dot(n, b), k3 = k * k * k; const V nb = n - b;
    // n x (nb x bd) = nb (n.bd) - bd (n.nb)
    const V acc = dot(n, bd) * nb - dot(n, nb) * bd;
    E = (q * (1 - dot(b, b)) / (k3 * R * R)) * nb + (q / (C * k3 * R)) * acc;   // (q/c) n x ((n - b) x bdot), bdot = a/c
    Bz = crs(n, E); };
  // field sampling: Catmull-Rom value and gradient of Re[psi e^{-i phi}] (stored psi' = psi e^{GAM t})
  auto cr = [](double t, double w[4], double d[4]) {
    const double t2 = t * t, t3 = t2 * t;
    w[0] = (-t3 + 2 * t2 - t) / 2; w[1] = (3 * t3 - 5 * t2 + 2) / 2; w[2] = (-3 * t3 + 4 * t2 + t) / 2; w[3] = (t3 - t2) / 2;
    d[0] = (-3 * t2 + 4 * t - 1) / 2; d[1] = (9 * t2 - 10 * t) / 2; d[2] = (-9 * t2 + 8 * t + 1) / 2; d[3] = (3 * t2 - 2 * t) / 2; };
  auto sample = [&](V r, double phi, double t, double& val, V& grad, double& amp) {
    const double gx = (r.x + BOX) / DX, gy = (r.y + BOX) / DX; const int ix = (int)std::floor(gx), iy = (int)std::floor(gy);
    val = 0; grad = {0, 0}; amp = 0; if (ix < 1 || iy < 1 || ix + 2 >= NG || iy + 2 >= NG) return;
    double wx[4], dx[4], wy[4], dy[4]; cr(gx - ix, wx, dx); cr(gy - iy, wy, dy);
    Cx f = 0, fx = 0, fy = 0;
    for (int a = 0; a < 4; ++a) for (int b = 0; b < 4; ++b) { const Cx z = F[(size_t)(ix - 1 + a) * NG + (iy - 1 + b)];
      f += wx[a] * wy[b] * z; fx += dx[a] * wy[b] * z; fy += wx[a] * dy[b] * z; }
    const Cx ph = std::exp(Cx(-GAM * t, -phi)); val = (f * ph).real(); grad = {(fx * ph).real() / DX, (fy * ph).real() / DX}; amp = std::abs(f) * std::exp(-GAM * t); };
  // derivatives
  struct D { V dx1, dx2, dp1, dp2; double dphi, om1, om2, act; };
  double tcur = 0;
  auto deriv = [&](const St& q, double t) {
    D d; const V u1 = vel(q.p1, M1), u2 = vel(q.p2, M2); d.dx1 = u1; d.dx2 = u2; d.om1 = d.om2 = 0;
    const V r = q.x2 - q.x1; const double rr = nrm(r);
    V F1, F2; double E1x = 0, E1y = 0, B1 = 0, E2x = 0, E2y = 0, B2 = 0;
    if (!MODE) { F2 = (-1 / (rr * rr * rr)) * r; F1 = (-1.0) * F2; }
    else { V E1, E2; lw(h2, -1.0, q.x1, t, E1, B1); lw(h1, +1.0, q.x2, t, E2, B2);   // field at 1 from 2 (charge -1) and at 2 from 1 (+1)
      F1 = (+1.0) * (E1 + (1 / C) * V{u1.y * B1, -u1.x * B1}); F2 = (-1.0) * (E2 + (1 / C) * V{u2.y * B2, -u2.x * B2});
      E1x = E1.x; E1y = E1.y; E2x = E2.x; E2y = E2.y; }
    double val, amp; V g; sample(r, q.phi, t, val, g, amp); const V Fg = KAPPA * g;
    d.dp1 = F1 - Fg; d.dp2 = F2 + Fg;
    d.act = dot(q.p1, u1) + dot(q.p2, u2);
    if (CLOCK == 3) d.dphi = d.act;
    else d.dphi = 0.5 * M1 * dot(u1, u1) + 0.5 * M2 * dot(u2, u2) + 1 / rr;
    if (spin) {   // Thomas-BMT about z: Omega = -(q/(m c)) [ (g/2 - 1 + 1/gamma) B - (g/2 - gamma/(gamma+1)) (beta x E) ]
      auto om = [&](double qq, double m, double g2, V u, double Ex, double Ey, double Bz) { const double b2 = dot(u, u) / (C * C), ga = 1 / std::sqrt(1 - b2);
        const double bxE = (u.x * Ey - u.y * Ex) / C; return -(qq / (m * C)) * ((g2 / 2 - 1 + 1 / ga) * Bz - (g2 / 2 - ga / (ga + 1)) * bxE); };
      d.om1 = om(+1.0, M1, G1, u1, E1x, E1y, B1); d.om2 = om(-1.0, M2, G2, u2, E2x, E2y, B2); d.dphi -= MS1 * d.om1 + MS2 * d.om2; }
    return d; };
  auto acc_of = [&](V p, V dp, double m) { if (!MODE) return (1 / m) * dp; const V u = vel(p, m); const double g = 1 / std::sqrt(1 - dot(u, u) / (C * C));
    return (1 / (g * m)) * (dp - (dot(dp, u) / (C * C)) * u); };
  auto push_hist = [&](const St& q, double t) { const D d = deriv(q, t);
    h1[hn % HR] = {t, q.x1, vel(q.p1, M1), acc_of(q.p1, d.dp1, M1)}; h2[hn % HR] = {t, q.x2, vel(q.p2, M2), acc_of(q.p2, d.dp2, M2)}; ++hn; };
  // seed history (static past: uniform motion backwards is not available; seed with constant state over a long past)
  if (MODE) {   // seed the past with the exact non-relativistic Kepler orbit (pericentre at t = 0)
    const double nm = std::sqrt(1 / (MU * A0 * A0 * A0)), tw = 2 * BOX / C + 1.0, sq = std::sqrt(1 - E0 * E0);
    for (double tp = -tw; tp < -1e-12; tp += 0.01) { double M = nm * tp, u = M; for (int it = 0; it < 80; ++it) u -= (u - E0 * std::sin(u) - M) / (1 - E0 * std::cos(u));
      const double ud = nm / (1 - E0 * std::cos(u)); const V r = {A0 * (std::cos(u) - E0), A0 * sq * std::sin(u)}, v = {-A0 * std::sin(u) * ud, A0 * sq * std::cos(u) * ud};
      const double rr = nrm(r); const V a = (-1 / (MU * rr * rr * rr)) * r;
      h1[hn % HR] = {tp, (-M2 / MT) * r, (-M2 / MT) * v, (-M2 / MT) * a}; h2[hn % HR] = {tp, (M1 / MT) * r, (M1 / MT) * v, (M1 / MT) * a}; ++hn; } }
  push_hist(s, 0.0);
  // writes FIFO
  struct Wr { V r; double phi, t, w; size_t k; };
  std::deque<Wr> fifo; const double gnorm = 1 / (std::sqrt(2 * M_PI) * W);
  auto deposit = [&](double t) {
    for (auto it = fifo.begin(); it != fifo.end();) {
      if (t < it->t + TAUEX) break;
      const int i0 = (int)std::lround((it->r.x + BOX) / DX), j0 = (int)std::lround((it->r.y + BOX) / DX); const Cx e = std::exp(Cx(0, it->phi));
      while (it->k < off.size()) { const int i = i0 + off[it->k].first, j = j0 + off[it->k].second;
        const double cx = -BOX + i * DX, cy = -BOX + j * DX, ddx = cx - it->r.x, ddy = cy - it->r.y, d2 = ddx * ddx + ddy * ddy;
        const double tarr = it->t + TAUEX + (WDELAY ? std::sqrt(d2) / C : 0.0); if (tarr > t) break;
        if (d2 <= 9 * W * W && i >= 0 && j >= 0 && i < NG && j < NG) F[(size_t)i * NG + j] += (it->w * std::exp(-d2 / (2 * W * W) + GAM * tarr)) * e;
        ++it->k; }
      if (it->k >= off.size()) it = fifo.erase(it); else ++it; } };
  std::printf("# 871 MODE %d CLOCK %d WDELAY %d M1 %.8f M2 %.8f C %.3f G1 %.10f G2 %.11f MS1 %+.1f MS2 %+.1f N0 %.4f E0 %.3f KAPPA %.6e GAM %.2e W %.3f DX %.3f TAUEX %.4f BOX %.1f TEND %.0f H %.4f\n",
              MODE, CLOCK, WDELAY, M1, M2, C, G1, G2, MS1, MS2, N0, E0, KAPPA, GAM, W, DX, TAUEX, BOX, TEND, H);
  std::printf("# k t J_tot A L dtheta-2pi E Phi1 Phi2 rmin rmax Rread amp\n");
  // cycle bookkeeping
  auto energy = [&](const St& q) { const V u1 = vel(q.p1, M1), u2 = vel(q.p2, M2), r = q.x2 - q.x1; const double rr = nrm(r);
    if (!MODE) return 0.5 * M1 * dot(u1, u1) + 0.5 * M2 * dot(u2, u2) - 1 / rr;
    const double g1 = 1 / std::sqrt(1 - dot(u1, u1) / (C * C)), g2 = 1 / std::sqrt(1 - dot(u2, u2) / (C * C)); const V n = (1 / rr) * r;
    return (g1 - 1) * M1 * C * C + (g2 - 1) * M2 * C * C - 1 / rr + (1 / (2 * C * C * rr)) * (dot(u1, u2) + dot(u1, n) * dot(u2, n)); };
  auto angmom = [&](const St& q) { return crs(q.x1, q.p1) + crs(q.x2, q.p2); };
  double Aacc = 0, P1 = 0, P2 = 0, Lacc = 0, Tacc = 0, rmin = 1e300, rmax = 0, rdlast = 0, Racc = 0, ampacc = 0; int kcyc = 0;
  double thUnw = 0, thPrev = std::atan2((s.x2 - s.x1).y, (s.x2 - s.x1).x), thCycStart = 0; double rpp = -1, rp1 = -1; St sp1 = s; double tp1 = 0;
  (void)rdlast;
  double t = 0; long step = 0;
  while (t < TEND) {
    const V rv = s.x2 - s.x1; const double rr = nrm(rv); const double he = H * std::min(1.0, std::pow(rr, 1.5));
    // RK4
    auto add = [&](const St& a, const D& d, double c) { St b = a; b.x1 = a.x1 + c * d.dx1; b.x2 = a.x2 + c * d.dx2; b.p1 = a.p1 + c * d.dp1; b.p2 = a.p2 + c * d.dp2; b.phi = a.phi + c * d.dphi; return b; };
    const D k1 = deriv(s, t), k2 = deriv(add(s, k1, he / 2), t + he / 2), k3 = deriv(add(s, k2, he / 2), t + he / 2), k4 = deriv(add(s, k3, he), t + he);
    St sn = s; sn.x1 = s.x1 + (he / 6) * (k1.dx1 + 2 * k2.dx1 + 2 * k3.dx1 + k4.dx1); sn.x2 = s.x2 + (he / 6) * (k1.dx2 + 2 * k2.dx2 + 2 * k3.dx2 + k4.dx2);
    sn.p1 = s.p1 + (he / 6) * (k1.dp1 + 2 * k2.dp1 + 2 * k3.dp1 + k4.dp1); sn.p2 = s.p2 + (he / 6) * (k1.dp2 + 2 * k2.dp2 + 2 * k3.dp2 + k4.dp2);
    sn.phi = s.phi + (he / 6) * (k1.dphi + 2 * k2.dphi + 2 * k3.dphi + k4.dphi);
    const double o1 = (k1.om1 + 2 * k2.om1 + 2 * k3.om1 + k4.om1) / 6, o2 = (k1.om2 + 2 * k2.om2 + 2 * k3.om2 + k4.om2) / 6;
    const double arate = (k1.act + 2 * k2.act + 2 * k3.act + k4.act) / 6;   // action rate sum p.v (any clock)
    // write (source at the start of the step)
    { const V u1 = vel(s.p1, M1), u2 = vel(s.p2, M2); fifo.push_back({rv, s.phi, t, nrm(u2 - u1) * he * gnorm, 0}); }
    double val, amp; V g; sample(rv, s.phi, t, val, g, amp); Racc += val * he; ampacc += amp * he;
    s = sn; t += he; ++step;
    if (MODE) push_hist(s, t);
    deposit(t);
    Aacc += arate * he; P1 += o1 * he; P2 += o2 * he; Lacc += angmom(s) * he; Tacc += he;
    const V rn = s.x2 - s.x1; const double rnn = nrm(rn); rmin = std::min(rmin, rnn); rmax = std::max(rmax, rnn);
    const double th = std::atan2(rn.y, rn.x); double dth = th - thPrev; if (dth > M_PI) dth -= 2 * M_PI; if (dth < -M_PI) dth += 2 * M_PI; thUnw += dth; thPrev = th;
    if (!std::isfinite(rnn) || rnn > BOX - 1) { std::printf("# STOP t %.4f r %.6e\n", t, rnn); break; }
    bool cyc = false; double tfrac = 0, thAt = thUnw;   // tfrac: time after the event within the last step(s)
    if (peri) { if (rp1 >= 0 && rpp >= 0 && rp1 < rpp && rp1 <= rnn) {   // minimum at the previous step
        const double a2 = (rpp - 2 * rp1 + rnn) / 2, b = (rnn - rpp) / 2; const double x = (a2 > 0) ? -b / (2 * a2) : 0.0;   // offset in steps from previous point
        tfrac = (1 - x) * he; thAt = thUnw - (1 - x) * dth; cyc = true; }
      rpp = rp1; rp1 = rnn; }
    else { const double target = thCycStart + 2 * M_PI * ((thUnw > thCycStart) ? 1 : -1);
      if ((thUnw - target) * (thUnw - dth - target) <= 0 && std::fabs(thUnw - thCycStart) >= 2 * M_PI - 1e-12) { tfrac = (thUnw - target) / dth * he; thAt = target; cyc = true; } }
    if (cyc) {
      const double Ac = Aacc - arate * tfrac, P1c = P1 - o1 * tfrac, P2c = P2 - o2 * tfrac, Tc = Tacc - tfrac, Lc = (Lacc - angmom(s) * tfrac) / Tc;
      const double dthc = thAt - thCycStart; const double J = Ac / (2 * M_PI) - Lc * (std::fabs(dthc) - 2 * M_PI) / (2 * M_PI);
      if (kcyc > 0) std::printf("%d %.4f %.9f %.9f %.9f %.6e %.9e %.6e %.6e %.6e %.6e %.4e %.4e\n", kcyc, t - tfrac, J, Ac / (2 * M_PI), Lc, std::fabs(dthc) - 2 * M_PI,
                                energy(s), P1c, P2c, rmin, rmax, Racc / Tc, ampacc / Tc);
      ++kcyc; Aacc = arate * tfrac; P1 = o1 * tfrac; P2 = o2 * tfrac; Lacc = angmom(s) * tfrac; Tacc = tfrac; thCycStart = thAt; rmin = 1e300; rmax = 0; Racc = ampacc = 0; }
    (void)sp1; (void)tp1;
  }
  std::printf("# END t %.4f steps %ld\n", t, step);
  return 0; }

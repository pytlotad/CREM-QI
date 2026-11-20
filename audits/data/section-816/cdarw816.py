# Audit 816: action-angle variables of the Darwin Lagrangian with exact kinetic energy, centre-of-mass frame (written before computing).
#   v1 = (m2/M) v, v2 = -(m1/M) v;  L = -m1 c^2/g1 - m2 c^2/g2 - kq/r - A (2 v_r^2 + v_t^2),  A = dar * kq m1 m2/(2 c^2 r M^2)
#   E = g1 m1 c^2 + g2 m2 c^2 + kq/r - A (2 v_r^2 + v_t^2);  p_t = (F - 2A) v_t, L_can = r p_t;  p_r = (F - 4A) v_r,
#   F = g1 m1 (m2/M)^2 + g2 m2 (m1/M)^2.   J_r(E, L) = (1/pi) int p_r dr;  Omega = -omega_r (dJ_r/dL + 1), omega_r = 1/(dJ_r/dE).
import mpmath as mp
mp.mp.dps = 50
HBAR = mp.mpf('6.62607015e-34')/(2*mp.pi)

class Pair:
    def __init__(s, m1, m2, kq, c, dar=1):
        s.m1, s.m2, s.kq, s.c, s.dar = mp.mpf(m1), mp.mpf(m2), mp.mpf(kq), mp.mpf(c), dar
        s.M = s.m1 + s.m2; s.mu = s.m1*s.m2/s.M; s.f1 = s.m2/s.M; s.f2 = s.m1/s.M
    def A(s, r): return s.dar*s.kq*s.m1*s.m2/(2*s.c**2*r*s.M**2)
    def gam(s, v2):
        c2 = s.c**2; return 1/mp.sqrt(1 - s.f1**2*v2/c2), 1/mp.sqrt(1 - s.f2**2*v2/c2)
    def F(s, v2):
        g1, g2 = s.gam(v2); return g1*s.m1*s.f1**2 + g2*s.m2*s.f2**2
    def energy(s, r, vr, vt):
        v2 = vr**2 + vt**2; g1, g2 = s.gam(v2)
        return g1*s.m1*s.c**2 + g2*s.m2*s.c**2 + s.kq/r - s.A(r)*(2*vr**2 + vt**2)
    def lcan(s, r, vr, vt): return r*(s.F(vr**2 + vt**2) - 2*s.A(r))*vt
    def e_rs(s, r, sq, L):       # energy at (r, v_r^2 = sq) on L_can = L (analytic in sq, also for sq < 0)
        vt = mp.findroot(lambda vt: r*(s.F(sq + vt**2) - 2*s.A(r))*vt - L, L/(s.mu*r))
        v2 = sq + vt**2; g1, g2 = s.gam(v2)
        return g1*s.m1*s.c**2 + g2*s.m2*s.c**2 + s.kq/r - s.A(r)*(2*sq + vt**2), vt
    def sq_of(s, r, E, L):
        s0 = 2*(E - s.M*s.c**2 - s.kq/r)/s.mu - (L/(s.mu*r))**2
        return mp.findroot(lambda q: s.e_rs(r, q, L)[0] - E, s0)
    def pr(s, r, E, L):
        q = s.sq_of(r, E, L); vt = s.e_rs(r, q, L)[1]
        return (s.F(q + vt**2) - 4*s.A(r))*mp.sqrt(q) if q > 0 else mp.mpf(0)
    def turning(s, E, L):
        En = E - s.M*s.c**2; a = s.kq/(2*En); e = mp.sqrt(max(mp.mpf(0), 1 + 2*En*L**2/(s.mu*s.kq**2)))
        g = lambda r: s.e_rs(r, 0, L)[0] - E; rc = L**2/(s.mu*(-s.kq))
        return (mp.findroot(g, (a*(1 - e)/2, rc), solver='anderson'), mp.findroot(g, (rc, 2*a*(1 + e)), solver='anderson'))
    def Jr(s, E, L):
        r1, r2 = s.turning(E, L)
        return mp.quad(lambda r: s.pr(r, E, L), [r1, r2])/mp.pi
    def omegas(s, E, L, rel=mp.mpf('1e-12')):
        hE = abs(E - s.M*s.c**2)*rel; hL = L*rel
        dE = (s.Jr(E + hE, L) - s.Jr(E - hE, L))/(2*hE); dL = (s.Jr(E, L + hL) - s.Jr(E, L - hL))/(2*hL)
        wr = 1/dE; return -wr*(dL + 1), wr
    def ecirc(s, L):             # J_r = 0: minimum of E(r, v_r = 0) at fixed L_can
        r0 = L**2/(s.mu*(-s.kq))
        rc = mp.findroot(lambda r: mp.diff(lambda x: s.e_rs(x, 0, L)[0], r), r0)
        return s.e_rs(rc, 0, L)[0]
    def etorus(s, n, L, E0=None):  # E with J_r(E, L) + L = n hbar
        if n*HBAR - L <= 0: return s.ecirc(L)
        E = E0 if E0 is not None else s.M*s.c**2 - s.mu*s.kq**2/(2*(n*HBAR)**2)
        for _ in range(30):
            f = s.Jr(E, L) + L - n*HBAR; hE = abs(E - s.M*s.c**2)*mp.mpf('1e-12')
            d = (s.Jr(E + hE, L) - s.Jr(E - hE, L))/(2*hE); dEn = f/d; E -= dEn
            if abs(dEn) < abs(E - s.M*s.c**2)*mp.mpf('1e-28'): break
        return E
    def start_invariants(s, st):
        x = [mp.mpf(float(q)) for q in st['x1'] - st['x2']]; v1 = [mp.mpf(float(q)) for q in st['v1']]
        v2 = [mp.mpf(float(q)) for q in st['v2']]; c2 = s.c**2
        r = mp.sqrt(sum(q*q for q in x)); n = [q/r for q in x]
        g1 = 1/mp.sqrt(1 - sum(q*q for q in v1)/c2); g2 = 1/mp.sqrt(1 - sum(q*q for q in v2)/c2)
        v12 = sum(a*b for a, b in zip(v1, v2)); n1 = sum(a*b for a, b in zip(v1, n)); n2 = sum(a*b for a, b in zip(v2, n))
        E = g1*s.m1*c2 + g2*s.m2*c2 + s.kq/r + s.dar*s.kq/(2*c2*r)*(v12 + n1*n2)
        v = [a - b for a, b in zip(v1, v2)]; vr = sum(a*b for a, b in zip(v, n)); vt = mp.sqrt(sum(q*q for q in v) - vr**2)
        return E, s.lcan(r, vr, vt), r, vr, vt

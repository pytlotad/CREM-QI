"""Audit 728 (offline): spin-orbit precession about J (H = H0 + A L.S), space reading of the orbital action phase alone and together with the spin's
geometric phase (S.J^ - |S|) phi_S.  Units mu = k = hbar = 1; z := J^.  See prereg728.txt."""
import math, os
import numpy as np
from scipy.integrate import solve_ivp
from scipy.spatial import cKDTree
HERE = os.path.dirname(os.path.abspath(__file__))
PSI_NODE = 0.3; BETA = 0.05; A_SO = 0.002; NS = 200; KS = (4000, 20000)
ALPHA_S = 1.0
NS_ = np.array([math.sin(ALPHA_S), 0.0, math.cos(ALPHA_S)])
out = []
def P(s): out.append(s); print(s, flush=True)


def cb(Jr, L, beta):
    Lam = math.sqrt(L*L + 2*beta); N = Jr + Lam; a = N*N; e = math.sqrt(max(0.0, 1 - (Lam/N)**2)); Om = 1/N**3
    return Lam, N, a, e, Om


def nu_of_u(u, e):
    k = np.floor(u/(2*np.pi)); w = u - 2*np.pi*k
    return 2*np.arctan2(np.sqrt(1 + e)*np.sin(w/2), np.sqrt(1 - e)*np.cos(w/2)) + 2*np.pi*k


def triangle(L, S, J):
    return (J*J - L*L - S*S)/2, (J*J + L*L - S*S)/(2*J), (J*J + S*S - L*L)/(2*J)   # L.S, L.J^, S.J^


class Torus:
    def __init__(s, Jr, L, S, J, spin_phase=True, A=A_SO, beta=BETA):
        s.Jr, s.L, s.S, s.J, s.sp, s.A, s.beta = Jr, L, S, J, spin_phase, A, beta
        s.LS, s.LJ, s.SJ = triangle(L, S, J)
        s.Lam, s.N, s.a, s.e, s.Om = cb(Jr, L, beta)
        s.ci = s.LJ/L; s.si = math.sqrt(max(0.0, 1 - s.ci**2)); s.Wp = A*J

    def state(s, u):
        e, a = s.e, s.a
        t = (u - e*np.sin(u))/s.Om
        nu = nu_of_u(u, e); psi = (s.L/s.Lam)*nu - PSI_NODE - s.A*s.L*t           # in-plane angle from node, with the -A|L| rotation
        r = a*(1 - e*np.cos(u)); rdot = a*e*np.sin(u)*s.Om/(1 - e*np.cos(u)); psidot = s.L/r**2 - s.A*s.L
        cp, sp = np.cos(psi), np.sin(psi)
        ub = np.stack([cp, sp*s.ci, sp*s.si]); eb = np.stack([-sp, cp*s.ci, cp*s.si])
        qb = r*ub; vb = rdot*ub + r*(s.L/r**2)*eb                                   # canonical p: H0 velocity, rotated with the plane
        c, sn = np.cos(s.Wp*t), np.sin(s.Wp*t)
        q = np.stack([c*qb[0] - sn*qb[1], sn*qb[0] + c*qb[1], qb[2]])
        p = np.stack([c*vb[0] - sn*vb[1], sn*vb[0] + c*vb[1], vb[2]])
        S0 = s.N*(u + e*np.sin(u)) - 2*s.beta*nu/s.Lam + s.A*s.LS*t
        phiS = 0.5*np.pi + s.Wp*t
        if s.sp: S0 = S0 + (s.SJ - s.S)*phiS
        nth = (np.floor((psi - 0.5*np.pi)/np.pi) - math.floor((-PSI_NODE - 0.5*np.pi)/np.pi)) if s.si > 0 else 0*u
        mas = np.floor(u/np.pi) + nth
        wind = np.stack([np.floor(u/(2*np.pi)), np.floor(psi/(2*np.pi)), np.floor(s.Wp*t/(2*np.pi))])
        return q, p, t, S0, mas, wind, phiS

    def crossings(s, K):
        us = []
        for k0 in range(0, K, 2000):
            k1 = min(K, k0 + 2000)
            u = np.linspace(2*np.pi*k0, 2*np.pi*k1, (k1 - k0)*NS + 1)
            y = (s.state(u)[0]*NS_[:, None]).sum(0)
            idx = np.where(np.sign(y[1:]) != np.sign(y[:-1]))[0]
            lo, hi = u[idx].copy(), u[idx + 1].copy(); ylo = y[idx]
            for _ in range(45):
                mid = 0.5*(lo + hi); ym = (s.state(mid)[0]*NS_[:, None]).sum(0); same = np.sign(ym) == np.sign(ylo)
                lo = np.where(same, mid, lo); hi = np.where(same, hi, mid); ylo = np.where(same, ym, ylo)
            us.append(0.5*(lo + hi))
        return np.concatenate(us)

    def B4(s, K):
        u = s.crossings(K)
        q, p, t, S0, mas, wind, phiS = s.state(u)
        rr = np.sqrt((q**2).sum(0)); rhat = q/rr
        th = np.arccos(np.clip(q[2]/rr, -1, 1)); rho = np.hypot(q[0], q[1])
        eth = np.stack([np.cos(th)*q[0]/rho, np.cos(th)*q[1]/rho, -np.sin(th)])
        pr = (p*rhat).sum(0); pth = rr*(p*eth).sum(0)
        ph = S0 - 0.5*np.pi*mas
        z = []; typ = []
        for b1 in (1, -1):
            for b2 in (1, -1):
                sel = np.where((np.sign(pr) == b1) & (np.sign(pth) == b2))[0]
                if len(sel) < 3: continue
                pts = q[:, sel].T; tree = cKDTree(pts)
                _, nn = tree.query(pts, k=2); j = nn[:, 1]; i = np.arange(len(sel))
                pairs = set((min(x, y), max(x, y)) for x, y in zip(i, j))
                A_, B_ = map(np.array, zip(*pairs)); Ai, Bi = sel[A_], sel[B_]
                dq = q[:, Bi] - q[:, Ai]; pm = 0.5*(p[:, Ai] + p[:, Bi])
                d = (ph[Bi] - ph[Ai]) - (pm*dq).sum(0)
                if s.sp:
                    dphi = phiS[Bi] - phiS[Ai]; d = d - (s.SJ - s.S)*((dphi + np.pi) % (2*np.pi) - np.pi)
                z.append(np.exp(1j*d)); dw = np.abs(wind[:, Bi] - wind[:, Ai]).astype(np.int64)
                typ.append(dw[0]*10**12 + dw[1]*10**6 + dw[2])
        z = np.concatenate(z); typ = np.concatenate(typ)
        vals = [z[typ == tt].mean().real for tt in np.unique(typ) if (typ == tt).sum() >= 5]
        return (min(vals) if vals else float('nan')), len(u), len(vals)


P('Audit 728: H = H0 + A L.S (A = %.3f, beta = %.2f); z = J^; section plane normal (sin 1, 0, cos 1); branches (p_r, p_theta); K = %s' % (A_SO, BETA, KS))
# ---- H728a full equations
P('\n(A) full equations (q, p, S), DOP853 rtol 1e-12, 10 radial periods')
okA = True
for (Jr, L, S, J) in ((0.5, 1.5, 1.5, 1.5), (1.5, 1.5, 1.5, 2.5), (0.5, 2.5, 1.5, 2.0)):
    LS, LJ, SJ = triangle(L, S, J); Lam, N, a, e, Om = cb(Jr, L, BETA)
    ci = LJ/L; si = math.sqrt(1 - ci*ci)
    Lv = L*np.array([0, -si, ci]); Sv = np.array([0, 0, J]) - Lv
    rp = a*(1 - e); vp = Lam/rp*0 + L/rp                                           # periapsis on the node line, speed L/r_p (p_r = 0)
    q0 = rp*np.array([1.0, 0, 0]); e2 = np.cross(Lv/L, q0/rp); p0 = vp*e2
    def f(t, y):
        q, p, Sx = y[0:3], y[3:6], y[6:9]; r = np.linalg.norm(q)
        grad = q/r**3 - 2*BETA*q/r**4
        Lq = np.cross(q, p)
        return np.concatenate([p + A_SO*np.cross(Sx, q), -grad + A_SO*np.cross(Sx, p), A_SO*np.cross(Lq, Sx)])
    T = 10*2*np.pi*N**3
    sol = solve_ivp(f, (0, T), np.concatenate([q0, p0, Sv]), method='DOP853', rtol=1e-12, atol=1e-14, dense_output=True)
    ts = np.linspace(0, T, 20001); Y = sol.sol(ts); q, p, Sx = Y[0:3], Y[3:6], Y[6:9]
    Lt = np.cross(q.T, p.T).T; Jt = Lt + Sx; Jn = np.linalg.norm(Jt, axis=0)
    dJ = abs(Jn/Jn[0] - 1).max(); Jh = Jt[:, 0]/Jn[0]
    dLJ = abs((Lt*Jh[:, None]).sum(0) - LJ).max(); dSJ = abs((Sx*Jh[:, None]).sum(0) - SJ).max()
    phiL = np.unwrap(np.arctan2(Lt[1], Lt[0])); rate = np.polyfit(ts, phiL, 1)[0]
    qd = np.stack([f(0, Y[:, k])[0:3] for k in range(0, 20001, 1000)], 1); pk = p[:, ::1000]; qk = q[:, ::1000]; Sk = Sx[:, ::1000]
    Tk = 0.5*(pk**2).sum(0); Lk = np.cross(qk.T, pk.T).T
    ident = abs((pk*qd).sum(0) - 2*Tk - A_SO*(Lk*Sk).sum(0)).max()
    ok = dJ < 1e-9 and dLJ < 1e-8 and dSJ < 1e-8 and abs(rate/(A_SO*J) - 1) < 1e-4 and ident < 1e-8; okA &= ok
    P('  (J_r %.1f, |L| %.1f, |S| %.1f, |J| %.1f): |J| drift %.1e, L.J^ dev %.1e, S.J^ dev %.1e, precession rate / A|J| = %.7f, p.qdot - 2T - A L.S %.1e' %
      (Jr, L, S, J, dJ, dLJ, dSJ, rate/(A_SO*J), ident))
P('H728a: %s' % okA)


def row(Jr, L, S, J, sp):
    v = [Torus(Jr, L, S, J, spin_phase=sp).B4(K) for K in KS]
    return [x[0] for x in v], v[-1][1], v[-1][2]


# ---- H728b orbital only, |S| = 3/2
P('\n(B) orbital phase only, |S| = 3/2, |L| = 3/2 (o-Ps 2^3P_J)')
okB = True
for nr in (0, 1):
    for J, kind in ((0.5, 'Q'), (1.5, 'Q'), (2.5, 'Q'), (2.0, 'int')):
        b, nc, nt = row(nr + 0.5, 1.5, 1.5, J, False); LS, LJ, SJ = triangle(1.5, 1.5, J)
        okB &= (max(b) < 0.5) if kind == 'Q' else (min(b) >= 0.99)
        P('  n_r %d |J| %.1f (%s; L.J^ %.2f): B4 %s (%d cr., %d types)' % (nr, J, kind, LJ, '/'.join('%.4f' % x for x in b), nc, nt))
P('H728b (orbital phase alone quantizes L.J^: rule-Q points < 0.5, L.J^ = 1 >= 0.99): %s' % okB)

# ---- H728c orbital + spin, |S| = 3/2
P('\n(C) orbital + spin geometric phase, |S| = 3/2')
okC = True
for nr in (0, 1):
    for J, kind in ((0.5, 'Q'), (1.5, 'Q'), (2.5, 'Q'), (1.0, 'int'), (2.0, 'int')):
        b, nc, nt = row(nr + 0.5, 1.5, 1.5, J, True); LS, LJ, SJ = triangle(1.5, 1.5, J)
        okC &= (min(b) >= 0.99) if kind == 'Q' else (max(b) < 0)
        P('  n_r %d |J| %.1f (%s; L.J^ %.2f, S.J^ %.2f): B4 %s (%d cr., %d types)' % (nr, J, kind, LJ, SJ, '/'.join('%.4f' % x for x in b), nc, nt))
P('H728c (rule Q |J| = J + 1/2 selected, integer |J| rejected): %s' % okC)

# ---- H728d orbital + spin, |S| = 1
P('\n(D) orbital + spin geometric phase, bare |S| = 1')
okD = True
for nr in (0, 1):
    for J, kind in ((1.0, 'int'), (2.0, 'int'), (1.5, 'half'), (2.5, 'half')):
        b, nc, nt = row(nr + 0.5, 1.5, 1.0, J, True); LS, LJ, SJ = triangle(1.5, 1.0, J)
        okD &= (min(b) >= 0.99) if kind == 'int' else (max(b) < 0)
        P('  n_r %d |J| %.1f (%s; L.J^ %.3f, S.J^ %.3f): B4 %s' % (nr, J, kind, LJ, SJ, '/'.join('%.4f' % x for x in b)))
P('H728d (|S| = 1: integer |J| selected, half-integer rejected): %s' % okD)

# ---- H728e random
P('\n(E) 100 random tori, orbital + spin phase, |S| = 3/2')
rng = np.random.default_rng(728); hit = []
for _ in range(100):
    Jr, L = rng.uniform(0.05, 2.0), rng.uniform(0.8, 3.0); J = rng.uniform(abs(L - 1.5) + 0.05, L + 1.5 - 0.05)
    b = [Torus(Jr, L, 1.5, J).B4(K)[0] for K in KS]
    if min(b) > 0.5: hit.append('(%.3f, %.3f, %.3f) %s' % (Jr, L, J, '/'.join('%.3f' % x for x in b)))
okE = len(hit) <= 2
P('  B4 > 0.5 at both K: %d / 100 %s' % (len(hit), '; '.join(hit)))
P('H728e (<= 2 %%): %s' % okE)
open(os.path.join(HERE, 'out', 'ana728.txt'), 'w').write('\n'.join(out) + '\n')

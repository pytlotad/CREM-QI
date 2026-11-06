"""Audit 727 (offline): L_z from the space reading of the action phase when the node precesses (H = H0 + omega L_z, linear orbital Zeeman).
Body frame = 726 orbit (Coulomb + beta/r^2, inclined plane, node fixed); lab = rotation by omega t about z.  Units mu = k = hbar = 1.  See prereg727.txt."""
import math, os
import numpy as np
from scipy.spatial import cKDTree
HERE = os.path.dirname(os.path.abspath(__file__))
PSI_NODE = 0.3; PHI0 = 0.7; BETA = 0.05; OMEGA = 0.0015; NS = 200
out = []
def P(s): out.append(s); print(s, flush=True)


def cb(Jr, L, beta):
    Lam = math.sqrt(L*L + 2*beta); N = Jr + Lam; a = N*N; e = math.sqrt(max(0.0, 1 - (Lam/N)**2)); Om = 1/N**3
    return Lam, N, a, e, Om


def nu_of_u(u, e):
    k = np.floor(u/(2*np.pi)); w = u - 2*np.pi*k
    return 2*np.arctan2(np.sqrt(1 + e)*np.sin(w/2), np.sqrt(1 - e)*np.cos(w/2)) + 2*np.pi*k


class Torus:
    def __init__(s, Jr, L, Lz, beta=BETA, omega=OMEGA):
        s.Jr, s.L, s.Lz, s.beta, s.om = Jr, L, Lz, beta, omega
        s.Lam, s.N, s.a, s.e, s.Om = cb(Jr, L, beta)
        s.ci = Lz/L; s.si = math.sqrt(max(0.0, 1 - s.ci**2))

    def state(s, u):
        """lab position, canonical momentum, time, S0 (body + omega Lz t), Maslov count, winding counters"""
        e, a = s.e, s.a
        t = (u - e*np.sin(u))/s.Om
        nu = nu_of_u(u, e); phi_in = (s.L/s.Lam)*nu; psi = phi_in - PSI_NODE
        r = a*(1 - e*np.cos(u)); rdot = a*e*np.sin(u)*s.Om/(1 - e*np.cos(u)); psidot = s.L/r**2
        cp, sp = np.cos(psi), np.sin(psi)
        ub = np.stack([cp, sp*s.ci, sp*s.si]); eb = np.stack([-sp, cp*s.ci, cp*s.si])
        qb = r*ub; vb = rdot*ub + r*psidot*eb
        c, sn = np.cos(s.om*t), np.sin(s.om*t)
        q = np.stack([c*qb[0] - sn*qb[1], sn*qb[0] + c*qb[1], qb[2]])
        p = np.stack([c*vb[0] - sn*vb[1], sn*vb[0] + c*vb[1], vb[2]])
        S0 = s.N*(u + e*np.sin(u)) - 2*s.beta*nu/s.Lam + s.om*s.Lz*t
        nth = (np.floor((phi_in - (PSI_NODE + 0.5*np.pi))/np.pi) + 1) if s.si > 0 else 0*u
        mas = np.floor(u/np.pi) + nth
        wind = np.stack([np.floor(u/(2*np.pi)), np.floor(phi_in/(2*np.pi)), np.floor(s.om*t/(2*np.pi))])
        return q, p, t, S0, mas, wind, r, rdot

    def yprime(s, u):
        q = s.state(u)[0]
        return NS_[0]*q[0] + NS_[1]*q[1] + NS_[2]*q[2], np.ones_like(q[0])

    def crossings(s, K):
        us = []
        for k0 in range(0, K, 2000):
            k1 = min(K, k0 + 2000)
            u = np.linspace(2*np.pi*k0, 2*np.pi*k1, (k1 - k0)*NS + 1)
            y, _ = s.yprime(u)
            idx = np.where(np.sign(y[1:]) != np.sign(y[:-1]))[0]
            lo, hi = u[idx].copy(), u[idx + 1].copy(); ylo = y[idx]
            for _ in range(45):
                mid = 0.5*(lo + hi); ym, _ = s.yprime(mid); same = np.sign(ym) == np.sign(ylo)
                lo = np.where(same, mid, lo); hi = np.where(same, hi, mid); ylo = np.where(same, ym, ylo)
            us.append(0.5*(lo + hi))
        return np.concatenate(us)

    def B4(s, K):
        """Addendum 1: section = plane with normal NS_ (not containing the z axis)."""
        u = s.crossings(K)
        q, p, t, S0, mas, wind, r, rdot = s.state(u)
        rr = np.sqrt((q**2).sum(0)); rhat = q/rr
        x1 = (q*E1_[:, None]).sum(0); x2 = (q*E2_[:, None]).sum(0)
        p1 = (p*E1_[:, None]).sum(0); p2 = (p*E2_[:, None]).sum(0)
        pr = (p*rhat).sum(0)
        th = np.arccos(np.clip(q[2]/rr, -1, 1)); rho = np.hypot(q[0], q[1])
        eth = np.stack([np.cos(th)*q[0]/rho, np.cos(th)*q[1]/rho, -np.sin(th)]); pth = rr*(p*eth).sum(0)
        ph = S0 - 0.5*np.pi*mas
        z = []; typ = []
        for b1 in (1, -1):
            for b2 in (1, -1):
                if True:
                    sel = np.where((np.sign(pr) == b1) & (np.sign(pth) == b2))[0]   # Addendum 2: branches (p_r, p_theta) as preregistered
                    if len(sel) < 3: continue
                    pts = np.stack([x1[sel], x2[sel]], 1); tree = cKDTree(pts)
                    _, nn = tree.query(pts, k=2); j = nn[:, 1]; i = np.arange(len(sel))
                    pairs = set((min(a, b), max(a, b)) for a, b in zip(i, j))
                    A_, Bn_ = map(np.array, zip(*pairs)); A, Bn = sel[A_], sel[Bn_]
                    d = (ph[Bn] - ph[A]) - 0.5*(p1[A] + p1[Bn])*(x1[Bn] - x1[A]) - 0.5*(p2[A] + p2[Bn])*(x2[Bn] - x2[A])
                    z.append(np.exp(1j*d)); dw = np.abs(wind[:, Bn] - wind[:, A]).astype(np.int64)
                    typ.append(dw[0]*10**12 + dw[1]*10**6 + dw[2])
        z = np.concatenate(z); typ = np.concatenate(typ)
        vals = [z[typ == tt].mean().real for tt in np.unique(typ) if (typ == tt).sum() >= 5]
        return (min(vals) if vals else float('nan')), len(u), len(vals)


ALPHA_S = 1.0
NS_ = np.array([math.sin(ALPHA_S), 0.0, math.cos(ALPHA_S)]); E1_ = np.array([math.cos(ALPHA_S), 0.0, -math.sin(ALPHA_S)]); E2_ = np.array([0.0, 1.0, 0.0])
KS = (4000, 20000)
P('Audit 727 (run 3, Addenda 1-2: section = plane with normal (sin 1, 0, cos 1), branches (sign p_r, sign p_theta)): L_z from the space reading with node precession, H = H0 + omega L_z; beta %.2f, omega %.4f, K = %s' % (BETA, OMEGA, KS))
# ---- H727a identities
P('\n(A) identities: int p.qdot dt (trapezoid, dense) vs S0_body + omega L_z t; theta caustics lab vs body')
okA = True
for (Jr, L, Lz) in ((0.5, 1.5, 1.0), (1.5, 2.5, -2.0), (0.5, 2.5, 0.0)):
    T = Torus(Jr, L, Lz); u = np.linspace(0, 2*np.pi*3, 600001)
    q, p, t, S0, mas, wind, r, rdot = T.state(u)
    qd = np.gradient(q, t, axis=1); integrand = (p*qd).sum(0)
    S_num = np.concatenate([[0], np.cumsum(0.5*(integrand[1:] + integrand[:-1])*np.diff(t))])
    rel = abs((S_num[-1] - (S0[-1] - S0[0]))/(S0[-1] - S0[0]))
    th_lab = np.arccos(np.clip(q[2]/np.sqrt((q**2).sum(0)), -1, 1)); d = np.diff(th_lab)
    c_lab = int(np.sum(np.sign(d[1:]) != np.sign(d[:-1])))
    Tb = Torus(Jr, L, Lz, omega=0.0); qb = Tb.state(u)[0]; thb = np.arccos(np.clip(qb[2]/np.sqrt((qb**2).sum(0)), -1, 1)); db = np.diff(thb)
    c_body = int(np.sum(np.sign(db[1:]) != np.sign(db[:-1])))
    okA &= rel < 1e-7 and c_lab == c_body
    P('  (J_r %.1f, L %.1f, L_z %+.1f): S rel. diff %.1e; theta caustics lab %d, body %d' % (Jr, L, Lz, rel, c_lab, c_body))
P('H727a (S to 1e-7, caustics equal): %s' % okA)

# ---- H727b
P('\n(B) Langer tori (n_r + 1/2, l + 1/2) with L_z = m, m + 1/2, m + 0.3 (omega = %.4f)' % OMEGA)
okB = True; Erep = []
for nr in (0, 1):
    for l in (1, 2):
        for m in range(-l, l + 1):
            Jr, L = nr + 0.5, l + 0.5
            row = []
            for lab, Lz in (('m', m), ('m+1/2', m + 0.5), ('m+0.3', m + 0.3)):
                if abs(Lz) >= L: row.append('%s n/a' % lab); continue
                v = [Torus(Jr, L, Lz).B4(K) for K in KS]; b = [x[0] for x in v]
                if lab == 'm': okB &= min(b) >= 0.99
                elif lab == 'm+1/2': okB &= max(b) < 0
                else: okB &= max(b) < 0.5
                row.append('%s %s (%d cr., %d types)' % (lab, '/'.join('%.4f' % x for x in b), v[-1][1], v[-1][2]))
            Lam, N, a, e, Om = cb(Jr, L, BETA); Erep.append((nr, l, m, -0.5/N**2 + OMEGA*m))
            P('  n_r %d l %d m %+d: %s' % (nr, l, m, ';  '.join(row)))
P('H727b (L_z = m >= 0.99; m + 1/2 < 0; m + 0.3 < 0.5; both K): %s' % okB)

# ---- H727c control omega = 0
P('\n(C) control omega = 0 (node fixed): Langer tori with L_z = m + 0.3')
okC = True
for nr, l, m in ((0, 1, 0), (0, 2, 1), (1, 1, -1), (1, 2, 0)):
    v = [Torus(nr + 0.5, l + 0.5, m + 0.3, omega=0.0).B4(K)[0] for K in KS]; okC &= min(v) >= 0.99
    P('  n_r %d l %d L_z %+.1f: B4 %s' % (nr, l, m + 0.3, '/'.join('%.4f' % x for x in v)))
P('H727c (L_z not selected without precession, B4 >= 0.99): %s' % okC)

# ---- H727d random
P('\n(D) 100 random tori (omega = %.4f)' % OMEGA)
rng = np.random.default_rng(727); hit = []
for _ in range(100):
    Jr, L = rng.uniform(0.05, 2.0), rng.uniform(0.6, 3.0); Lz = L*rng.uniform(-0.95, 0.95)
    v = [Torus(Jr, L, Lz).B4(K)[0] for K in KS]
    if min(v) > 0.5: hit.append('(%.3f, %.3f, %+.3f) %s' % (Jr, L, Lz, '/'.join('%.3f' % x for x in v)))
okD = len(hit) <= 2
P('  B4 > 0.5 at both K: %d / 100 %s' % (len(hit), '; '.join(hit)))
P('H727d (<= 2 %%): %s' % okD)
P('\n(E) energies E = E0 + omega m at the quantized tori (normal orbital Zeeman):')
for nr, l, m, E in Erep: P('  n_r %d l %d m %+d: E %.12f' % (nr, l, m, E))
open(os.path.join(HERE, 'out', 'ana727.txt'), 'w').write('\n'.join(out) + '\n')

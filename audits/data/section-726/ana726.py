"""Audit 726 (offline): the space reading of the action phase S0/hbar (723) on the 3D-embedded orbit of Coulomb + beta/r^2: Maslov -pi/2 at every
radial AND polar (theta) turning point.  Planar orbit inclined by i (cos i = L_z/L), node at in-plane angle PSI_NODE from periapsis.
Units mu = k = hbar = 1.  See prereg726.txt."""
import math, os
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
PHI0 = 1.0; PSI_NODE = 0.3
out = []
def P(s): out.append(s); print(s, flush=True)


def cb(Jr, L, beta):
    Lam = math.sqrt(L*L + 2*beta); N = Jr + Lam; a = N*N; e = math.sqrt(max(0.0, 1 - (Lam/N)**2)); Om = 1/N**3
    return Lam, N, a, e, Om


def nu_of_u(u, e):
    k = np.floor(u/(2*np.pi)); w = u - 2*np.pi*k
    return 2*np.arctan2(np.sqrt(1 + e)*np.sin(w/2), np.sqrt(1 - e)*np.cos(w/2)) + 2*np.pi*k


def u_of_nu(nu, e):
    k = np.floor(nu/(2*np.pi)); w = nu - 2*np.pi*k
    return 2*np.arctan2(np.sqrt(1 - e)*np.sin(w/2), np.sqrt(1 + e)*np.cos(w/2)) + 2*np.pi*k


def S0_cb(u, Jr, L, beta):
    Lam, N, a, e, Om = cb(Jr, L, beta)
    return N*(u + e*np.sin(u)) - 2*beta*nu_of_u(u, e)/Lam


def n_theta_turn(phi, cosi):
    """theta turning points passed in (0, phi]: extrema of u_z = sin(psi) sin i, psi = phi - PSI_NODE, i.e. cos psi = 0; none if i = 0."""
    if abs(cosi) >= 1.0: return np.zeros_like(phi)
    t0 = PSI_NODE + 0.5*np.pi
    return np.floor((phi - t0)/np.pi) + 1


def B4(Jr, L, beta, K, cosi):
    Lam, N, a, e, Om = cb(Jr, L, beta)
    nutot = 2*np.pi*K
    m = np.arange(0, int((L/Lam)*nutot/(2*np.pi)) + 2)
    nu = (PHI0 + 2*np.pi*m)*Lam/L; keep = nu < nutot; nu = nu[keep]; m = m[keep]
    u = u_of_nu(nu, e); kr = np.floor(u/(2*np.pi)); phi = (L/Lam)*nu
    r = a*(1 - e*np.cos(u)); pr = a*e*np.sin(u)*Om/(1 - e*np.cos(u))
    ph = S0_cb(u, Jr, L, beta) - 0.5*np.pi*(np.floor(u/np.pi) + n_theta_turn(phi, cosi))
    z = []; typ = []
    for br in (1, -1):
        sel = (pr*br > 0)
        if sel.sum() < 3: continue
        o = np.argsort(r[sel]); rr = r[sel][o]; pp = pr[sel][o]; hh = ph[sel][o]; mm = m[sel][o]; kk = kr[sel][o]
        d = np.diff(hh) - 0.5*(pp[1:] + pp[:-1])*np.diff(rr)
        z.append(np.exp(1j*d)); typ.append(np.abs(np.diff(mm))*100000 + np.abs(np.diff(kk)))
    z = np.concatenate(z); typ = np.concatenate(typ)
    return min(z[typ == t].mean().real for t in np.unique(typ) if (typ == t).sum() >= 5)


KS = (400, 4000, 40000)
P('Audit 726: action phase S0/hbar, space reading on the 3D-embedded orbit (Coulomb + beta/r^2), Maslov -pi/2 at r and theta turning points; '
  'B4 = min over winding types of Re<exp(i d)>; K = %s radial periods' % (KS,))
# ---- H726a geometry
P('\n(A) theta turning points from dense sampling of theta(t) vs predicted 2 per in-plane revolution (0 if equatorial)')
okA = True
for (Jr, L, beta, cosi) in ((0.5, 1.5, 0.05, 0.5), (1.5, 2.5, 0.05, 0.9), (0.3, 1.2, 0.2, 0.5), (2.0, 0.8, 0.2, 0.9), (0.5, 1.5, 0.05, 1.0)):
    Lam, N, a, e, Om = cb(Jr, L, beta)
    u = np.linspace(0, 2*np.pi*5, 400001); nu = nu_of_u(u, e); phi = (L/Lam)*nu
    si = math.sqrt(max(0.0, 1 - cosi*cosi)); psi = phi - PSI_NODE
    uz = np.sin(psi)*si; theta = np.arccos(np.clip(uz, -1, 1))
    dth = np.diff(theta); cnt = int(np.sum(np.sign(dth[1:]) != np.sign(dth[:-1]))) if si > 0 else 0
    pred = int(n_theta_turn(np.array([phi[-1]]), cosi)[0])
    rev = phi[-1]/(2*np.pi); okA &= (cnt == pred)
    P('  J_r %.1f L %.1f beta %.2f cos i %.1f: 5 radial periods, %.4f revolutions -> theta turning points sampled %d, predicted %d' % (Jr, L, beta, cosi, rev, cnt, pred))
P('H726a (sampled = predicted, 2 per revolution, 0 equatorial): %s' % okA)

# ---- H726b Langer vs 2D points
P('\n(B) Langer points (n_r + 1/2, l + 1/2) and 2D points (n_r + 1/2, m), cos i = 0.5')
okB = True
for beta in (0.05, 0.2):
    for nr in range(3):
        for l in range(3):
            Jr = nr + 0.5; L = l + 0.5; v = [B4(Jr, L, beta, K, 0.5) for K in KS]; okB &= min(v) >= 0.999
            Lam, N, a, e, Om = cb(Jr, L, beta); Nqm = nr + 0.5 + math.sqrt((l + 0.5)**2 + 2*beta)
            P('  beta %.2f Langer (%.1f, %.1f): B4 %s;  E %.12f = 3D QM %.12f (N %.10f)' % (beta, Jr, L, ' / '.join('%.4f' % x for x in v), -0.5/N**2, -0.5/Nqm**2, Nqm))
        for mm in (1, 2, 3):
            v = [B4(nr + 0.5, float(mm), beta, K, 0.5) for K in KS]; okB &= max(v) < 0
            P('  beta %.2f 2D     (%.1f, %d  ): B4 %s' % (beta, nr + 0.5, mm, ' / '.join('%.4f' % x for x in v)))
P('H726b (Langer B4 >= 0.999 at all K; 2D points B4 < 0 at all K): %s' % okB)

# ---- H726c selectivity
P('\n(C) 200 random tori, cos i = 0.5')
rng = np.random.default_rng(726); hit = []; n = 0
for beta in (0.05, 0.2):
    for _ in range(100):
        Jr, L = rng.uniform(0.05, 3.0), rng.uniform(0.6, 3.4); n += 1
        v = [B4(Jr, L, beta, K, 0.5) for K in KS]
        if min(v) > 0.5: hit.append('(%.3f, %.3f, beta %.2f) %s' % (Jr, L, beta, '/'.join('%.3f' % x for x in v)))
okC = len(hit) <= 0.02*n
P('  B4 > 0.5 at all K: %d / %d %s' % (len(hit), n, '; '.join(hit)))
P('H726c (<= 2 %%): %s' % okC)

# ---- H726d equatorial
P('\n(D) equatorial orbit (cos i = 1, no theta caustic)')
okD = True
for beta in (0.05, 0.2):
    for nr in range(2):
        v2 = [[B4(nr + 0.5, float(mm), beta, K, 1.0) for K in KS] for mm in (1, 2)]
        vl = [[B4(nr + 0.5, l + 0.5, beta, K, 1.0) for K in KS] for l in (1, 2)]
        okD &= all(min(v) >= 0.999 for v in v2) and all(max(v) < 0 for v in vl)
        P('  beta %.2f n_r %d: 2D points L = 1, 2: %s;  Langer L = 1.5, 2.5: %s' % (beta, nr, ' | '.join('/'.join('%.4f' % x for x in v) for v in v2),
                                                                                 ' | '.join('/'.join('%.4f' % x for x in v) for v in vl)))
P('H726d (equatorial: 2D points >= 0.999, Langer < 0): %s' % okD)
open(os.path.join(HERE, 'out', 'ana726.txt'), 'w').write('\n'.join(out) + '\n')

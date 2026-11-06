"""Audit 723 (offline): the action-phase clock S0/hbar read in TIME (one point, long-memory field: amplitude A) and in SPACE (same point, same branch,
successive orbits, Maslov -pi/2 per turning point: coherence B) on Coulomb + beta/r^2 (non-closed orbits), with Kepler and oscillator controls.
Units mu = k = hbar = omega = 1, planar motion.  See prereg723.txt."""
import math, os
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import minimize_scalar
HERE = os.path.dirname(os.path.abspath(__file__))
K = 400; PHI0 = 1.0
out = []
def P(s): out.append(s); print(s, flush=True)


def wrap(x): return (x + math.pi) % (2*math.pi) - math.pi


# ---------- analytic Coulomb + beta/r^2 torus (beta = 0: Kepler)
def cb(Jr, L, beta):
    Lam = math.sqrt(L*L + 2*beta); N = Jr + Lam; a = N*N; e = math.sqrt(max(0.0, 1 - (Lam/N)**2)); Om = 1/N**3
    return Lam, N, a, e, Om


def nu_of_u(u, e):   # continuous true anomaly
    k = np.floor(u/(2*np.pi)); w = u - 2*np.pi*k
    return 2*np.arctan2(np.sqrt(1 + e)*np.sin(w/2), np.sqrt(1 - e)*np.cos(w/2)) + 2*np.pi*k


def u_of_nu(nu, e):
    k = np.floor(nu/(2*np.pi)); w = nu - 2*np.pi*k
    return 2*np.arctan2(np.sqrt(1 - e)*np.sin(w/2), np.sqrt(1 + e)*np.cos(w/2)) + 2*np.pi*k


def S0_cb(u, Jr, L, beta):
    Lam, N, a, e, Om = cb(Jr, L, beta)
    return N*(u + e*np.sin(u)) - 2*beta*nu_of_u(u, e)/Lam


def A_time(Jr, L, beta):
    Lam, N, a, e, Om = cb(Jr, L, beta)
    Mn = (np.arange(2048) + 0.5)*2*np.pi/2048; u = Mn.copy()
    for _ in range(50): u = u - (u - e*np.sin(u) - Mn)/(1 - e*np.cos(u))
    c0 = np.mean(np.exp(-1j*S0_cb(u, Jr, L, beta)))
    dS = 2*np.pi*(Jr + L*L/Lam)
    s = np.sum(np.exp(-1j*dS*np.arange(K)))
    return abs(c0)*abs(s)/K


def B_space(Jr, L, beta, kind='min'):
    """kind 'pool': run-1 estimator (all neighbour pairs pooled); 'min': min over winding types (k_r, k_phi) with >= 5 pairs (Addendum 1)."""
    Lam, N, a, e, Om = cb(Jr, L, beta)
    nutot = 2*np.pi*K
    m = np.arange(0, int((L/Lam)*nutot/(2*np.pi)) + 2)
    nu = (PHI0 + 2*np.pi*m)*Lam/L; keep = nu < nutot; nu = nu[keep]; m = m[keep]
    u = u_of_nu(nu, e); kr = np.floor(u/(2*np.pi))
    r = a*(1 - e*np.cos(u)); pr = a*e*np.sin(u)*Om/(1 - e*np.cos(u))
    ph = S0_cb(u, Jr, L, beta) - 0.5*np.pi*np.floor(u/np.pi)
    z = []; typ = []
    for br in (1, -1):
        sel = (pr*br > 0)
        if sel.sum() < 3: continue
        o = np.argsort(r[sel]); rr = r[sel][o]; pp = pr[sel][o]; hh = ph[sel][o]; mm = m[sel][o]; kk = kr[sel][o]
        d = np.diff(hh) - 0.5*(pp[1:] + pp[:-1])*np.diff(rr)
        z.append(np.exp(1j*d)); typ.append(np.abs(np.diff(mm))*100000 + np.abs(np.diff(kk)))
    z = np.concatenate(z); typ = np.concatenate(typ)
    if kind == 'pool': return abs(z.mean())
    if kind == 'abs': vals = [abs(z[typ == t].mean()) for t in np.unique(typ) if (typ == t).sum() >= 5]   # run 2 (B3)
    else: vals = [z[typ == t].mean().real for t in np.unique(typ) if (typ == t).sum() >= 5]           # B4 (Addendum 2)
    return min(vals)


# ---------- generic one-period DOP853 (control and oscillator)
def period(Vp, E, L, r0, tmax):
    def f(t, y):
        r, pr = y[0], y[1]
        return [pr, L*L/r**3 - Vp(r), L/r**2, pr*pr + L*L/r**2, 0.5*(pr*pr + L*L/r**2)]
    def ev(t, y): return y[1]
    ev.terminal = True; ev.direction = -1
    s = solve_ivp(f, (0, tmax), [r0, 0, 0, 0, 0], method='DOP853', rtol=1e-12, atol=1e-14, events=ev)
    y = s.y_events[0][-1]; half = s.t_events[0][-1]
    Tr = 2*half; return {'Om_r': 2*math.pi/Tr, 'Om_phi': 2*y[2]/Tr, 'Jr': y[3]/math.pi - L*2*y[2]/(2*math.pi)*0, 'S0': 2*y[3], 'T': y[4]/half,
                         'Jr_true': (2*y[3] - L*2*y[2])/(2*math.pi)}


P('Audit 723: action-phase clock S0/hbar -- time reading A vs space reading B (mu = k = hbar = omega = 1; planar; K = %d radial periods)' % K)
# ---------- H723a
P('\n(A) identity 2<T> = J_r Omega_r + L Omega_phi (DOP853, one period) and the analytic S0 per period')
worst = 0.0; n = 0
for beta in (0.0, 0.05, 0.2):
    for Jr in (0.3, 0.5, 1.2, 2.0):
        for L in (0.7, 1.0, 2.0):
            Lam, N, a, e, Om = cb(Jr, L, beta); E = -1/(2*N*N); rp = a*(1 - e)
            p = period(lambda r: 1/r**2 - 2*beta/r**3, E, L, rp, 20*math.pi*N**3)
            lhs = 2*p['T']; rhs = p['Jr_true']*p['Om_r'] + L*p['Om_phi']
            dS_an = 2*math.pi*(Jr + L*L/Lam)
            w = max(abs(lhs/rhs - 1), abs(p['Jr_true'] - Jr), abs(p['S0'] - dS_an)/dS_an, abs(p['Om_r']/Om - 1)); worst = max(worst, w); n += 1
            if n % 6 == 1:
                P('  beta %.2f J_r %.1f L %.1f: 2<T> %.12f  J_r Om_r + L Om_phi %.12f  J_r(int) - J_r %+.1e  S0/period vs analytic %+.1e' %
                  (beta, Jr, L, lhs, rhs, p['Jr_true'] - Jr, (p['S0'] - dS_an)/dS_an))
for E, L in ((1.0, 0.3), (2.5, 1.5), (4.0, 0.4)):
    rp = math.sqrt(E - math.sqrt(E*E - L*L)); p = period(lambda r: r, E, L, rp, 10*math.pi)
    lhs = 2*p['T']; rhs = p['Jr_true']*p['Om_r'] + L*p['Om_phi']; worst = max(worst, abs(lhs/rhs - 1), abs(p['Om_phi']/p['Om_r'] - 0.5)); n += 1
    P('  oscillator E %.1f L %.1f: 2<T> %.12f  J_r Om_r + L Om_phi %.12f  Om_phi/Om_r %.12f  J_r %.10f (= (E - L)/2 %.10f)' %
      (E, L, lhs, rhs, p['Om_phi']/p['Om_r'], p['Jr_true'], (E - L)/2))
P('  %d tori, max relative deviation %.1e' % (n, worst))
# 722 roots
t722 = open(os.path.join(HERE, '../section-722/out/ana722.txt')).read().splitlines()
w722 = 0.0
for line in t722:
    if line.strip().startswith('beta') and 'clock N' in line:
        w = line.split(); beta = float(w[1]); j = int(w[3]); L = float(w[5].rstrip(':')); N = float(w[8])
        Lam = math.sqrt(L*L + 2*beta); Jr = N - Lam; v = Jr + L*L/Lam; w722 = max(w722, abs(v - j))
        P('  722 root beta %.2f j %d L %.1f: N %.10f -> J_r + L^2/Lambda = %.10f (j = %d)' % (beta, j, L, N, v, j))
P('  722 oscillator roots E = 2j: J_r + L/2 = (E - L)/2 + L/2 = E/2 = j identically')
okA = worst < 1e-9 and w722 < 1e-6
P('H723a (identity to 1e-9; 722 roots on J_r + L Om_phi/Om_r = j to 1e-6): %s (%.1e, %.1e)' % (okA, worst, w722))


def scan(fun, xs):
    return np.array([fun(x) for x in xs])


def maxima(xs, ys, thr):
    idx = [i for i in range(1, len(xs) - 1) if ys[i] > ys[i - 1] and ys[i] >= ys[i + 1] and ys[i] > thr]
    return idx


def refine(fun, x0, h):
    r = minimize_scalar(lambda x: -fun(x), bounds=(x0 - h, x0 + h), method='bounded', options={'xatol': 1e-10})
    return r.x, -r.fun


# ---------- H723b: time reading
P('\n(B) time reading A = |<exp(-i S0/hbar)>_t| (no Maslov: the clock has none), beta = 0.05')
okB = True; beta = 0.05
JR = np.round(np.arange(0.0005, 3.0, 0.0005), 4)
for L in (1.0, 2.0):
    Lam = math.sqrt(L*L + 2*beta)
    ys = scan(lambda x: A_time(x, L, beta), JR); mx = maxima(JR, ys, 0.5)
    curve = [j - L*L/Lam for j in range(1, 6) if 0 < j - L*L/Lam < 3]
    res = []
    for i in mx:
        x, v = refine(lambda q: A_time(q, L, beta), JR[i], 0.001); dc = min(abs(x - c) for c in curve) if curve else 9
        okB &= dc < 0.002; res.append('%.6f (A %.3f, curve %+.1e)' % (x, v, dc))
    ebk = [(nr + 0.5, A_time(nr + 0.5, L, beta)) for nr in range(3)]
    okB &= all(v < 0.2 for _, v in ebk)
    P('  L %.0f: maxima A > 0.5 at J_r = %s; curve J_r = j - L^2/Lambda: %s' % (L, ', '.join(res), ', '.join('%.6f' % c for c in curve)))
    P('        A at EBK points J_r = 0.5, 1.5, 2.5: %s' % ', '.join('%.4f' % v for _, v in ebk))
P('H723b (time maxima on the R0 curve, EBK points A < 0.2): %s' % okB)

# ---------- H723c: space reading
P('\n(C) space reading B4 = min over winding types of Re<exp(i d)> (Addendum 2, POST HOC; run 1 pooled |.|: out/ana723_run1.txt; run 2 min |.|: out/ana723_run2.txt) (same point on a ray phi = %.1f, same branch, Maslov -pi/2 per turning point)' % PHI0)
okC = True; LS = np.round(np.arange(0.6, 3.4, 0.0001), 4); JRf = np.round(np.arange(0.0001, 3.0, 0.0001), 4)
for beta in (0.05, 0.2):
    P('  beta %.2f' % beta)
    for L in (1.0, 2.0, 3.0):
        Lam = math.sqrt(L*L + 2*beta)
        ys = scan(lambda x: B_space(x, L, beta), JRf); mx = maxima(JRf, ys, 0.9); res = []
        for i in mx:
            x, v = refine(lambda q: B_space(q, L, beta), JRf[i], 0.0001); d = min(abs(x - (k + 0.5)) for k in range(4))
            okC &= d < 0.002; N = x + Lam; Nq = round(x - 0.5) + 0.5 + Lam
            res.append('%.7f (B %.4f, |dJ| %.1e, E %.10f vs 2D QM %.10f, rel %.1e)' % (x, v, d, -0.5/N**2, -0.5/Nq**2, abs(Nq**2/N**2 - 1)))
            okC &= abs(Nq**2/N**2 - 1) < 1e-5
        r0 = [j - L*L/Lam for j in range(1, 6) if 0 < j - L*L/Lam < 3]
        br0 = [B_space(c, L, beta) for c in r0]; okC &= all(v < 0.5 for v in br0)
        P('    L %.0f, J_r scan: maxima B > 0.9 at %s' % (L, '; '.join(res) if res else 'none'))
        P('          B at the R0 curve (time-reading resonances) J_r = %s: %s' % (', '.join('%.4f' % c for c in r0), ', '.join('%.4f' % v for v in br0)))
    for Jr in (0.5, 1.5):
        ys = scan(lambda x: B_space(Jr, x, beta), LS); mx = maxima(LS, ys, 0.9); res = []
        for i in mx:
            x, v = refine(lambda q: B_space(Jr, q, beta), LS[i], 0.0001); d = abs(x - round(x)); okC &= d < 0.002
            res.append('%.7f (B %.4f, |dL| %.1e)' % (x, v, d))
        P('    J_r %.1f, L scan [0.6, 3.4]: maxima B > 0.9 at %s' % (Jr, '; '.join(res) if res else 'none'))
rng = np.random.default_rng(723); hi = 0; tot = 0
for beta in (0.05, 0.2):
    for _ in range(100):
        v = B_space(rng.uniform(0.05, 3.0), rng.uniform(0.6, 3.4), beta); tot += 1; hi += v > 0.5
P('  200 random tori: B4 > 0.5 in %d (%.1f %%)' % (hi, 100*hi/tot))
P('  EBK points B4: ' + ', '.join('beta %.2f (%.1f, %d) %.4f' % (b, nr + 0.5, mm, B_space(nr + 0.5, mm, b)) for b in (0.05, 0.2) for nr in range(3) for mm in (1, 2, 3)))
P('H723c\'\' (B4 space maxima at EBK J_r = n_r + 1/2, L = m to 0.002; E = 2D QM to 1e-5; B < 0.5 on the R0 curve): %s' % okC)

# ---------- H723d: closed-orbit controls
P('\n(D) closed orbits')
okD = True
vals = []
for L in (0.7, 1.0, 1.3, 2.2):
    for n in (1, 2, 3):
        Jr = n + 0.5 - L
        if Jr <= 0: continue
        v = B_space(Jr, L, 0.0); vals.append(v); okD &= v > 0.9
off = [B_space(1.2, 1.0, 0.0), B_space(0.9, 1.3, 0.0)]
P('  Kepler on the line J_r + L = n + 1/2 (L = 0.7, 1.0, 1.3, 2.2; n = 1-3): B min %.4f max %.4f; off the line (J_r + L = 2.2, 2.2): %s' %
  (min(vals), max(vals), ', '.join('%.4f' % v for v in off)))


def B_osc(E, L):
    # oscillator: x = A cos t, y = B sin t, A B = L, (A^2 + B^2)/2 = E; S0(t) = int 2T dt; analytic
    s = math.sqrt(E*E - L*L); A = math.sqrt(E + s); Bb = math.sqrt(E - s)
    # ray phi = PHI0: tan(phi) = (Bb/A) tan t -> t_c = atan((A/Bb) tan PHI0) + pi k
    t0 = math.atan((A/Bb)*math.tan(PHI0)); tc = t0 + np.pi*np.arange(0, 2*K)
    x = A*np.cos(tc); y = Bb*np.sin(tc); keep = (np.cos(np.arctan2(y, x) - PHI0) > 0); tc = tc[keep]
    r = np.sqrt((A*np.cos(tc))**2 + (Bb*np.sin(tc))**2)
    pr = ((Bb*Bb - A*A)*np.sin(tc)*np.cos(tc))/r
    S0 = E*tc + 0.5*(Bb*Bb - A*A)*np.sin(2*tc)*0.5*(-1)*(-1)      # int (A^2 sin^2 + B^2 cos^2) dt = E t + (B^2 - A^2) sin 2t / 4
    S0 = E*tc + (Bb*Bb - A*A)*np.sin(2*tc)/4
    # radial turning points: d r^2/dt = (B^2 - A^2) sin 2t = 0 -> t = k pi/2 (t > 0)
    nturn = np.floor(tc/(np.pi/2))
    ph = S0 - 0.5*np.pi*nturn
    z = []
    for br in (1, -1):
        sel = pr*br > 0
        if sel.sum() < 2: continue
        o = np.argsort(r[sel]); rr = r[sel][o]; pp = pr[sel][o]; hh = ph[sel][o]
        if len(rr) < 2: continue
        d = np.diff(hh) - 0.5*(pp[1:] + pp[:-1])*np.diff(rr); z.append(np.exp(1j*d))
    z = np.concatenate(z); return abs(z.mean())


ov = [B_osc(E, L) for E in (2.0, 3.0, 4.0) for L in (0.5, 1.3, 1.9) if L < E]
oo = [B_osc(E, L) for E in (2.5, 3.3) for L in (0.5, 1.3)]
okD &= min(ov) > 0.9
P('  oscillator at E = 2, 3, 4 (any L): B min %.4f; at E = 2.5, 3.3: %s  (2D QM E = (2 n_r + |m| + 1): integers)' % (min(ov), ', '.join('%.4f' % v for v in oo)))
P('  (on a closed orbit returns coincide: the same point recurs with the same (k_r, k_phi) -- one condition only, as in QM degeneracy)')
P('H723d (closed orbits: one condition only -- Kepler line, oscillator E integer): %s' % okD)
open(os.path.join(HERE, 'out', 'ana723.txt'), 'w').write('\n'.join(out) + '\n')

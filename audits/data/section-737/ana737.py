"""Audit 737 (offline): second-order (hbar^2) correction of the EBK wave v2 -- iterate the uniform oscillator mapping with the Schwarzian term
D[sigma] = sigma'^{1/2} (sigma'^{-1/2})'' (order 0 = Langer, D = -1/4r^2).  E = -1/2n^2 fixed.  See prereg737.txt.
Run 2 (after run 1, see out/run1/): high-precision mapping -- phase integral in singularity-free variables (theta between the turning points, t^2 beyond
them, ln r near 0; cumulative Simpson), exact inversion of the oscillator action (sigma = -/+ s0 cos eps, s0 cosh psi), Gamma without cancellation,
D in closed form from h = ln(c Q/Gamma)."""
import math, os, time
import numpy as np
from scipy.special import genlaguerre, gammaln
from scipy.integrate import cumulative_simpson
from scipy.interpolate import CubicSpline
from scipy.optimize import brentq
HERE = os.path.dirname(os.path.abspath(__file__))
out = []
def P(s): out.append(s); print(s, flush=True)
ALPHA = 1/137.035999177; TAU_AU = 2.4188843265864e-17; MU_H = 1/(1 + 1/1836.152673426)
NQ = 200001          # quadrature points per region
WIN = 2e-3           # fill window around a turning point (relative to r_t)


def grid(nmax, fine=1):
    Rmax = 3.0*nmax*nmax + 60
    return np.concatenate([np.geomspace(1e-7, 1e-2, 3000*fine)[:-1], np.linspace(1e-2, Rmax, 1500000*fine + 1)])


def u_qm(n, l, r):
    rho = 2*r/n
    lognorm = 0.5*(3*math.log(2/n) + gammaln(n - l) - math.log(2*n) - gammaln(n + l + 1))
    with np.errstate(divide='ignore', invalid='ignore', over='ignore'):
        u = r*np.exp(lognorm - rho/2 + l*np.log(rho))*genlaguerre(n - l - 1, 2*l + 1)(rho)
    u = np.nan_to_num(u)
    return u/math.sqrt(np.trapezoid(u*u, r))


def hermite_fn(k, x):
    h0 = np.pi**-0.25*np.exp(-x*x/2)
    if k == 0: return h0
    h1 = math.sqrt(2)*x*h0
    for j in range(1, k):
        h0, h1 = h1, math.sqrt(2/(j + 1))*x*h1 - math.sqrt(j/(j + 1))*h0
    return h1


def solve_cubic_like(T, f, fp, x0):
    """solve f(x) = T (f ~ (4/3) x^3 near 0, monotone) by Newton on f^{1/3}; T >= 0 array."""
    x = x0.copy(); T3 = np.cbrt(T)
    for _ in range(40):
        fx = f(x); c = np.cbrt(np.maximum(fx, 1e-300))
        d = (c - T3)/(fp(x)/(3*c*c))
        d = np.where(fx <= 0, -x/2, d)
        x = np.maximum(x - d, x*0.1)
        if np.max(np.abs(d)) < 1e-15: break
    return x


def inv_allowed(F, s0):
    """F in [0, pi s0^2/2] (action from r1): returns sigma and Gamma = s0^2 - sigma^2 = s0^2 sin^2 eps."""
    T = 4*F/(s0*s0); up = T > math.pi; Ts = np.where(up, 2*math.pi - T, T); Ts = np.clip(Ts, 0, math.pi)
    e0 = np.clip(np.cbrt(0.75*Ts), 1e-12, math.pi/2)
    e = solve_cubic_like(Ts, lambda x: 2*x - np.sin(2*x), lambda x: 4*np.sin(x)**2, e0)
    sig = np.where(up, s0*np.cos(e), -s0*np.cos(e))
    return sig, s0*s0*np.sin(e)**2


def inv_forbidden(Gv, s0):
    """Gv >= 0 (action from the turning point outward): returns |sigma| = s0 cosh psi and Gamma = -s0^2 sinh^2 psi."""
    T = 4*Gv/(s0*s0)
    p0 = np.where(T < 1, np.cbrt(0.75*T), 0.5*np.arcsinh(T/2 + 1) + 0.3); p0 = np.maximum(p0, 1e-12)
    p = solve_cubic_like(T, lambda x: np.sinh(2*x) - 2*x, lambda x: 4*np.sinh(x)**2, p0)
    return s0*np.cosh(p), -s0*s0*np.sinh(p)**2


def action(Qf, r1, r2, r):
    """C(r) - C(r1) with sign convention: inside >= 0 rising to Phi at r2; outside: action measured outward from the nearer turning point (>= 0)."""
    th = np.linspace(0, math.pi, NQ); rr = r1 + (r2 - r1)*(1 - np.cos(th))/2
    ig = np.sqrt(np.maximum(Qf(rr), 0))*(r2 - r1)/2*np.sin(th)
    Cin = cumulative_simpson(ig, x=th, initial=0); Phi = Cin[-1]; spin = CubicSpline(th, Cin)
    to = np.linspace(0, math.sqrt(r[-1] - r2 + 1), NQ); ro = r2 + to*to
    Co = cumulative_simpson(np.sqrt(np.maximum(-Qf(ro), 0))*2*to, x=to, initial=0); spo = CubicSpline(to, Co)
    ti = np.linspace(0, math.sqrt(r1/2), NQ); ri = r1 - ti*ti
    Ci = cumulative_simpson(np.sqrt(np.maximum(-Qf(ri), 0))*2*ti, x=ti, initial=0); spi = CubicSpline(ti, Ci)
    xl = np.linspace(math.log(r1/2), math.log(r[0]) - 0.1, NQ); rl = np.exp(xl)
    Cl = Ci[-1] + cumulative_simpson(np.sqrt(np.maximum(-Qf(rl), 0))*rl, x=-xl, initial=0); spl = CubicSpline(-xl, Cl)
    C = np.empty_like(r)
    a = r < r1; b = (r >= r1) & (r <= r2); c = r > r2
    C[b] = spin(np.arccos(np.clip(1 - 2*(r[b] - r1)/(r2 - r1), -1, 1)))
    C[c] = spo(np.sqrt(r[c] - r2))
    aa = a & (r >= r1/2); ab = r < r1/2
    C[aa] = spi(np.sqrt(r1 - r[aa])); C[ab] = spl(-np.log(r[ab]))
    return C, Phi, (a, b, c)


def mapping(Qf, r, nr, guess, glob=False):
    s0 = math.sqrt(2*nr + 1)
    def root(g, lo, hi):
        x = np.linspace(lo, hi, 4001); v = Qf(x); i = np.where(np.sign(v[1:]) != np.sign(v[:-1]))[0]
        j = i[np.argmin(np.abs(x[i] - g))]
        return brentq(Qf, x[j], x[j + 1], xtol=1e-15, rtol=1e-15)
    mid = 0.5*(guess[0] + guess[1])
    r1 = root(guess[0], guess[0]*0.5, mid); r2 = root(guess[1], mid, guess[1]*1.5)
    C, Phi, (a, b, c) = action(Qf, r1, r2, r)
    scale = (math.pi*s0*s0/2)/Phi
    sig = np.empty_like(r); Gam = np.empty_like(r)
    sig[b], Gam[b] = inv_allowed(C[b]*scale, s0)
    fo = scale if glob else 1.0                     # glob: one phase scale everywhere (smooth sigma'); else 733 recipe (inside only)
    x, g = inv_forbidden(C[a]*fo, s0); sig[a] = -x; Gam[a] = g
    x, g = inv_forbidden(C[c]*fo, s0); sig[c] = x; Gam[c] = g
    cc = np.where(b, scale*scale, fo*fo)
    return sig, Gam, cc, Phi/(math.pi*(nr + 0.5)), (r1, r2)


def wave(sig, r, nr):
    ds = np.gradient(sig, r)
    with np.errstate(divide='ignore', invalid='ignore'):
        u = hermite_fn(nr, sig)/np.sqrt(np.abs(ds))
    u = np.nan_to_num(u, nan=0.0, posinf=0.0)
    return u/math.sqrt(np.trapezoid(u*u, r))


def smooth_fill(r, y, tp, win):
    """replace y near each turning point by a cubic fit from |x| in [win, 6 win] r_t, blended in with a C-infinity taper (1 inside win, 0 beyond 4 win)."""
    y = y.copy()
    for rt in tp:
        w = win*rt; x = r - rt; ax = np.abs(x)
        fit = (ax >= w) & (ax < 6*w); near = ax < 4*w
        cf = np.polyfit(x[fit]/w, y[fit], 3); pv = np.polyval(cf, x[near]/w)
        t = np.clip((ax[near] - w)/(3*w), 0, 1)
        with np.errstate(divide='ignore', over='ignore'):
            b = np.where(t <= 0, 1.0, np.where(t >= 1, 0.0, 1/(1 + np.exp(np.clip(1/(1 - t) - 1/np.maximum(t, 1e-300), -700, 700)))))
        yn = y[near]; good = np.isfinite(yn)
        y[near] = np.where(good, (1 - b)*np.where(good, yn, 0) + b*pv, pv)
    return y


def schwarz(sig, Gam, cc, r, Q, Qp, Qpp, tp, win=WIN):
    """D = h'^2/16 - h''/4, h = ln(c Q/Gamma), sigma' = sqrt(c Q/Gamma); 0/0 at the turning points filled by a cubic fit from [win, 4 win] r_t."""
    with np.errstate(divide='ignore', invalid='ignore'):
        sp = np.sqrt(np.abs(cc*Q/Gam))
    if Qp is None:                                  # k >= 1: Q_k known only on the grid -> D = f''/f, f = sigma'^{-1/2} (sigma' analytic from Q/Gamma)
        sp = smooth_fill(r, sp, tp, win)            # sigma' through the turning points
        f = sp**-0.5
        D = np.gradient(np.gradient(f, r), r)/f
        D = smooth_fill(r, D, tp, win)
        D[:2] = D[2]; D[-2:] = D[-3]
        return D, sp
    with np.errstate(divide='ignore', invalid='ignore'):
        hp = Qp/Q + 2*sig*sp/Gam
        spp = 0.5*hp*sp
        hpp = Qpp/Q - (Qp/Q)**2 + 2*sp*sp/Gam + 2*sig*spp/Gam + 4*sig*sig*sp*sp/(Gam*Gam)
        D = hp*hp/16 - hpp/4
    return smooth_fill(r, D, tp, win), sp


def u_iter(n, l, r, k, check=None, glob=False):
    """k = 0: 733 (Langer); k >= 1: k Schwarzian iterations.  Returns u and the phase checks int sqrt(Q_k)/pi(n_r+1/2)."""
    E = -0.5/n**2; nr = n - l - 1; L2 = (l + 0.5)**2
    QL = lambda x: 2*E + 2/x - L2/x**2
    QLp = lambda x: -2/x**2 + 2*L2/x**3
    QLpp = lambda x: 4/x**3 - 6*L2/x**4
    disc = math.sqrt(1 + 2*E*L2); tp = ((1 - disc)/(-2*E), (1 + disc)/(-2*E))
    Qf, Q, Qp, Qpp = QL, QL(r), QLp(r), QLpp(r)
    sig, Gam, cc, ph, tp = mapping(Qf, r, nr, tp, glob)
    phs = [ph]
    for it in range(k):
        D, sp = schwarz(sig, Gam, cc, r, Q, Qp, Qpp, tp)
        if check is not None and it == 0:
            f = sp**-0.5; mid = 0.5*(tp[0] + tp[1]); m = np.abs(r - mid) < 0.1*(tp[1] - tp[0])
            Dn = np.gradient(np.gradient(f, r), r)/f
            check.append((np.max(np.abs(Dn[m][5:-5] - D[m][5:-5])), np.max(np.abs(D[m])), -4e-8*np.interp(1e-4, r, D)))
        eta = CubicSpline(r, 0.25 + r*r*D)          # Q_{k+1} = Q_L + eta/r^2 (eta = 1/4 + r^2 D)
        Qf = (lambda QL_, et: (lambda x: QL_(x) + et(x)/x**2))(QL, eta)
        Q = QL(r) + (0.25 + r*r*D)/r**2; Qp = Qpp = None
        sig, Gam, cc, ph, tp = mapping(Qf, r, nr, tp, glob)
        phs.append(ph)
    return wave(sig, r, nr), phs


_cache = {}
def ebk(n, l, r, k, key, glob=False):
    kk = (n, l, k, key, glob)
    if kk not in _cache:
        u, ph = u_iter(n, l, r, k, glob=glob)
        if np.trapezoid(u*u_qm(n, l, r), r) < 0: u = -u
        _cache[kk] = (u, ph)
    return _cache[kk]


def radial(ni, li, nf, lf, k, fine=1, glob=False):
    r = grid(max(ni, nf), fine); key = (max(ni, nf), fine)
    q = np.trapezoid(u_qm(ni, li, r)*u_qm(nf, lf, r)*r, r)
    ui, pi_ = ebk(ni, li, r, k, key, glob); uf, pf = ebk(nf, lf, r, k, key, glob)
    return q, np.trapezoid(ui*uf*r, r), pi_, pf


def A(ni, li, nf, lf, R, mu=MU_H):
    om = 0.5/nf**2 - 0.5/ni**2
    return (4/3)*ALPHA**3*om**3*(max(li, lf)/(2*li + 1))*R*R*mu/TAU_AU


if __name__ == '__main__':
    t0 = time.time()
    P('Audit 737 (run 2, high-precision mapping): EBK wave v2 + one Schwarzian (hbar^2) iteration vs exact QM; E = -1/2n^2')
    P('variants of the phase scale at k = 1: "in" = 733 recipe (scale between the turning points only; preregistered), "gl" = one scale everywhere')
    P('(POST HOC; smooth sigma\'); at k = 0 both are 733 exactly (scale = 1).  Iterations k >= 2: numerically ill-conditioned here (see prereg addendum).')
    TR = [((2, 1), (1, 0)), ((3, 1), (1, 0)), ((3, 1), (2, 0)), ((3, 0), (2, 1)), ((3, 2), (2, 1)), ((4, 1), (1, 0)), ((4, 3), (3, 2))]
    r733 = [0.99624, 0.99551, 0.99778, 1.00524, 0.99892, 0.99524, 0.99948]
    nm = lambda i, f: '%d%s -> %d%s' % (i[0], 'spdfghiklmnoqrtuvwxyz'[i[1]], f[0], 'spdfghiklmnoqrtuvwxyz'[f[1]])

    P('\n(A) seven transitions: R/R_QM at k = 0 (must equal 733) and k = 1 (in, gl)')
    res = {}; okE1 = True
    for (i, f), o in zip(TR, r733):
        q, e0, _, _ = radial(*i, *f, 0); _, e1, pa, pb = radial(*i, *f, 1); _, eg, _, _ = radial(*i, *f, 1, glob=True)
        res[(i, f)] = (e0/q, e1/q, eg/q); okE1 &= abs(e0/q - o) < 1.5e-5
        P('  %s: k0 %.6f (733 %.5f)  k1-in %.6f  k1-gl %.6f   |.-1|: %.1e -> %.1e / %.1e' % (nm(i, f), e0/q, o, e1/q, eg/q, abs(e0/q - 1), abs(e1/q - 1), abs(eg/q - 1)))
    better = sum(abs(res[t][1] - 1) < abs(res[t][0] - 1) for t in TR); betg = sum(abs(res[t][2] - 1) < abs(res[t][0] - 1) for t in TR)
    okA = better >= 6 and abs(res[TR[0]][1] - 1) < 1e-3
    P('H737a (k=1-in better in >= 6/7 and 2p->1s |.-1| < 1e-3): %s (%d/7, 2p->1s %.2e); gl: %d/7, 2p->1s %.2e' % (okA, better, abs(res[TR[0]][1] - 1), betg, abs(res[TR[0]][2] - 1)))

    P('\n(D) quantization at order 1: int sqrt(Q_1) between its turning points / pi(n_r+1/2); closed-form D vs finite differences mid-orbit; -4 r^2 D at r = 1e-4')
    for (n, l) in [(1, 0), (2, 0), (2, 1), (3, 0), (3, 1), (3, 2), (4, 1), (4, 3), (8, 1), (8, 7)]:
        chk = []
        _, ph = u_iter(n, l, grid(n), 1, chk)
        P('  %d%s: phase k0 %.7f  k1 %.7f (shift %+.2e);  |D_fd - D| %.1e (max|D| %.1e);  -4r^2 D(1e-4) %.5f' % (n, 'spdfghiklmnoqrtuvwxyz'[l], ph[0], ph[1], ph[1] - 1, *chk[0]))

    P('\n(B) np lifetimes, all E1 channels down (n\'s, n\'d), n = 2..8')
    rat = {}
    for n in range(2, 9):
        Aq = 0.0; Ae = [0.0]*3
        for nf in range(1, n):
            for lf in (0, 2):
                if lf > nf - 1: continue
                for j, (k, g) in enumerate(((0, False), (1, False), (1, True))):
                    q, e, _, _ = radial(n, 1, nf, lf, k, glob=g)
                    if j == 0: Aq += A(n, 1, nf, lf, q)
                    Ae[j] += A(n, 1, nf, lf, e)
        rat[n] = [Aq/x for x in Ae]
        P('  %dp: tau QM %.5f ns;  k0 %.5f (x QM %.5f)  k1-in %.5f (%.5f)  k1-gl %.5f (%.5f)' % (n, 1e9/Aq, 1e9/Ae[0], rat[n][0], 1e9/Ae[1], rat[n][1], 1e9/Ae[2], rat[n][2]))
        if n == 2: tq, tk = 1e9/Aq, [1e9/x for x in Ae]
        _cache.clear()
    okB = all(-0.002 <= rat[n][1] - 1 <= 0.002 for n in range(2, 9)); okBg = all(-0.002 <= rat[n][2] - 1 <= 0.002 for n in range(2, 9))
    P('H737b (k=1-in: tau/QM - 1 in [-0.002, 0.002], n = 2..8): %s;  gl: %s' % (okB, okBg))
    P('  tau(2p) hydrogen vs measured 1.600(4) ns: QM %.5f (%+.2f sigma); k0 %.5f (%+.2f); k1-in %.5f (%+.2f); k1-gl %.5f (%+.2f)'
      % (tq, (tq - 1.6)/0.004, tk[0], (tk[0] - 1.6)/0.004, tk[1], (tk[1] - 1.6)/0.004, tk[2], (tk[2] - 1.6)/0.004))
    P('  Ps 2P (same element, a_Ps, omega/2; p-Ps and o-Ps identical): QM %.4f ns; k0 %.4f; k1-in %.4f; k1-gl %.4f' % (2*tq*MU_H, 2*tk[0]*MU_H, 2*tk[1]*MU_H, 2*tk[2]*MU_H))

    P('\n(C\') POST HOC replacement of H737c: order-1 gain vs n -- circular n -> n-1 and np -> (n-1)s, |R/R_QM - 1| at k0 and k1-in')
    for n in (5, 10, 20):
        for (i, f) in (((n, n - 1), (n - 1, n - 2)), ((n, 1), (n - 1, 0))):
            q, e0, _, _ = radial(*i, *f, 0); _, e1, _, _ = radial(*i, *f, 1)
            P('  %s: k0 %.2e  k1 %.2e' % (nm(i, f), abs(e0/q - 1), abs(e1/q - 1)))
        _cache.clear()

    P('\n(E) grid 2x finer, k = 1 (in, gl)')
    okE2 = True
    for (i, f) in TR[:3]:
        _cache.clear()
        q, e1, _, _ = radial(*i, *f, 1, fine=2); _, eg, _, _ = radial(*i, *f, 1, fine=2, glob=True)
        d = (abs(e1/q - res[(i, f)][1]), abs(eg/q - res[(i, f)][2])); okE2 &= max(d) < 1e-5
        P('  %s: fine in %.7f gl %.7f; diff %.1e %.1e' % (nm(i, f), e1/q, eg/q, *d))
    P('H737e (order 0 = 733 to 1e-5 and grid 2x diff at k=1 < 1e-5): %s (%s, %s)' % (okE1 and okE2, okE1, okE2))
    P('\nruntime %.0f s' % (time.time() - t0))
    open(os.path.join(HERE, 'out', 'ana737.txt'), 'w').write('\n'.join(out) + '\n')

"""Audit 733 (offline): EBK wave v2 -- uniform two-turning-point approximation (mapping onto the harmonic oscillator, Hermite functions), Langer,
E = -1/2n^2 -- E1 radial elements and np lifetimes vs exact QM, the engine (n-dependence) and measurement.  See prereg733.txt."""
import math, os
import numpy as np
from scipy.special import genlaguerre, gammaln
HERE = os.path.dirname(os.path.abspath(__file__))
out = []
def P(s): out.append(s); print(s, flush=True)
ALPHA = 1/137.035999177; TAU_AU = 2.4188843265864e-17; MU_H = 1/(1 + 1/1836.152673426)


def grid(nmax):
    Rmax = 3.0*nmax*nmax + 60
    return np.concatenate([np.geomspace(1e-7, 1e-2, 3000)[:-1], np.linspace(1e-2, Rmax, 1500001)])


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


def u_ebk2(n, l, r):
    E = -0.5/n**2; L2 = (l + 0.5)**2; nr = n - l - 1; s0 = math.sqrt(2*nr + 1)
    Q = 2*E + 2/r - L2/r**2
    disc = math.sqrt(1 + 2*E*L2); r1 = (1 - disc)/(-2*E); r2 = (1 + disc)/(-2*E)
    sq = np.sqrt(np.abs(Q)); C = np.concatenate([[0.0], np.cumsum(0.5*(sq[1:] + sq[:-1])*np.diff(r))])
    C1, C2 = np.interp(r1, r, C), np.interp(r2, r, C)
    Phi_tot = C2 - C1                                   # should be pi (n_r + 1/2)
    sig = np.empty_like(r)
    # tables
    st = np.linspace(-s0, s0, 200001); F = 0.5*(st*np.sqrt(np.maximum(s0*s0 - st*st, 0)) + s0*s0*np.arcsin(st/s0)) + math.pi*s0*s0/4
    sg = np.concatenate([[s0], s0 + np.geomspace(1e-9, 200.0, 200000)]); G = 0.5*(sg*np.sqrt(sg*sg - s0*s0) - s0*s0*np.arccosh(sg/s0))
    a = r < r1; b = (r >= r1) & (r <= r2); c = r > r2
    sig[b] = np.interp((C[b] - C1)*(math.pi*s0*s0/2)/Phi_tot, F, st)
    sig[a] = -np.interp(C1 - C[a], G, sg)
    sig[c] = np.interp(C[c] - C2, G, sg)
    ds = np.gradient(sig, r)
    with np.errstate(divide='ignore', invalid='ignore'):
        u = hermite_fn(nr, sig)/np.sqrt(np.abs(ds))
    u = np.nan_to_num(u, nan=0.0, posinf=0.0)
    u = u/math.sqrt(np.trapezoid(u*u, r))
    return u, Phi_tot/(math.pi*(nr + 0.5))


def radial(ni, li, nf, lf):
    r = grid(max(ni, nf))
    q = np.trapezoid(u_qm(ni, li, r)*u_qm(nf, lf, r)*r, r)
    ui, ci = u_ebk2(ni, li, r); uf, cf = u_ebk2(nf, lf, r)
    if np.trapezoid(ui*u_qm(ni, li, r), r) < 0: ui = -ui
    if np.trapezoid(uf*u_qm(nf, lf, r), r) < 0: uf = -uf
    e = np.trapezoid(ui*uf*r, r)
    return q, e, ci, cf


def A(ni, li, nf, lf, R, mu=MU_H):
    om = 0.5/nf**2 - 0.5/ni**2
    return (4/3)*ALPHA**3*om**3*(max(li, lf)/(2*li + 1))*R*R*mu/TAU_AU


P('Audit 733: EBK wave v2 (uniform two-turning-point, Hermite mapping, Langer, E = -1/2n^2) vs exact QM')
P('\n(A) the seven transitions of 732 (732 ratios: 0.9826 0.9683 1.0075 0.9474 0.9819 0.9587 0.9848)')
r732 = {((2, 1), (1, 0)): 0.98260, ((3, 1), (1, 0)): 0.96826, ((3, 1), (2, 0)): 1.00749, ((3, 0), (2, 1)): 0.94741, ((3, 2), (2, 1)): 0.98191,
        ((4, 1), (1, 0)): 0.95870, ((4, 3), (3, 2)): 0.98482}
better = 0
for (i, f), old in r732.items():
    q, e, ci, cf = radial(*i, *f)
    if e*q < 0: e = -e
    better += abs(e/q - 1) < abs(old - 1)
    P('  %d%s -> %d%s: R_QM %+.6f  R_v2 %+.6f  ratio %.5f  (732: %.5f)  [EBK phase check %.6f %.6f]' % (i[0], 'spdf'[i[1]], f[0], 'spdf'[f[1]], q, e, e/q, old, ci, cf))
    if (i, f) == ((2, 1), (1, 0)): q21, e21 = q, e
t_v2 = 1e9/A(2, 1, 1, 0, e21); t_qm = 1e9/A(2, 1, 1, 0, q21)
P('  tau(2p) hydrogen: v2 %.5f ns, QM %.5f, engine 1.6747, measured 1.600(4) -> v2 - measured %+.4f ns (%+.1f sigma)' % (t_v2, t_qm, t_v2 - 1.600, (t_v2 - 1.600)/0.004))
P('  Ps 2P (p-Ps and o-Ps): v2 %.4f ns, QM %.4f ns' % (2*1e9/A(2, 1, 1, 0, e21, 1.0), 2*1e9/A(2, 1, 1, 0, q21, 1.0)))
okA = better >= 6
P('H733a (v2 better than 732 in >= 6/7): %s (%d/7)' % (okA, better))

P('\n(B) np lifetimes, all E1 channels down (n\'s, n\'d), n = 2..8; engine ratio from Model.md / Artykul (612, 585, 638)')
eng = {2: 1.6747/1.5962, 3: 5.5848/5.2715, 4: 13.060/12.30, 5: 25.25/23.78, 6: 1.0615, 7: 1.0615, 8: 1.0615}
rat = {}
for n in range(2, 9):
    Aq = Ae = 0.0
    for nf in range(1, n):
        for lf in (0, 2):
            if lf > nf - 1: continue
            q, e, _, _ = radial(n, 1, nf, lf)
            if e*q < 0: e = -e
            Aq += A(n, 1, nf, lf, q); Ae += A(n, 1, nf, lf, e)
    rat[n] = Aq/Ae
    P('  %dp: tau QM %.4f ns, v2 %.4f ns, v2/QM %.5f; engine/QM %.4f%s' % (n, 1e9/Aq, 1e9/Ae, Aq/Ae, eng[n], ' (n = 6..8: 1.061-1.062, 612)' if n >= 6 else ''))
okB = abs(rat[8] - 1) < 0.5*abs(rat[2] - 1)
P('H733b (|v2/QM - 1| at n = 8 < half of n = 2): %s (%.5f vs %.5f)' % (okB, rat[8] - 1, rat[2] - 1))
okC = okB and all(abs(rat[n] - 1) < 0.02 for n in range(4, 9))
P('H733c (v2 within 2 %% of QM for n >= 4 while the engine stays at 1.06 -> 732d reading wrong): %s' % okC)

P('\n(D) large n: circular n -> n-1 (l = n-1 -> n-2) and np -> (n-1)s')
for n in (5, 10, 20, 30):
    q, e, _, _ = radial(n, n - 1, n - 1, n - 2)
    q2, e2, _, _ = radial(n, 1, n - 1, 0)
    if e*q < 0: e = -e
    if e2*q2 < 0: e2 = -e2
    P('  n %2d: circular ratio %.6f (|.-1| %.1e);  np -> (n-1)s ratio %.6f (|.-1| %.1e)' % (n, e/q, abs(e/q - 1), e2/q2, abs(e2/q2 - 1)))
open(os.path.join(HERE, 'out', 'ana733.txt'), 'w').write('\n'.join(out) + '\n')

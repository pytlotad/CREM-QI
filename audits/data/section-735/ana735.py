"""Audit 735 (offline): E2 rates from the EBK wave (r^2 element) and M1 rates from the spin wave (Schwinger lengths, midpoint projection) vs QM and
measurement.  See prereg735.txt."""
import math, os, sys
import numpy as np
from scipy.special import sph_harm_y
from sympy.physics.wigner import wigner_6j
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '../section-733'))
HERE = os.path.dirname(os.path.abspath(__file__))
import importlib.util
spec = importlib.util.spec_from_file_location('a733', os.path.join(HERE, '../section-733/ana733.py'))
src = open(os.path.join(HERE, '../section-733/ana733.py')).read().split("P('Audit 733")[0]
g733 = {'__file__': os.path.join(HERE, '../section-733/ana733.py')}; exec(src, g733)
grid, u_qm, u_ebk2 = g733['grid'], g733['u_qm'], g733['u_ebk2']
out = []
def P(s): out.append(s); print(s, flush=True)
ALPHA = 1/137.035999177; TAU_AU = 2.4188843265864e-17; ME_MP = 1/1836.152673426; MU_H = 1/(1 + ME_MP)
KAPPA2 = (1 - ME_MP)/(1 + ME_MP)       # |kappa_2|/e for p + e- (= (m_p - m_e)/M)

# ---- angular sums on the sphere
th = np.arccos(np.polynomial.legendre.leggauss(64)[0]); wth = np.polynomial.legendre.leggauss(64)[1]
ph = np.linspace(0, 2*np.pi, 128, endpoint=False); TH, PH = np.meshgrid(th, ph, indexing='ij')
W = np.outer(wth, np.full(128, 2*np.pi/128))
n = np.stack([np.sin(TH)*np.cos(PH), np.sin(TH)*np.sin(PH), np.cos(TH)])
def Y(l, m): return sph_harm_y(l, m, TH, PH)
def ang_E1(li, lf):
    s = 0.0
    for mf in range(-lf, lf + 1):
        for i in range(3): s += abs(np.sum(np.conj(Y(lf, mf))*n[i]*Y(li, 0)*W))**2
    return s
def ang_E2(li, lf):
    s = 0.0
    for mf in range(-lf, lf + 1):
        for i in range(3):
            for j in range(3):
                Q = 3*n[i]*n[j] - (1.0 if i == j else 0.0)
                s += abs(np.sum(np.conj(Y(lf, mf))*Q*Y(li, 0)*W))**2
    return s

P('Audit 735 (offline): E2 from the EBK wave, M1 from the spin wave')
P('\nformula check: E1 angular sum vs l_>/(2l_i + 1): ' + ', '.join('%d->%d %.10f/%.10f' % (a, b, ang_E1(a, b), max(a, b)/(2*a + 1)) for a, b in ((1, 0), (2, 1), (0, 1), (3, 2))))
# classical limit check of A_E2 = k w^5 sum|Q_fi|^2/(90 hbar c^5): circular orbit, large n, l = n-1 -> n-2: QM vs classical power of the rotating quadrupole
def rad(ni, li, nf, lf, power):
    r = grid(max(ni, nf))
    uq_i, uq_f = u_qm(ni, li, r), u_qm(nf, lf, r)
    ue_i, _ = u_ebk2(ni, li, r); ue_f, _ = u_ebk2(nf, lf, r)
    if np.trapezoid(ue_i*uq_i, r) < 0: ue_i = -ue_i
    if np.trapezoid(ue_f*uq_f, r) < 0: ue_f = -ue_f
    q = np.trapezoid(uq_i*uq_f*r**power, r); e = np.trapezoid(ue_i*ue_f*r**power, r)
    if e*q < 0: e = -e
    return q, e
def A_E2(ni, li, nf, lf, R2):    # s^-1, hydrogen, reduced mass, kappa_2
    om = 0.5/nf**2 - 0.5/ni**2
    return KAPPA2**2*ALPHA**5*om**5*ang_E2(li, lf)*R2*R2/90*MU_H/TAU_AU
def A_E1(ni, li, nf, lf, R):
    om = 0.5/nf**2 - 0.5/ni**2
    return (4/3)*ALPHA**3*om**3*(max(li, lf)/(2*li + 1))*R*R*MU_H/TAU_AU
# classical-limit check (corrected after run 1: charge quadrupole P = (8/5) k e^2 a^4 w0^6/c^5, not the gravitational 32/5; midpoint orbit n_m):
# circular n -> n-2 (E2, Delta l = -2) and n -> n-1 (E1), QM count / classical count, n = 10, 20, 40 -- must tend to 1 as O(1/n)
for nC in (10, 20, 40):
    q1, _ = rad(nC, nC - 1, nC - 1, nC - 2, 1); om1 = 0.5/(nC - 1)**2 - 0.5/nC**2
    A1 = (4/3)*ALPHA**3*om1**3*(nC - 1)/(2*nC - 1)*q1*q1; nm = nC - 0.5; A1cl = (2/3)*ALPHA**3*nm**4/nm**9
    q2, _ = rad(nC, nC - 1, nC - 2, nC - 3, 2); om2 = 0.5/(nC - 2)**2 - 0.5/nC**2
    A2 = ALPHA**5*om2**5*ang_E2(nC - 1, nC - 3)*q2*q2/90; nm2 = nC - 1; A2cl = (8/5)*ALPHA**5*nm2**8/nm2**18/(2/nm2**3)
    P('classical-limit check circular n = %d: E1 QM/classical %.4f, E2 QM/classical %.4f' % (nC, A1/A1cl, A2/A2cl))

P('\n(E2) hydrogen, R2 = int u_i u_f r^2 dr (a.u.), A_E2 with kappa_2 = %.6f e, reduced mass' % KAPPA2)
okA = True; e2tab = {}
for (ni, li, nf, lf) in ((3, 2, 1, 0), (4, 2, 1, 0), (4, 2, 2, 0), (4, 3, 2, 1), (5, 4, 3, 2), (4, 2, 3, 2)):
    q, e = rad(ni, li, nf, lf, 2); okA &= abs(e/q - 1) < 0.01
    aq, ae = A_E2(ni, li, nf, lf, q), A_E2(ni, li, nf, lf, e); e2tab[(ni, li, nf, lf)] = (aq, ae)
    P('  %d%s -> %d%s: R2_QM %+.5f  R2_EBK %+.5f  ratio %.5f;  A_E2 QM %.4e  EBK %.4e s^-1  (angular %.5f)' %
      (ni, 'spdfg'[li], nf, 'spdfg'[lf], q, e, e/q, aq, ae, ang_E2(li, lf)))
P('H735a (|R2_EBK/R2_QM - 1| < 1 %% for six transitions): %s' % okA)
q1, e1 = rad(3, 2, 2, 1, 1); A3d_E1 = A_E1(3, 2, 2, 1, q1)
P('H735b (report): 3d: A_E1(3d->2p) QM %.4e s^-1; A_E2(3d->1s) QM %.4e, EBK %.4e -> E2 share %.3e (QM), engine (396) 6.5e-6' %
  (A3d_E1, e2tab[(3, 2, 1, 0)][0], e2tab[(3, 2, 1, 0)][1], e2tab[(3, 2, 1, 0)][0]/A3d_E1))

# ---- M1
P('\n(M1) single spin: |<m-1|S_-|m>|^2 = S(S+1) - m(m-1) vs semiclassical (S + 1/2)^2 - (m - 1/2)^2')
okC = True
for S in (0.5, 1.0, 1.5):
    m = S; vals = []
    while m > -S + 1e-9:
        qm = S*(S + 1) - m*(m - 1); sc = (S + 0.5)**2 - (m - 0.5)**2; okC &= abs(qm - sc) < 1e-12; vals.append('m %+.1f: %.4f/%.4f' % (m, qm, sc)); m -= 1
    P('  S = %.1f: %s' % (S, '; '.join(vals)))
P('H735c (exact): %s' % okC)
mu0 = 4e-7*math.pi*(1 + 5.5e-10); hbar = 1.054571817646e-34; c = 299792458.0; muB = 9.2740100657e-24
def A_M1(omega, g, X): return mu0*omega**3/(3*math.pi*hbar*c**3)*(g*muB)**2*X
def X_qm(I, s, F, Ff):   # sum_{M'} |<Ff M'| s |F M>|^2
    return float((2*Ff + 1)*wigner_6j(s, Ff, I, F, s, 1)**2*s*(s + 1)*(2*s + 1))
def X_sc(I, s, F):       # 1/2 s_perp^2, triangle (|I| = I + 1/2, |s| = s + 1/2, F_m = F)
    Il, sl, Fm = I + 0.5, s + 0.5, F
    sF = (Fm*Fm + sl*sl - Il*Il)/(2*Fm); return 0.5*(sl*sl - sF*sF)
from sympy import Rational as Rr
# hydrogen 21 cm
om21 = 2*math.pi*1420.405751768e6
xq, xs = X_qm(Rr(1, 2), Rr(1, 2), 1, 0), X_sc(0.5, 0.5, 1)
okE = abs(xs/xq - 1) < 0.05
P('\n(M1) hydrogen 1s F = 1 -> 0 (21 cm): X_QM %.5f, X_spin-wave %.5f (ratio %.4f); A QM %.4e s^-1 (lit. 2.87e-15), spin wave %.4e' %
  (xq, xs, xs/xq, A_M1(om21, 2.00231930436092, xq), A_M1(om21, 2.00231930436092, xs)))
P('H735e (within 5 %%): %s' % okE)
# Bi 209
Z = 83; za = Z/137.035999177; g_bound = (2/3)*(1 + 2*math.sqrt(1 - za*za))
omBi = 2*math.pi*c/243.87e-9
xq, xs = X_qm(Rr(9, 2), Rr(1, 2), 5, 4), X_sc(4.5, 0.5, 5)
okD = abs(xs/xq - 1) < 0.05
tq, ts = 1e3/A_M1(omBi, g_bound, xq), 1e3/A_M1(omBi, g_bound, xs)
P('(M1) 209Bi82+ 1s F = 5 -> 4 (243.87 nm, 5.0840 eV): X_QM %.5f, X_spin-wave %.5f (ratio %.4f); bound-electron g (Dirac) %.4f;' % (xq, xs, xs/xq, g_bound))
P('     tau QM (this formula, no QED / Bohr-Weisskopf) %.4f ms, spin wave %.4f ms; measured 0.351(16) ms (Klaft 1994) -> QM %+.1f sigma, spin wave %+.1f sigma' %
  (tq, ts, (tq - 0.351)/0.016, (ts - 0.351)/0.016))
P('H735d (within 5 %% of QM): %s' % okD)
# Ps o-Ps -> p-Ps
xqPs = 1.0    # |<00| s1z - s2z |10>|^2, sum over M' (M = 0); same for every M of the triplet
Sm = 1.0; s1F = Sm/2; s1perp2 = 1.0 - s1F**2; xsPs = 0.5*(2*math.sqrt(s1perp2))**2
omPs = 2*math.pi*203.394e9
P('(M1) Ps 1^3S1 -> 1^1S0 (203.39 GHz, moment g muB (s1 - s2)): X_QM %.4f, spin wave (midpoint |S| = 1, |s| = 1) %.4f (ratio %.3f); A QM %.3e s^-1 vs annihilation 7.04e6' %
  (xqPs, xsPs, xsPs/xqPs, A_M1(omPs, 2.00231930436092, xqPs)))
open(os.path.join(HERE, 'out', 'ana735.txt'), 'w').write('\n'.join(out) + '\n')

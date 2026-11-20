# Audit 821: is the ordering term (S.L)/2 of the 521 projection rule the Weyl symmetrization in the Schwinger polarization?
# Prescription W: <O> = Weyl symbol of O averaged on the EBK torus with rule-Q lengths |L| = l+1/2, |S| = S+1/2, |J| = J+1/2.
# Usage: python3 ana821.py -> out/ana821.txt.  Written after the prereg.
import os, math, itertools
import numpy as np, sympy as sp
from scipy.special import sph_harm_y
pr = lambda *a: print(*a, flush=True); D = os.path.dirname(os.path.abspath(__file__))
# ---------------- Moyal product (hbar = 1, z <-> a, zb <-> a^+)
a, ab, b, bb, c, cb, d, db = sp.symbols('a ab b bb c cb d db')
PAIRS = [(a, ab), (b, bb), (c, cb), (d, db)]
def star(f, g):
    out = sp.Integer(0); terms = [(sp.expand(f), sp.expand(g), sp.Integer(1))]; n = 0
    while terms:
        out += sum(cf*fp*gp for fp, gp, cf in terms)/sp.factorial(n)*sp.Rational(1, 2)**n
        new = []
        for fp, gp, cf in terms:
            for z, zb in PAIRS:
                for s, (x, y) in ((1, (z, zb)), (-1, (zb, z))):
                    f1, g1 = sp.diff(fp, x), sp.diff(gp, y)
                    if f1 != 0 and g1 != 0: new.append((f1, g1, cf*s))
        terms = new; n += 1
    return sp.expand(out)
def vec(z1, z1b, z2, z2b):  # Schwinger: K = (1/2) v^+ sigma v
    return [sp.Rational(1, 2)*(z1b*z2 + z2b*z1), sp.Rational(1, 2)*(-sp.I*z1b*z2 + sp.I*z2b*z1), sp.Rational(1, 2)*(z1b*z1 - z2b*z2)]
Lv, Sv = vec(a, ab, b, bb), vec(c, cb, d, db)
L2c = sp.expand(sum(x**2 for x in Lv)); S2c = sp.expand(sum(x**2 for x in Sv)); SLc = sp.expand(sum(x*y for x, y in zip(Sv, Lv)))
# ---------------- H821a
chk = []
chk.append(sp.expand(star(a, ab) - star(ab, a)) == 1)
chk.append(sp.expand(star(ab, a) - (a*ab - sp.Rational(1, 2))) == 0)
symS = [[sp.expand((star(Sv[i], Sv[j]) + star(Sv[j], Sv[i]))/2) for j in range(3)] for i in range(3)]
chk.append(all(sp.expand(symS[i][j] - Sv[i]*Sv[j] + (sp.Rational(1, 8) if i == j else 0)) == 0 for i in range(3) for j in range(3)))
chk.append(sp.expand(sum(star(x, x) for x in Sv) - (S2c - sp.Rational(3, 8))) == 0)                     # P3
SLsq = star(SLc, SLc)
chk.append(sp.expand(SLsq - (SLc**2 - SLc/2 - (S2c + L2c)/8 + sp.Rational(3, 64))) == 0)
symL = [[sp.expand((star(Lv[i], Lv[j]) + star(Lv[j], Lv[i]))/2) for j in range(3)] for i in range(3)]
Osym = sp.expand(sum(star(symS[i][j], symL[i][j]) for i in range(3) for j in range(3)))
chk.append(sp.expand(Osym - (SLc**2 - (S2c + L2c)/8 + sp.Rational(3, 64))) == 0)                         # P1 symbol
okA = all(chk)
pr(f"H821a (algebra Moyala = symbole z przewidywania): {'TAK' if okA else 'NIE'}  {chk}")
# ---------------- QM: spin and orbit matrices, r^ r^ by quadrature on the sphere
def jmats(j):
    m = np.arange(j, -j - 1, -1); n = len(m); Jz = np.diag(m).astype(complex); Jp = np.zeros((n, n), complex)
    for k in range(1, n): Jp[k - 1, k] = math.sqrt(j*(j + 1) - m[k]*(m[k] + 1))
    return [(Jp + Jp.T.conj())/2, (Jp - Jp.T.conj())/(2j), Jz]
def rr(l):  # <l m| r^_i r^_j |l m'> by Gauss-Legendre x uniform phi quadrature (exact for these polynomials)
    xg, wg = np.polynomial.legendre.leggauss(4*l + 8); ph = np.linspace(0, 2*np.pi, 8*l + 16, endpoint=False); wp = 2*np.pi/len(ph)
    th = np.arccos(xg); TH, PH = np.meshgrid(th, ph, indexing='ij'); W = np.outer(wg, np.full(len(ph), wp))
    ms = np.arange(l, -l - 1, -1); Y = np.array([sph_harm_y(l, m, TH, PH) for m in ms])
    n = [np.sin(TH)*np.cos(PH), np.sin(TH)*np.sin(PH), np.cos(TH)]
    return [[np.einsum('agh,bgh,gh->ab', Y.conj(), Y, n[i]*n[j]*W) for j in range(3)] for i in range(3)]
def qm(l, S):  # {J: <(S.r^)^2>} from exact matrices
    Sm, Lm, R = jmats(S), jmats(l), rr(l); nS, nL = Sm[0].shape[0], Lm[0].shape[0]; I_S, I_L = np.eye(nS), np.eye(nL)
    O = sum(np.kron(Sm[i] @ Sm[j], R[i][j]) for i in range(3) for j in range(3))
    J2 = sum((np.kron(Sm[i], I_L) + np.kron(I_S, Lm[i])) @ (np.kron(Sm[i], I_L) + np.kron(I_S, Lm[i])) for i in range(3))
    w, V = np.linalg.eigh(J2); out = {}
    for Jv in np.unique(np.round((-1 + np.sqrt(1 + 4*w))/2*2)/2):
        P = V[:, np.abs((-1 + np.sqrt(1 + 4*w))/2 - Jv) < 1e-6]; out[float(Jv)] = float(np.real(np.trace(P.conj().T @ O @ P))/P.shape[1])
    return out
def rule521(l, S, J):
    x = (J*(J + 1) - l*(l + 1) - S*(S + 1))/2; Ll = l + 0.5
    return S*(S + 1)/3 - (x*(x + 0.5) - S*(S + 1)*l*(l + 1)/3)/(2*(Ll - 1)*(Ll + 1))
def classical(l, S, J):
    x = (J*(J + 1) - l*(l + 1) - S*(S + 1))/2
    return S*(S + 1)/3 - (x*x - S*(S + 1)*l*(l + 1)/3)/(2*l*(l + 1))
def W(l, S, J):  # Weyl symbol of sym(S_iS_j) (S_iS_j - d_ij/8) times orbit average (d_ij - L^_iL^_j)/2 on the rule-Q torus
    Sl, Ll, Jl = S + 0.5, l + 0.5, J + 0.5; xc = (Jl**2 - Ll**2 - Sl**2)/2
    return (Sl**2 - xc**2/Ll**2 - 0.25)/2
def W_O(l, S, J):
    x = (J*(J + 1) - l*(l + 1) - S*(S + 1))/2; xc = ((J + .5)**2 - (l + .5)**2 - (S + .5)**2)/2
    return xc**2 - ((S + .5)**2 + (l + .5)**2)/8 + 3/64, x*x + x/2, x
# ---------------- H821q, H821w, H821o
okQ, okW, okO, maxW = True, True, True, 0.0
pr('  o-Ps (S = 1): <S12> = 6 <(S.r^)^2> - 4;  kolumny QM / 521 / W / klasyczny (Born, bez zamian)')
for l in range(1, 7):
    q = qm(l, 1)
    for J in (l - 1, l, l + 1):
        Q, R5, Wv, Cl = q[float(J)], rule521(l, 1, J), W(l, 1, J), classical(l, 1, J)
        okQ &= abs(Q - R5) <= 1e-12; e = abs(Wv - Q); maxW = max(maxW, e); okW &= e <= 1e-12
        pr(f'    l = {l} J = {J}: <S12> {6*Q - 4:+.6f} / {6*R5 - 4:+.6f} / {6*Wv - 4:+.6f} / {6*Cl - 4:+.6f}   (W - QM w <S12>: {6*(Wv - Q):+.4f})')
        wo, qo, x = W_O(l, 1, J); okO &= abs((wo - qo) - (-3*x/4 - (l*(l + 1) + 2)/8)) <= 1e-12
pr('  wodor (S = 1/2): (S.r^)^2 elektronu, F nie wchodzi; kolumny QM / 521 / W')
for l in range(1, 7):
    q = qm(l, 0.5)
    for j in (l - 0.5, l + 0.5):
        Q, R5, Wv = q[j], rule521(l, 0.5, j), W(l, 0.5, j); okQ &= abs(Q - R5) <= 1e-12; e = abs(Wv - Q); maxW = max(maxW, e); okW &= e <= 1e-12
        for F in (j - 0.5, j + 0.5):
            pr(f'    l = {l} j = {j} F = {F:g} ({"PHEN=1" if F == 0 else "PHEN=2" if F == 1 else "F > 1"}): {Q:.6f} / {R5:.6f} / {Wv:.6f}')
        wo, qo, x = W_O(l, 0.5, j); okO &= abs((wo - qo) - (-3*x/4 - (l*(l + 1) + 0.75)/8)) <= 1e-12
pr(f"H821q (QM z kwadratury = regula 521 do 1e-12, S = 1 i 1/2, l = 1..6): {'TAK' if okQ else 'NIE'}")
pr(f"H821w (GLOWNA: przepis W = QM do 1e-12): {'TAK' if okW else 'NIE'}  (max |W - QM| w <(S.r^)^2> = {maxW:.4f})")
pr(f"H821o (W(O) - QM = -3x/4 - [l(l+1) + S(S+1)]/8 dokladnie): {'TAK' if okO else 'NIE'}")
pr('  RAPORT p-Ps (S = 0, singlet): QM (S.r^)^2 = 0; W: ' + ', '.join(f'l = {l}: {W(l, 0, l):+.6f}' for l in range(1, 7)))
# ---------------- RAPORT: multiplet 2^3P o-Ps with <S12> from W
import re
E = {m.group(1): float(m.group(2)) for m in re.finditer(r'(2\^\dP\d|2\^3S1) (\S+)', open(os.path.join(D, '..', 'section-820', 'out', 'p640_820.txt')).read())}
X = {0: -2, 1: -1, 2: 1}; Tq = {J: 6*rule521(1, 1, J) - 4 for J in (0, 1, 2)}
A = np.array([[1, X[J], Tq[J]] for J in (0, 1, 2)], float); e0, aa, bb_ = np.linalg.solve(A, [E[f'2^3P{J}'] for J in (0, 1, 2)])
pr(f'  RAPORT 2^3P silnika: E_c {e0:.4f}, a {aa:.4f}, b {bb_:.4f} MHz (z poziomow p640, 820)')
HEV, C = 4.135667696e-15, 299792458.0; MEAS = {0: (18501.02, 0.61, 'Gurung 2020'), 1: (13012.42, 1.67, 'Hagena 1993'), 2: (8627.94, 0.95, 'Sheldon 2023')}
for J in (0, 1, 2):
    nu = E['2^3S1'] - E[f'2^3P{J}']; nuW = nu - bb_*((6*W(1, 1, J) - 4) - Tq[J]); m, s, src = MEAS[J]
    pr(f'    nu{J}: silnik {nu:.2f} MHz; z W {nuW:.2f} MHz = {nuW*1e6*HEV*1e6:.4f} ueV = {C/(nuW*1e6)*100:.4f} cm; pomiar {m}({s}) [{src}]: '
       f'silnik {(nu - m)/s:+.1f} sigma, W {(nuW - m)/s:+.1f} sigma')

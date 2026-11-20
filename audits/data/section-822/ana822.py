# Audit 822: Born product over neighbouring transitions l -> l +- 1 for the spin components: <(S.r^)^2> = sum_l' A^2(l -> l').
# Usage: python3 ana822.py -> out/ana822.txt.  Written after the prereg.
import os, re, math
import numpy as np
from scipy.special import sph_harm_y
pr = lambda *a: print(*a, flush=True); D = os.path.dirname(os.path.abspath(__file__))
def jmats(j):
    m = np.arange(j, -j - 1, -1); n = len(m); Jz = np.diag(m).astype(complex); Jp = np.zeros((n, n), complex)
    for k in range(1, n): Jp[k - 1, k] = math.sqrt(j*(j + 1) - m[k]*(m[k] + 1))
    return [(Jp + Jp.T.conj())/2, (Jp - Jp.T.conj())/(2j), Jz]
def rmat(lp, l):  # <lp m'| r^_i | l m> by quadrature
    N = 2*(l + lp) + 8; xg, wg = np.polynomial.legendre.leggauss(N); ph = np.linspace(0, 2*np.pi, 2*N, endpoint=False)
    TH, PH = np.meshgrid(np.arccos(xg), ph, indexing='ij'); W = np.outer(wg, np.full(len(ph), 2*np.pi/len(ph)))
    Yp = np.array([sph_harm_y(lp, m, TH, PH) for m in np.arange(lp, -lp - 1, -1)]); Y = np.array([sph_harm_y(l, m, TH, PH) for m in np.arange(l, -l - 1, -1)])
    n = [np.sin(TH)*np.cos(PH), np.sin(TH)*np.sin(PH), np.cos(TH)]
    return [np.einsum('agh,bgh,gh->ab', Yp.conj(), Y, n[i]*W) for i in range(3)]
def projJ(S, l):
    Sm, Lm = jmats(S), jmats(l); IS, IL = np.eye(Sm[0].shape[0]), np.eye(Lm[0].shape[0])
    J2 = sum((np.kron(Sm[i], IL) + np.kron(IS, Lm[i])) @ (np.kron(Sm[i], IL) + np.kron(IS, Lm[i])) for i in range(3))
    w, V = np.linalg.eigh(J2); Jv = (-1 + np.sqrt(1 + 4*w))/2
    return {round(float(J)*2)/2: V[:, np.abs(Jv - J) < 1e-6] for J in np.unique(np.round(Jv*2)/2)}
def qm_amp(S, l, lp, J):  # (1/(2J+1)) sum_M |<lp S J M| S.r^ |l S J M>|^2
    Sm = jmats(S); R = rmat(lp, l); T = sum(np.kron(Sm[i], R[i]) for i in range(3))
    P, Pp = projJ(S, l), projJ(S, lp)
    if J not in Pp: return 0.0
    B = Pp[J].conj().T @ T @ P[J]; return float(np.real(np.trace(B.conj().T @ B)))/P[J].shape[1]
def allowed(S, lp, J): return lp >= 0 and abs(lp - S) <= J <= lp + S
def rule521(l, S, J):
    x = (J*(J + 1) - l*(l + 1) - S*(S + 1))/2; Ll = l + 0.5
    return S*(S + 1)/3 - (x*(x + 0.5) - S*(S + 1)*l*(l + 1)/3)/(2*(Ll - 1)*(Ll + 1))
def amp(rule, S, l, lp, J):
    lg = max(l, lp)
    if rule == 'B1':
        S2, L2, Lam2 = (S + .5)**2, lg**2, (J + .5)**2; SL = (Lam2 - L2 - S2)/2; return (S2 - SL**2/L2)/4
    if rule == 'B3':
        S2, L2 = S*(S + 1), (l + .5)*(lp + .5); x = lambda q: (J*(J + 1) - q*(q + 1) - S*(S + 1))/2; return (S2 - x(l)*x(lp)/L2)/4
    if rule == 'B5':
        S2, L2, Lam2 = (S + .5)**2, (lg - 1)*(lg + 1), (J + .5)**2; SL = (Lam2 - L2 - S2)/2
        return (S2 - (SL**2/L2 if L2 != 0 else (0.0 if abs(SL) < 1e-12 else float('inf'))))/4
def total(rule, S, l, J): return sum(amp(rule, S, l, lp, J) for lp in (l - 1, l + 1) if allowed(S, lp, J))
okQ = okH = okP = True; okB = {'B1': True, 'B3': True}
pr('  wodor (S = 1/2): (S.r^)^2 -- QM / 521 / B5 / B1 / B3;  amplitudy B5 vs QM po l\'')
for l in range(1, 7):
    for j in (l - .5, l + .5):
        A = {lp: qm_amp(.5, l, lp, j) for lp in (l - 1, l + 1) if allowed(.5, lp, j)}; Q = sum(A.values())
        okQ &= abs(Q - rule521(l, .5, j)) <= 1e-12; v5 = total('B5', .5, l, j); okH &= abs(v5 - Q) <= 1e-12
        for r in okB: okB[r] &= abs(total(r, .5, l, j) - Q) <= 1e-12
        for F in (j - .5, j + .5):
            pr(f'    l = {l} j = {j} F = {F:g} ({"PHEN=1" if F == 0 else "PHEN=2" if F == 1 else "F > 1"}): {Q:.6f} / {rule521(l, .5, j):.6f} / {v5:.6f} / '
               f'{total("B1", .5, l, j):.6f} / {total("B3", .5, l, j):.6f}   ' + ' '.join(f"l'={lp}: {amp('B5', .5, l, lp, j):.4f} vs {a:.4f}" for lp, a in A.items()))
pr('  o-Ps (S = 1): <S12> = 6<(S.r^)^2> - 4 -- QM / 521 / B5 / B1 / B3;  amplitudy (S.r^)^2 B5 vs QM po l\'')
for l in range(1, 7):
    for J in (l - 1, l, l + 1):
        A = {lp: qm_amp(1, l, lp, J) for lp in (l - 1, l + 1) if allowed(1, lp, J)}; Q = sum(A.values())
        okQ &= abs(Q - rule521(l, 1, J)) <= 1e-12; v5 = total('B5', 1, l, J); okP &= abs(v5 - Q) <= 1e-12
        for r in okB: okB[r] &= abs(total(r, 1, l, J) - Q) <= 1e-12
        pr(f'    l = {l} J = {J}: {6*Q - 4:+.6f} / {6*rule521(l, 1, J) - 4:+.6f} / {6*v5 - 4:+.6f} / {6*total("B1", 1, l, J) - 4:+.6f} / {6*total("B3", 1, l, J) - 4:+.6f}   '
           + ' '.join(f"l'={lp}: {amp('B5', 1, l, lp, J):.5f} vs {a:.5f}" for lp, a in A.items()))
pr(f"H822q (QM suma po l' = 521 do 1e-12): {'TAK' if okQ else 'NIE'}")
pr(f"H822h (wodor: B5 = QM, F = 0 i F = 1): {'TAK' if okH else 'NIE'}")
pr(f"H822p (GLOWNA, o-Ps: B5 = QM): {'TAK' if okP else 'NIE'}")
pr(f"H822b (B1 lub B3 = QM): {'TAK' if any(okB.values()) else 'NIE'}  {okB}")
pr('  RAPORT p-Ps (S = 0): QM 0; ' + ', '.join(f'l = {l}: B5 {total("B5", 0, l, l):+.5f}, B1 {total("B1", 0, l, l):+.5f}, B3 {total("B3", 0, l, l):+.5f}' for l in (1, 2, 3)))
E = {m.group(1): float(m.group(2)) for m in re.finditer(r'(2\^\dP\d|2\^3S1) (\S+)', open(os.path.join(D, '..', 'section-820', 'out', 'p640_820.txt')).read())}
X = {0: -2, 1: -1, 2: 1}; Tq = {J: 6*rule521(1, 1, J) - 4 for J in (0, 1, 2)}
e0, aa, bb = np.linalg.solve(np.array([[1, X[J], Tq[J]] for J in (0, 1, 2)], float), [E[f'2^3P{J}'] for J in (0, 1, 2)])
HEV, C = 4.135667696e-15, 299792458.0; MEAS = {0: (18501.02, 0.61, 'Gurung 2020'), 1: (13012.42, 1.67, 'Hagena 1993'), 2: (8627.94, 0.95, 'Sheldon 2023')}
for J in (0, 1, 2):
    nu = E['2^3S1'] - E[f'2^3P{J}']; nB = nu - bb*((6*total('B5', 1, 1, J) - 4) - Tq[J]); m, s, src = MEAS[J]
    pr(f'  RAPORT nu{J}: silnik {nu:.2f}, z B5 {nB:.2f} MHz = {nB*1e6*HEV*1e6:.4f} ueV = {C/(nB*1e6)*100:.4f} cm; pomiar {m}({s}) [{src}]: silnik {(nu - m)/s:+.1f} sigma, B5 {(nB - m)/s:+.1f} sigma')

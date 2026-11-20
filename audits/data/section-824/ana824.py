# Audit 824: triplet <-> singlet amplitude from the singlet torus (rules S1-S4).  Usage: python3 ana824.py -> out/ana824.txt.
import os, re, math
import numpy as np
from scipy.special import sph_harm_y
pr = lambda *a: print(*a, flush=True); D = os.path.dirname(os.path.abspath(__file__))
def jmats(j):
    m = np.arange(j, -j - 1, -1); n = len(m); Jz = np.diag(m).astype(complex); Jp = np.zeros((n, n), complex)
    for k in range(1, n): Jp[k - 1, k] = math.sqrt(j*(j + 1) - m[k]*(m[k] + 1))
    return [(Jp + Jp.T.conj())/2, (Jp - Jp.T.conj())/(2j), Jz]
def rmat(lp, l):
    N = 2*(l + lp) + 8; xg, wg = np.polynomial.legendre.leggauss(N); ph = np.linspace(0, 2*np.pi, 2*N, endpoint=False)
    TH, PH = np.meshgrid(np.arccos(xg), ph, indexing='ij'); W = np.outer(wg, np.full(len(ph), 2*np.pi/len(ph)))
    Yp = np.array([sph_harm_y(lp, m, TH, PH) for m in np.arange(lp, -lp - 1, -1)]); Y = np.array([sph_harm_y(l, m, TH, PH) for m in np.arange(l, -l - 1, -1)])
    n = [np.sin(TH)*np.cos(PH), np.sin(TH)*np.sin(PH), np.cos(TH)]
    return [np.einsum('agh,bgh,gh->ab', Yp.conj(), Y, n[i]*W) for i in range(3)]
h = jmats(0.5); I2 = np.eye(2); s1 = [np.kron(x, I2) for x in h]; s2 = [np.kron(I2, x) for x in h]
Sv = [a + b for a, b in zip(s1, s2)]; uv = [(a - b)/2 for a, b in zip(s1, s2)]
def proj(l):
    Lm = jmats(l); IL = np.eye(Lm[0].shape[0]); I4 = np.eye(4); Jt = [np.kron(Sv[i], IL) + np.kron(I4, Lm[i]) for i in range(3)]
    J2 = sum(x @ x for x in Jt); S2 = sum(np.kron(Sv[i] @ Sv[i], IL) for i in range(3)); w, V = np.linalg.eigh(J2 + 0.01*S2); out = {}
    for k in range(len(w)):
        v = V[:, k]; Sq = float(np.real(v.conj() @ S2 @ v)); Jq = float(np.real(v.conj() @ J2 @ v))
        out.setdefault((round((-1 + math.sqrt(1 + 4*Sq))/2), round((-1 + math.sqrt(1 + 4*Jq))/2)), []).append(v)
    return {k: np.array(v).T for k, v in out.items()}
PJ = {l: proj(l) for l in range(0, 8)}
def A2qm(l, S, J, lp, Sp):
    if (Sp, J) not in PJ[lp] or (S, J) not in PJ[l]: return 0.0
    R = rmat(lp, l); U = sum(np.kron(uv[i], R[i]) for i in range(3)); B = PJ[lp][(Sp, J)].conj().T @ U @ PJ[l][(S, J)]
    return float(np.real(np.trace(B.conj().T @ B)))/PJ[l][(S, J)].shape[1]
def rule(name, lp, up):   # singlet l' -> triplet l = l' + up (up = +1 / -1), J = l'
    Lp = lp + .5; c = -1/(4*Lp); sgn = -1 if up == 1 else 1   # (Delta l = +1) ~ (1 - c), (-1) ~ (1 + c)
    if name == 'S1': return (15/16)*(1 + sgn*c)**2/16
    if name == 'S2': return (15/16)*(1 + sgn*c)**2/8
    if name == 'S3': return 3/32
    if name == 'S4': return (1/8)*(1 + 2*sgn*c)
def exists(l, S, J): return l >= 0 and abs(l - S) <= J <= l + S
RULES = ('S1', 'S2', 'S3', 'S4'); ok = {r: True for r in RULES}; okQ = okD = True
pr("  singlet p-Ps (l', 0, J = l') -> tryplet (l' +- 1, 1, l'): A^2 QM / docelowe l_>/(8L') / S1 / S2 / S3 / S4")
for lp in range(1, 7):
    tot = 0.0
    for up in (1, -1):
        l = lp + up
        if not exists(l, 1, lp): continue
        q = A2qm(lp, 0, lp, l, 1); tot += q; tgt = max(l, lp)/(8*(lp + .5)); okQ &= abs(q - tgt) <= 1e-12
        v = {r: rule(r, lp, up) for r in RULES}
        for r in RULES: ok[r] &= abs(v[r] - q) <= 1e-12
        c = -1/(4*(lp + .5)); sgn = -1 if up == 1 else 1
        okD &= abs((v['S2'] - v['S4']) - ((15/16)/8*c*c + (15/16 - 1)/8*(1 + 2*sgn*c))) <= 1e-12
        pr(f"    l' = {lp} -> l = {l}: {q:.6f} / {tgt:.6f} / " + ' / '.join(f'{v[r]:.6f}' for r in RULES))
    okQ &= abs(tot - .25) <= 1e-12
    pr(f"    l' = {lp}: suma QM {tot:.6f}; S1 {sum(rule('S1', lp, u) for u in (1, -1)):.6f}, S2 {sum(rule('S2', lp, u) for u in (1, -1)):.6f}, S3 {sum(rule('S3', lp, u) for u in (1, -1)):.6f}, S4 {sum(rule('S4', lp, u) for u in (1, -1)):.6f}")
pr(f"H824q (QM = l_>/(8L'), suma 1/4): {'TAK' if okQ else 'NIE'}")
pr(f"H824s (GLOWNA: S1, S2 lub S3 = QM): {'TAK' if any(ok[r] for r in ('S1', 'S2', 'S3')) else 'NIE'}  {dict((r, ok[r]) for r in ('S1', 'S2', 'S3'))}")
pr(f"H824r (S4 = QM, tozsamosc): {'TAK' if ok['S4'] else 'NIE'}")
pr(f"H824d (S2 - S4 = (|u|^2/8) c^2 + (|u|^2 - 1)/8 (1 -+ 2c)): {'TAK' if okD else 'NIE'}")
def S12(r, l, J):   # triplet (l, 1, J): (S.r^)^2 = 1 - 4 A^2(triplet -> singlet l' = J) = 1 - 4 A^2(singlet J -> triplet l)
    if J == l or J < 0: return 2.0
    A2 = rule(r, J, l - J); return 6*(1 - 4*A2) - 4
pr('  RAPORT o-Ps <S12> (J = l - 1 / l / l + 1): QM / S1 / S2 / S3 / S4')
QMv = lambda l, J: {l - 1: -2*(l + 1)/(2*l - 1), l: 2.0, l + 1: -2*l/(2*l + 3)}[J]
for l in range(1, 7):
    pr(f'    l = {l}: ' + '  '.join(f'J={J}: {QMv(l, J):+.4f} / ' + ' / '.join(f'{S12(r, l, J):+.4f}' for r in RULES) for J in (l - 1, l, l + 1)))
E = {m.group(1): float(m.group(2)) for m in re.finditer(r'(2\^\dP\d|2\^3S1) (\S+)', open(os.path.join(D, '..', 'section-820', 'out', 'p640_820.txt')).read())}
X = {0: -2, 1: -1, 2: 1}; Tq = {J: QMv(1, J) for J in (0, 1, 2)}
e0, aa, bb = np.linalg.solve(np.array([[1, X[J], Tq[J]] for J in (0, 1, 2)], float), [E[f'2^3P{J}'] for J in (0, 1, 2)])
HEV, C = 4.135667696e-15, 299792458.0; MEAS = {0: (18501.02, 0.61, 'Gurung 2020'), 1: (13012.42, 1.67, 'Hagena 1993'), 2: (8627.94, 0.95, 'Sheldon 2023')}
for r in RULES:
    for J in (0, 1, 2):
        nu = E['2^3S1'] - E[f'2^3P{J}']; nB = nu - bb*(S12(r, 1, J) - Tq[J]); m, s, src = MEAS[J]
        pr(f'  RAPORT {r} nu{J}: {nB:.2f} MHz = {nB*1e6*HEV*1e6:.4f} ueV = {C/(nB*1e6)*100:.4f} cm; pomiar {m}({s}) [{src}]: {(nB - m)/s:+.1f} sigma (silnik {(nu - m)/s:+.1f})')

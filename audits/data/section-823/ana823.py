# Audit 823: triplet as two spins 1/2 -- (S.r^)^2 = 1 - 4 (u.r^)^2, u = (s1 - s2)/2; Born transitions triplet <-> singlet for u.r^.
# Usage: python3 ana823.py -> out/ana823.txt.  Written after the prereg.
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
h = jmats(0.5); I2 = np.eye(2)
s1 = [np.kron(x, I2) for x in h]; s2 = [np.kron(I2, x) for x in h]; Sv = [a + b for a, b in zip(s1, s2)]; uv = [(a - b)/2 for a, b in zip(s1, s2)]
def proj(l):  # {(S, J): columns} in spin(4) x orbit(2l+1)
    Lm = jmats(l); IL = np.eye(Lm[0].shape[0]); I4 = np.eye(4)
    Jt = [np.kron(Sv[i], IL) + np.kron(I4, Lm[i]) for i in range(3)]; J2 = sum(x @ x for x in Jt); S2 = sum(np.kron(Sv[i] @ Sv[i], IL) for i in range(3))
    w, V = np.linalg.eigh(J2 + 0.01*S2 + 0*1j); out = {}
    for k in range(len(w)):
        v = V[:, k]; Sq = float(np.real(v.conj() @ S2 @ v)); Jq = float(np.real(v.conj() @ J2 @ v))
        S = round((-1 + math.sqrt(1 + 4*Sq))/2); J = round((-1 + math.sqrt(1 + 4*Jq))/2); out.setdefault((S, J), []).append(v)
    return {k: np.array(v).T for k, v in out.items()}
PJ = {l: proj(l) for l in range(0, 8)}
def A2qm(l, S, J, lp, Sp):
    if (Sp, J) not in PJ[lp] or (S, J) not in PJ[l]: return 0.0
    R = rmat(lp, l); U = sum(np.kron(uv[i], R[i]) for i in range(3)); B = PJ[lp][(Sp, J)].conj().T @ U @ PJ[l][(S, J)]
    return float(np.real(np.trace(B.conj().T @ B)))/PJ[l][(S, J)].shape[1]
def Sr2qm(l, S, J):
    P = PJ[l][(S, J)]; tot = 0.0   # (S.r^)^2 within l = sum over intermediate l' = l +- 1 of |(S.r^)|^2
    for lp in (l - 1, l + 1):
        if lp < 0: continue
        R1 = rmat(lp, l); T = sum(np.kron(Sv[i], R1[i]) for i in range(3)); B = T @ P; tot += float(np.real(np.trace(B.conj().T @ B)))/P.shape[1]
    return tot
def ur2qm(l, S, J):
    P = PJ[l][(S, J)]; tot = 0.0
    for lp in (l - 1, l + 1):
        if lp < 0: continue
        R1 = rmat(lp, l); U = sum(np.kron(uv[i], R1[i]) for i in range(3)); B = U @ P; tot += float(np.real(np.trace(B.conj().T @ B)))/P.shape[1]
    return tot
def rule521(l, S, J):
    x = (J*(J + 1) - l*(l + 1) - S*(S + 1))/2; Ll = l + 0.5
    return S*(S + 1)/3 - (x*(x + 0.5) - S*(S + 1)*l*(l + 1)/3)/(2*(Ll - 1)*(Ll + 1))
def exists(l, S, J): return l >= 0 and S in (0, 1) and abs(l - S) <= J <= l + S
def A2cl(var, l, S, J, lp, Sp):
    if not exists(lp, Sp, J): return 0.0
    dl, dS = lp - l, Sp - S; eta = 1 if (dl, dS) in ((1, -1), (-1, 1)) else -1; lg = max(l, lp)
    if var == 'V1': s2, Sb2, L2, Lam2 = 1.0, 1.0, float(lg*lg), (J + .5)**2
    else: s2, Sb2, L2, Lam2 = .75, .75, float((lg - 1)*(lg + 1)), (J + .5)**2
    u2 = s2 - Sb2/4; SL = (Lam2 - L2 - Sb2)/2
    if L2 == 0: return float('nan')
    c = SL/math.sqrt(Sb2*L2); return u2*(1 + eta*c)**2/16
def ur2cl(var, l, S, J): return sum(A2cl(var, l, S, J, lp, 1 - S) for lp in (l - 1, l + 1))
okI = okS = okA = True; okV = {'V1': True, 'V2': True}
pr("  tryplet o-Ps (S = 1): <S12> QM / 521 / V1 / V2;  (u.r^)^2 QM / V1 / V2;  amplitudy QM vs docelowe l_>/(8(l'+1/2))")
for l in range(1, 7):
    for J in (l - 1, l, l + 1):
        Q = Sr2qm(l, 1, J); U2 = ur2qm(l, 1, J); okI &= abs(Q - (1 - 4*U2)) <= 1e-12 and abs(Q - rule521(l, 1, J)) <= 1e-12
        amps = {lp: A2qm(l, 1, J, lp, 0) for lp in (l - 1, l + 1) if lp >= 0}
        okS &= all((a > 1e-12) == (lp == J) for lp, a in amps.items()) and (J != l or U2 <= 1e-12)
        for lp, a in amps.items():
            if lp == J: okA &= abs(a - max(l, lp)/(8*(lp + .5))) <= 1e-12
        v = {k: ur2cl(k, l, 1, J) for k in okV}
        for k in okV: okV[k] &= np.isfinite(v[k]) and abs(v[k] - U2) <= 1e-12
        pr(f'    l = {l} J = {J}: <S12> {6*Q - 4:+.6f} / {6*rule521(l, 1, J) - 4:+.6f} / {6*(1 - 4*v["V1"]) - 4:+.6f} / {6*(1 - 4*v["V2"]) - 4:+.6f};  '
           f'(u.r^)^2 {U2:.6f} / {v["V1"]:.6f} / {v["V2"]:.6f};  ' + ' '.join(f"l'={lp}: {a:.6f} (doc. {max(l, lp)/(8*(lp + .5)) if lp == J else 0:.6f})" for lp, a in amps.items()))
pr("  singlet p-Ps (S = 0, J = l): (u.r^)^2 QM / V1 / V2 (QM 1/4);  amplitudy do trypletow l' = l +- 1")
for l in range(1, 7):
    U2 = ur2qm(l, 0, l); Q = Sr2qm(l, 0, l); okI &= abs(U2 - 0.25) <= 1e-12 and abs(Q) <= 1e-12
    amps = {lp: A2qm(l, 0, l, lp, 1) for lp in (l - 1, l + 1)}
    for lp, a in amps.items(): okA &= abs(a - max(l, lp)/(8*(l + .5))) <= 1e-12   # reverse of triplet(lp, J=l) -> singlet(l)
    v = {k: ur2cl(k, l, 0, l) for k in okV}
    for k in okV: okV[k] &= np.isfinite(v[k]) and abs(v[k] - U2) <= 1e-12
    pr(f'    l = {l}: {U2:.6f} / {v["V1"]:.6f} / {v["V2"]:.6f};  ' + ' '.join(f"l'={lp}: {a:.6f}" for lp, a in amps.items()))
pr(f"H823i (tozsamosc (S.r^)^2 = 1 - 4(u.r^)^2 i QM = 521, tryplet i singlet): {'TAK' if okI else 'NIE'}")
pr(f"H823s (regula wyboru: tylko singlet l' = J; J = l -> (u.r^)^2 = 0, <S12> = +2): {'TAK' if okS else 'NIE'}")
pr(f"H823a (amplitudy QM = l_>/(8(l'+1/2)) [singlet: odwrotne przejscie]): {'TAK' if okA else 'NIE'}")
pr(f"H823v (GLOWNA: V1 lub V2 = QM, tryplet i singlet): {'TAK' if any(okV.values()) else 'NIE'}  {okV}")
E = {m.group(1): float(m.group(2)) for m in re.finditer(r'(2\^\dP\d|2\^3S1) (\S+)', open(os.path.join(D, '..', 'section-820', 'out', 'p640_820.txt')).read())}
X = {0: -2, 1: -1, 2: 1}; Tq = {J: 6*rule521(1, 1, J) - 4 for J in (0, 1, 2)}
e0, aa, bb = np.linalg.solve(np.array([[1, X[J], Tq[J]] for J in (0, 1, 2)], float), [E[f'2^3P{J}'] for J in (0, 1, 2)])
HEV, C = 4.135667696e-15, 299792458.0; MEAS = {0: (18501.02, 0.61, 'Gurung 2020'), 1: (13012.42, 1.67, 'Hagena 1993'), 2: (8627.94, 0.95, 'Sheldon 2023')}
for var in ('V1', 'V2'):
    for J in (0, 1, 2):
        T = 6*(1 - 4*ur2cl(var, 1, 1, J)) - 4; nu = E['2^3S1'] - E[f'2^3P{J}']; m, s, src = MEAS[J]
        if not np.isfinite(T): pr(f'  RAPORT {var} nu{J}: <S12> nieokreslone (L^2 = 0)'); continue
        nB = nu - bb*(T - Tq[J])
        pr(f'  RAPORT {var} nu{J}: <S12> {T:+.4f}; {nB:.2f} MHz = {nB*1e6*HEV*1e6:.4f} ueV = {C/(nB*1e6)*100:.4f} cm; pomiar {m}({s}) [{src}]: {(nB - m)/s:+.1f} sigma (silnik {(nu - m)/s:+.1f})')

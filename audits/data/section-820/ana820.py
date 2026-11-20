# Audit 820: decomposition of the Ps n = 2 residual (engine - measurement) into a common shift and J structure.
# Usage: python3 ana820.py -> out/ana820.txt.  Written after the prereg, before reading any fit result.
import os, re, math
import numpy as np
from scipy import stats
D = os.path.dirname(os.path.abspath(__file__)); pr = lambda *a: print(*a, flush=True)
HEV, C = 4.135667696e-15, 299792458.0
u3 = lambda f: f'{f:+.3f} MHz = {f*1e6*HEV*1e6:+.5f} ueV' + (f' = {C/abs(f*1e6)*100:.1f} cm' if f else '')
E = {m.group(1): float(m.group(2)) for m in re.finditer(r'(2\^\dP\d|2\^3S1) (\S+)', open(f'{D}/out/p640_820.txt').read())}
nu = {J: E['2^3S1'] - E[f'2^3P{J}'] for J in (0, 1, 2)}; nuS = E['2^3S1'] - E['2^1P1']
# ---- H820t
ref = {0: 18496.12, 1: 13010.90, 2: 8625.26}
okT = all(abs(nu[J] - ref[J]) <= 0.01 for J in ref) and abs(nuS - 11183.90) <= 0.01
pr(f"H820t (silnik = 18 496,12 / 13 010,90 / 8 625,26 / 11 183,90 MHz do 0,01): {'TAK' if okT else 'NIE'}  " +
   ' '.join(f'nu{J} {nu[J]:.4f}' for J in nu) + f' singlet {nuS:.4f}')
for J in nu: pr(f'  nu{J} silnika {nu[J]:.4f} MHz = {nu[J]*1e6*HEV*1e6:.4f} ueV = {C/(nu[J]*1e6)*100:.4f} cm')
pr(f'  2^3S1 - 2^1P1 silnika {nuS:.4f} MHz = {nuS*1e6*HEV*1e6:.4f} ueV = {C/(nuS*1e6)*100:.4f} cm')
# ---- data
M = [(0, 18501.02, 0.61, 'Gurung 2020'), (0, 18499.65, math.hypot(1.20, 4.00), 'Hagena 1993'), (1, 13012.42, math.hypot(0.67, 1.54), 'Hagena 1993'),
     (2, 8627.94, 0.95, 'Sheldon 2023'), (2, 8624.38, math.hypot(0.54, 1.40), 'Hagena 1993')]
LEY = (11180.0, math.hypot(5, 4))
X = {0: -2.0, 1: -1.0, 2: 1.0}; T = {0: -4.0, 1: 2.0, 2: -0.4}
TH = {0: 18498.25, 1: 13012.41, 2: 8626.71}; THS = 11185.37
def fit(r, s, Js, cols):
    A = np.array([[{'c': 1.0, 'a': -X[J], 'b': -T[J]}[k] for k in cols] for J in Js]); W = np.diag(1/np.array(s)**2)
    cov = np.linalg.inv(A.T @ W @ A); p = cov @ A.T @ W @ np.array(r); res = np.array(r) - A @ p
    return p, np.sqrt(np.diag(cov)), float(res @ W @ res), res
r = [nu[J] - m for J, m, s, src in M]; s = [x[2] for x in M]; Js = [x[0] for x in M]
pr('  reszty silnik - pomiar:'); [pr(f'    nu{J} {src}: {ri:+.3f} +/- {si:.3f} MHz ({ri/si:+.2f} sigma)') for (J, m, si, src), ri in zip(M, r)]
MODELS = {'A': ['c'], 'B': ['c', 'a'], 'C': ['c', 'b'], 'D': ['c', 'a', 'b']}; FIT = {}
for k, cols in MODELS.items():
    p, e, chi, res = fit(r, s, Js, cols); dof = len(r) - len(cols); FIT[k] = (p, e, chi, dof)
    pr(f'  model {k}: chi2 {chi:.3f} / {dof} st. sw., p = {stats.chi2.sf(chi, dof):.4f}; ' + ', '.join(f'{c} = {pv:+.3f} +/- {ev:.3f} MHz' for c, pv, ev in zip(cols, p, e)))
pA = stats.chi2.sf(FIT['A'][2], FIT['A'][3]); okA = pA >= 0.05
pr(f"H820a (GLOWNA: model A wystarcza, p >= 0,05): {'TAK' if okA else 'NIE'} (p = {pA:.4f})")
dchi = FIT['A'][2] - FIT['D'][2]; pJ = stats.chi2.sf(dchi, 2); okJ = pJ < 0.05
pr(f"H820j (struktura J istotna, spadek chi2 A -> D na 2 st. sw. p < 0,05): {'TAK' if okJ else 'NIE'} (spadek {dchi:.3f}, p = {pJ:.4f})")
cA, eA = FIT['A'][0][0], FIT['A'][1][0]; rS = nuS - LEY[0]; sig = math.hypot(LEY[1], eA); okS = abs(rS - cA) <= 2*sig
pr(f"H820s (p-Ps 2^3S1 - 2^1P1: reszta {rS:+.3f} +/- {LEY[1]:.3f} MHz vs c_A {cA:+.3f} +/- {eA:.3f}; roznica {(rS - cA)/sig:+.2f} sigma): {'TAK' if okS else 'NIE'}")
d0 = (18501.02 - 18499.65)/math.hypot(0.61, M[1][2]); d2 = (8627.94 - 8624.38)/math.hypot(0.95, M[4][2]); okP = abs(d0) <= 2 and abs(d2) <= 2
pr(f"H820p (spojnosc pomiarow: nu0 Gurung - Hagena {d0:+.2f} sigma, nu2 Sheldon - Hagena {d2:+.2f} sigma, prog 2): {'TAK' if okP else 'NIE'}")
pr('  RAPORT parametry modelu A i D w trzech jednostkach:')
pr(f'    A: c = {u3(cA)} (+/- {eA:.3f} MHz)')
for c, pv, ev in zip(MODELS['D'], FIT['D'][0], FIT['D'][1]): pr(f'    D: {c} = {u3(pv)} (+/- {ev:.3f} MHz, {pv/ev:+.2f} sigma)')
# ---- POST HOC-labelled split via theory (Karshenboim 2004), same bases
pr('  RAPORT (podzial przez teorie Karshenboim 2004 -- nie odniesienie):')
rt = [nu[J] - TH[J] for J in (0, 1, 2)]; pt = np.linalg.solve(np.array([[1.0, X[J], T[J]] for J in (0, 1, 2)]) * np.array([1, -1, -1]), rt)
pr('    silnik - teoria: ' + ' / '.join(f'{x:+.3f}' for x in rt) + f' MHz; singlet {nuS - THS:+.3f}; dokladnie c, a, b = ' + ', '.join(f'{x:+.3f}' for x in pt))
rm = [(TH[J] - m, si, J, src) for J, m, si, src in M]
pm, em, chim, resm = fit([x[0] for x in rm], s, Js, ['c']); pr(f'    teoria - pomiar: ' + '; '.join(f'nu{J} {src} {x:+.3f} ({x/si:+.2f} sigma)' for x, si, J, src in rm) +
   f'; model A teorii: c = {pm[0]:+.3f} +/- {em[0]:.3f}, chi2 {chim:.2f} / 4 (p = {stats.chi2.sf(chim, 4):.4f}); singlet {THS - LEY[0]:+.2f}')

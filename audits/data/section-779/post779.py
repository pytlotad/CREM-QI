# Audit 779, POST HOC (does not change verdicts) -> out/post779.txt
# (1) reconstruction of P_D with the true window: predict776 integrates 300 * Tr0 (not 300 torus periods T); K_eff = 300 Tr0 / T, Fejer kernel with real K_eff;
#     tori with |Phi_closed - Phi_num| > 1e-3 (near-radial orbits, beta 0.05) excluded.
# (2) gradient and Hessian of Phi_n at the own EBK torus (finite differences of the closed form); local quadratic estimate of Phi at the shifted ends of 774/775.
import os, sys, numpy as np, io, contextlib
D = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, D)
with contextlib.redirect_stdout(io.StringIO()): import ana779 as A9
from phase779 import Phi, wrap, fejer, Nmode
pr = lambda *a: print(*a, flush=True)
S, pd, am = A9.S, A9.pd, A9.am
ok = {k: float(np.max(np.abs(wrap(h['Pc'] - h['Pn'])))) < 1e-3 for k, h in S.items()}
for c, lab in ((-1.5707963267948966, 'CAUST -pi/2'), (0.0, 'CAUST 0')):
    R = [k for k, h in S.items() if abs(h['caust'] - c) < 1e-5 and ok[k]]; hit = 0
    for k in R:
        h = S[k]; N0 = 0.5 + np.sqrt(1 + 2 * h['beta']); Keff = 300 * N0 ** 3 / h['N'] ** 3
        F = h['I2'] * fejer(h['Pc'], Keff); hit += am(h, F) == pd(k)[0]
    pr(f"(1) rekonstrukcja z K_eff = 300 Tr0 / T, {lab}: {hit}/{len(R)} (wykluczone torusy z bledem fazy: {sum(abs(h['caust'] - c) < 1e-5 and not ok[k] for k, h in S.items())})")
pr("(2) gradient i hesjan Phi_n na wlasnym torusie EBK (CAUST -pi/2), roznice skonczone eps 1e-4:")
eps = 1e-4
for beta in (0.05, 0.2, 0.8):
    for nr, m in ((0, 0), (0, 1), (0, 2), (1, -1), (1, 2), (2, 0), (2, -2), (0, 3), (3, 2)):
        J0, L0 = nr + 0.5, m; f = lambda j, l: float(Phi(j, l, nr, m, beta, -np.pi / 2))
        gJ = (f(J0 + eps, L0) - f(J0 - eps, L0)) / (2 * eps); gL = (f(J0, L0 + eps) - f(J0, L0 - eps)) / (2 * eps)
        e2 = 1e-3; H = np.array([[f(J0 + e2, L0) - 2 * f(J0, L0) + f(J0 - e2, L0), (f(J0 + e2, L0 + e2) - f(J0 + e2, L0 - e2) - f(J0 - e2, L0 + e2) + f(J0 - e2, L0 - e2)) / 4],
                                 [0, f(J0, L0 + e2) - 2 * f(J0, L0) + f(J0, L0 - e2)]]) / e2 ** 2; H[1, 0] = H[0, 1]
        ev = np.linalg.eigvalsh(H)
        pr(f"   beta {beta} {(nr, m)}: Phi {f(J0, L0):+.1e}, grad ({gJ:+.1e}, {gL:+.1e}), hesjan [[{H[0,0]:+.3f} {H[0,1]:+.3f}] [{H[1,1]:+.3f}]], wartosci wlasne {ev[0]:+.3f} {ev[1]:+.3f} -> {'okreslony' if ev[0] * ev[1] > 0 else 'siodlo'}")
pr("(3) konce z (iv): Phi zwyciezcy wobec 'typowej' wartosci przy gradiencie 2 pi (|Phi| ~ 2 pi d) -- stosunek")
A = os.path.join(D, '..')
for arm, pat, beta, ids in (('774G', 'section-774/out/runs/G_{}.txt', 0.2, range(1, 41)), ('775B05', 'section-775/out/runs/B05_{}.txt', 0.05, range(1, 21)), ('775B80', 'section-775/out/runs/B80_{}.txt', 0.8, range(1, 21))):
    r = []
    for k in ids:
        q = np.loadtxt(os.path.join(A, pat.format(f"{k:02d}")), comments='#', ndmin=2); j, l = q[-20:, 2:4].mean(0); nr, m = A9.near3(j, l)[0]; d = np.hypot(j - nr - 0.5, l - m)
        if 0.01 < d < 0.15 and q[-1, 5] >= 0.9: r.append((abs(float(Phi(j, l, nr, m, beta, -np.pi / 2))), d))
    if r: pr(f"   {arm}: {len(r)} koncow z odl. 0,01-0,15; mediana |Phi| {np.median([x[0] for x in r]):.4f}, mediana odl. {np.median([x[1] for x in r]):.4f}, mediana |Phi| / (2 pi d) {np.median([x[0] / (2 * np.pi * x[1]) for x in r]):.4f}")

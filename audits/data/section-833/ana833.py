# Audit 833: static Thomas comparison.  Usage: python3 ana833.py -> out/ana833.txt.  Written before the runs.  nu_of as ana830/832.
import os, math, glob
import numpy as np
D = os.path.dirname(os.path.abspath(__file__)); R = os.path.join(D, 'out', 'runs')
HEV, C = 4.135667696e-15, 299792458.0
pr = lambda *a: print(*a, flush=True)
units = lambda f: f'{f/1e6:+.4f} MHz = {f*HEV*1e6:+.5f} ueV' + (f' (|lambda| {C/abs(f)*100:.4g} cm)' if f else '')
def gauss(N):
    x, w = np.polynomial.legendre.leggauss(N); u = 0.75 + 0.25*x[::-1]; return 1/u, 0.25*w[::-1]
NODES = dict(zip([f'u5{c}' for c in 'abcde'], zip(*gauss(5))))
nu_of = lambda Om: sum(w*L**2*Om[k] for k, (L, w) in NODES.items())/(2*math.pi)
DT = {'ops': 2.9540e6, 'H0': 9.6536e6, 'H1': -1.6983e6}   # 832 delta_T, Hz
SYS = {'ops': 'o-Ps (PHEN=2)', 'H0': 'wodor F = 0 (PHEN=1)', 'H1': 'wodor F = 1 (PHEN=2)', 'pps': 'p-Ps (PHEN=1)'}
KEYS = ('C1', 'C2', 'L1', 'L2', 'Lb', 'N1', 'N2')
res = {}; okfiles = True
for sk in SYS:
    for nd in NODES:
        p = f'{R}/{sk}_{nd}.out'
        if not os.path.exists(p): okfiles = False; continue
        an = None; rows = []
        for ln in open(p):
            w = ln.split()
            if w[0] == 'SYS': an = (float(w[w.index('an1') + 1]), float(w[w.index('an2') + 1]))
            if w[0] == 'OM': rows.append({k: float(w[w.index(k) + 1]) for k in KEYS})
        okfiles &= len(rows) == 8 and an is not None
        res[(sk, nd)] = (an, rows)
pr(f'pliki kompletne (20 x 8 azymutow): {okfiles}')
def nu(sk, f):   # f(row, an) -> Omega; mean, min, max over azimuth (same azimuth index on all nodes)
    vals = [nu_of({nd: f(res[(sk, nd)][1][i], res[(sk, nd)][0]) for nd in NODES}) for i in range(8)]
    return float(np.mean(vals)), min(vals), max(vals)
V = {}
for sk, lab in SYS.items():
    V[sk] = {
        'C': nu(sk, lambda r, a: r['C1'] + r['C2']), 'L': nu(sk, lambda r, a: r['Lb']), 'Lsum': nu(sk, lambda r, a: r['L1'] + r['L2']),
        'N': nu(sk, lambda r, a: r['N1'] + r['N2']), 'A': nu(sk, lambda r, a: a[0] + a[1]),
        'C-L': nu(sk, lambda r, a: r['C1'] + r['C2'] - r['Lb']),
        'C1-L1': nu(sk, lambda r, a: r['C1'] - r['L1']), 'C2-L2': nu(sk, lambda r, a: r['C2'] - r['L2']),
        'C1': nu(sk, lambda r, a: r['C1']), 'C2': nu(sk, lambda r, a: r['C2']), 'L1': nu(sk, lambda r, a: r['L1']), 'L2': nu(sk, lambda r, a: r['L2']),
        'N1': nu(sk, lambda r, a: r['N1']), 'N2': nu(sk, lambda r, a: r['N2'])}
    pr(f'  {lab}:')
    for k, (m, lo, hi) in V[sk].items(): pr(f'    nu_{k:6s} {units(m)}; azymut {lo/1e6:+.4f} .. {hi/1e6:+.4f} MHz')
    if sk in DT: pr(f'    delta_T (832) {units(DT[sk])}; (nu_C - nu_L) / delta_T {V[sk]["C-L"][0]/DT[sk]:+.3f}')
    pr('    wezly Omega [rad/s] C, L, N, analit.: ' + '; '.join(
        f'{nd} {np.mean([r["C1"]+r["C2"] for r in res[(sk,nd)][1]]):+.3e} {np.mean([r["Lb"] for r in res[(sk,nd)][1]]):+.3e} '
        f'{np.mean([r["N1"]+r["N2"] for r in res[(sk,nd)][1]]):+.3e} {sum(res[(sk,nd)][0]):+.3e}' for nd in NODES))
M3 = ('ops', 'H0', 'H1')
rel = lambda x, y: abs(x - y) <= 0.05*abs(y) if y else abs(x) == 0.0
pr(f"H833k (nu_L = analityczne do 5 %): {'TAK' if all(rel(V[s]['L'][0], V[s]['A'][0]) for s in M3) else 'NIE'}  "
   + ', '.join(f"{s} {V[s]['L'][0]/V[s]['A'][0]:.4f}" for s in M3 if V[s]['A'][0]))
pr(f"H833r (|nu_C - nu_N| <= 0,05 |nu_C|): {'TAK' if all(rel(V[s]['N'][0], V[s]['C'][0]) for s in M3) else 'NIE'}  "
   + ', '.join(f"{s} {V[s]['N'][0]/V[s]['C'][0]:.4f}" for s in M3 if V[s]['C'][0]))
pr(f"H833s (GLOWNA: |(nu_C - nu_L) - delta_T| <= 0,25 |delta_T|): "
   f"{'TAK' if all(abs(V[s]['C-L'][0] - DT[s]) <= 0.25*abs(DT[s]) for s in M3) else 'NIE'}")
# hydrogen: first particle = proton
pr(f"H833p (wodor: |nu_C1 - nu_L1| proton < 0,1 MHz): {'TAK' if all(abs(V[s]['C1-L1'][0]) < 0.1e6 for s in ('H0', 'H1')) else 'NIE'}  "
   + ', '.join(f"{s} {V[s]['C1-L1'][0]/1e6:+.5f} MHz" for s in ('H0', 'H1')))

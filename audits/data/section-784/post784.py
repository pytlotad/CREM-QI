# Audit 784 POST HOC (does not change verdicts) -> out/post784.txt
# (1) Excited pure modes that did not descend: do the dipole-connected lower modes grow (censoring) or not?  Weights from #W (13 lowest modes), ln w slopes in the last third.
# (2) R arms: does w_max grow (selection under way) or stay; R_p of mixtures vs pure modes; weight of (0,0) in time.
import os, numpy as np
D = os.path.dirname(os.path.abspath(__file__)); O = os.path.join(D, 'out', 'runs')
pr = lambda *a: print(*a, flush=True)
LN = [(0,0),(0,1),(0,-1),(1,0),(0,2),(0,-2),(1,1),(1,-1),(2,0),(1,2),(1,-2),(2,1),(2,-1)]
BETA = 0.2; GREF = 0.0003; Nm = lambda t: t[0] + 0.5 + np.sqrt(t[1] ** 2 + 2 * BETA); En = lambda t: -0.5 / Nm(t) ** 2
def W(fn):
    k, w = [], []
    for line in open(fn):
        if line.startswith('#W'): p = line.split(); k.append(int(p[1])); w.append([float(x) for x in p[2:]])
    return np.array(k), np.array(w)
pr("(1) czyste mody wzbudzone: wagi nizszych modow z |dm| = 1 (z 13 najnizszych) -- w(5000), w(15000), w(30000), nachylenie ln w w okresach 20000-30000 [/okres, / Gamma_ref]")
for p, md in (('1e6', (1, 1)), ('1e6', (1, -2)), ('1e6', (2, 0)), ('1e8', (1, -2)), ('1e8', (2, 0)), ('1e8', (1, 1))):
    k, w = W(os.path.join(O, f"N{p}_{md[0]}_{md[1]}.txt")); i5, i15, i30 = [np.searchsorted(k, x) for x in (5000, 15000, 30000)]
    sel = (k >= 20000)
    pr(f"  N{p} {md}: w(start) {w[i5, LN.index(md)]:.4f} / {w[i15, LN.index(md)]:.4f} / {w[-1, LN.index(md)]:.4f}")
    for t in LN:
        if abs(t[1] - md[1]) != 1 or En(t) >= En(md): continue
        y = np.log(np.maximum(w[sel, LN.index(t)], 1e-300)); sl = np.polyfit(k[sel], y, 1)[0]
        pr(f"      -> {str(t):8s} w {w[i5, LN.index(t)]:.2e} / {w[i15, LN.index(t)]:.2e} / {w[-1, LN.index(t)]:.2e}   nachylenie {sl:+.2e} ({sl / GREF:+.2f} Gamma_ref)")
    rest = 1 - w[-1].sum(); pr(f"      waga poza 13 najnizszymi na koncu: {rest:.3f}")
pr("(2) R: w_max (srednia z okresow 1-1000 / 7000-8000 / 14000-15000), w(0,0) na koncu (z kolumn: dominujacy), R_p mediana")
for p in ('6', '8'):
    a = []
    for i in range(1, 21):
        q = np.loadtxt(os.path.join(O, f"R{p}_{i:02d}.txt"), comments='#', ndmin=2)
        a.append((q[:1000, 5].mean(), q[7000:8000, 5].mean(), q[14000:, 5].mean(), np.median(q[:, 10]), (int(q[-1, 6]), int(q[-1, 7]))))
    a = np.array([x[:4] for x in a]); pr(f"  R{p}: mediana w_max {np.median(a[:, 0]):.3f} -> {np.median(a[:, 1]):.3f} -> {np.median(a[:, 2]):.3f}; w_max rosnie (koniec > poczatek + 0,05) {int(np.sum(a[:, 2] > a[:, 0] + 0.05))}/20, "
                      f"maleje {int(np.sum(a[:, 2] < a[:, 0] - 0.05))}/20; mediana R_p {np.median(a[:, 3]):.2e} / okres (czyste mody N{'1e6' if p == '6' else '1e8'}: ~1e-3 / ~1e-4)")

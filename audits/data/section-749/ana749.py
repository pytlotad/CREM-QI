# Audit 749: analysis -> out/ana749.txt
import os, numpy as np
D = os.path.dirname(os.path.abspath(__file__)); O = os.path.join(D, 'out', 'runs'); pr = lambda *a: print(*a, flush=True)
cm = lambda a: np.angle(np.mean(np.exp(1j * a)))
def run(f):
    d = np.loadtxt(os.path.join(O, f), comments='#'); dump = []
    for line in open(os.path.join(O, f)):
        if line.startswith('# DUMP ') and not line.startswith('# DUMP n_r'):
            a = line.split()[2:]; dump.append((int(a[0]), int(a[1]), float(a[2]), float(a[3]), float(a[4])))
    return cm(d[-20:, 4]), dump
res = {}
for f in sorted(os.listdir(O)):
    s, dump = run(f); res[f] = (s, dump); pr(f"{f}: poslizg {s:+.4f} rad (|s| - pi = {abs(s) - np.pi:+.4f})")
    if dump:
        top = sorted(dump, key=lambda r: -r[3])[:8]
        pr("    najwieksze udzialy |c_mem|^2: " + '; '.join(f"({r[0]},{r[1]:+d}) {r[3]:.4f} [w psi(x_p) {r[4]:.3f}]" for r in top))
A = res['A_full_Jr1.5_L0.0.txt'][1]; d = {(r[0], r[1]): r[3] for r in A}
asym = max(abs(d[(n, m)] - d[(n, -m)]) / max(d[(n, m)], 1e-300) for (n, m) in d if m > 0)
pr(f"H749a: max wzgl. roznica udzialow +m / -m = {asym:.2e} -> {'TAK' if asym < 1e-6 else 'NIE'}")
sA, sB, sC = res['A_full_Jr1.5_L0.0.txt'][0], res['B_m0_Jr1.5_L0.0.txt'][0], res['C_mpos_Jr1.5_L0.0.txt'][0]
pr(f"H749b: |s_C - s_A| = {abs(np.angle(np.exp(1j*(sC-sA)))):.4f} -> {'TAK' if abs(np.angle(np.exp(1j*(sC-sA)))) < 0.1 else 'NIE'}")
pr(f"H749c: B (tylko m = 0): | |s| - pi | = {abs(abs(sB) - np.pi):.4f} -> {'TAK' if abs(abs(sB) - np.pi) < 0.1 else 'NIE'}")
own = d[(1, 0)]; dD = {(r[0], r[1]): r[3] for r in res['D_full_Jr0.5_L1.0.txt'][1]}
pr(f"H749d: udzial wlasnego modu A (1,0) {own:.4f}, D (0,1) {dD[(0, 1)]:.4f} -> {'TAK' if own < 0.5 and dD[(0, 1)] > 0.5 else 'NIE'}")

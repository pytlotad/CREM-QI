# Audit 750: analysis -> out/ana750.txt
import os, numpy as np
D = os.path.dirname(os.path.abspath(__file__)); O = os.path.join(D, 'out', 'runs'); pr = lambda *a: print(*a, flush=True)
ld = lambda f: np.loadtxt(os.path.join(O, f), comments='#', ndmin=2)
c = ld('ctrl745.txt')[:, :8]; r = np.loadtxt(os.path.join(D, '..', 'section-745', 'out', 'runs', 'I_k0_Jr0.5_L1.0_eps1e-3.txt'))
pr(f"kontrola: lock750 bez nowych zmiennych = 745 I_k0_Jr0.5_L1.0_eps1e-3 (kolumny 1-8): max |roznica| {np.abs(c - r).max():.1e}")
def ebk(d):
    J, L = d[:, 2], d[:, 3]; n = np.maximum(np.round(J - 0.5), 0); m = np.round(L)
    return n, m, np.hypot(J - n - 0.5, L - m)
def longest(mask):
    best = cur = 0
    for v in mask: cur = cur + 1 if v else 0; best = max(best, cur)
    return best
res = {}
for f in ['R0', 'R0p', 'R1', 'R2', 'R2p', 'R3']:
    d = ld(f + '.txt'); n, m, dE = ebk(d); res[f] = d
    agree = np.mean((d[:, 8] == n) & (np.abs(d[:, 9]) == np.abs(m)))
    slip = np.mean(np.abs(np.abs(d[:, 4]) - np.pi) < 0.3)
    t_in = int((dE < 0.05).sum()); plat = longest(dE < 0.05)
    pr(f"=== {f}: okresy {len(d)}; J_r {d[0,2]:.3f} -> {d[-1,2]:.3f}, L {d[0,3]:.3f} -> {d[-1,3]:.3f}")
    pr(f"  mod dominujacy = (n_r*, +-m*): {agree:.3f}; poslizg w pi +- 0,3: {slip:.3f}; okresy z d_EBK < 0,05: {t_in} (najdluzszy ciag {plat})")
    sw = [i for i in range(1, len(d)) if (d[i, 8], d[i, 9]) != (d[i - 1, 8], d[i - 1, 9])]
    pr(f"  zmian modu dominujacego: {len(sw)}; pierwsze 12: " + '; '.join(f"okr {d[i,0]:.0f} J_r {d[i,2]:.3f} L {d[i,3]:.3f} ({d[i-1,8]:.0f},{d[i-1,9]:+.0f})->({d[i,8]:.0f},{d[i,9]:+.0f}) udzial {d[i-1,10]:.2f}->{d[i,10]:.2f}" for i in sw[:12]))
    for a in range(0, len(d), max(1, len(d) // 12)):
        pr(f"    okr {d[a,0]:5.0f} J_r {d[a,2]:.4f} L {d[a,3]:.4f} poslizg {d[a,4]:+.3f} mod ({d[a,8]:.0f},{d[a,9]:+.0f}) {d[a,10]:.3f} d_EBK {dE[a]:.3f}")
    res[f] = (d, agree, slip, t_in, plat)
pr(f"H750a (R0): {'TAK' if res['R0'][1] >= 0.7 else 'NIE'} ({res['R0'][1]:.3f});  H750b (R0): {'TAK' if res['R0'][2] >= 0.7 else 'NIE'} ({res['R0'][2]:.3f})")
for f, ref in [('R1', 'R0'), ('R2', 'R0'), ('R2p', 'R0p'), ('R3', 'R0')]:
    ok = res[f][3] >= 3 * max(res[ref][3], 1) and res[f][4] >= 100
    pr(f"H750c {f}: czas d<0,05 {res[f][3]} wobec {ref} {res[ref][3]}, plateau {res[f][4]} -> {'TAK' if ok else 'NIE'}")

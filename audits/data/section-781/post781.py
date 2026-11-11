# Audit 781 POST HOC (does not change verdicts): (1) classical Larmor power <a^2> on EBK tori relative to (0,1), eccentricity;
# (2) period at which m = 0, n_r >= 1 ends switch to (0,0) in R1 / R10 / P; (3) deficit and w_max of non-single-mode ends.
import os, numpy as np
D = os.path.dirname(os.path.abspath(__file__)); O = os.path.join(D, 'out', 'runs'); B = 0.2
pr = lambda *a: print(*a, flush=True)
def a2(nr, m, n=200000):
    Lam = np.sqrt(m * m + 2 * B); N = nr + 0.5 + Lam; e = np.sqrt(max(0, 1 - (Lam / N) ** 2))
    M = (np.arange(n) + 0.5) / n * 2 * np.pi; u = M.copy()
    for _ in range(60): u -= (u - e * np.sin(u) - M) / (1 - e * np.cos(u))
    r = N * N * (1 - e * np.cos(u)); return e, np.mean((1 / r**2 - 2 * B / r**3) ** 2), N
e01, ref, _ = a2(0, 1)
pr("(1) Larmor <a^2> na torusach EBK (beta 0,2) wzgledem (0,1):  torus  e  N  <a^2>/<a^2>(0,1)  moc / (E - E00)")
E = lambda N: -0.5 / N**2; _, _, N00 = a2(0, 0)
for t in [(0, 0), (1, 0), (2, 0), (0, 1), (1, 1), (0, 2), (1, 2), (2, 2), (0, 3)]:
    e, v, N = a2(*t); pr(f"  {str(t):8s} {e:.4f}  {N:.3f}  {v / ref:10.3f}   {'-' if t == (0, 0) else f'{v / ref / ((E(N) - E(N00)) / (E(a2(0, 1)[2]) - E(N00))):.3f}'}")
pr("(2) przejscie na (0,0): okres pierwszego prowadzenia (0,0) i poczatek ostatniej serii")
for f in sorted(os.listdir(O)):
    if not f.startswith(('R', 'P')): continue
    q = np.loadtxt(os.path.join(O, f), comments='#', ndmin=2); dm = q[:, 6:8].astype(int)
    z = np.where((dm[:, 0] == 0) & (dm[:, 1] == 0))[0]
    if len(z) and (dm[-1] == 0).all():
        k = len(q)
        while k > 1 and (dm[k - 2] == 0).all(): k -= 1
        pr(f"  {f[:-4]:10s} start dominujacy {tuple(int(x) for x in dm[0])}  pierwsze (0,0) @{z[0] + 1}  ostatnia seria od @{k}  w_max konca {q[-1, 5]:.3f}  J,L konca {q[-20:, 2].mean():.4f} {q[-20:, 3].mean():+.4f}")

# Audit 783 POST HOC (does not change verdicts) -> out/post783.txt
# (1) Gamma_n / Gamma_ref (s = 30) for all modes, recomputed in numpy from the guide783 definition (dipole radial integrals, sum over lower |dm| = 1); check vs #G.
# (2) Width selection: in arms WD and W (P) and in R (WD, s = 1, vs 774 G) -- does a change of dominant mode go to a NARROWER mode (smaller Gamma_n)?
# (3) Arm D: |dm| of every settled step; field norm P growth (Gamma_00 = 0, and in D all widths 0).
import os, numpy as np
from scipy.special import eval_genlaguerre
D = os.path.dirname(os.path.abspath(__file__)); A = os.path.join(D, '..'); O = os.path.join(D, 'out', 'runs')
pr = lambda *a: print(*a, flush=True)
BETA = 0.2; NR, MM = 12, 6
G = 200000; s = (np.arange(G) + 0.5) * np.sqrt(3000.0) / G; r = s * s; wgt = 2 * s * (np.sqrt(3000.0) / G) * r      # r dr
modes = [(nr, m) for m in range(-MM, MM + 1) for nr in range(NR + 1)]
Nm = lambda nr, m: nr + 0.5 + np.sqrt(m * m + 2 * BETA); E = {q: -0.5 / Nm(*q) ** 2 for q in modes}
R = {}
for nr, m in modes:
    nu = np.sqrt(m * m + 2 * BETA); N = Nm(nr, m); rho = 2 * r / N
    f = np.exp(nu * np.log(rho) - rho / 2) * eval_genlaguerre(nr, 2 * nu, rho); R[(nr, m)] = f / np.sqrt(np.sum(wgt * f * f))
I = lambda a, b: np.sum(wgt * r * R[a] * R[b])
g = lambda a, b: abs(E[a] - E[b]) ** 3 * I(a, b) ** 2
ref = g((0, 1), (0, 0)); S = 30.0
Gn = {a: S * sum(g(a, b) for b in modes if abs(a[1] - b[1]) == 1 and E[b] < E[a]) / ref for a in modes}
pr("(1) Gamma_n / Gamma_ref przy s = 30 (numpy; kontrola wobec #G z PWD_0_1)")
Gt = {(int(l.split()[1]), int(l.split()[2])): float(l.split()[4]) for l in open(os.path.join(O, 'PWD_0_1.txt')) if l.startswith('#G')}
pr(f"  max |numpy - #G| / #G (13 modow, bez (0,0)): {max(abs(Gn[k] - v) / v for k, v in Gt.items() if v > 0):.2e}")
for nr in range(4): pr("  n_r=%d: " % nr + '  '.join(f"m{m:+d} {Gn[(nr, m)]:6.2f}" for m in range(0, 5)))
def dom(fn):
    q = np.loadtxt(fn, comments='#', ndmin=2); return q, [(int(a), int(b)) for a, b in q[:, 6:8]]
def settled(dm, L=200):
    seq = []; i = 0
    while i < len(dm):
        j = i
        while j < len(dm) and dm[j] == dm[i]: j += 1
        if j - i >= L and (not seq or seq[-1] != dm[i]): seq.append(dm[i])
        i = j
    return seq
pr("(2) selekcja szerokoscia: kazdy krok miedzy stanami osiadlymi (>= 200 okresow) -- Gamma_n przed -> po (s = 30 skala; w R s = 1, te same proporcje)")
PM = ['0_1', '0_2', '1_1', '1_-2', '2_0']; nar = tot = 0
for arm in ('WD', 'W'):
    for p in PM:
        q, dm = dom(os.path.join(O, f"P{arm}_{p}.txt")); sq = settled(dm)
        st = ' -> '.join(f"{t}[{Gn[t]:.2f}]" for t in sq)
        for a, b in zip(sq, sq[1:]): tot += 1; nar += Gn[b] < Gn[a]
        pr(f"  {arm:2s} {p:5s}: {st}   E: {' -> '.join(f'{E[t]:.4f}' for t in sq)}")
pr(f"  P (WD, W): krokow do wezszego modu {nar}/{tot}")
nr2 = tot2 = 0; up = 0
for i in range(1, 21):
    _, dg = dom(os.path.join(A, 'section-774', 'out', 'runs', f"G_{i:02d}.txt")); _, d3 = dom(os.path.join(O, f"R_{i:02d}.txt"))
    a, b = dg[-1], d3[-1]
    if a != b:
        tot2 += 1; nr2 += Gn[b] < Gn[a]; up += E[b] > E[a]
        pr(f"  R_{i:02d}: 774G {a}[{Gn[a]:.2f}] -> 783 {b}[{Gn[b]:.2f}]  dE {E[b] - E[a]:+.4f}  osiadle 783: {' -> '.join(map(str, settled(d3)))}")
pr(f"  R: koniec 783 inny niz 774G {tot2}/20; wezszy {nr2}/{tot2}; wyzsza energia {up}/{tot2}")
pr("(3) ramie D (bez szerokosci): |dm| krokow miedzy stanami osiadlymi; wzrost normy pola P")
for p in PM + ['0_0']:
    q, dm = dom(os.path.join(O, f"PD_{p}.txt")); sq = settled(dm)
    pr(f"  D {p:5s}: {' -> '.join(map(str, sq))}  |dm| {[abs(b[1] - a[1]) for a, b in zip(sq, sq[1:])]}  P(1) {q[0, 4]:.3e}  P(koniec) {q[-1, 4]:.3e}")
for arm in ('WD', 'W'):
    q, _ = dom(os.path.join(O, f"P{arm}_0_0.txt")); pr(f"  {arm} 0_0: P(1) {q[0, 4]:.3e}  P(10000) {q[9999, 4]:.3e}  P(koniec) {q[-1, 4]:.3e}  (Gamma_00 = 0: naped bez tlumienia)")

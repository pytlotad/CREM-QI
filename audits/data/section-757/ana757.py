# Audit 757: work of the engine-type coupling on the orbit (fixed points of the lock in the field picture) -> out/ana757.txt
import os, numpy as np
D = os.path.dirname(os.path.abspath(__file__)); O = os.path.join(D, 'out', 'runs'); pr = lambda *a: print(*a, flush=True)
ld = lambda f: np.loadtxt(os.path.join(O, f), comments='#', ndmin=2)
c = ld('ctrl744_k0_Jr0.50.txt')[:, :7]; r = np.loadtxt(os.path.join(D, '..', 'section-744', 'out', 'runs', 'k0_Jr0.50.txt'), comments='#')
pr(f"kontrola: lock757 kolumny 1-7 = 744 k0_Jr0.50 : max |roznica| {np.abs(c - r).max():.1e}")
def scan(kind, kt, xs):
    W0, W9, sl, n = [], [], [], []
    for x in xs:
        d = ld(f"{kind}_k{kt}_{'Jr' if kind == 'J' else 'L'}{x}.txt"); n.append(len(d))
        W0.append(d[-20:, 7].mean()); W9.append(d[-20:, 8].mean()); sl.append(np.angle(np.mean(np.exp(1j * d[-20:, 4]))))
    return np.array(W0), np.array(W9), np.array(sl), n
def zeros(xs, W):
    z = []
    for i in range(len(xs) - 1):
        if W[i] == 0 or W[i] * W[i + 1] < 0: z.append((xs[i] - W[i] * (xs[i + 1] - xs[i]) / (W[i + 1] - W[i]), np.sign(W[i + 1] - W[i])))
    return z
Js = [f"{x:.2f}" for x in np.arange(0.20, 1.301, 0.05)]; Ls = [f"{x:.2f}" for x in np.arange(0.60, 1.401, 0.05)]
xJ = np.array([float(x) for x in Js]); xL = np.array([float(x) for x in Ls])
res = {}
for kt in ['10', '0']:
    for kind, xs, xv in [('J', Js, xJ), ('L', Ls, xL)]:
        W0, W9, sl, n = scan(kind, kt, xs); res[(kt, kind)] = (W0, W9)
        pr(f"=== kappa T {kt}, skan {'J_r (L = 1)' if kind == 'J' else 'L (J_r = 0,5)'}; okresy min {min(n)} ===")
        m = max(np.abs(W0).max(), np.abs(W9).max())
        pr("   x      W0/max     W90/max    poslizg")
        for x, a, b, s in zip(xs, W0, W9, sl): pr(f"   {x}  {a/m:+.4f}   {b/m:+.4f}   {s:+.3f}")
        for name, W in [('W0', W0), ('W90', W9)]:
            pr(f"  zera {name}: " + '; '.join(f"{z:.4f} ({'rosnie' if s > 0 else 'maleje'})" for z, s in zeros(xv, W)) + f";  max|{name}|/max {np.abs(W).max()/m:.3f}")
def check_a(kt):
    W0, W9 = res[(kt, 'J')]; W = W0 if np.abs(W0).max() >= np.abs(W9).max() else W9; nm = 'W0' if W is W0 else 'W90'
    z = zeros(xJ, W); a = [p for p in z if abs(p[0] - 0.5) < 0.02]; b = [p for p in z if abs(p[0] - 1.0) < 0.05]
    return nm, W, z, bool(a and b and a[0][1] != b[0][1])
nm, W, z, okA = check_a('10'); pr(f"H757a ({nm}): {'TAK' if okA else 'NIE'} (zera {[round(p[0],4) for p in z]})")
G = 0.1; q = np.exp(-G)   # memory per radial period (Gamma T_r = 0,1 at EBK)
F = lambda th: q * np.sin(th) / (1 - 2 * q * np.cos(th) + q * q)
best = max((abs(np.corrcoef(W, f(2 * np.pi * (xJ - 0.5) + ph))[0, 1]), ph, fn) for fn, f in [('F', F), ('sin', np.sin)] for ph in [0, np.pi / 2, np.pi, -np.pi / 2])
pr(f"H757b: max |r| = {best[0]:.4f} ({best[2]}, przesuniecie fazy {best[1]:+.3f}) -> {'TAK' if best[0] > 0.9 else 'NIE'}")
W0, W9 = res[('10', 'L')]; WL = W0 if np.abs(W0).max() >= np.abs(W9).max() else W9; zl = zeros(xL, WL)
pr(f"H757c: zera W(L): {[round(p[0],4) for p in zl]} -> {'TAK' if any(abs(p[0]-1) < 0.02 for p in zl) else 'NIE'}")
W0, W9 = res[('0', 'J')]; W = W0 if np.abs(W0).max() >= np.abs(W9).max() else W9; z0 = zeros(xJ, W)
pr(f"H757d (kappa 0): zera {[round(p[0],4) for p in z0]} -> {'TAK' if not any(abs(p[0]-0.5) < 0.02 for p in z0) else 'NIE'}")

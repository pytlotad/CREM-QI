# Audit 769: caustic counter in the source phase -- verdicts H769a-f, h -> out/ana769.txt
import os, filecmp, numpy as np
D = os.path.dirname(os.path.abspath(__file__)); O = os.path.join(D, 'out', 'runs'); pr = lambda *a: print(*a, flush=True)
ld = lambda f: np.loadtxt(os.path.join(O, f), comments='#', ndmin=2)
same = filecmp.cmp(os.path.join(O, 'ctrl757_Jr0.50.txt'), os.path.join(D, '..', 'section-757', 'out', 'runs', 'J_k0_Jr0.50.txt'), shallow=False)
pr(f"H769a (CAUST 0 = 757 J_k0_Jr0.50 bajt w bajt): {'TAK' if same else 'NIE'}")
d = ld('J_m2_Jr0.50.txt'); s = np.abs(np.abs(d[-20:, 4]) - np.pi).max()
pr(f"H769b (|poslizg| = pi do 1e-6): max odchylenie {s:.2e} -> {'TAK' if s < 1e-6 else 'NIE'}")
Js = [f"{x:.2f}" for x in np.arange(0.20, 1.701, 0.05)]; Ls = [f"{x:.2f}" for x in np.arange(0.60, 2.401, 0.05)]
xJ = np.array([float(x) for x in Js]); xL = np.array([float(x) for x in Ls]); AR = {'m2': '-pi/2', 'z': '0', 'p2': '+pi/2', 'm4': '-pi/4'}
def scan(kind, a):
    W0, W9, P = [], [], []
    for x in (Js if kind == 'J' else Ls):
        q = ld(f"{kind}_{a}_{'Jr' if kind == 'J' else 'L'}{x}.txt"); assert len(q) == 100, (kind, a, x, len(q))
        W0.append(q[-20:, 7].mean()); W9.append(q[-20:, 8].mean()); P.append(q[-20:, 5].mean())
    return np.array(W0), np.array(W9), np.array(P)
def zeros(xs, W):
    return [(xs[i] - W[i] * (xs[i + 1] - xs[i]) / (W[i + 1] - W[i]), np.sign(W[i + 1] - W[i])) for i in range(len(xs) - 1) if W[i] == 0 or W[i] * W[i + 1] < 0]
near = lambda z, x, t: [p for p in z if abs(p[0] - x) < t]
def cond(WJ, WL):
    zj, zl = zeros(xJ, WJ), zeros(xL, WL); a, b, c = near(zj, 0.5, 0.02), near(zj, 1.5, 0.02), near(zj, 1.0, 0.05)
    ok = bool(a and b and c and a[0][1] == b[0][1] and any(p[1] != a[0][1] for p in c) and near(zl, 1.0, 0.02) and near(zl, 2.0, 0.02))
    return ok, zj, zl
R = {}; C = {}
for a in AR:
    J0, J9, JP = scan('J', a); L0, L9, LP = scan('L', a); R[a] = (J0, J9, L0, L9, JP, LP); m = max(np.abs(J0).max(), np.abs(J9).max()); pr(f"=== ramie CAUST {AR[a]} ===")
    pr("  J_r    W0/max   W90/max   <P>");
    for x, u, v, p in zip(Js, J0, J9, JP): pr(f"  {x}  {u/m:+.4f}  {v/m:+.4f}  {p:.3e}")
    m2 = max(np.abs(L0).max(), np.abs(L9).max()); pr("  L      W0/max   W90/max   <P>")
    for x, u, v, p in zip(Ls, L0, L9, LP): pr(f"  {x}  {u/m2:+.4f}  {v/m2:+.4f}  {p:.3e}")
    for nm, WJ, WL in [('chi 0', J0, L0), ('chi pi/2', J9, L9)]:
        ok, zj, zl = cond(WJ, WL); C[a, nm] = ok
        pr(f"  {nm}: warunek H769c {'TAK' if ok else 'NIE'};  zera J_r: " + ', '.join(f"{z:.3f}{'+' if s > 0 else '-'}" for z, s in zj) + ";  zera L: " + ', '.join(f"{z:.3f}{'+' if s > 0 else '-'}" for z, s in zl))
Hc = C['m2', 'chi 0'] or C['m2', 'chi pi/2']; pr(f"H769c (ramie -pi/2): {'TAK' if Hc else 'NIE'}")
q = np.exp(-0.1); F = lambda th: q * np.sin(th) / (1 - 2 * q * np.cos(th) + q * q); th = 2 * np.pi * (xJ - 0.5)
b = max((abs(np.corrcoef(W, f(th + p))[0, 1]), nm, fn, p) for W, nm in [(R['m2'][0], 'W0'), (R['m2'][1], 'W90')] for fn, f in [('F', F), ('sin', np.sin)] for p in [0, np.pi / 2, np.pi, -np.pi / 2])
pr(f"H769d: max |r| = {b[0]:.4f} ({b[1]}, {b[2]}, przesuniecie {b[3]:+.3f}) -> {'TAK' if b[0] > 0.9 else 'NIE'}")
e = not (C['z', 'chi 0'] or C['z', 'chi pi/2'] or C['m4', 'chi 0'] or C['m4', 'chi pi/2'])
pr(f"H769e (warunek H769c NIE w ramionach 0 i -pi/4): {'TAK' if e else 'NIE'};  ramie +pi/2 (raport): chi 0 {C['p2', 'chi 0']}, chi pi/2 {C['p2', 'chi pi/2']}")
ST = [("0.35", "1.00"), ("0.65", "1.00"), ("0.50", "0.85"), ("0.50", "1.15"), ("0.40", "0.90"), ("0.60", "1.10"), ("0.40", "1.10"), ("0.60", "0.90")]; EB = ("0.50", "1.00")
def end(a, c, s, ep, S):
    q = ld(f"Z_{a}_c{c}_s{s}_e{ep}_Jr{S[0]}_L{S[1]}.txt"); n = {'1e-5': 800, '1e-4': 300}[ep]
    if len(q) < n: return np.inf, len(q), np.nan, np.nan
    j, l = q[-20:, 2].mean(), q[-20:, 3].mean(); return np.hypot(j - 0.5, l - 1), len(q), j, l
pr("=== zatrzask: liczba startow poza EBK (z 8) z koncem < 0,02 od (0,5; 1); odleglosc startu EBK; ramie -pi/2 / ramie 0 ===")
f = False
for c in '01':
    for s in '+-':
        for ep in ['1e-5', '1e-4']:
            cnt = {}; eb = {}; unb = {}; dd = {}
            for a in ['m2', 'z']:
                v = [end(a, c, s, ep, S) for S in ST]; cnt[a] = sum(x[0] < 0.02 for x in v); eb[a] = end(a, c, s, ep, EB)[0]; unb[a] = sum(not np.isfinite(x[0]) for x in v); dd[a] = np.median([x[0] for x in v])
            g = cnt['m2'] >= 6 and eb['m2'] < 0.02 and cnt['z'] <= 2; f |= g
            pr(f"  chi {'pi/2' if c == '1' else '0'}  znak {s}  |eps| {ep}:  {cnt['m2']}/8 (EBK {eb['m2']:.4f}, mediana {dd['m2']:.3f}, unbound {unb['m2']})  /  {cnt['z']}/8 (EBK {eb['z']:.4f}, mediana {dd['z']:.3f}, unbound {unb['z']})  {'*' if g else ''}")
pr(f"H769f: {'TAK' if f else 'NIE'}")
pr("H769g: " + ("WARUNKOWA -- potrzebne biegi dodatkowe" if f else "nie dotyczy (H769f NIE)"))
J0, J9, L0, L9 = R['m2'][:4]; i = np.argmin(abs(xJ - 0.5)); chi = -np.arctan2(J0[i], J9[i]); m = np.abs(np.cos(chi) * J0 + np.sin(chi) * J9).max()
w = lambda a, b: np.cos(chi) * a + np.sin(chi) * b
k = np.argmin(abs(xJ - 1.5)); l1 = np.argmin(abs(xL - 1.0)); l2 = np.argmin(abs(xL - 2.0))
pr(f"H769h (raport): chi* = {chi:+.4f} rad; W_chi*/max: (0,5; 1) {w(J0, J9)[i]/m:+.4f}; (1,5; 1) {w(J0, J9)[k]/m:+.4f}; L-skan (0,5; 1) {w(L0, L9)[l1]/np.abs(w(L0, L9)).max():+.4f}; (0,5; 2) {w(L0, L9)[l2]/np.abs(w(L0, L9)).max():+.4f}")

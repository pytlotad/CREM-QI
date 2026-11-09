# Audit 760: branches of the source slip and fixed points of the work -> out/ana760.txt
import os, numpy as np
D = os.path.dirname(os.path.abspath(__file__)); O = os.path.join(D, 'out', 'runs'); pr = lambda *a: print(*a, flush=True)
xs = [f"{x:.3f}" for x in np.arange(0.200, 1.4001, 0.025)]; xJ = np.array([float(x) for x in xs])
def branch(s):
    if -1.4 <= s <= -0.6: return 'B'
    if abs(np.angle(np.exp(1j * (s + np.pi)))) <= 0.4: return 'A'
    return '-'
okA = True; okB = False; summary = {}
for tag, name in [('k1', 'kappa T 1'), ('k3', 'kappa T 3'), ('sol', 'soliton (10; 0,1)')]:
    W0, W9, S = [], [], []
    for x in xs:
        d = np.loadtxt(os.path.join(O, f'{tag}_Jr{x}.txt'), comments='#', ndmin=2)
        W0.append(d[-20:, 7].mean()); W9.append(d[-20:, 8].mean()); S.append(np.angle(np.mean(np.exp(1j * d[-20:, 4]))))
    W0, W9, S = map(np.array, (W0, W9, S)); br = [branch(s) for s in S]; m = max(np.abs(W0).max(), np.abs(W9).max())
    pr(f"=== {name} ===")
    pr("   J_r    poslizg  galaz  W0/max   W90/max")
    for x, s, b, a, c in zip(xs, S, br, W0, W9): pr(f"   {x}  {s:+.3f}   {b}    {a/m:+.4f}  {c/m:+.4f}")
    sB = S[[b == 'B' for b in br]]
    if len(sB): pr(f"  galaz B: {len(sB)} punktow, poslizg srednio {sB.mean():+.3f}, zakres {sB.min():+.3f} ... {sB.max():+.3f} (-pi/3 = {-np.pi/3:+.3f})")
    summary[tag] = sB
    zB = []
    for nm, W in [('W0', W0), ('W90', W9)]:
        for i in range(len(xs) - 1):
            if W[i] * W[i + 1] < 0 and br[i] == br[i + 1] and br[i] in 'AB':
                z = xJ[i] - W[i] * (xJ[i + 1] - xJ[i]) / (W[i + 1] - W[i]); s = S[i] + (S[i + 1] - S[i]) * (z - xJ[i]) / (xJ[i + 1] - xJ[i])
                c = [k - s / (2 * np.pi) for k in range(-1, 3)] + [k + 0.5 - s / (2 * np.pi) for k in range(-1, 3)] if br[i] == 'B' else [k * 0.5 for k in range(0, 4)]
                dmin = min(abs(z - v) for v in c); ok = dmin < 0.03; okA &= ok
                slope = np.sign(W[i + 1] - W[i])
                pr(f"  zero {nm} na galezi {br[i]}: J_r {z:.4f}, poslizg {s:+.3f}, przewidywane najblizsze {min(c, key=lambda v: abs(z-v)):.4f} (odl. {dmin:.4f}) "
                   f"{'zgodne' if ok else 'NIE'}; nachylenie {'+' if slope > 0 else '-'} (stabilne dla eps {'<' if slope > 0 else '>'} 0)")
                if br[i] == 'B': zB.append((nm, z, slope))
    for nm, z, sl in zB:
        if any(n2 == nm and abs(z2 - z - 1) < 0.05 and s2 == sl for n2, z2, s2 in zB): okB = True; pr(f"  powtorzenie z okresem 1: {nm} {z:.4f} i ~{z+1:.4f}")
allB = np.concatenate([v for v in summary.values() if len(v)])
okC = all(len(v) == 0 or np.ptp(v) <= 0.3 for v in summary.values())
pr(f"H760a (zera = samouzgodnienie): {'TAK' if okA else 'NIE'}")
pr(f"H760b (zero na galezi B powtarza sie z okresem 1): {'TAK' if okB else 'NIE'}")
pr(f"H760c (poslizg galezi B staly w +-0,15): {'TAK' if okC else 'NIE'}")
pr(f"H760d (srednia galezi B = -pi/3 do 0,1): {allB.mean():+.3f} -> {'TAK' if abs(allB.mean() + np.pi/3) < 0.1 else 'NIE'}")

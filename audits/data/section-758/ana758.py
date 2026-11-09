# Audit 758: kappa scan -- Maslov take-over vs mismatch sensing -> out/ana758.txt
import os, numpy as np
D = os.path.dirname(os.path.abspath(__file__)); pr = lambda *a: print(*a, flush=True)
O = os.path.join(D, 'out', 'runs'); O7 = os.path.join(D, '..', 'section-757', 'out', 'runs')
Js = [f"{x:.2f}" for x in np.arange(0.20, 1.301, 0.05)]; xJ = np.array([float(x) for x in Js])
def load(kt):
    base = O7 if kt in ('0', '10') else O
    W0, W9, sl = [], [], []
    for x in Js:
        d = np.loadtxt(os.path.join(base, f'J_k{kt}_Jr{x}.txt'), comments='#', ndmin=2)
        W0.append(d[-20:, 7].mean()); W9.append(d[-20:, 8].mean()); sl.append(np.angle(np.mean(np.exp(1j * d[-20:, 4]))))
    return np.array(W0), np.array(W9), np.array(sl)
def zeros(W):
    return [(xJ[i] - W[i] * (xJ[i + 1] - xJ[i]) / (W[i + 1] - W[i]), np.sign(W[i + 1] - W[i])) for i in range(len(W) - 1) if W[i] * W[i + 1] < 0]
i5, i4, i6 = Js.index('0.50'), Js.index('0.40'), Js.index('0.60'); okA = False; okC = True; rows = []
for kt in ['0', '0.03', '0.1', '0.3', '1', '3', '10', '30']:
    W0, W9, sl = load(kt); s = sl[i5]; mas = abs(s + np.pi) < 0.3
    feel = {}
    for nm, W in [('W0', W0), ('W90', W9)]:
        z = [p for p in zeros(W) if abs(p[0] - 0.5) < 0.05]; feel[nm] = (bool(z) and W[i4] * W[i6] < 0, z)
    f = feel['W0'][0] or feel['W90'][0]; okA |= mas and f
    if float(kt) >= 1: okC &= not f
    m = max(np.abs(W0).max(), np.abs(W9).max())
    pr(f"kappa T {kt:>4}: poslizg EBK {s:+.3f} (|s + pi| {abs(s+np.pi):.3f}, Maslow {'TAK' if mas else 'nie'}); "
       f"W0/max przy 0,40/0,50/0,60: {W0[i4]/m:+.3f} {W0[i5]/m:+.3f} {W0[i6]/m:+.3f}; W90: {W9[i4]/m:+.3f} {W9[i5]/m:+.3f} {W9[i6]/m:+.3f}; "
       f"zero przy EBK: {'TAK' if f else 'nie'}")
    pr(f"      zera W0: {[round(p[0],3) for p in zeros(W0)]}; zera W90: {[round(p[0],3) for p in zeros(W9)]}; poslizg na skanie: " + ' '.join(f"{v:+.2f}" for v in sl[::2]))
    rows.append((float(kt), s))
pr(f"H758a: {'TAK' if okA else 'NIE'};  H758c: {'TAK' if okC else 'NIE'}")
k = [r[0] for r in rows]; s = [r[1] for r in rows]
for a, b in zip(rows[:-1], rows[1:]):
    if (a[1] + np.pi / 2) * (b[1] + np.pi / 2) < 0: pr(f"H758b: poslizg przechodzi przez -pi/2 miedzy kappa T {a[0]} ({a[1]:+.3f}) i {b[0]} ({b[1]:+.3f})")

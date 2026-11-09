# Audit 759: soliton source -- Maslov take-over vs mismatch sensing -> out/ana759.txt
import os, numpy as np
D = os.path.dirname(os.path.abspath(__file__)); O = os.path.join(D, 'out', 'runs'); pr = lambda *a: print(*a, flush=True)
Js = ['0.30', '0.35', '0.40', '0.45', '0.50', '0.55', '0.60', '0.65', '0.70', '0.75', '0.80', '1.00']; xJ = np.array([float(x) for x in Js])
ref758 = {'3': -0.904, '10': -3.103, '30': -3.127}
i4, i5, i6 = Js.index('0.40'), Js.index('0.50'), Js.index('0.60'); okA = False; okB = True
for ke in ['3', '10', '30']:
    for g in ['0.03', '0.1', '0.3', '1', '10']:
        W0, W9, sl, am, n = [], [], [], [], []
        for x in Js:
            d = np.loadtxt(os.path.join(O, f'S_ke{ke}_g{g}_Jr{x}.txt'), comments='#', ndmin=2); n.append(len(d))
            W0.append(d[-20:, 7].mean()); W9.append(d[-20:, 8].mean()); sl.append(np.angle(np.mean(np.exp(1j * d[-20:, 4])))); am.append(d[-1, 10])
        W0, W9, sl = map(np.array, (W0, W9, sl)); s = sl[i5]; mas = abs(s + np.pi) < 0.3
        def feel(W):
            if W[i4] * W[i6] >= 0: return None
            for i in range(i4, i6):
                if W[i] * W[i + 1] < 0:
                    z = xJ[i] - W[i] * (xJ[i + 1] - xJ[i]) / (W[i + 1] - W[i])
                    if abs(z - 0.5) < 0.05: return z
            return None
        f0, f9 = feel(W0), feel(W9); f = f0 is not None or f9 is not None; okA |= mas and f
        if g == '10': okB &= abs(s - ref758[ke]) < 0.15
        m = max(np.abs(W0).max(), np.abs(W9).max())
        pr(f"kappa_eff T {ke:>2}, gamma T {g:>4}: poslizg EBK {s:+.3f} ({'Maslow' if mas else '-'}); W0/max 0,40/0,50/0,60 {W0[i4]/m:+.3f} {W0[i5]/m:+.3f} {W0[i6]/m:+.3f}, "
           f"W90 {W9[i4]/m:+.3f} {W9[i5]/m:+.3f} {W9[i6]/m:+.3f}; zero przy EBK: {'W0 %.3f' % f0 if f0 is not None else ''} {'W90 %.3f' % f9 if f9 is not None else ''}{'' if f else 'nie'}; "
           f"max|a-1| {max(am):.3f}; okresy min {min(n)}")
        pr("      poslizg na skanie: " + ' '.join(f"{v:+.2f}" for v in sl))
pr(f"H759a: {'TAK' if okA else 'NIE'};  H759b (gamma T 10 = 758 do 0,15): {'TAK' if okB else 'NIE'}")

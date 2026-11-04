# Audit 701 POST HOC: order p and limit x0 of n_clock - n_E from the step densities 512 / 2048 / 8192 (x(h) = x0 + a h^p, h ~ 1/steps),
# per run; the mean of H and Hbar against Sommerfeld; the H - Hbar difference against the fine-structure interval.
# usage: python3 post701.py <701 runs dir with s512/ s2048/ s8192/>
import sys, os, math
import numpy as np
from ana701 import parse, derived, ALPHA, PAIRS
R = sys.argv[1]; out = []; x0s = {}
INTERVAL = 10969.045e6*4.135667696e-15    # eV: 2P3/2 - 2P1/2 = 9911.200 (Hagley & Pipkin 1994) + 1057.845 MHz (Lundeen & Pipkin 1981)
for s in ('0', '1'):
    for tag, name in PAIRS:
        key = f'{tag}_B0.01_A2_{s}.err'; v = {}
        for st in ('512', '2048', '8192'):
            rows = derived(parse(f'{R}/s{st}/{key}')); v[st] = np.mean([r['ncl'] - r['nE'] for r in rows])
            k = np.mean([r['X'] for r in rows]); n = np.mean([r['nE'] for r in rows])
        r1, r2 = v['512'] - v['2048'], v['2048'] - v['8192']
        p = math.log(r1/r2)/math.log(4.0) if r1*r2 > 0 else float('nan')
        x0 = v['8192'] - r2/(4.0**p - 1.0)
        x0s[(s, name)] = (x0, p, ALPHA**2*(1/(2*k) - 3/(8*n)))
        out.append(f'{name:12s} seed {s}: n_clock - n_E 512 {v["512"]:+.4e} 2048 {v["2048"]:+.4e} 8192 {v["8192"]:+.4e} | p {p:.3f} x0 {x0:+.4e}')
for s in ('0', '1'):
    (xa, pa, sa), (xh, ph, sh) = x0s[(s, 'antihydrogen')], x0s[(s, 'hydrogen')]
    mean, pred = 0.5*(xa + xh), 0.5*(sa + sh); ok = all(min(abs(q - 1), abs(q - 2)) < 0.2 for q in (pa, ph))
    verdict = ('zgodna' if abs(mean/pred - 1) <= 0.2 else 'niezgodna') if ok else 'nierozstrzygnieta'
    d = xh - xa; e = d*2*13.598/8
    out.append(f'seed {s}: mean (H, Hbar) x0 {mean:+.4e} vs Sommerfeld {pred:+.4e} (ratio {mean/pred:.3f}) -> {verdict} | H - Hbar {d:+.4e} -> {e*1e6:.2f} ueV'
               f' = {e/4.135667696e-15/1e9:.3f} GHz = {e/INTERVAL:.3f} x measured 2P fine-structure interval ({INTERVAL*1e6:.2f} ueV)')
text = '\n'.join(out); print(text); os.makedirs('out', exist_ok=True); open('out/post701.txt', 'w').write(text + '\n')

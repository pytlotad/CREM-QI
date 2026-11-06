"""Audit 723, POST HOC diagnostic (no criterion): does the B4 selectivity sharpen with K?  Needle points from run 3 and random tori."""
import numpy as np, sys, os
HERE = os.path.dirname(os.path.abspath(__file__))
src = open(os.path.join(HERE, 'ana723.py')).read().split("P('Audit 723")[0]
g = {'__file__': os.path.join(HERE, 'ana723.py')}; exec(src, g)
out = []
pts = [(0.0123288, 1.0, 0.05), (0.0342466, 1.0, 0.05), (0.2909836, 2.0, 0.05), (0.6174891, 0.5, 0.05)]
pts[3] = (0.5, 0.6174891, 0.05)
for K in (400, 4000, 40000):
    g['K'] = K
    vals = ['(%.4f, %.4f) %.4f' % (j, l, g['B_space'](j, l, b)) for j, l, b in pts]
    ebk = min(g['B_space'](nr + 0.5, m, b) for b in (0.05, 0.2) for nr in range(3) for m in (1, 2, 3))
    rng = np.random.default_rng(723); hi = 0; n = 0
    for b in (0.05, 0.2):
        for _ in range(100):
            v = g['B_space'](rng.uniform(0.05, 3.0), rng.uniform(0.6, 3.4), b); n += 1; hi += v > 0.5
    s = 'K %6d: run-3 needle points B4 %s; EBK min %.4f; random tori B4 > 0.5: %d / %d' % (K, ', '.join(vals), ebk, hi, n)
    out.append(s); print(s, flush=True)
open(os.path.join(HERE, 'out', 'diag723.txt'), 'w').write('\n'.join(out) + '\n')
# stability across K: a torus counts only if B4 > 0.5 at every K
rng = np.random.default_rng(723); allk = 0; n = 0; keep = []
for b in (0.05, 0.2):
    for _ in range(100):
        j, l = rng.uniform(0.05, 3.0), rng.uniform(0.6, 3.4); n += 1; ok = True; vs = []
        for K in (400, 4000, 40000):
            g['K'] = K; v = g['B_space'](j, l, b); vs.append(v); ok &= v > 0.5
        if ok: allk += 1; keep.append('(%.3f, %.3f, beta %.2f) %s' % (j, l, b, '/'.join('%.3f' % v for v in vs)))
s = 'random tori with B4 > 0.5 at all K = 400, 4000, 40000: %d / %d %s' % (allk, n, '; '.join(keep))
out.append(s); print(s)
open(os.path.join(HERE, 'out', 'diag723.txt'), 'w').write('\n'.join(out) + '\n')

# Audit 687, POST HOC (after ana687.py): (1) H687e pooled -- every cell has z < 0, so the combined z = sum z / sqrt(8) and the overall j = 3/2
# frequency; (2) nj sc - QM = (v1 - v3) + (v3 - QM) by H687a, so the PHEN 2 residual of nj is the photon-time-state term (A) of 684/686.
# usage: python3 post687.py <687 runs dir> <686 runs dir> <685 runs dir> (the ana687 arguments; ana687 runs on import, rewriting out/ana687.txt identically)
import sys, glob, math
from scipy.stats import norm
from ana687 import parse, PAIRS, R6
zs = []; c = e = 0.0; out = []
for tag, name in PAIRS:
    for B in ('0.01', '0.1'):
        for P in ('1', '2'):
            rows = [r for r in (parse(f) for f in sorted(glob.glob(f'{R6}/born/{tag}_B{B}_A{P}_*.err'))) if 'jcl' in r]
            p3 = [r['kept'] if r['jcl'] == 1.5 else 1 - r['kept'] for r in rows]
            cc = sum(r['jcl'] == 1.5 for r in rows); ee = sum(p3)
            zs.append((cc - ee)/math.sqrt(sum(x*(1 - x) for x in p3))); c += cc; e += ee; n = len(rows)
z = sum(zs)/math.sqrt(len(zs))
out.append(f'H687e pooled: {len(zs)} cells, all z < 0: {all(x < 0 for x in zs)}; combined z {z:.2f} (two-sided p {2*norm.sf(abs(z)):.2g}),'
           f' sum z^2 {sum(x*x for x in zs):.1f} on {len(zs)}; classical j = 3/2 {c:.0f}/1024 = {c/1024:.3f}, Born {e/1024:.3f}')
text = '\n'.join(out); print(text)
open('out/post687.txt', 'w').write(text + '\n')

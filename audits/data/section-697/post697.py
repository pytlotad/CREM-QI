# Audit 697 POST HOC: paired change of the path gap sc - QM (inst minus 693 poi, same seeds) and of the sampled 1S label distribution.
# usage: python3 post697.py <697 runs dir with on/> <693 runs dir with poi/>
import sys, glob, os, math
import numpy as np
from scipy.stats import chi2
from ana697 import parse, hot, PAIRS
R, R3 = sys.argv[1:3]
out = ['paired gap change (inst - 693 poi) per cell: mean per 1S state, Hotelling on 3 states (chi^2(3)); label change rate; cross-check gap pattern vs (A)']
for tag, name in PAIRS:
    T = 0.0
    for P in ('1', '2'):
        fs = sorted(glob.glob(f'{R}/on/{tag}_B1.0329_A{P}_*.err'))
        a = [parse(f) for f in fs]; b = [parse(f'{R3}/poi/{os.path.basename(f)}') for f in fs]
        dG = np.array([[x - y for x, y in zip(p['sc'], q['sc'])] for p, q in zip(a, b)])
        mu, t2 = hot(dG); T += t2
        moved = sum(max(abs(x - y) for x, y in zip(p['sc'], q['sc'])) > 1e-6 for p, q in zip(a, b))
        Ai = np.array([[x - y for x, y in zip(p['v3'], p['qm'])] for p in a]).mean(0)
        G = np.array([[x - y for x, y in zip(p['sc'], p['qm'])] for p in a]).mean(0)
        out.append(f'{name:12s} PHEN {P} n {len(fs)} d(sc) ' + ' '.join(f'{x:+.4f}' for x in mu) + f' T2 {t2:.2f} p {chi2.sf(t2, 3):.3g}'
                   f' | sc changed in {moved}/{len(fs)} | corr(mean gap, mean (A)) over 4 states {np.corrcoef(G, Ai)[0, 1]:+.2f}')
    out.append(f'{name:12s} pooled paired chi^2(6) {T:.1f} p {chi2.sf(T, 6):.3g}')
text = '\n'.join(out); print(text)
os.makedirs('out', exist_ok=True); open('out/post697.txt', 'w').write(text + '\n')

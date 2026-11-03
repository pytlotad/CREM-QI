# Audit 693 POST HOC: (A) of the earliest photons (tau < 0.15 ns, a few Larmor periods at 1 T) and the hydrogen PHEN 2 gap against 690.
# usage: python3 post693.py <693 runs dir with poi/>
import sys, glob, math
import numpy as np
R = sys.argv[1]
PAIRS = (('antiproton', 'antihydrogen'), ('proton', 'hydrogen'))


def vec(t, key, n=4):
    i = t.index(key); return [float(x) for x in t[i+1:i+1+n]]


def parse(path):
    d = {}
    for line in open(path):
        t = line.split()
        if line.startswith('PAIR QM1S684'): d['v3'] = vec(t, 'v3')
        elif line.startswith('PAIR QM1S path'): d['qm'] = vec(t, 'qm')
        elif line.startswith('END'): d['tau'] = float.fromhex(t[t.index('tau')+1])*1e9
    return d


out = []
for tag, name in PAIRS:
    for P in ('1', '2'):
        rows = [parse(f) for f in sorted(glob.glob(f'{R}/poi/{tag}_B1.0329_A{P}_*.err'))]
        A = np.array([[a - b for a, b in zip(r['v3'], r['qm'])] for r in rows]); tau = np.array([r['tau'] for r in rows])
        s = int(np.argmax(np.abs(A.mean(0)))); e = tau < 0.15
        out.append(f'{name:12s} PHEN {P} min tau {tau.min()*1e3:.1f} ps; tau < 0.15 ns n {e.sum()}: (A)_{s} mean {A[e, s].mean():+.4f} +- {A[e, s].std(ddof=1)/math.sqrt(e.sum()):.4f}'
                   f' rms |A| {np.sqrt((A[e]**2).sum(1).mean()):.4f} | all: (A)_{s} {A[:, s].mean():+.4f} +- {A[:, s].std(ddof=1)/math.sqrt(len(A)):.4f} rms |A| {np.sqrt((A**2).sum(1).mean()):.4f}')
text = '\n'.join(out); print(text); open('out/post693.txt', 'w').write(text + '\n')

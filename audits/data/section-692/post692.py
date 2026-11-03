# Audit 692 POST HOC: the photons with the largest step sensitivity (max |d sc|, |d dE|) per step, and the size of (A) itself at 0.05.
# usage: python3 post692.py <692 runs dir> <686 runs dir> <690 runs dir>
import sys, glob, os
import numpy as np
sys.argv = sys.argv[:4]
from ana692 import parse, ref, PAIRS, R
out = []
for v in ('0.025', '0.1'):
    rows = []
    for f in sorted(glob.glob(f'{R}/s{v}/*.err')):
        p = parse(f); q = parse(ref(os.path.basename(f).split('_')[1][1:], os.path.basename(f)))
        rows.append((max(abs(x - y) for x, y in zip(p['sc'], q['sc'])), abs(p['dE'] - q['dE']), os.path.basename(f), p['born'], q['born']))
    rows.sort(reverse=True)
    out.append(f'step {v}: photons with max |d sc| > 1e-3: {sum(r[0] > 1e-3 for r in rows)}/{len(rows)}; top 3:')
    out += [f'  {n} max|d sc| {a:.2e} |d dE| {b:.3f} ueV born {c}/{d}' for a, b, n, c, d in rows[:3]]
for tag, name in PAIRS:
    for B in ('0.01', '1.0329'):
        for P in ('1', '2'):
            A = np.array([[x - y for x, y in zip(d['v3'], d['qm'])] for d in (parse(ref(B, os.path.basename(f))) for f in sorted(glob.glob(f'{R}/s0.1/{tag}_B{B}_A{P}_*.err')))])
            out.append(f'(A) at 0.05, {name:12s} B {B:>6} T PHEN {P}: mean ' + ' '.join(f'{x:+.4f}' for x in A.mean(0)))
text = '\n'.join(out); print(text); open('out/post692.txt', 'w').write(text + '\n')

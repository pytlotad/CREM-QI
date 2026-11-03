# Audit 689, POST HOC (after ana689.py): at PHEN 2 the classical j follows p3 more sharply than a Born draw -- frequency of jcl = 3/2 in bins
# of p3 (pooled over pairs and fields, per PHEN) against the Born diagonal P(jcl = 3/2 | p3) = p3.
# usage: python3 post689.py <689 runs dir> <686 runs dir> (the ana689 arguments; ana689 runs on import, rewriting out/ana689.txt identically)
import glob, math
from ana689 import parse, PAIRS, R
EDGES = (0.0, 0.2, 0.4, 0.6, 0.8, 0.9, 1.0001)
out = ['P(jcl = 3/2 | p3 bin) pooled over pairs and fields; Born draw = mean p3 of the bin; z binomial']
for P in ('1', '2'):
    rows = [r for tag, _ in PAIRS for B in ('0.01', '0.1') for r in (parse(f) for f in sorted(glob.glob(f'{R}/born/{tag}_B{B}_A{P}_*.err'))) if 'p3' in r]
    for lo, hi in zip(EDGES, EDGES[1:]):
        b = [r for r in rows if lo <= r['p3'] < hi]
        if not b: continue
        f = sum(r['jcl'] == 1.5 for r in b)/len(b); m = sum(r['p3'] for r in b)/len(b)
        v = sum(r['p3']*(1 - r['p3']) for r in b)
        z = (sum(r['jcl'] == 1.5 for r in b) - sum(r['p3'] for r in b))/math.sqrt(v) if v > 0 else float('nan')
        out.append(f'PHEN {P} p3 [{lo:.1f}, {min(hi, 1):.1f}) n {len(b):3d} jcl 3/2 {f:.3f} Born {m:.3f} z {z:+.2f}')
text = '\n'.join(out); print(text)
open('out/post689.txt', 'w').write(text + '\n')
# (2) cost of a Born draw of j at the photon: expected fraction of photons whose drawn j differs from the classical j, mean I (1 - p3) + (1 - I) p3.
out2 = ['expected fraction of photons with a Born-drawn j != classical j (each would leave the 2P fine-structure gap, ~45 ueV, in the ledger)']
for P in ('1', '2'):
    for tag, name in PAIRS:
        rows = [r for B in ('0.01', '0.1') for r in (parse(f) for f in sorted(glob.glob(f'{R}/born/{tag}_B{B}_A{P}_*.err'))) if 'p3' in r]
        mis = sum(r['p3'] if r['jcl'] == 0.5 else 1 - r['p3'] for r in rows)/len(rows)
        out2.append(f'PHEN {P} {name:12s} n {len(rows)} mismatch {mis:.3f}')
text2 = '\n'.join(out2); print(text2)
open('out/post689.txt', 'a').write(text2 + '\n')

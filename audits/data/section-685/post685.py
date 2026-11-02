# Audit 685, POST HOC: (1) pooled z^2 of the 1S gap over all cells, born vs default; (2) the ledger import together with the booked dE of the clipped
# L' (dE is not in the 666 ledger): rms of imp + dE and imp - dE (sign convention of dE relative to the ledger not fixed in advance).
# usage: python3 post685.py <685 runs dir> <684 runs dir>
import sys, glob, math
sys.argv += []
from ana685 import parse, R, R4, FIELDS  # noqa: E402  (ana685 writes out/ana685.txt again, unchanged)
out = []
tot = {'born': [0.0, 0], 'def': [0.0, 0]}
for B in FIELDS:
    for P in ('1', '2'):
        fs = sorted(glob.glob(f'{R}/born/B{B}_A{P}_*.err'))
        for name, rows in (('born', [parse(f) for f in fs]), ('def', [parse(f.replace(R + '/born/', R4 + '/def/')) for f in fs])):
            rows = [r for r in rows if 'sc' in r]; n = len(rows)
            for i in range(4):
                d = [r['sc'][i] - r['qm'][i] for r in rows]; mu = sum(d)/n
                s = math.sqrt(sum((x - mu)**2 for x in d)/(n - 1)/n)
                if s > 0: tot[name][0] += (mu/s)**2; tot[name][1] += 1
        rows = [parse(f) for f in fs]
        both = [(r['imp'], r.get('dE', 0.0)) for r in rows if 'imp' in r]
        dflt = [parse(f.replace(R + '/born/', R4 + '/def/'))['imp'] for f in fs]
        rp = math.sqrt(sum((a + b)**2 for a, b in both)/len(both)); rm = math.sqrt(sum((a - b)**2 for a, b in both)/len(both))
        rd = math.sqrt(sum(x*x for x in dflt)/len(dflt))
        out.append(f'B {B:>6} T PHEN {P}: import rms def {rd:.2f}, born imp+dE {rp:.2f} ({rp/rd:.2f}x), imp-dE {rm:.2f} ({rm/rd:.2f}x) ueV')
for k, (c, n) in tot.items():
    out.append(f'pooled 1S gap {k}: sum z^2 {c:.1f} for {n} cells (4 states per photon sum to 1, so ~{3*n//4} independent)')
text = '\n'.join(out); print(text); open('out/post685.txt', 'w').write(text + '\n')

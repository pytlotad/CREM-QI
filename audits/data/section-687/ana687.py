# Audit 687: the Born-Delta-m path without the 668/670 j rule (CREM_NO_PAIR_J_PROJECTION=1, "nj") against QM and against the 686 born path
# (same seeds); the j rule diagnosis on the 686 born files: classical j = 3/2 frequency against the Born P(j = 3/2) of the photon-time state.
# usage: python3 ana687.py <687 runs dir with nj/ ctl/> <686 runs dir with born/> <685 runs dir with born/>
import sys, glob, math, os
import numpy as np
from scipy.stats import chi2
R, R6, R5 = sys.argv[1:4]
UEV = 1.602176634e-25
PAIRS = (('antiproton', 'antihydrogen'), ('proton', 'hydrogen'))


def vec(t, key, n=4):
    i = t.index(key); return [float(x) for x in t[i+1:i+1+n]]


def parse(path):
    d = {}
    for line in open(path):
        t = line.split()
        if line.startswith('PAIR QM1S684'):
            d['line684'] = line; d['v1'] = vec(t, 'v1'); d['v3'] = vec(t, 'v3')
        elif line.startswith('PAIR QM1S path'):
            d['sc'], d['qm'] = vec(t, 'sc'), vec(t, 'qm')
        elif line.startswith('PAIR eigen ledger'):
            e = {t[i]: float(t[i+1]) for i in range(3, len(t)-1, 2)}
            d['imp'] = (e['B'] + e['S2P'] - e['S1S'])/UEV
        elif line.startswith('PAIR j projection jcl'):
            d['jcl'] = float(t[t.index('jcl')+1]); d['kept'] = float(t[t.index('kept')+1])
    return d


def rms(xs): return math.sqrt(sum(x*x for x in xs)/len(xs)) if xs else float('nan')


def hot(D):
    D = np.array(D); n = len(D); mu = D.mean(0)
    t2 = float(n*mu[:3] @ np.linalg.solve(np.cov(D[:, :3], rowvar=False), mu[:3])); return n, mu, t2


out = ['nj = born + CREM_NO_PAIR_J_PROJECTION=1; born = 686 files, same seeds; 1S states in energy order; Hotelling on 3 states ~ chi^2(3)']
res = {}; same684 = tot684 = 0; maxv1 = 0.0; okD = True
for tag, name in PAIRS:
    pooled = 0.0; perP = {'1': 0.0, '2': 0.0}
    for B in ('0.01', '0.1'):
        for P in ('1', '2'):
            fs = sorted(glob.glob(f'{R}/nj/{tag}_B{B}_A{P}_*.err'))
            nj = [parse(f) for f in fs]; bo = [parse(f.replace(R + '/nj/', R6 + '/born/')) for f in fs]
            for a, b in zip(nj, bo):
                if 'line684' in a or 'line684' in b:
                    tot684 += 1; same684 += a.get('line684') == b.get('line684')
                if 'sc' in a and 'v1' in b: maxv1 = max(maxv1, max(abs(x - y) for x, y in zip(a['sc'], b['v1'])))
            n, mu, t2 = hot([[x - y for x, y in zip(r['sc'], r['qm'])] for r in nj if 'sc' in r])
            nb, mub, t2b = hot([[x - y for x, y in zip(r['sc'], r['qm'])] for r in bo if 'sc' in r])
            pooled += t2; perP[P] += t2
            inj = rms([r['imp'] for r in nj if 'imp' in r]); ib = rms([r['imp'] for r in bo if 'imp' in r])
            okD &= inj <= 2*ib
            res[(tag, B, P)] = (t2, inj/ib)
            out.append(f'{name:12s} B {B:>4} T PHEN {P} n {n} nj diff ' + ' '.join(f'{x:+.4f}' for x in mu)
                       + f' max {max(abs(mu)):.4f} T2 {t2:6.2f} p {chi2.sf(t2, 3):.3g} | born max {max(abs(mub)):.4f} T2 {t2b:6.2f} p {chi2.sf(t2b, 3):.3g}'
                       + f' | import rms nj {inj:.2f} born {ib:.2f} ueV ({inj/ib:.2f}x)')
    res[tag] = (pooled, perP)
    out.append(f'{name:12s} nj pooled chi^2(12) {pooled:.1f} p {chi2.sf(pooled, 12):.3g}; PHEN 1 chi^2(6) {perP["1"]:.1f} p {chi2.sf(perP["1"], 6):.3g};'
               f' PHEN 2 chi^2(6) {perP["2"]:.1f} p {chi2.sf(perP["2"], 6):.3g}')
ctl = sorted(glob.glob(f'{R}/ctl/*.err'))
cs = sum(open(f, 'rb').read() == open(f.replace(R + '/ctl/', R5 + '/born/'), 'rb').read() for f in ctl)
out.append(f'H687a: QM1S684 lines nj = born 686 {same684}/{tot684}; max |sc(nj) - v1(born)| {maxv1:.2e}; 1.0329 T nj = born 685 byte for byte {cs}/{len(ctl)}:'
           f' {same684 == tot684 and maxv1 <= 1e-9 and cs == len(ctl)}')
out.append('H687b (nj pooled p > 0.05 both pairs): ' + str(all(chi2.sf(res[t][0], 12) > 0.05 for t, _ in PAIRS))
           + '  ' + ', '.join(f'{n} p {chi2.sf(res[t][0], 12):.3g}' for t, n in PAIRS))
a2 = res['antiproton'][1]['2']
out.append(f'H687c (antihydrogen PHEN 2 nj chi^2(6) < 19.1 and p > 0.05): {a2 < 19.1 and chi2.sf(a2, 6) > 0.05}  chi^2(6) {a2:.1f} p {chi2.sf(a2, 6):.3g}')
out.append(f'H687d (import nj <= 2 x born every cell; False = the j rule is needed for energy): {okD}  ratios '
           + ', '.join(f'{k[0]} {k[1]} P{k[2]} {v[1]:.2f}' for k, v in res.items() if len(k) == 3))
# H687e on the 686 born files
out.append('H687e: classical j = 3/2 frequency against mean Born P(j = 3/2) of the photon-time state (686 born files)')
for tag, name in PAIRS:
    for B in ('0.01', '0.1'):
        for P in ('1', '2'):
            rows = [r for r in (parse(f) for f in sorted(glob.glob(f'{R6}/born/{tag}_B{B}_A{P}_*.err'))) if 'jcl' in r]
            p3 = [r['kept'] if r['jcl'] == 1.5 else 1 - r['kept'] for r in rows]
            c = sum(r['jcl'] == 1.5 for r in rows); e = sum(p3); v = sum(x*(1 - x) for x in p3)
            out.append(f'{name:12s} B {B:>4} T PHEN {P} n {len(rows)} jcl 3/2 {c} Born {e:.1f} z {(c - e)/math.sqrt(v):+.2f}')
text = '\n'.join(out); print(text)
os.makedirs('out', exist_ok=True); open('out/ana687.txt', 'w').write(text + '\n')

# Audit 692: step convergence of the spin-orbit transport (CREM_SPIN_SUBSTEP 0.025 / 0.1 against the default 0.05), paired by seed.
# usage: python3 ana692.py <692 runs dir with s0.025/ s0.1/ s0.05/> <686 runs dir with born/> <690 runs dir with born/>
import sys, glob, math, os
import numpy as np
from scipy.stats import chi2
R, R6, R9 = sys.argv[1:4]
UEV = 1.602176634e-25
PAIRS = (('antiproton', 'antihydrogen'), ('proton', 'hydrogen'))


def vec(t, key, n=4):
    i = t.index(key); return [float(x) for x in t[i+1:i+1+n]]


def parse(path):
    d = {}
    for line in open(path):
        t = line.split()
        if line.startswith('PAIR QM1S684'): d['v3'] = vec(t, 'v3')
        elif line.startswith('PAIR QM1S path'): d['sc'], d['qm'] = vec(t, 'sc'), vec(t, 'qm')
        elif line.startswith('PAIR BORNDM'): d['dE'] = float(t[t.index('dE')+1])/UEV; d['born'] = t[t.index('born')+1]
        elif line.startswith('PAIR j projection jcl'): d['jcl'] = t[t.index('jcl')+1]
        elif line.startswith('PAIR eigen ledger'):
            e = {t[i]: float(t[i+1]) for i in range(3, len(t)-1, 2)}; d['imp'] = (e['B'] + e['S2P'] - e['S1S'])/UEV
        elif line.startswith('END'): d['end'] = (t[t.index('stop')+1], t[t.index('photons')+1])
    return d


def ref(B, name): return (R6 if B == '0.01' else R9) + '/born/' + name


def paired(D):
    D = np.array(D); n = len(D); mu = D.mean(0); S = np.cov(D[:, :3], rowvar=False)
    t2 = float(n*mu[:3] @ np.linalg.solve(S, mu[:3])) if np.linalg.det(S) > 0 else float('nan')
    return mu, t2


def rms(xs): return math.sqrt(sum(x*x for x in xs)/len(xs))


out = []
ctl = sorted(glob.glob(f'{R}/s0.05/*.err')); same = sum(open(f, 'rb').read() == open(ref('1.0329', os.path.basename(f)), 'rb').read() for f in ctl)
out.append(f'H692a: CREM_SPIN_SUBSTEP=0.05 = default byte for byte {same}/{len(ctl)}: {same == len(ctl) == 16}')
okB = okC = okD = okE = True; diffs = {}
out.append('paired (step - 0.05): (A) = v3 - QM and gap sc - QM mean differences per 1S state, Hotelling p; max |d sc| per photon; dE / import rms ratio; discrete changes')
for v in ('0.025', '0.1'):
    for tag, name in PAIRS:
        for B in ('0.01', '1.0329'):
            for P in ('1', '2'):
                fs = sorted(glob.glob(f'{R}/s{v}/{tag}_B{B}_A{P}_*.err'))
                a = [parse(f) for f in fs]; b = [parse(ref(B, os.path.basename(f))) for f in fs]
                dA = [[x - y - (z - w) for x, y, z, w in zip(p['v3'], p['qm'], q['v3'], q['qm'])] for p, q in zip(a, b)]
                dG = [[x - y - (z - w) for x, y, z, w in zip(p['sc'], p['qm'], q['sc'], q['qm'])] for p, q in zip(a, b)]
                mA, tA = paired(dA); mG, tG = paired(dG)
                pA = chi2.sf(tA, 3) if not math.isnan(tA) else 1.0; pG = chi2.sf(tG, 3) if not math.isnan(tG) else 1.0
                mx = max(max(abs(x) for x in r) for r in dG)
                rde = rms([p['dE'] for p in a])/rms([q['dE'] for q in b]); rim = rms([p['imp'] for p in a])/rms([q['imp'] for q in b])
                disc = sum(p.get('born') != q.get('born') or p.get('jcl') != q.get('jcl') or p['end'] != q['end'] for p, q in zip(a, b))
                okB &= max(abs(mA)) < 0.0065 and pA > 0.05; okC &= max(abs(mG)) < 0.0065 and pG > 0.05
                okD &= abs(rde - 1) < 0.05 and abs(rim - 1) < 0.05; okE &= disc <= 2
                diffs[(v, tag, B, P)] = mG
                out.append(f'{v:5s} {name:12s} B {B:>6} T PHEN {P} n {len(fs)} (A) ' + ' '.join(f'{x:+.1e}' for x in mA) + f' p {pA:.2g} | gap '
                           + ' '.join(f'{x:+.1e}' for x in mG) + f' p {pG:.2g} | max|d sc| {mx:.1e} | dE {rde:.4f} imp {rim:.4f} | discrete {disc}')
out.append(f'H692b ((A) |mean d| < 0.0065 and p > 0.05 every cell): {okB}')
out.append(f'H692c (gap |mean d| < 0.0065 and p > 0.05 every cell): {okC}')
out.append(f'H692d (dE and import rms within 5 %): {okD}')
out.append(f'H692e (discrete changes <= 2 of 32 per cell): {okE}')
rat = []
for tag, _ in PAIRS:
    for B in ('0.01', '1.0329'):
        for P in ('1', '2'):
            x = np.abs(diffs[('0.1', tag, B, P)]).max(); y = np.abs(diffs[('0.025', tag, B, P)]).max()
            if y > 0: rat.append(x/y)
out.append('report: max |mean gap difference| ratio (0.1 - 0.05)/(0.025 - 0.05) per cell: ' + ' '.join(f'{r:.1f}' for r in rat)
           + ' (2 for first order, 4 for second order in the step)')
text = '\n'.join(out); print(text)
os.makedirs('out', exist_ok=True); open('out/ana692.txt', 'w').write(text + '\n')

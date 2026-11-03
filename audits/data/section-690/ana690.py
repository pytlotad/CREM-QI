# Audit 690: 1S distribution of the Born-Delta-m path (default since 688) and the old basin path (CREM_NO_PAIR_BORN_DELTA_M=1) against QM at 1.0329 T,
# 128 photons per cell, hydrogen and antihydrogen.  Hotelling T^2 on the first three 1S states per cell ~ chi^2(3), a pair's 2 cells ~ chi^2(6).
# usage: python3 ana690.py <690 runs dir with def/ born/ prov/> <684 runs dir> <685 runs dir>
import sys, glob, math, os
import numpy as np
from scipy.stats import chi2
R, R4, R5 = sys.argv[1:4]
UEV = 1.602176634e-25
PAIRS = (('antiproton', 'antihydrogen'), ('proton', 'hydrogen'))


def vec(t, key, n=4):
    i = t.index(key); return [float(x) for x in t[i+1:i+1+n]]


def parse(path):
    d = {}
    for line in open(path):
        t = line.split()
        if line.startswith('PAIR BORNDM'):
            d['dJz'] = float(t[t.index('dJz')+1]); d['dE'] = float(t[t.index('dE')+1])/UEV
            d['chg'] = t[t.index('born')+1] != t[t.index('classical')+1]
        elif line.startswith('PAIR QM1S path'):
            d['sc'], d['qm'] = vec(t, 'sc'), vec(t, 'qm'); d['m'] = [int(x) for x in vec(t, 'm')]
        elif line.startswith('PAIR eigen ledger'):
            e = {t[i]: float(t[i+1]) for i in range(3, len(t)-1, 2)}
            d['imp'] = (e['B'] + e['S2P'] - e['S1S'])/UEV
    return d


def rms(xs): return math.sqrt(sum(x*x for x in xs)/len(xs)) if xs else float('nan')


out = ['per cell: n, mean P_sc - P_QM (4 states, energy order), max |.| (2 sigma), Hotelling T^2 (3 states) and p; import rms [ueV]']
res = {}
for tag, name in PAIRS:
    for path in ('def', 'born'):
        T = 0.0
        for B in ('1.0329',):
            for P in ('1', '2'):
                rows = [parse(f) for f in sorted(glob.glob(f'{R}/{path}/{tag}_B{B}_A{P}_*.err'))]
                ph = [r for r in rows if 'sc' in r]; n = len(ph)
                D = np.array([[r['sc'][i] - r['qm'][i] for i in range(4)] for r in ph])
                mu = D.mean(0); se2 = 2*D.std(0, ddof=1)/math.sqrt(n)
                X = D[:, :3]; S = np.cov(X, rowvar=False)
                t2 = float(n*mu[:3] @ np.linalg.solve(S, mu[:3])); T += t2
                i = int(np.argmax(np.abs(mu)))
                imp = rms([r['imp'] for r in rows if 'imp' in r])
                extra = ''
                if path == 'born':
                    dj = [abs(r['dJz']) for r in rows if 'dJz' in r]; de = [r['dE'] for r in rows if 'dE' in r]
                    both = [(r['imp'], r['dE']) for r in rows if 'imp' in r and 'dE' in r]
                    extra = (f'  changed {sum(r.get("chg", False) for r in rows)}/{len(dj)} mean|dJz| {sum(dj)/len(dj):.3f} hbar rms dE {rms(de):.3f}'
                             f'  imp+dE {rms([a + b for a, b in both]):.2f} imp-dE {rms([a - b for a, b in both]):.2f}')
                res[(tag, path, B, P)] = dict(n=n, mx=abs(mu[i]), s=se2[i], t2=t2, p=float(chi2.sf(t2, 3)), imp=imp)
                out.append(f'{name:12s} {path:4s} B {B:>4} T PHEN {P} n {n:3d} m {ph[0]["m"]} diff ' + ' '.join(f'{x:+.4f}' for x in mu)
                           + f'  max {abs(mu[i]):.4f} ({se2[i]:.4f})  T2 {t2:6.2f} p {chi2.sf(t2, 3):.3g}  import {imp:.2f}' + extra)
        res[(tag, path)] = (T, float(chi2.sf(T, 6)))
        out.append(f'{name:12s} {path:4s} pooled chi^2(6) {T:.1f}  p {chi2.sf(T, 6):.3g}')
# H690a provenance
same = tot = 0
for f in sorted(glob.glob(f'{R}/prov/*.err')):
    b = os.path.basename(f); path, rest = b.split('_', 1)
    ref = (R5 + '/born/' if path == 'born' else R4 + '/def/') + rest
    tot += 1; same += open(f, 'rb').read() == open(ref, 'rb').read()
out.append(f'H690a: seeds 0-3 rerun = 685 born / 684 def byte for byte: {same}/{tot}')
out.append('H690b (born pooled p > 0.05 both pairs): ' + str(all(res[(t, 'born')][1] > 0.05 for t, _ in PAIRS))
           + '  ' + ', '.join(f'{n} p {res[(t, "born")][1]:.3g}' for t, n in PAIRS))
out.append('H690c (def pooled p < 0.01 both pairs): ' + str(all(res[(t, 'def')][1] < 0.01 for t, _ in PAIRS))
           + '  ' + ', '.join(f'{n} p {res[(t, "def")][1]:.3g}' for t, n in PAIRS))
cells = [k for k in res if len(k) == 4 and k[1] == 'born']
out.append('H690d (import rms born <= 1.2 x def every cell): ' + str(all(res[k]['imp'] <= 1.2*res[(k[0], 'def', k[2], k[3])]['imp'] for k in cells))
           + '  ratios ' + ', '.join(f'{k[0]} P{k[3]} {res[k]["imp"]/res[(k[0], "def", k[2], k[3])]["imp"]:.2f}' for k in cells))
text = '\n'.join(out); print(text)
os.makedirs('out', exist_ok=True); open('out/ana690.txt', 'w').write(text + '\n')

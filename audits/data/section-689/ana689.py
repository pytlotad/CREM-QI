# Audit 689 (print only): the 668/670 j rule at the photon -- V1 split by the j label, the classical j against the Born P(j = 3/2) p3.
# usage: python3 ana689.py <689 runs dir with born/ reg/ validation.txt> <686 runs dir with born/>
import sys, glob, math, os
import numpy as np
R, R6 = sys.argv[1:3]
PAIRS = (('antiproton', 'antihydrogen'), ('proton', 'hydrogen'))


def vec(t, key, n=4):
    i = t.index(key); return [float(x) for x in t[i+1:i+1+n]]


def parse(path):
    d = {}
    for line in open(path):
        t = line.split()
        if line.startswith('PAIR QM1S689'):
            d['jcl'] = float(t[t.index('jcl')+1]); d['p3'] = float(t[t.index('p3')+1]); d['q3'] = float(t[t.index('q3')+1])
            d['v12'] = vec(t, 'v12'); d['v32'] = vec(t, 'v32')
        elif line.startswith('PAIR j projection jcl'):
            d['jrule'] = float(t[t.index('jcl')+1]); d['kept'] = float(t[t.index('kept')+1])
        elif line.startswith('PAIR QM1S684'):
            d['v1'] = vec(t, 'v1')
        elif line.startswith('PAIR QM1S path'):
            d['sc'], d['qm'] = vec(t, 'sc'), vec(t, 'qm')
    return d


out = []
# H689a
same = tot = 0
for f in sorted(glob.glob(f'{R}/born/*.err')):
    a = ''.join(l for l in open(f) if not l.startswith('PAIR QM1S689'))
    tot += 1; same += a == open(f.replace(R + '/born/', R6 + '/born/')).read()
reg = 0
for tag in ('1', '2', 'L2', '5Ps', '5H', '6'):
    reg += all(open(f'{R}/reg/old/{tag}/{x}', 'rb').read() == open(f'{R}/reg/new/{tag}/{x}', 'rb').read() for x in ('out.txt', 'err.txt'))
val = ' '.join(l.strip() for l in open(f'{R}/validation.txt') if 'Enforced checks' in l or 'EXIT' in l)
out.append(f'H689a: files without QM1S689 = born 686 {same}/{tot}; regression dedcb03f 6 experiments identical {reg}/6; validation {val}')
# H689b, c, d
maxid = 0.0; nid = 0; jmis = 0; out_c = []; out_d = []; okC = True; dd = []
for tag, name in PAIRS:
    for B in ('0.01', '0.1'):
        for P in ('1', '2'):
            rows = [r for r in (parse(f) for f in sorted(glob.glob(f'{R}/born/{tag}_B{B}_A{P}_*.err'))) if 'p3' in r and 'sc' in r]
            I = np.array([r['jcl'] == 1.5 for r in rows], float); p3 = np.array([r['p3'] for r in rows])
            D = np.array([[a - b for a, b in zip(r['v32'], r['v12'])] for r in rows])
            G = np.array([[a - b for a, b in zip(r['sc'], r['v1'])] for r in rows])
            for r, i, p, dr, g in zip(rows, I, p3, D, G):
                jmis += r['jrule'] != r['jcl']
                if (i == 1 and p == 0) or (i == 0 and p == 1): continue
                nid += 1; maxid = max(maxid, float(np.max(np.abs((i - p)*dr - g))))
            n = len(rows); x = I - p3
            freq = x.mean()*D.mean(0); corr = (x[:, None]*D).mean(0) - freq; gm = G.mean(0); se = 2*G.std(0, ddof=1)/math.sqrt(n)
            s = int(np.argmax(np.abs(gm)))
            if P == '2': okC &= abs(corr[s]) > abs(freq[s])
            out_c.append(f'{name:12s} B {B:>4} T PHEN {P} n {n} sc - v1 ' + ' '.join(f'{v:+.4f}' for v in gm) + f' (2 sigma {se[s]:.4f})'
                         + f' | state {s}: frequency {freq[s]:+.4f} correlation {corr[s]:+.4f}')
            m1 = p3[I == 1].mean(); m0 = p3[I == 0].mean(); sd = math.sqrt(p3[I == 1].var(ddof=1)/sum(I) + p3[I == 0].var(ddof=1)/sum(1 - I))
            m = p3.mean(); born = p3.var()/(m*(1 - m)); dd.append((m1 - m0, sd, born))
            out_d.append(f'{name:12s} B {B:>4} T PHEN {P} E[p3|jcl 3/2] {m1:.4f} (n {int(sum(I))}) E[p3|jcl 1/2] {m0:.4f} diff {m1 - m0:+.4f} z {(m1 - m0)/sd:+.2f}'
                         f' | Born draw {born:.4f} | mean p3 {m:.4f} q3 {np.mean([r["q3"] for r in rows]):.4f}')
out.append(f'H689b: sc - v1 = (I - p3)(v32 - v12) in {nid} photons, max |residual| {maxid:.2e}; jcl of QM1S689 != jcl of the j rule line in {jmis}: {maxid <= 1e-8 and jmis == 0}')
out.append('H689c: mean sc - v1 = frequency (mean(I - p3) mean D) + correlation (cov(I - p3, D)), state of the largest |mean|')
out += out_c
out.append(f'H689c (PHEN 2: |correlation| > |frequency| in every cell): {okC}')
out.append('H689d: E[p3 | jcl] difference against 0 and against a Born draw of j (Var p3 / (m (1 - m)))')
out += out_d
w = sum(1/s**2 for _, s, _ in dd); md = sum(d/s**2 for d, s, _ in dd)/w
out.append(f'H689d pooled (inverse variance): diff {md:+.4f} +- {1/math.sqrt(w):.4f} (z {md*math.sqrt(w):+.2f}); Born draw mean {np.mean([b for _, _, b in dd]):.4f}')
text = '\n'.join(out); print(text)
os.makedirs('out', exist_ok=True); open('out/ana689.txt', 'w').write(text + '\n')

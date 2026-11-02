# Audit 683, POST HOC (after ana683.py): (1) the 666 ledger on the formation path sets Hbar = E2P of the label (B = 0), so the Born part moves into S2P;
# compare the whole 2P side B + S2P - Zsb and the whole import instead.  (2) H683c pooled: z of every (field, PHEN, state) cell and its chi^2.
# usage: python3 post683.py <runs dir with def/ form/>
import sys, glob, math
R = sys.argv[1]
UEV = 1.602176634e-25
FIELDS = ('0.01', '0.1', '1.0329', '3')


def parse(path):
    qm = sc = led = None
    for line in open(path):
        t = line.split()
        if line.startswith('PAIR QM1S'):
            i = t.index('qm'); qm = [float(x) for x in t[i+1:i+5]]
            i = t.index('sc'); sc = [float(x) for x in t[i+1:i+5]]
        elif line.startswith('PAIR eigen ledger'):
            d = {t[i]: float(t[i+1]) for i in range(3, len(t)-1, 2)}
            led = ((d['B'] + d['S2P'] - d['Zsb'])/UEV, (d['S1S'] - d['Zsa'])/UEV, (d['B'] + d['S2P'] - d['S1S'])/UEV)
    return qm, sc, led


def rms(xs): return math.sqrt(sum(x*x for x in xs)/len(xs)) if xs else float('nan')


out = ['rms over photons [ueV]: 2P side B + S2P* (S2P* = S2P - Zsb), 1S side S1S*, whole import B + S2P - S1S; def / form']
z = []
for B in FIELDS:
    for P in ('1', '2'):
        line = f'B {B:>6} T PHEN {P}'
        for path in ('def', 'form'):
            rows = [parse(f) for f in sorted(glob.glob(f'{R}/{path}/B{B}_A{P}_*.err'))]
            led = [r[2] for r in rows if r[2]]
            line += f'  {path:4s} n {len(led):2d} 2P {rms([l[0] for l in led]):7.2f} 1S {rms([l[1] for l in led]):6.2f} imp {rms([l[2] for l in led]):7.2f}'
            if path == 'form':
                rr = [r for r in rows if r[0]]; n = len(rr)
                for i in range(4):
                    d = [r[1][i] - r[0][i] for r in rr]; mu = sum(d)/n
                    s = math.sqrt(sum((x - mu)**2 for x in d)/(n - 1)/n)
                    if s > 0: z.append((B, P, i, mu/s))
        out.append(line)
chi2 = sum(v[3]**2 for v in z)
out.append(f'H683c pooled (formation path): {len(z)} cells, |z| > 2 in {sum(abs(v[3]) > 2 for v in z)} (expected {0.0455*len(z):.1f}),'
           f' chi^2 {chi2:.1f} for {len(z)} dof (cells correlated: the 4 states of a photon sum to 1)')
out.append('  |z| > 2: ' + ', '.join(f'{b} T P{p} state {i} z {v:+.2f}' for b, p, i, v in z if abs(v) > 2))
text = '\n'.join(out); print(text)
open('out/post683.txt', 'w').write(text + '\n')

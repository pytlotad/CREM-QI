# Audit 684: the default path's 1S distribution decomposed -- V0 (default) - V1 = j rule (C), V1 - V2 = 1S from the m_L block (D),
# V2 - V3 = weights in the m_L block (B), V3 - V4 = photon-time state instead of the formation state (A); photon block iL against Born pB.
# usage: python3 ana684.py <runs dir> <683 runs dir>
import sys, glob, math, os
R, R3 = sys.argv[1], sys.argv[2]
FIELDS = ('0.01', '0.1', '1.0329')


def vec(t, key, n=4):
    i = t.index(key); return [float(x) for x in t[i+1:i+1+n]]


def parse(path):
    d = {}
    for line in open(path):
        t = line.split()
        if line.startswith('PAIR QM1S684'):
            d['v1'], d['v2'], d['v3'] = vec(t, 'v1'), vec(t, 'v2'), vec(t, 'v3')
            d['iL'] = int(t[t.index('iL')+1]); d['pB'] = vec(t, 'pB', 3)
        elif line.startswith('PAIR QM1S path'):
            d['v0'], d['v4'] = vec(t, 'sc'), vec(t, 'qm'); d['m'] = [int(x) for x in vec(t, 'm')]
    return d


def keep(path):  # lines compared with 683
    return [l for l in open(path) if l.startswith('END') or l.startswith('PAIR eigen ledger') or l.startswith('PAIR QM1S path')]


out = ['mean over photons, 1S states in energy order; components max_i |mean difference| (2 sigma of the mean of that difference)']
same = tot = 0
okB = True; okC = True; chi2 = 0.0; dof = 0; blocks = []
for B in FIELDS:
    for P in ('1', '2'):
        files = sorted(glob.glob(f'{R}/def/B{B}_A{P}_*.err'))
        for f in files:
            tot += 1; same += keep(f) == keep(f.replace(R + '/def/', R3 + '/def/'))
        rows = [d for d in (parse(f) for f in files) if 'v1' in d]; n = len(rows)
        mean = {k: [sum(r[k][i] for r in rows)/n for i in range(4)] for k in ('v0', 'v1', 'v2', 'v3', 'v4')}

        def comp(a, b):
            d = [[r[a][i] - r[b][i] for r in rows] for i in range(4)]
            mu = [sum(x)/n for x in d]
            sd = [2*math.sqrt(sum((x - m)**2 for x in xs)/(n - 1)/n) for xs, m in zip(d, mu)]
            i = max(range(4), key=lambda j: abs(mu[j])); return mu, abs(mu[i]), sd[i]
        out.append(f'B {B:>6} T PHEN {P} n {n} m {rows[0]["m"]}')
        for k in ('v0', 'v1', 'v2', 'v3', 'v4'):
            out.append(f'    {k.upper()} ' + ' '.join(f'{x:.4f}' for x in mean[k]))
        tot_gap = comp('v0', 'v4')
        parts = {}
        for name, a, b in (('total V0-V4', 'v0', 'v4'), ('(C) j rule V0-V1', 'v0', 'v1'), ('(D) 1S block V1-V2', 'v1', 'v2'),
                           ('(B) block weights V2-V3', 'v2', 'v3'), ('(A) state V3-V4', 'v3', 'v4'), ('(B)+(D) V1-V3', 'v1', 'v3')):
            mu, mx, s = comp(a, b); parts[name] = mx
            out.append(f'    {name:24s} ' + ' '.join(f'{x:+.4f}' for x in mu) + f'   max {mx:.4f} ({s:.4f})')
        if B in ('0.01', '0.1') and P == '1':
            okB &= parts['(B)+(D) V1-V3'] >= 0.6*parts['total V0-V4']
        okC &= parts['(A) state V3-V4'] <= 0.02
        cnt = [sum(r['iL'] == i for r in rows) for i in range(3)]
        exp = [sum(r['pB'][i] for r in rows) for i in range(3)]
        var = [sum(r['pB'][i]*(1 - r['pB'][i]) for r in rows) for i in range(3)]
        for i in range(3):
            if var[i] > 0: chi2 += (cnt[i] - exp[i])**2/var[i]; dof += 1
        out.append('    photon block m_L +1/0/-1: count ' + ' '.join(map(str, cnt)) + '  Born ' + ' '.join(f'{x:.1f}' for x in exp))
p = math.exp(-chi2/2) * sum((chi2/2)**k/math.factorial(k) for k in range(dof//2)) if dof % 2 == 0 else float('nan')
out.append(f'H684a: END, ledger and QM1S path lines = 683 def: {same}/{tot}')
ps = sum(open(f, 'rb').read() == open(f.replace('/on_', '/off_'), 'rb').read() for f in sorted(glob.glob(f'{R}/ps/on_*')))
out.append(f'H684a: Ps with CREM_PAIR_QM1S=1 = without, byte for byte: {ps}/{len(glob.glob(R + "/ps/on_*"))} files')
out.append(f'H684b ((B)+(D) >= 0.6 of the gap at 0.01 / 0.1 T PHEN 1): {okB}')
out.append(f'H684c ((A) <= 0.02 at every field and PHEN): {okC}')
out.append(f'H684d (photon block = Born): chi^2 {chi2:.1f} for {dof} cells (blocks of a photon sum to 1), p {p:.3g}: {p > 0.05}')
text = '\n'.join(out); print(text)
os.makedirs('out', exist_ok=True); open('out/ana684.txt', 'w').write(text + '\n')

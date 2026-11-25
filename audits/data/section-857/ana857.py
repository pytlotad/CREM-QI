# Audit 857: dipole analogue of Peters as the integrator check.  Usage: LC_ALL=C python3 ana857.py -> out/ana857.txt.  Written before any result.
import os, re, math
SP = '/tmp/claude-1000/-home-teddy-Projekty-CREM-QI/4fa83d5c-8878-4cfd-898e-1a9b53a01ab2/scratchpad'
R = os.environ.get('R857', SP + '/r857'); pr = lambda *a: print(*a, flush=True)
rd = lambda v, t, x: open(f'{R}/reg/{v}/{t}/{x}', 'rb').read()
same = lambda a, b, t: all(rd(a, t, x) == rd(b, t, x) for x in ('out.txt', 'err.txt', 'rc.txt'))
SETS = ('1', '2', 'L2', '5Ps', '5H', '6'); NEW = (b'  closed form, eccentric ', b'  mean / eccentric form  ')
for a, b in (('old', 'new'), ('old', 'new_legacy')):
    pr(f'  {a} vs {b}: ' + ', '.join(f"{t} {'SAME' if same(a, b, t) else 'DIFF'}" for t in SETS))
pr(f"H857o (new_legacy = old, 6 zestawow): {'TAK' if all(same('old', 'new_legacy', t) for t in SETS) else 'NIE'}")
strip = lambda t: b'\n'.join(l for l in rd('new', t, 'out.txt').split(b'\n') if not l.startswith(NEW))
z1 = all(same('old', 'new', t) for t in ('5Ps', '5H', '6'))
z2 = all(strip(t) == rd('old', t, 'out.txt') and rd('old', t, x) == rd('new', t, x) for t in ('1', '2', 'L2') for x in ('err.txt', 'rc.txt'))
added = {t: sum(1 for l in rd('new', t, 'out.txt').split(b'\n') if l.startswith(NEW)) for t in ('1', '2', 'L2')}
pr(f"H857z (5Ps, 5H, 6 bajt w bajt; 1, 2, L2 rozne tylko nowymi wierszami): {'TAK' if z1 and z2 else 'NIE'}  (z1 {z1}, z2 {z2}, nowe wiersze {added})")
for v in ('old', 'new'):
    s = open(f'{R}/validation_{v}.txt').read().strip().split('\n'); pr(f'  walidacja {v}: ' + ' | '.join(s[-4:]))
pr('H857v: ocena recznie z wierszy powyzej')
def info(t):
    s = rd('new', t, 'out.txt').decode()
    g = lambda pat: (re.search(pat, s).groups() if re.search(pat, s) else None)
    L = []
    for ln in rd('new', t, 'err.txt').decode(errors='replace').split('\n'):
        if ln.startswith('TRANSIT855'):
            w = ln.split(); L.append({w[i]: w[i + 1] for i in range(1, len(w) - 1, 2)})
    return dict(mean=g(r'mean of completed runs ([\d.e+-]+) \+/- ([\d.e+-]+) (\w+)'), circ=g(r'closed form\s+([\d.e+-]+) (\w+)  \(mu'),
                rc_=g(r'mean / closed form\s+([\d.e+-]+)'), ecc=g(r'closed form, eccentric ([\d.e+-]+) (\w+)'), re_=g(r'mean / eccentric form\s+([\d.e+-]+)'),
                life=g(r'lifetime at the final state ([\d.]+ \w+)'), rc=rd('new', t, 'rc.txt').decode().strip(), tr=L)
okA, okM, res = True, True, {}
for p, lab in (('1', 'p-Ps (PHEN=1)'), ('2', 'o-Ps (PHEN=2)')):
    for n in ('1', '2'):
        for kind in ('D', 'C'):
            I = info(f'{kind}_p{p}_n{n}'); res[(kind, p, n)] = I
            e0 = [float(d['e0']) for d in I['tr'] if 'e0' in d]; rt = [float(d['tlab']) / float(d['eref']) for d in I['tr'] if d.get('decision') == 'used']
            pr(f"  {kind} {lab} n = {n}: mean {I['mean']}, circular {I['circ']} (ratio {I['rc_']}), eccentric {I['ecc']} (ratio {I['re_']}), rc {I['rc']}, life {I['life']}")
            if e0: pr(f"    e0 min {min(e0):.6f} max {max(e0):.6f}; per-seed tlab/eref min {min(rt):.4f} max {max(rt):.4f} mean {sum(rt) / len(rt):.4f} (N {len(rt)})")
            if kind == 'C':
                if not (I['rc_'] and I['re_'] and abs(float(I['re_'][0]) / float(I['rc_'][0]) - 1) <= 1e-5 and e0 and max(e0) <= 1e-6): okA = False
            else:
                want = math.sqrt(1 - (int(n) - 0.5)**2 / int(n)**2)
                if not (e0 and all(abs(x - want) <= 1e-3 for x in e0)): okA = False
                if not (I['re_'] and abs(float(I['re_'][0]) - 1) <= 0.06): okM = False
pr(f"H857a (C: eccentric = circular, e0 = 0; D: e0 = sqrt(1 - (n - 1/2)^2/n^2) do 1e-3): {'TAK' if okA else 'NIE'}")
pr(f"H857m (GLOWNA; D: |mean / eccentric - 1| <= 0,06 w p-Ps i o-Ps, n = 1 i 2): {'TAK' if okM else 'NIE'}")
for t in ('1', '2', 'L2'):
    pr(f"  zgodnosc D z REG new: zestaw {t} out.txt = D: {rd('new', t, 'out.txt') == rd('new', {'1': 'D_p1_n1', '2': 'D_p2_n1', 'L2': 'D_p2_n2'}[t], 'out.txt')}")

# Audit 856: transit from the prepared start.  Usage: LC_ALL=C python3 ana856.py -> out/ana856.txt.  Written before any result.
import os, re
SP = '/tmp/claude-1000/-home-teddy-Projekty-CREM-QI/4fa83d5c-8878-4cfd-898e-1a9b53a01ab2/scratchpad'
R = os.environ.get('R856', SP + '/r856'); pr = lambda *a: print(*a, flush=True)
rd = lambda v, t, x: open(f'{R}/reg/{v}/{t}/{x}', 'rb').read()
same = lambda a, b, t: all(rd(a, t, x) == rd(b, t, x) for x in ('out.txt', 'err.txt', 'rc.txt'))
SETS = ('1', '2', 'L2', '5Ps', '5H', '6')
for a, b in (('old', 'new'), ('old', 'new_legacy')):
    pr(f'  {a} vs {b}: ' + ', '.join(f"{t} {'SAME' if same(a, b, t) else 'DIFF'} (rc {rd(a, t, 'rc.txt').decode().strip()}/{rd(b, t, 'rc.txt').decode().strip()})" for t in SETS))
pr(f"H856o (new_legacy = old, 6 zestawow): {'TAK' if all(same('old', 'new_legacy', t) for t in SETS) else 'NIE'}")
pr(f"H856z (new = old w 1, 2, 5Ps, 5H, 6): {'TAK' if all(same('old', 'new', t) for t in SETS if t != 'L2') else 'NIE'}")
for v in ('old', 'new'):
    s = open(f'{R}/validation_{v}.txt').read().strip().split('\n')
    pr(f'  walidacja {v}: ' + ' | '.join(s[-4:]))
pr('H856v: ocena recznie z wierszy powyzej (72/72 i EXIT 0 w new)')
def coll(v, t):
    s = rd(v, t, 'out.txt').decode()
    u = re.search(r'over (\d+) usable trajectories \((\d+) observed, (\d+) right-censored\)\.\s+Completion fraction: ([\d.]+)%', s)
    m = re.search(r'mean of completed runs ([\d.e+-]+) \+/- ([\d.e+-]+) (\w+)', s)
    c = re.search(r'closed form\s+([\d.e+-]+) (\w+)', s); r = re.search(r'mean / closed form\s+([\d.e+-]+)', s)
    lt = re.search(r'lifetime at the final state ([\d.]+ \w+)', s); et = re.search(r'E\[T\] from the preparation \(cascade \+ annihilation\) ([\d.]+ \w+)', s)
    hd = re.search(r'Collapse time = transit from (.*?) to 0.005', s)
    return dict(usable=u.groups() if u else None, mean=m.groups() if m else None, closed=c.groups() if c else None, ratio=float(r.group(1)) if r else None,
                life=lt.group(1) if lt else None, et=et.group(1) if et else None, head=hd.group(1) if hd else None,
                rc=rd(v, t, 'rc.txt').decode().strip(), wall=rd(v, t, 'wall.txt').decode().strip())
for v in ('old', 'new'):
    pr(f'  L2 {v}: {coll(v, "L2")}')
n = coll('new', 'L2')
m = bool(n['usable'] and n['usable'][0] == '16' and n['usable'][1] == '16' and n['rc'] == '0' and n['mean'] and n['ratio'] is not None and 0.25 <= n['ratio'] <= 0.32)
pr(f"H856m (GLOWNA; L2 new usable 16 / observed 16, kod 0, srednia / wzor w [0,25; 0,32]): {'TAK' if m else 'NIE'}")
for t in ('1', '2'): pr(f'  zestaw {t} new: {coll("new", t)}')
ok = True
for p, lab in (('1', 'p-Ps (PHEN=1)'), ('2', 'o-Ps (PHEN=2)')):
    for nn in ('1', '2'):
        c = coll('new', f'C_p{p}_n{nn}'); ok = ok and c['ratio'] is not None and 0.8 <= c['ratio'] <= 1.2
        es = []
        for ln in rd('new', f'C_p{p}_n{nn}', 'err.txt').decode(errors='replace').split('\n'):
            if ln.startswith('TRANSIT855'):
                w = ln.split(); d = {w[i]: w[i + 1] for i in range(1, len(w) - 1, 2)}
                es.append((d['decision'], round(float(d['a0ratio']), 6), round(float(d['tlab']) / float(d['ref_start']), 4)))
        pr(f'  C {lab} n = {nn} (--circular-start): ratio {c["ratio"]}, mean {c["mean"]}, closed {c["closed"]}, usable {c["usable"]}, rc {c["rc"]}, head {c["head"]}')
        pr(f'    TRANSIT855 {es}')
pr(f"H856c (start kolowy: srednia / wzor w [0,8; 1,2], p-Ps i o-Ps, n = 1 i 2): {'TAK' if ok else 'NIE'}")

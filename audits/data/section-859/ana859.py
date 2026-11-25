# Audit 859: jump cap at large eccentricity.  Usage: LC_ALL=C python3 ana859.py -> out/ana859.txt.  Written before any result.
import os, re
SP = '/tmp/claude-1000/-home-teddy-Projekty-CREM-QI/4fa83d5c-8878-4cfd-898e-1a9b53a01ab2/scratchpad'
R = os.environ.get('R859', SP + '/r859'); pr = lambda *a: print(*a, flush=True)
rd = lambda v, t, x: open(f'{R}/reg/{v}/{t}/{x}', 'rb').read()
same = lambda a, b, t: all(rd(a, t, x) == rd(b, t, x) for x in ('out.txt', 'err.txt', 'rc.txt'))
SETS = ('1', '2', 'L2', '5Ps', '5H', '6')
for a, b in (('old', 'new'), ('new', 'new_on')):
    pr(f'  {a} vs {b}: ' + ', '.join(f"{t} {'SAME' if same(a, b, t) else 'DIFF'} (rc {rd(a, t, 'rc.txt').decode().strip()}/{rd(b, t, 'rc.txt').decode().strip()},"
                                    f" wall {rd(a, t, 'wall.txt').decode().strip()}/{rd(b, t, 'wall.txt').decode().strip()} s)" for t in SETS))
pr(f"H859k (new = old, 6 zestawow): {'TAK' if all(same('old', 'new', t) for t in SETS) else 'NIE'}")
for v in ('new', 'new_on'):
    s = open(f'{R}/validation_{v}.txt').read().strip().split('\n'); pr(f'  walidacja {v}: ' + ' | '.join(s[-4:]))
pr('H859v: ocena recznie z wierszy powyzej (new)')
def info(v, t):
    s = rd(v, t, 'out.txt').decode()
    g = lambda pat: (re.search(pat, s).groups() if re.search(pat, s) else None)
    T, O = [], {}
    for ln in rd(v, t, 'err.txt').decode(errors='replace').split('\n'):
        w = ln.split(); d = {w[i]: w[i + 1] for i in range(1, len(w) - 1, 2)}
        if ln.startswith('TRANSIT855'): T.append(d)
        elif ln.startswith('ORBIT858'): O.setdefault(d['seed'], []).append(d)
    lab = s[s.find('Laboratory observables'):] if 'Laboratory observables' in s else ''
    return dict(mean=g(r'mean of completed runs ([\d.e+-]+) \+/- ([\d.e+-]+) (\w+)'), re_=g(r'mean / eccentric form\s+([\d.e+-]+)'),
                rc_=g(r'mean / closed form\s+([\d.e+-]+)'), usable=g(r'over (\d+) usable trajectories \((\d+) observed, (\d+) right-censored\)'),
                life=re.findall(r'lifetime at the final state [^\n]*', s), et=re.findall(r'E\[T\] from the preparation[^\n]*', s), lab=lab,
                rc=rd(v, t, 'rc.txt').decode().strip(), wall=rd(v, t, 'wall.txt').decode().strip(), tr=T, orb=O)
# H859p / H859w
okP = True
for t in ('1', '2', 'L2'):
    a, b = info('new', t), info('new_on', t)
    eq = a['life'] == b['life'] and a['et'] == b['et'] and a['lab'] == b['lab']
    if t == 'L2' and not eq: okP = False
    pr(f"  zestaw {t}: OFF usable {a['usable']} rc {a['rc']} | ON usable {b['usable']} rc {b['rc']}; zycie/E[T]/Laboratory ON = OFF: {eq}")
    pr(f"    OFF {a['life']} {a['et']}; ON {b['life']} {b['et']}; mean OFF {a['mean']} ON {b['mean']}; ratio ecc OFF {a['re_']} ON {b['re_']}")
okP = okP and all(same('new', 'new_on', t) for t in ('5Ps', '5H', '6'))
pr(f"H859p (L2 zycie / E[T] / Laboratory ON = OFF; 5Ps, 5H, 6 ON = OFF bajt w bajt): {'TAK' if okP else 'NIE'}")
w = [info('new_on', t) for t in ('1', '2')]
pr(f"H859w (ON, budzet 90 s, zestawy 1 i 2: przejscie cenzurowane w wiekszosci ziaren, kod 2): {'TAK' if all(x['rc'] == '2' for x in w) else 'NIE'}  "
   f"(kody {[x['rc'] for x in w]}, decyzje {[sorted(set(d.get('decision', '?') for d in x['tr'])) for x in w]}, outcome {[sorted(set(d.get('outcome', '?') for d in x['tr'])) for x in w]})")
# H859c
REF857 = {'1': 0.998765, '2': 1.00004}
okC = True
for p in ('1', '2'):
    I = info('new_on', f'C_p{p}_n2'); r = float(I['rc_'][0]) if I['rc_'] else float('nan')
    if not abs(r - REF857[p]) <= 1e-6: okC = False
    pr(f"  C p{p} n2 ON: mean {I['mean']} mean/closed {r} (857 {REF857[p]}), rc {I['rc']}")
pr(f"H859c (start kolowy ON = 857 do 1e-6): {'TAK' if okC else 'NIE'}")
# H859m, H859e, H859s
REF = {('1', '2'): 0.91545, ('2', '2'): 0.915939, ('1', '1'): 0.831682, ('2', '1'): 0.890971}
SPREAD857 = {'1': 0.336, '2': 0.403}
okM, okE, okS = True, True, True
for p, lab in (('1', 'p-Ps (PHEN=1)'), ('2', 'o-Ps (PHEN=2)')):
    for n in ('2', '1'):
        I = info('new_on', f'D_p{p}_n{n}')
        rt = [float(d['tlab']) / float(d['eref']) for d in I['tr'] if d.get('decision') == 'used']
        re_ = float(I['re_'][0]) if I['re_'] else float('nan')
        walls = [float(d['wall']) for d in I['tr'] if 'wall' in d]
        pr(f"  D {lab} n = {n} ON: mean {I['mean']}, mean/eccentric {re_} (OFF 857 {REF[(p, n)]}), usable {I['usable']}, rc {I['rc']}, wall {I['wall']} s, life {I['life']}")
        if rt: pr(f"    per-seed tlab/eref min {min(rt):.4f} max {max(rt):.4f} rozrzut {max(rt) - min(rt):.4f} (N {len(rt)}); wall przejscia {min(walls):.0f}-{max(walls):.0f} s; outcome {sorted(set(d['outcome'] for d in I['tr']))}")
        if n == '2':
            if not abs(re_ - 1) <= 0.02: okM = False
            dev, nck = [], []
            for s_, L in I['orb'].items():
                nck.append(len(L)); dev += [abs(float(d['e1']) / float(d['einv']) - 1) for d in L if float(d['einv']) >= 0.05]
            if not (dev and max(dev) <= 0.05): okE = False
            pr(f"    ORBIT858: ziaren {len(I['orb'])}, punktow {min(nck) if nck else 0}-{max(nck) if nck else 0}, max |e/einv - 1| {max(dev) if dev else float('nan'):.4f}")
        else:
            if not (rt and len(rt) == 8 and max(rt) - min(rt) <= SPREAD857[p] / 2 and abs(re_ - 1) <= 0.05): okS = False
pr(f"H859m (GLOWNA; n = 2 ON |mean/eccentric - 1| <= 0,02, p-Ps i o-Ps): {'TAK' if okM else 'NIE'}")
pr(f"H859e (n = 2 ON |e/einv - 1| <= 0,05 przy einv >= 0,05): {'TAK' if okE else 'NIE'}")
pr(f"H859s (n = 1 ON: rozrzut <= polowa 857 i |mean/eccentric - 1| <= 0,05): {'TAK' if okS else 'NIE'}")

# Audit 860: cheaper variant -- eccentricity-aware envelope shape (same jumps).  Usage: LC_ALL=C python3 ana860.py -> out/ana860.txt.  Written before any result.
import os, re
SP = '/tmp/claude-1000/-home-teddy-Projekty-CREM-QI/4fa83d5c-8878-4cfd-898e-1a9b53a01ab2/scratchpad'
R = os.environ.get('R860', SP + '/r860'); pr = lambda *a: print(*a, flush=True)
rd = lambda v, t, x: open(f'{R}/reg/{v}/{t}/{x}', 'rb').read()
same = lambda a, b, t: all(rd(a, t, x) == rd(b, t, x) for x in ('out.txt', 'err.txt', 'rc.txt'))
SETS = ('1', '2', 'L2', '5Ps', '5H', '6')
for a, b in (('old', 'new'), ('new', 'new_on')):
    pr(f'  {a} vs {b}: ' + ', '.join(f"{t} {'SAME' if same(a, b, t) else 'DIFF'} (rc {rd(a, t, 'rc.txt').decode().strip()}/{rd(b, t, 'rc.txt').decode().strip()},"
                                    f" wall {rd(a, t, 'wall.txt').decode().strip()}/{rd(b, t, 'wall.txt').decode().strip()} s)" for t in SETS))
pr(f"H860k (new = old, 6 zestawow): {'TAK' if all(same('old', 'new', t) for t in SETS) else 'NIE'}")
for v in ('new', 'new_on'):
    s = open(f'{R}/validation_{v}.txt').read().strip().split('\n'); pr(f'  walidacja {v}: ' + ' | '.join(s[-4:]))
pr('H860v: ocena recznie z wierszy powyzej (new)')
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
# H860p / H860w
okP = True
for t in ('1', '2', 'L2'):
    a, b = info('new', t), info('new_on', t)
    eq = a['life'] == b['life'] and a['et'] == b['et'] and a['lab'] == b['lab']
    if t == 'L2' and not eq: okP = False
    pr(f"  zestaw {t}: OFF usable {a['usable']} rc {a['rc']} | ON usable {b['usable']} rc {b['rc']}; zycie/E[T]/Laboratory ON = OFF: {eq}")
    pr(f"    OFF {a['life']} {a['et']}; ON {b['life']} {b['et']}; mean OFF {a['mean']} ON {b['mean']}; ratio ecc OFF {a['re_']} ON {b['re_']}")
okP = okP and all(same('new', 'new_on', t) for t in ('5Ps', '5H', '6'))
pr(f"H860p (L2 zycie / E[T] / Laboratory ON = OFF; 5Ps, 5H, 6 ON = OFF bajt w bajt): {'TAK' if okP else 'NIE'}")
# H860c
REF857 = {'1': 0.998765, '2': 1.00004}
okC = True
for p in ('1', '2'):
    I = info('new_on', f'C_p{p}_n2'); r = float(I['rc_'][0]) if I['rc_'] else float('nan')
    if not abs(r - REF857[p]) <= 1e-3: okC = False
    pr(f"  C p{p} n2 ON: mean {I['mean']} mean/closed {r} (857 {REF857[p]}), rc {I['rc']}")
pr(f"H860c (start kolowy ON = 857 do 1e-3): {'TAK' if okC else 'NIE'}")
# H860m, H860e, H860s
REF = {('1', '2'): 0.91545, ('2', '2'): 0.915939, ('1', '1'): 0.831682, ('2', '1'): 0.890971}
SPREAD857 = {'1': 0.336, '2': 0.403}
okM, okE, okS, okT = True, True, True, True
WALL857 = {('1', '1'): 251, ('1', '2'): 326, ('2', '1'): 256, ('2', '2'): 477}
M860 = {}
import re as _re
t859 = open('../section-859/out/ana859.txt').read() if os.path.exists('../section-859/out/ana859.txt') else ''
for p, lab in (('1', 'p-Ps (PHEN=1)'), ('2', 'o-Ps (PHEN=2)')):
    for n in ('2', '1'):
        I = info('new_on', f'D_p{p}_n{n}')
        rt = [float(d['tlab']) / float(d['eref']) for d in I['tr'] if d.get('decision') == 'used']
        re_ = float(I['re_'][0]) if I['re_'] else float('nan')
        walls = [float(d['wall']) for d in I['tr'] if 'wall' in d]; M860[(p, n)] = re_
        if not (float(I['wall']) <= 1.2 * WALL857[(p, n)] and I['rc'] == '0' and I['usable'] and I['usable'][0] == '16'): okT = False
        pr(f"  D {lab} n = {n} ON: mean {I['mean']}, mean/eccentric {re_} (OFF 857 {REF[(p, n)]}), usable {I['usable']}, rc {I['rc']}, wall {I['wall']} s, life {I['life']}")
        if rt: pr(f"    per-seed tlab/eref min {min(rt):.4f} max {max(rt):.4f} rozrzut {max(rt) - min(rt):.4f} (N {len(rt)}); wall przejscia {min(walls):.0f}-{max(walls):.0f} s; outcome {sorted(set(d['outcome'] for d in I['tr']))}")
        if n == '2':
            if not abs(re_ - 1) <= 0.02: okM = False
            dev, nck = [], []
            for s_, L in I['orb'].items():
                nck.append(len(L)); dev += [abs(float(d['e1']) / float(d['einv']) - 1) for d in L if float(d['einv']) >= 0.05]
            if not (dev and max(dev) <= 0.02): okE = False
            pr(f"    ORBIT858: ziaren {len(I['orb'])}, punktow {min(nck) if nck else 0}-{max(nck) if nck else 0}, max |e/einv - 1| {max(dev) if dev else float('nan'):.4f}")
        else:
            if not abs(re_ - 1) <= 0.05: okS = False
pr(f"H860m (GLOWNA; n = 2 ON |mean/eccentric - 1| <= 0,02, p-Ps i o-Ps): {'TAK' if okM else 'NIE'}")
pr(f"H860e (n = 2 ON |e/einv - 1| <= 0,02 przy einv >= 0,05): {'TAK' if okE else 'NIE'}")
pr(f"H860s (n = 1 ON: |mean/eccentric - 1| <= 0,05): {'TAK' if okS else 'NIE'}")
pr(f"H860t (zegar D ON <= 1,2 x 857, kod 0, usable 16 -- takze n = 1 przy budzecie 90 s): {'TAK' if okT else 'NIE'}")
okX, rep = True, []
for p, lab in (('1', 'p-Ps (PHEN=1)'), ('2', 'o-Ps (PHEN=2)')):
    m = _re.search(r'D ' + _re.escape(lab) + r' n = 2 ON: .*?mean/eccentric ([\d.e+-]+)', t859)
    v859 = float(m.group(1)) if m else float('nan'); rep.append(f'{lab}: 860 {M860[(p, "2")]} / 859 {v859}')
    if not abs(M860[(p, '2')] / v859 - 1) <= 0.01: okX = False
pr(f"H860x (n = 2: |860 / 859 - 1| <= 0,01): {'TAK' if okX else 'NIE'}  ({'; '.join(rep)})")

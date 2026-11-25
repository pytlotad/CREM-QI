# Audit 855: why the collapse time is missing at n = 2.  Usage: LC_ALL=C python3 ana855.py -> out/ana855.txt.  Written before any result.
import os, re
SP = '/tmp/claude-1000/-home-teddy-Projekty-CREM-QI/4fa83d5c-8878-4cfd-898e-1a9b53a01ab2/scratchpad'
REG = os.environ.get('R855REG', SP + '/r855/reg'); pr = lambda *a: print(*a, flush=True)
rd = lambda t, x: open(f'{REG}/{t}/{x}', 'rb').read()
def lines(t):
    out = []
    for ln in rd(t, 'err.txt').decode(errors='replace').split('\n'):
        if ln.startswith('TRANSIT855'):
            w = ln.split(); d = {w[i]: w[i + 1] for i in range(1, len(w) - 1, 2)}; out.append(d)
    return out
k = all(rd('K_L2', x) == rd('OFF_L2', x) for x in ('out.txt', 'err.txt', 'rc.txt'))
o2 = rd('ON_p2_n2', 'out.txt') == rd('OFF_L2', 'out.txt')
cnt = {t: len(lines(t)) for t in ('ON_p1_n1', 'ON_p1_n2', 'ON_p2_n1', 'ON_p2_n2')}
off = len(lines('OFF_L2')) + len(lines('K_L2'))
pr(f"H855k (K = OFF bajt w bajt z kodem wyjscia; ON o-Ps n = 2 out = OFF; 8 wierszy w kazdym ON, 0 w OFF): "
   f"{'TAK' if k and o2 and all(v == 8 for v in cnt.values()) and off == 0 else 'NIE'}  (K {k}; out {o2}; {cnt}; OFF {off})")
for t in ('K_L2', 'OFF_L2', 'ON_p1_n1', 'ON_p1_n2', 'ON_p2_n1', 'ON_p2_n2'): pr(f'  rc {t} = {rd(t, "rc.txt").decode().strip()}')
OUT = {'0': 'ReachedCutoff', '1': 'ObservationLimit', '2': 'NumericalFailure'}
STOP = {'0': 'None', '1': 'Compton', '2': 'Retardation', '3': 'Floor', '4': 'ChannelClosed'}
res = {}
for p, lab in (('1', 'p-Ps (PHEN=1)'), ('2', 'o-Ps (PHEN=2)')):
    for n in ('1', '2'):
        L = lines(f'ON_p{p}_n{n}'); res[(p, n)] = L
        pr(f'  {lab} n = {n}:')
        for d in L:
            tl, rf = float(d['tlab']), float(d['ref_start'])
            pr(f"    seed {d['seed']:>20} a0ratio {float(d['a0ratio']):.6f} {OUT.get(d['outcome'], d['outcome']):>16} stop {STOP.get(d['stop'], d['stop']):>13}"
               f" tlab {tl:.4e} s ref_start {rf:.4e} s tlab/ref {tl / rf if rf else float('nan'):.4f} aterm {float(d['aterm']):.4e} wall {float(d['wall']):8.1f} s  {d['decision']}")
        m = re.search(r'over (\d+) usable trajectories \((\d+) observed, (\d+) right-censored\)', rd(f'ON_p{p}_n{n}', 'out.txt').decode())
        pr(f'    out.txt: usable {m.group(1)}, observed {m.group(2)}, censored {m.group(3)}' if m else '    out.txt: brak wiersza usable')
ok = lambda p, n, f: sum(1 for d in res[(p, n)] if f(d))
l = all(ok(p, '1', lambda d: d['decision'] == 'used' and abs(float(d['a0ratio']) - 1) <= 1e-6) == 8 for p in '12')
pr(f"H855l (n = 1: used 8/8, a0ratio = 1 do 1e-6): {'TAK' if l else 'NIE'}")
a = all(ok(p, '2', lambda d: d['decision'] == 'rejected_start' and 3.9 <= float(d['a0ratio']) <= 4.1) == 8 for p in '12')
pr(f"H855a (GLOWNA; n = 2: rejected_start 8/8, a0ratio w [3,9; 4,1]): {'TAK' if a else 'NIE'}  (p-Ps {ok('1', '2', lambda d: d['decision'] == 'rejected_start')}, o-Ps {ok('2', '2', lambda d: d['decision'] == 'rejected_start')})")
o = all(ok(p, '2', lambda d: d['outcome'] == '0') >= 7 for p in '12')
pr(f"H855o (n = 2: ReachedCutoff >= 7/8): {'TAK' if o else 'NIE'}  (p-Ps {ok('1', '2', lambda d: d['outcome'] == '0')}, o-Ps {ok('2', '2', lambda d: d['outcome'] == '0')})")
for p in '12':
    for n in '12':
        L = res[(p, n)]; w = sum(float(d['wall']) for d in L)
        pr(f"  RAPORT phen {p} n {n}: suma zegara sciany biegu ciaglego {w:.0f} s; srednia tlab/ref {sum(float(d['tlab']) / float(d['ref_start']) for d in L) / max(1, len(L)):.4f}")

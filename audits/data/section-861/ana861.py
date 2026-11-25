# Audit 861: eccentric jump shape (860) made default.  Usage: LC_ALL=C python3 ana861.py -> out/ana861.txt.  Written before any result.
import os, re
SP = '/tmp/claude-1000/-home-teddy-Projekty-CREM-QI/4fa83d5c-8878-4cfd-898e-1a9b53a01ab2/scratchpad'
R = os.environ.get('R861', SP + '/r861'); pr = lambda *a: print(*a, flush=True)
rd = lambda base, v, t, x: open(f'{base}/reg/{v}/{t}/{x}', 'rb').read()
NEW, REF = R, R + '/ref860'
SETS = ('1', '2', 'L2', '5Ps', '5H', '6')
def same(v, rv, t, files=('out.txt', 'err.txt', 'rc.txt')):
    return all(rd(NEW, v, t, x) == rd(REF, rv, t, x) for x in files)
for v, rv, rule in (('new_legacy', 'new', 'H861l'), ('new', 'new_on', 'H861n')):
    res = {t: same(v, rv, t) for t in SETS}
    pr(f'  {v} vs 860 {rv}: ' + ', '.join(f"{t} {'SAME' if res[t] else 'DIFF'} (rc {rd(NEW, v, t, 'rc.txt').decode().strip()}/{rd(REF, rv, t, 'rc.txt').decode().strip()},"
                                          f" wall {rd(NEW, v, t, 'wall.txt').decode().strip()}/{rd(REF, rv, t, 'wall.txt').decode().strip()} s)" for t in SETS))
    if rule == 'H861l': okL = all(res.values())
    else: okN = all(res.values())
pr(f"H861l (new_legacy = 860 new bajt w bajt, 6 zestawow): {'TAK' if okL else 'NIE'}")
# D p-Ps n = 2: out and rc byte for byte; err with the wall field removed
strip = lambda b: re.sub(rb' wall [^ \n]+', b' wall X', b)
dOut = rd(NEW, 'new', 'D_p1_n2', 'out.txt') == rd(REF, 'new_on', 'D_p1_n2', 'out.txt') and rd(NEW, 'new', 'D_p1_n2', 'rc.txt') == rd(REF, 'new_on', 'D_p1_n2', 'rc.txt')
eA = sorted(strip(rd(NEW, 'new', 'D_p1_n2', 'err.txt')).split(b'\n')); eB = sorted(strip(rd(REF, 'new_on', 'D_p1_n2', 'err.txt')).split(b'\n'))
dErr = eA == eB
pr(f'  D p-Ps n = 2: out + rc = 860 {dOut}; err (bez pola wall, wiersze posortowane) = 860 {dErr}')
pr(f"H861n (GLOWNA; new = 860 new_on bajt w bajt w 6 zestawach + D p-Ps n = 2): {'TAK' if okN and dOut and dErr else 'NIE'}")
for v in ('new', 'new_legacy'):
    s = open(f'{R}/validation_{v}.txt').read().strip().split('\n'); pr(f'  walidacja {v}: ' + ' | '.join(s[-4:]))
pr('H861v: ocena recznie z wierszy powyzej (new)')
for t, lab in (('1', 'p-Ps (PHEN=1) n = 1'), ('2', 'o-Ps (PHEN=2) n = 1'), ('L2', 'o-Ps (PHEN=2) n = 2'), ('D_p1_n2', 'p-Ps (PHEN=1) n = 2')):
    s = rd(NEW, 'new', t, 'out.txt').decode()
    g = lambda pat: (re.search(pat, s).group(1) if re.search(pat, s) else '?')
    pr(f"  RAPORT {lab}: czas zapadania {g(r'mean of completed runs ([^\n]+)')}; / wzor eliptyczny {g(r'mean / eccentric form\s+([\d.e+-]+)')}; "
       f"{g(r'(lifetime at the final state [^;]+)')}")

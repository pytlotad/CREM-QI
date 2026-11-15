# Audit 796: two-body target of the periapsis start (canonical L = nominal, E_T = E_S - nu mu c^2 alpha^4/(8 n^4)); CREM_ONEBODY_START=1 restores 795.
# usage: python3 ana796.py   (regression / validation copied from the scratchpad into out/; 795 G unpacked into $R795; reg files in $R796REG)
import sys, os, math, re
import numpy as np
D = os.path.dirname(os.path.abspath(__file__)); O = os.path.join(D, 'out', 'runs')
sys.path.insert(0, os.path.join(D, '..', 'section-699')); import ana699
SP = '/tmp/claude-1000/-home-teddy-Projekty-CREM-QI/4fa83d5c-8878-4cfd-898e-1a9b53a01ab2/scratchpad'
R795 = os.environ.get('R795', SP + '/r795x/runs'); REG = os.environ.get('R796REG', SP + '/r796reg') + '/reg'
HEV, C = 4.135667696e-15, 299792458.0; FS = 10969.045e6
pr = lambda *a: print(*a, flush=True)
units = lambda f: f'{f/1e6:.2f} MHz = {f*HEV*1e6:.4f} ueV = {C/f*100:.4f} cm'
TORI = [(3, 1), (3, 1.5), (3, 2), (3, 2.5), (3, 3), (2, 1), (2, 1.5), (2, 2)]
lab = lambda L: '%g' % L
def omega(path):
    m = ana699.parse(path)['inc']['mechB0']
    return np.linalg.norm(m['s21'] - m['s20'])/(m['dt']*np.linalg.norm(np.cross(m['L0'], m['s20'])))/(2*math.pi)
# --- H796o, H796v
reg = [l for l in open(os.path.join(D, 'out', 'regression796.txt')).read().split('\n') if l]
pr('\n'.join(reg)); okO = all(l.endswith('SAME') for l in reg if l.startswith('old vs new_off'))
pr(f"H796o (old = new_off bajt w bajt, 6 zestawow): {'TAK' if okO else 'NIE'}")
V = {}
for v in ('old', 'new'):
    t = open(os.path.join(D, 'out', f'validation_{v}.txt')).read()
    m = re.search(r'Enforced checks:\s*(\d+)/(\d+)', t); e = re.search(r'EXIT (\d+)', t)
    V[v] = (int(m.group(1)), int(m.group(2)), int(e.group(1))) if m and e else None; pr(f'  walidacja {v}: {V[v]}')
okV = None not in V.values() and V['new'][2] == 0 and V['new'][0] == V['new'][1] and V['new'][:2] == V['old'][:2]
pr(f"H796v (walidacja new EXIT 0, wszystkie kontrole, liczba jak old): {'TAK' if okV else 'NIE'}")
# --- t796
rows = []
for ln in open(os.path.join(D, 'out', 't796.txt')):
    parts = [p.split() for p in ln.split('|')]; h = parts[0]
    r = dict(pair=h[0], dip=int(h[2]), n=float(h[4]), L=float(h[6]))
    for p in parts[1:]:
        md = p[0].rstrip(':'); r[md] = {p[i]: float(p[i + 1]) for i in range(1, len(p) - 1, 2)}
    rows.append(r)
okK = all(abs(r['nowy']['Lkan/L-1']) <= 1e-9 for r in rows)
okE = all(abs(r['nowy']['Eorb/E_T-1']) <= 1e-9 for r in rows if r['L'] != r['n'])
okC = all(abs(r['nowy']['Eorb/E_T-1']) <= 2e-9 for r in rows if r['L'] == r['n'])
okS = all(abs(r['795']['Lkan/L-1']/r['795']['pred'] - 1) <= 0.01 for r in rows)
pr(f'=== t796: {len(rows)} linii ===')
for r in rows:
    if r['dip'] != 0: continue
    a, b, k = r['nowy'], r['795'], r['Kepler']
    pr(f"  {r['pair'][:6]} n {r['n']:g} L {r['L']:3g}: nowy E/E_T-1 {a['Eorb/E_T-1']:+.2e} Lkan-1 {a['Lkan/L-1']:+.2e} Lmech-1 {a['Lmech/L-1']:+.2e} | 795 E/E_T-1 {b['Eorb/E_T-1']:+.2e}"
       f" Lkan-1 {b['Lkan/L-1']:+.3e} (pred {b['pred']:.3e}) | Kepler E/E_T-1 {k['Eorb/E_T-1']:+.2e}")
pr(f"  maks.: |E/E_T-1| mimosrodowe {max(abs(r['nowy']['Eorb/E_T-1']) for r in rows if r['L'] != r['n']):.1e}, kolowe {max(abs(r['nowy']['Eorb/E_T-1']) for r in rows if r['L'] == r['n']):.1e};"
   f" |Lkan/L-1| {max(abs(r['nowy']['Lkan/L-1']) for r in rows):.1e}; 795 Lkan/pred-1 maks. {max(abs(r['795']['Lkan/L-1']/r['795']['pred'] - 1) for r in rows):.1e}")
pr(f"H796k (|L_kan/L_nom - 1| <= 1e-9, 48 linii): {'TAK' if okK else 'NIE'}")
pr(f"H796e (mimosrodowe |E_orb/E_T - 1| <= 1e-9): {'TAK' if okE else 'NIE'}")
pr(f"H796c (kolowe |E_orb/E_T - 1| <= 2e-9): {'TAK' if okC else 'NIE'}")
pr(f"H796s (795: L_kan/L - 1 = nu alfa^2 a/r_p do 1 %): {'TAK' if okS else 'NIE'}")
# --- G
okG = True; mx = 0.0; W = {}
pr('=== G (ziarno 5, F = 1, start dwucialowy) wobec 795 G ===')
for t in ('proton', 'antiproton'):
    W[t] = {k: omega(f'{O}/G/{t}_n{k[0]}L{lab(k[1])}_5.err') for k in TORI}
    q = {k: W[t][k]/omega(f'{R795}/G/{t}_n{k[0]}L{lab(k[1])}_5.err') - 1 for k in TORI}
    mx = max(mx, max(abs(x) for x in q.values())); okG &= all(abs(x) <= 1e-5 for x in q.values())
    pr(f'  {t}: ' + ', '.join(f'n{k[0]}L{lab(k[1])} {q[k]:+.1e}' for k in TORI))
    for n, l, st in ((2, 1, '2p'), (3, 1, '3p'), (3, 2, '3d')):
        od = l*W[t][(n, l)] - (l + 1)*W[t][(n, l + 1)]; qq = od/W[t][(n, l + 0.5)]/((l + 0.5)**2/(l*(l + 1)))
        pr(f'    {st}: {units(od)} | {qq:.5f} x regula' + (f' | {od/FS:.5f} pomiaru 10 969,045(15) MHz (Hagley i Pipkin 1994)' if st == '2p' else ''))
pr(f"H796g (Omega G / 795 G - 1 w 1e-5, maks. {mx:.1e}): {'TAK' if okG else 'NIE'}")
# --- H796z and report
def grab(v, s):
    t = open(f'{REG}/{v}/{s}/out.txt').read()
    m = re.search(r'mean of completed runs\s+([\d.e+-]+) \+/- ([\d.e+-]+) (\w+)', t); a = re.search(r'lifetime at the final state ([\d.]+ \w+)', t)
    e = re.search(r'E\[T\] from the preparation \(cascade \+ annihilation\) ([\d.]+ \w+)', t)
    return (float(m.group(1)), float(m.group(2)), m.group(3)) if m else None, a.group(1) if a else None, e.group(1) if e else None
okZ = True; cls = []
for s, ch in (('1', 'p-Ps'), ('2', 'o-Ps')):
    o, oa, oe = grab('old', s); n, na, ne = grab('new', s); d = n[0]/o[0] - 1; sg = math.hypot(o[1], n[1])/o[0]
    c = 'gladka' if abs(d) <= 0.005 else ('chaotyczna' if abs(d) <= 0.10 else 'poza obu')
    cls.append(c); pr(f'  {ch} (zestaw {s}): zapadanie 795 {o[0]:.6g} {o[2]} -> 796 {n[0]:.6g} {n[2]} ({d*100:+.2f} %, sigma roznicy {sg*100:.1f} %) -> {c};'
                      f' anihilacja {oa} -> {na}; E[T] {oe} -> {ne}')
pr(f"H796z (odczyt zagadki zapadania): p-Ps {cls[0]}, o-Ps {cls[1]}")
for s, ch in (('L2', 'o-Ps n = 2'), ('5H', 'wodor 5H'), ('5Ps', 'Ps 5'), ('6', 'Ps 6')):
    o = open(f'{REG}/old/{s}/out.txt').read().split('\n'); n = open(f'{REG}/new/{s}/out.txt').read().split('\n')
    d = [(a, b) for a, b in zip(o, n) if a != b]; pr(f'  {ch}: {len(d)} linii rozne' + ''.join(f'\n    old: {a.strip()[:140]}\n    new: {b.strip()[:140]}' for a, b in d[:6]))
    _, a, e = grab('new', s) if s in ('L2',) else (None, None, None)
    if a: pr(f'    anihilacja {a}, E[T] {e}')
pr('  o-Ps: czynnik z wodoru (Hagena i in., PRL 71, 2887 (1993)); p-Ps bez czesci spinowej; F = 0 tozsamy w tym ustawieniu (795), nie biegany.')

# Audit 826: photon target on the whole rung.  Usage: python3 ana826.py -> out/ana826.txt.  Written during the runs, before reading them.
# Needs ../section-825/runs/runs825.tar.xz (byte check of new_off, H826o).
import os, re, glob, tarfile, io
from collections import defaultdict
D = os.path.dirname(os.path.abspath(__file__)); O = os.path.join(D, 'out', 'runs'); pr = lambda *a: print(*a, flush=True)
kv = lambda s: dict(re.findall(r'([A-Za-z_][A-Za-z_0-9]*)=(\S+)', s))
SYS = {'H0': 'wodor F = 0 (PHEN=1)', 'H1': 'wodor F = 1 (PHEN=2)', 'pPs': 'p-Ps (PHEN=1)', 'oPs': 'o-Ps (PHEN=2)'}
TOL = {'H0': 5e-4, 'H1': 5e-4, 'pPs': 1e-5, 'oPs': 1e-5}
def events(path):
    ev = []; cur = None
    for ln in open(path):
        if ln.startswith('LADDER E_rel'):
            d = kv(ln); cur = dict(n=float(d['n']), photon=float(d['photon'])); ev.append(cur)
        elif ln.startswith('RULES') and cur is not None:
            d = kv(ln); cur['nr'] = float(d['n']); cur['nf'] = int(d['n_final'])
    return ev
okR = okK = okL = True; rows = defaultdict(list); refl = defaultdict(list); reach = defaultdict(list)
runs = sorted(glob.glob(f'{O}/new/*.err')); pr(f'  biegi new: {len(runs)}')
for p in runs:
    name = os.path.basename(p)[:-4]; sysk, st, seed = name.split('_'); ev = events(p)
    nref = sum(1 for e in ev if e['photon'] == 0.0 and e.get('nr', e['n']) >= 1.5); okR &= nref == 0
    got = any(e['photon'] > 0 and e.get('nf') == 1 for e in ev); okK &= got
    refl[(sysk, st)].append(nref); reach[(sysk, st)].append(got)
    for i, e in enumerate(ev):
        if e['photon'] > 0 and i + 1 < len(ev) and 'nf' in e:
            nx = ev[i + 1]; dev = nx.get('nr', nx['n']) - e['nf']
            okL &= abs(dev) <= TOL[sysk]; rows[(sysk, f"{e['nf'] + round(e.get('nr', e['n']) - e['nf'])}->{e['nf']}")].append(dev)
for sysk, lab in SYS.items():
    pr(f'  {lab} (prog |dn| {TOL[sysk]:.0e}):')
    for key in sorted(k for k in rows if k[0] == sysk):
        v = rows[key]; pr(f'    {key[1]}: N {len(v)}; n(nastepny) - n_final {min(v):+.2e} .. {max(v):+.2e}; max |dn| {max(abs(x) for x in v):.2e}')
    for st in ('4f', '4d', '4p', '3d'):
        if (sysk, st) in refl: pr(f'    start {st}: odrzucone {refl[(sysk, st)]}; do 1s {reach[(sysk, st)]}')
# H826o: new_off vs 825, byte for byte
tf = tarfile.open(os.path.join(D, '..', 'section-825', 'runs', 'runs825.tar.xz')); okO = True; noff = 0
for p in sorted(glob.glob(f'{O}/off/*')):
    b = os.path.basename(p); m = tf.extractfile(f'runs/{b}'); same = m is not None and m.read() == open(p, 'rb').read()
    okO &= same; noff += 1
    if not same: pr(f'    H826o ROZNICA: {b}')
okO &= noff == 24
pr(f'  new_off: {noff} plikow (12 biegow x out, err) porownanych z 825')
reg = open(os.path.join(D, 'out', 'regression826.txt')).read().splitlines() if os.path.exists(os.path.join(D, 'out', 'regression826.txt')) else []
okB = len(reg) == 12 and all(l.endswith('SAME') for l in reg)
for v in ('old', 'new'):
    f = os.path.join(D, 'out', f'validation_{v}.txt'); t = open(f).read() if os.path.exists(f) else ''
    ok = 'Enforced checks:     72/72' in t and 'validation verdict: PASS' in t; okB &= ok; pr(f'  walidacja {v}: {"72/72 PASS" if ok else "BRAK / NIE"}')
pr('  regresja: ' + '; '.join(reg))
pf = os.path.join(D, 'out', 'probes826.txt'); pr('  sondy: ' + ('; '.join(open(pf).read().splitlines()) if os.path.exists(pf) else 'brak'))
pr(f"H826b (produkcja bajt w bajt, walidacja 72/72 old i new): {'TAK' if okB else 'NIE'}")
pr(f"H826o (new_off = 825 bajt w bajt, 12 biegow 4f): {'TAK' if okO else 'NIE'}")
pr(f"H826r (GLOWNA: 0 odrzucen i wszystkie do 1s): {'TAK' if okR and okK else 'NIE'} (odrzucenia {'0' if okR else '>0'}, do 1s {'wszystkie' if okK else 'nie wszystkie'})")
pr(f"H826l (|n_nastepny - n_final| <= 5e-4 wodor, 1e-5 Ps): {'TAK' if okL else 'NIE'}")

# Audit 825: landing on the torus after each photon on the mechanical path.  Usage: python3 ana825.py -> out/ana825.txt.
# Written during the runs, before reading them; POPRAWKA 2 (16:31): key regex with digits, open/closed window split.
import os, re, glob, math
from collections import defaultdict
D = os.path.dirname(os.path.abspath(__file__)); O = os.path.join(D, 'out', 'runs'); pr = lambda *a: print(*a, flush=True)
kv = lambda s: dict(re.findall(r'([A-Za-z_][A-Za-z_0-9]*)=(\S+)', s))
HEV, C, SCALE = 4.135667696e-15, 299792458.0, 8e5
def parse(path):
    ev = []; cur = None
    for ln in open(path):
        if ln.startswith('LADDER E_rel'):
            d = kv(ln); cur = dict(E=float(d['E_rel']), n=float(d['n']), photon=float(d['photon'])); ev.append(cur)
        elif ln.startswith('LADDER   e=') and cur is not None:
            cur['k'] = float(kv(ln)['harmonic'])
        elif ln.startswith('RULES') and cur is not None:
            d = kv(ln); cur['nr'] = float(d['n']); cur['L'] = float(d['L']); cur['nf'] = int(d['n_final']); cur['down'] = float(d['downE1'])
        elif ln.startswith('SPLIT t=') and cur is not None:
            cur['jrem'] = float(kv(ln)['Jrem'])
        elif ln.startswith('WINDOW closed'):
            d = kv(ln.replace('unpaid J', 'unpaidJ').replace(' L=', ' Lun='))
            w = [e for e in ev if e['photon'] > 0 and 'win' not in e]
            if cur is not None: cur['closed_after'] = True
            if w: w[-1]['win'] = dict(n=float(d['n']), E_open=float(d['E_open']), E_close=float(d['E_close']), paid=float(d['paid']), queued=float(d['queued']), unpaidJ=float(d['unpaidJ']), t=float(d['t']))
        elif ln.startswith('SPLIT closed'):
            d = kv(ln.replace('[', ' ').replace(']', ' ').replace('kin ', 'kin=').replace('coul ', 'coul=').replace('dd ', 'dd=').replace('cd ', 'cd=').replace('darwin ', 'darwin=').replace('constr ', 'constr='))
            w = [e for e in ev if e['photon'] > 0 and 'split' not in e]
            if w: w[-1]['split'] = {k: float(d[k]) for k in ('photon', 'rotation', 'dynamics', 'kin', 'coul', 'dd', 'cd', 'darwin', 'constr')}
    return ev
runs = sorted(glob.glob(f'{O}/*.err')); SYS = {'H0': 'wodor F = 0 (PHEN=1)', 'H1': 'wodor F = 1 (PHEN=2)', 'pPs': 'p-Ps (PHEN=1)', 'oPs': 'o-Ps (PHEN=2)'}
okC = okW = okL = okR = True; agg = defaultdict(list); refuse = defaultdict(list); reach = defaultdict(list); worstC = 0.0
for p in runs:
    name = os.path.basename(p)[:-4]; sysk, st, seed = name.split('_'); ev = parse(p)
    if not ev: pr(f'  {name}: brak fotonow'); continue
    B = -ev[0]['E']*ev[0]['n']**2   # from 6-digit n: refined below with RULES n
    if 'nr' in ev[0]: B = -ev[0]['E']*ev[0]['nr']**2
    nrefused = 0; t_first_ref = None; last_n = None
    for i, e in enumerate(ev):
        n = e.get('nr', e['n'])
        refused = e['photon'] == 0.0 and n >= 1.5
        predicted = n >= 1.5 and n < 2 - 1e-4 and e.get('k', 1) >= 1
        okR &= (refused == predicted) if n >= 1.5 else True
        if refused: nrefused += 1
        if e['photon'] > 0:
            nt = n - min(e.get('k', 1), math.floor(n - 1 + 1e-4))
            nxt = next((x for x in ev[i + 1:]), None)
            dn_next = (nxt.get('nr', nxt['n']) - nt) if nxt else float('nan')
            dn_next_int = (nxt.get('nr', nxt['n']) - e.get('nf', round(nt))) if nxt else float('nan')
            w = e.get('win'); s = e.get('split')
            if w:
                nclose = math.sqrt(B/(-w['E_close'])); dn_close = nclose - nt
                okW &= abs(dn_close) <= 1e-6
                if s:
                    lhs = w['E_close'] - w['E_open']; rhs = -w['paid'] + s['photon'] + s['rotation'] + s['dynamics']
                    rel = abs(lhs - rhs)/max(abs(lhs), 1e-300); worstC = max(worstC, rel); okC &= rel <= 1e-3
                    conv = lambda dE: dE*nt**3/(2*B)
                    comp = {k: conv(s[k]) for k in ('photon', 'rotation', 'dynamics', 'kin', 'coul', 'darwin', 'constr')}
                    comp['unpaid'] = conv(w['paid'] - w['queued'])
                else: comp = {}
            else: dn_close = float('nan'); comp = {}
            if nxt: okL &= dn_next >= -1e-4
            agg[(sysk, f"{round(n)}->{round(nt)}")].append((dn_close, dn_next, dn_next_int, comp, name, e.get('closed_after', False), e.get('jrem', float('nan'))))
    reach[(sysk, st)].append(min(x.get('nr', x['n']) for x in ev) < 1.5 or any(x['photon'] > 0 and round(x.get('nr', x['n']) - x.get('k', 1)) == 1 for x in ev))
    refuse[(sysk, st)].append(nrefused)
pr(f'  biegi: {len(runs)}')
for sysk, lab in SYS.items():
    pr(f'  {lab}:')
    for key in sorted(k for k in agg if k[0] == sysk):
        v = agg[key]; dc = [x[0] for x in v if not math.isnan(x[0])]; dnx = [x[1] for x in v if not math.isnan(x[1])]
        comps = defaultdict(list)
        for x in v:
            for k2, val in x[3].items(): comps[k2].append(val)
        cs = ', '.join(f'{k2} {sum(vv)/len(vv):+.2e}' for k2, vv in comps.items())
        pr(f'    {key[1]}: N {len(v)}; dn_close {min(dc) if dc else float("nan"):+.2e} .. {max(dc) if dc else float("nan"):+.2e}; dn_next {min(dnx) if dnx else float("nan"):+.2e} .. {max(dnx) if dnx else float("nan"):+.2e}; srednie skladniki dn: {cs}')
        for flag, lab2 in ((True, 'okno zamkniete'), (False, 'okno otwarte')):
            dd = [x[1] for x in v if x[5] == flag and not math.isnan(x[1])]; jr = [x[6] for x in v if x[5] == flag and not math.isnan(x[6])]
            if dd: pr(f'      {lab2} przed nastepnym fotonem: N {len(dd)}, dn_next {min(dd):+.2e} .. {max(dd):+.2e} (mediana {sorted(dd)[len(dd)//2]:+.2e}); ostatnie Jrem {min(jr) if jr else float("nan"):.3g} .. {max(jr) if jr else float("nan"):.3g}')
    for st in ('4f', '4d', '4p', '3d'):
        if (sysk, st) in refuse: pr(f'    start {st}: odrzucone fotony {refuse[(sysk, st)]}; doszlo do 1s (foton do n = 1) {reach[(sysk, st)]}')
pr(f"H825c (zamkniecie SPLIT do 1e-3): {'TAK' if okC else 'NIE'} (najgorsze {worstC:.2e})")
pr(f"H825w (okno zamyka sie na celu, |dn_close| <= 1e-6): {'TAK' if okW else 'NIE'}")
pr(f"H825l (GLOWNA: dn_next >= -1e-4 w kazdym fotonie): {'TAK' if okL else 'NIE'}")
pr(f"H825r (odrzucenie <=> n < 2 - 1e-4): {'TAK' if okR else 'NIE'}")

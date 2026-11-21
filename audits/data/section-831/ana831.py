# Audit 831: J4L3 payout weight with a slower window (WS: eps 0.002) vs W0 (826), V3 (827, eps 0.002 uniform), WJL (829, eps 0.01 J4L3).
# Usage: python3 ana831.py -> out/ana831.txt.  Written before the runs; analyse() copied from ana829.py.
import os, re, tarfile, io
from collections import defaultdict
D = os.path.dirname(os.path.abspath(__file__)); O = os.path.join(D, 'out', 'runs'); pr = lambda *a: print(*a, flush=True)
kv = lambda s: dict(re.findall(r'([A-Za-z_][A-Za-z_0-9]*)=(\S+)', s))
SYS = {'H0': 'wodor F = 0 (PHEN=1)', 'H1': 'wodor F = 1 (PHEN=2)', 'pPs': 'p-Ps (PHEN=1)', 'oPs': 'o-Ps (PHEN=2)'}
VAR = {'W0': 'eps 0,01 rownomiernie (826)', 'V3': 'eps 0,002 rownomiernie (827)', 'WJL': 'eps 0,01 J4L3 (829)', 'WS': 'eps 0,002 J4L3 (831)'}
ARC = {'W0': ('section-826/runs/runs826.tar.xz', 'runs/new'), 'V3': ('section-827/runs/runs827.tar.xz', 'runs/V3'), 'WJL': ('section-829/runs/runs829.tar.xz', 'runs/WJL')}
def analyse(lines):
    ev = []; cur = None; splits = []
    for ln in lines:
        if ln.startswith('LADDER E_rel'):
            d = kv(ln); cur = dict(n=float(d['n']), photon=float(d['photon'])); ev.append(cur)
        elif ln.startswith('RULES') and cur is not None:
            d = kv(ln); cur['nr'] = float(d['n']); cur['nf'] = int(d['n_final'])
        elif ln.startswith('SPLIT closed'): splits.append(ln.strip())
    out = defaultdict(list)
    for i, e in enumerate(ev):
        if e['photon'] > 0 and i + 1 < len(ev) and 'nf' in e:
            nx = ev[i + 1]; out[f"{round(e['nr'])}->{e['nf']}"].append(nx.get('nr', nx['n']) - e['nf'])
    nref = sum(1 for e in ev if e['photon'] == 0.0 and e.get('nr', e['n']) >= 1.5)
    got = any(e['photon'] > 0 and e.get('nf') == 1 for e in ev)
    return out, nref, got, splits
data = {}
for v, (a, sub) in ARC.items():
    tf = tarfile.open(os.path.join(D, '..', a))
    for sysk in SYS:
        for i in range(3): data[(v, sysk, i)] = analyse(io.TextIOWrapper(tf.extractfile(f'{sub}/{sysk}_4f_{i}.err')).read().splitlines())
for sysk in SYS:
    for i in range(3):
        p = f'{O}/WS/{sysk}_4f_{i}.err'
        if os.path.exists(p): data[('WS', sysk, i)] = analyse(open(p).read().splitlines())
Dmax = {}; okR = True
for sysk, lab in SYS.items():
    pr(f'  {lab}:')
    for v, vl in VAR.items():
        for tr in ('4->3', '3->2'):
            vals = [x for i in range(3) if (v, sysk, i) in data for x in data[(v, sysk, i)][0].get(tr, [])]
            if not vals: pr(f'    {v} {vl} {tr}: brak'); continue
            Dmax[(v, sysk, tr)] = max(abs(x) for x in vals)
            pr(f'    {v} {vl} {tr}: N {len(vals)}; dn ' + ' '.join(f'{x:+.2e}' for x in vals) + f'; max |dn| {Dmax[(v, sysk, tr)]:.2e}')
        refs = [data[(v, sysk, i)][1] for i in range(3) if (v, sysk, i) in data]; got = [data[(v, sysk, i)][2] for i in range(3) if (v, sysk, i) in data]
        pr(f'    {v}: odrzucenia {refs}; do 1s {got}')
        if v == 'WS': okR &= len(refs) == 3 and all(r == 0 for r in refs) and all(got)
        if sysk in ('pPs', 'oPs') and v in ('W0', 'WJL', 'WS'):
            for i in range(3):
                if (v, sysk, i) in data:
                    for s in data[(v, sysk, i)][3]: pr(f'      SPLIT {v} ziarno {i} {s[:230]}')
for sysk in SYS:
    for tr in ('4->3', '3->2'):
        a = Dmax.get(('W0', sysk, tr))
        if a: pr(f'  {SYS[sysk]} {tr} stosunek Vx/W0: ' + ', '.join(f'{v} {Dmax[(v, sysk, tr)]/a:.3g}' for v in ('V3', 'WJL', 'WS') if (v, sysk, tr) in Dmax))
g = lambda v, k: Dmax.get((v, k, '4->3'))
okA = all(g('WS', k) is not None and g('WS', k) <= g('W0', k)/10 for k in ('H0', 'H1'))
okV = all(g('WS', k) is not None and g('WS', k) <= g('V3', k)/2 for k in ('H0', 'H1'))
okP = all(g('WS', k) is not None and g('WS', k) <= g('W0', k) for k in ('pPs', 'oPs'))
pr(f"H831a (GLOWNA: WS <= W0/10, wodor F = 0 i F = 1): {'TAK' if okA else 'NIE'}")
pr(f"H831v (WS <= V3/2, wodor F = 0 i F = 1): {'TAK' if okV else 'NIE'}")
pr(f"H831p (WS <= W0, p-Ps i o-Ps): {'TAK' if okP else 'NIE'}")
pr(f"H831r (0 odrzucen, 12 kaskad do 1s): {'TAK' if okR else 'NIE'}")

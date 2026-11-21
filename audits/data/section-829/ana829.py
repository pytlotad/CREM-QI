# Audit 829: window payout weight WJL (J r^-4, L r^-3) and WF (phi) vs uniform W0.  Usage: python3 ana829.py -> out/ana829.txt.
# Written before the runs (from ana827.py).  W0 = 4f runs of 826 (../section-826/runs/runs826.tar.xz); file check W0 here = 826 byte for byte.
import os, re, glob, tarfile, io, statistics
from collections import defaultdict
D = os.path.dirname(os.path.abspath(__file__)); O = os.path.join(D, 'out', 'runs'); pr = lambda *a: print(*a, flush=True)
kv = lambda s: dict(re.findall(r'([A-Za-z_][A-Za-z_0-9]*)=(\S+)', s))
SYS = {'H0': 'wodor F = 0 (PHEN=1)', 'H1': 'wodor F = 1 (PHEN=2)', 'pPs': 'p-Ps (PHEN=1)', 'oPs': 'o-Ps (PHEN=2)'}
VAR = {'V0': 'W0 rownomiernie (826)', 'WJL': 'J r^-4, L r^-3', 'WF': 'rownomiernie w kacie'}
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
tf = tarfile.open(os.path.join(D, '..', 'section-826', 'runs', 'runs826.tar.xz'))
for sysk in SYS:
    for i in range(3):
        data[('V0', sysk, i)] = analyse(io.TextIOWrapper(tf.extractfile(f'runs/new/{sysk}_4f_{i}.err')).read().splitlines())
        for v in ('WJL', 'WF'):
            p = f'{O}/{v}/{sysk}_4f_{i}.err'
            if os.path.exists(p): data[(v, sysk, i)] = analyse(open(p).read().splitlines())
Dmax = {}
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
        for i in range(3):
            if (v, sysk, i) in data:
                for s in data[(v, sysk, i)][3]: pr(f'      ziarno {i} {s[:200]}')
def hyp(v, cond):
    ok = True
    for sysk in ('H0', 'H1'):
        a = Dmax.get(('V0', sysk, '4->3')); b = Dmax.get((v, sysk, '4->3'))
        ok &= a is not None and b is not None and cond(a, b)
    return ok
for sysk in ('H0', 'H1', 'pPs', 'oPs'):
    a = Dmax.get(('V0', sysk, '4->3'))
    pr(f'  {SYS[sysk]} 4->3 stosunek Vx/V0: ' + ', '.join(f'{v} {Dmax[(v, sysk, "4->3")]/a:.3g}' for v in ('WJL', 'WF') if a and (v, sysk, '4->3') in Dmax))
okB = all(open(f'{O}/W0/{b}', 'rb').read() == tf.extractfile(f'runs/new/{b}').read() for b in ('H1_4f_0.err', 'H1_4f_0.out', 'oPs_4f_0.err', 'oPs_4f_0.out'))
pr(f"H829b (bez przelacznika = 826 bajt w bajt, 2 biegi): {'TAK' if okB else 'NIE'}")
pr(f"H829j (GLOWNA: WJL <= W0/4, wodor F = 0 i F = 1): {'TAK' if hyp('WJL', lambda a, b: b <= a/4) else 'NIE'}")
pr(f"H829f (WF <= W0/4): {'TAK' if hyp('WF', lambda a, b: b <= a/4) else 'NIE'}")
okP = all(Dmax.get(('WJL', k, '4->3'), 9) <= 1.5*Dmax[('V0', k, '4->3')] for k in ('pPs', 'oPs'))
pr(f"H829p (p-Ps i o-Ps WJL <= 1,5 W0): {'TAK' if okP else 'NIE'}")

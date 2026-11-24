# Audit 849: payout-window drift vs window fraction eps.  Usage: LC_ALL=C python3 ana849.py -> out/ana849.txt.
# Written during the runs, before reading them.  analyse() from ana827.py.  H / Ps eps 0.01 = 826 runs/new, eps 0.002 = 827 runs/V3.
import os, re, tarfile, io, math, filecmp
from collections import defaultdict
D = os.path.dirname(os.path.abspath(__file__)); O = os.path.join(D, 'out', 'runs'); pr = lambda *a: print(*a, flush=True)
kv = lambda s: dict(re.findall(r'([A-Za-z_][A-Za-z_0-9]*)=(\S+)', s))
SYS = {'H0': 'wodor F = 0 (PHEN=1)', 'H1': 'wodor F = 1 (PHEN=2)', 'pPs': 'p-Ps (PHEN=1)', 'oPs': 'o-Ps (PHEN=2)',
       'M0': 'mion F = 0 (PHEN=1)', 'M1': 'mion F = 1 (PHEN=2)'}
EPS = [0.01, 0.005, 0.002, 0.001]
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
A = os.path.join(D, '..')
t826 = tarfile.open(os.path.join(A, 'section-826', 'runs', 'runs826.tar.xz')); t827 = tarfile.open(os.path.join(A, 'section-827', 'runs', 'runs827.tar.xz'))
rd = lambda tf, m: tf.extractfile(m).read()
# H849k
k1 = open(f'{O}/K/H0_4f_0_e0.01.err', 'rb').read() == rd(t826, 'runs/new/H0_4f_0.err')
k2 = open(f'{O}/K/H0_4f_0_e0.002.err', 'rb').read() == rd(t827, 'runs/V3/H0_4f_0.err')
hk = k1 and k2
pr(f"H849k (K1 = 826, K2 = 827 V3 bajt w bajt): {'TAK' if hk else 'NIE'}  (K1 {k1}, K2 {k2})")
data = {}
for s in SYS:
    for e in EPS:
        for i in range(3):
            if s[0] in 'Hpo' and e in (0.01, 0.002):
                if not hk: continue
                m = f'runs/new/{s}_4f_{i}.err' if e == 0.01 else f'runs/V3/{s}_4f_{i}.err'
                data[(s, e, i)] = analyse(rd(t826 if e == 0.01 else t827, m).decode().splitlines())
            else:
                p = f'{O}/e{e}/{s}_4f_{i}.err'
                if os.path.exists(p): data[(s, e, i)] = analyse(open(p).read().splitlines())
def nnls2(xs, ys):   # min sum (y - a x - b x^2)^2, a, b >= 0
    best = None
    S = lambda a, b: sum((y - a*x - b*x*x)**2 for x, y in zip(xs, ys))
    sxx = sum(x*x for x in xs); sx3 = sum(x**3 for x in xs); sx4 = sum(x**4 for x in xs); sxy = sum(x*y for x, y in zip(xs, ys)); sx2y = sum(x*x*y for x, y in zip(xs, ys))
    det = sxx*sx4 - sx3*sx3; cands = []
    if det:
        a = (sxy*sx4 - sx2y*sx3)/det; b = (sxx*sx2y - sx3*sxy)/det
        if a >= 0 and b >= 0: cands.append((a, b))
    cands += [(max(sxy/sxx, 0), 0.0), (0.0, max(sx2y/sx4, 0))]
    return min(cands, key=lambda ab: S(*ab))
P = {}
for s, lab in SYS.items():
    pr(f'  {lab}:')
    Ds = {}
    for e in EPS:
        for tr in ('4->3', '3->2'):
            vals = [x for i in range(3) if (s, e, i) in data for x in data[(s, e, i)][0].get(tr, [])]
            if not vals: pr(f'    eps {e} {tr}: brak'); continue
            m = max(abs(x) for x in vals)
            if tr == '4->3': Ds[e] = m
            pr(f'    eps {e} {tr}: N {len(vals)}; dn ' + ' '.join(f'{x:+.2e}' for x in vals) + f'; max |dn| {m:.3e}')
        refs = [data[(s, e, i)][1] for i in range(3) if (s, e, i) in data]; got = [data[(s, e, i)][2] for i in range(3) if (s, e, i) in data]
        pr(f'    eps {e}: odrzucenia {refs}; do 1s {got}')
        for i in range(3):
            if (s, e, i) in data:
                for sp in data[(s, e, i)][3]: pr(f'      ziarno {i} {sp[:200]}')
    if len(Ds) >= 2:
        xs = [math.log(e) for e in Ds]; ys = [math.log(Ds[e]) for e in Ds]; n = len(xs); mx = sum(xs)/n; my = sum(ys)/n
        p = sum((x-mx)*(y-my) for x, y in zip(xs, ys))/sum((x-mx)**2 for x in xs); P[s] = p
        es = sorted(Ds); a, b = nnls2(es, [Ds[e] for e in es]); fit = a*0.01 + b*1e-4
        loc = ', '.join(f'{es[j+1]}/{es[j]}: {math.log(Ds[es[j+1]]/Ds[es[j]])/math.log(es[j+1]/es[j]):.2f}' for j in range(len(es)-1))
        pr(f'    4->3: p = {p:.3f} ({len(Ds)} pkt); lokalnie {loc}; NNLS a {a:.3e}, b {b:.3e}; f_lin(0,01) {a*0.01/fit if fit else float("nan"):.3f}')
rng = lambda s, lo, hi: s in P and lo <= P[s] <= hi
hp = rng('pPs', 1.7, 2.3) and rng('oPs', 1.7, 2.3)
hh = 'H0' in P and 'H1' in P and P['H0'] <= 1.3 and P['H1'] <= 1.3
hm = all(k in P for k in ('H0', 'H1', 'M0', 'M1', 'pPs', 'oPs')) and P['H0'] - 0.1 <= P['M0'] <= P['pPs'] + 0.1 and P['H1'] - 0.1 <= P['M1'] <= P['oPs'] + 0.1
pr(f"H849p (p-Ps i o-Ps p w [1,7; 2,3]): {'TAK' if hp else 'NIE'}")
pr(f"H849h (GLOWNA: wodor p <= 1,3 w F = 0 i F = 1): {'TAK' if hh else 'NIE'}")
pr(f"H849m (p_H - 0,1 <= p_M <= p_Ps + 0,1; F = 0 z p-Ps, F = 1 z o-Ps): {'TAK' if hm else 'NIE'}")

# Audit 850: source of the payout-window drift floor (print-only DRIFT850).  Usage: LC_ALL=C python3 ana850.py -> out/ana850.txt.
# Written during the runs, before reading them.
import os, re, tarfile, math
D = os.path.dirname(os.path.abspath(__file__)); O = os.path.join(D, 'out', 'runs'); pr = lambda *a: print(*a, flush=True)
kv = lambda s: dict(re.findall(r'([A-Za-z_][A-Za-z_0-9]*)=(\S+)', s))
SYS = {'H0': 'wodor F = 0 (PHEN=1)', 'H1': 'wodor F = 1 (PHEN=2)', 'M0': 'mion F = 0 (PHEN=1)', 'M1': 'mion F = 1 (PHEN=2)',
       'pPs': 'p-Ps (PHEN=1)', 'oPs': 'o-Ps (PHEN=2)'}
PARTS = ['kin_wzgl', 'kin_sm', 'coul', 'dd', 'cd(w)', 'darwin', 'constr']
SEEDS = range(10)
# H850k
k = open(f'{O}/K/H0_4f_0.err', 'rb').read()
t849 = tarfile.open(os.path.join(D, '..', 'section-849', 'runs', 'runs849.tar.xz'))
ref = t849.extractfile('runs/e0.001/H0_4f_0.err').read()
lines_ok = all(b'DRIFT850' in open(f'{O}/W/{s}_4f_{i}.err', 'rb').read() for s in SYS for i in SEEDS if os.path.exists(f'{O}/W/{s}_4f_{i}.err'))
nW = sum(os.path.exists(f'{O}/W/{s}_4f_{i}.err') for s in SYS for i in SEEDS)
hk = (k == ref) and b'DRIFT850' not in k and lines_ok
pr(f"H850k (K = 849 bajt w bajt; DRIFT850 we wszystkich {nW} biegach W, brak w K): {'TAK' if hk else 'NIE'}  (K {k == ref}, wiersze {lines_ok})")
def parse(path):
    dec = []; cur = None
    for ln in open(path):
        if ln.startswith('DRIFT850'):
            t = ln.split(); d = {t[j]: t[j+1] for j in range(1, 31, 2)}
            cur = dict(nr=None, nf=None, D={k2: float(v) for k2, v in d.items() if k2 not in ('started',)}, started=int(d['started']))
            i = t.index('Pp'); cur['Pp'] = [float(x) for x in t[i+1:i+8]]; i = t.index('Dp'); cur['Dp'] = [float(x) for x in t[i+1:i+8]]
            dec.append(cur)
        elif ln.startswith('RULES') and cur is not None and cur['nr'] is None:
            d = kv(ln); cur['nr'] = float(d['n']); cur['nf'] = int(d['n_final'])
    return dec
res = {}; worst = 0.0; nphot = 0
for s, lab in SYS.items():
    rows = []
    for i in SEEDS:
        p = f'{O}/W/{s}_4f_{i}.err'
        if not os.path.exists(p): continue
        dec = parse(p)
        for a, b in zip(dec, dec[1:]):
            if a['nf'] is None or b['nr'] is None or round(a['nr']) != 4 or a['nf'] != 3: continue
            dn = b['nr'] - a['nf']; d = b['D']; n = d['n_orb']; B = -d['E'] * n * n; f = n**3 / (2 * B)
            c = {x: d[x] * f for x in ('P', 'R', 'Dy', 'X')}
            cl = dn - sum(c.values()); worst = max(worst, abs(cl)); nphot += 1
            rows.append(dict(seed=i, dn=dn, c=c, cl=cl, Pp=[x * f for x in b['Pp']], Dp=[x * f for x in b['Dp']], dt=d['dt'], incr=d['incr']))
    res[s] = rows
    if not rows: pr(f'  {lab}: brak fotonow 4 -> 3'); continue
    rms = lambda v: math.sqrt(sum(x * x for x in v) / len(v))
    R0 = rms([r['dn'] for r in rows])
    pr(f'  {lab}: fotony 4 -> 3 N {len(rows)}; rms dn {R0:.3e}; max |dn| {max(abs(r["dn"]) for r in rows):.3e}')
    for r in rows:
        pr(f'    ziarno {r["seed"]}: dn {r["dn"]:+.3e} = P {r["c"]["P"]:+.3e} R {r["c"]["R"]:+.3e} Dy {r["c"]["Dy"]:+.3e} X {r["c"]["X"]:+.3e}'
           f' (domkn. {r["cl"]:+.1e}; przyrosty {r["incr"]:.0f}, dt {r["dt"]:.2e} s)')
    u = {x: (rms([r['c'][x] for r in rows]) / R0 if R0 else float('nan')) for x in ('P', 'R', 'Dy', 'X')}
    pr('    udzialy u_c: ' + ', '.join(f'{x} {v:.3f}' for x, v in u.items()))
    pr('    czesci P (rms, w n): ' + ', '.join(f'{PARTS[j]} {rms([r["Pp"][j] for r in rows]):.2e}' for j in range(7)))
    pr('    czesci Dy (rms, w n): ' + ', '.join(f'{PARTS[j]} {rms([r["Dp"][j] for r in rows]):.2e}' for j in range(7)))
    res[s] = (rows, u)
pr(f"H850z (domkniecie |dn - sum dn_c| <= 1e-6 dla kazdego z {nphot} fotonow 4 -> 3): {'TAK' if worst <= 1e-6 else 'NIE'}  (max {worst:.2e})")
hm = all(isinstance(res.get(s), tuple) and res[s][1]['P'] >= 0.7 for s in ('H0', 'H1', 'M0', 'M1'))
pr(f"H850m (GLOWNA: u_P >= 0,7 w H0, H1, M0, M1): {'TAK' if hm else 'NIE'}  (" + ', '.join(f"{s} {res[s][1]['P']:.3f}" for s in ('H0', 'H1', 'M0', 'M1') if isinstance(res.get(s), tuple)) + ')')

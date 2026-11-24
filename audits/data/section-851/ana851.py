# Audit 851: CREM_PHOTON_INCREMENT_STABLE test.  Usage: LC_ALL=C python3 ana851.py -> out/ana851.txt.
# Written during the runs, after the static probe was read (22:21), before reading any W run.  Parse from ana850.py.
import os, re, tarfile, math, statistics as st
D = os.path.dirname(os.path.abspath(__file__)); O = os.path.join(D, 'out', 'runs'); pr = lambda *a: print(*a, flush=True)
kv = lambda s: dict(re.findall(r'([A-Za-z_][A-Za-z_0-9]*)=(\S+)', s))
SYS = {'H0': 'wodor F = 0 (PHEN=1)', 'H1': 'wodor F = 1 (PHEN=2)', 'M0': 'mion F = 0 (PHEN=1)', 'M1': 'mion F = 1 (PHEN=2)',
       'pPs': 'p-Ps (PHEN=1)', 'oPs': 'o-Ps (PHEN=2)'}
SEEDS = range(10)
# 850 values written in prereg851
MED850 = {'H0': -7.16e-9, 'H1': -7.14e-9, 'M0': -2.20e-9, 'M1': -2.13e-9}
RMS850 = {'H0': 5.43e-5, 'H1': 5.23e-5, 'M0': 1.254e-5, 'M1': 1.231e-5, 'pPs': 9.85e-7, 'oPs': 9.93e-7}
rms = lambda v: math.sqrt(sum(x * x for x in v) / len(v))
# H851k
k = open(f'{O}/K/H0_4f_0.err', 'rb').read()
ref = tarfile.open(os.path.join(D, '..', 'section-850', 'runs', 'runs850.tar.xz')).extractfile('runs/W/H0_4f_0.err').read()
miss = []; nrul = 0
for s in SYS:
    for i in SEEDS:
        t = open(f'{O}/W/{s}_4f_{i}.err', 'rb').read()
        if b'\nRULES' in t or t.startswith(b'RULES'):
            nrul += 1
            if b'PAIR PHOTONSTABLE851 on' not in t: miss.append(f'{s}_{i}')
hk = k == ref and b'PHOTONSTABLE851' not in k and not miss
pr(f"H851k (K = 850 W/H0_4f_0 bajt w bajt; PHOTONSTABLE851 we wszystkich {nrul} biegach W z RULES, brak w K): {'TAK' if hk else 'NIE'}"
   f"  (K {k == ref}, brak wiersza w {miss})")
# static
S = {}
for f in sorted(os.listdir(f'{O}/S')):
    for ln in open(f'{O}/S/{f}'):
        if ln.startswith('S851'):
            t = ln.split(); S[(t[2], int(t[4]))] = dict(mean=float(t[16]), rms=float(t[19]), neg=int(t[25]), N=int(t[8]))
            pr('  ' + ln.strip())
s0 = {'H': S[('proton,electron', 0)]['mean'], 'M': S[('antimuon,electron', 0)]['mean']}
okS0 = all(v < 0 and 1/3 <= v / r <= 3 for v, r in ((s0['H'], -7.2e-9), (s0['M'], -2.2e-9)))
pr(f"H851s0 (sonda bez przelacznika: srednie P / przyrost ujemne, w czynniku 3 od 850 -7,2e-9 / -2,2e-9): {'TAK' if okS0 else 'NIE'}"
   f"  (wodor {s0['H']:+.3e}, mion {s0['M']:+.3e})")
s1 = {'H': S[('proton,electron', 1)]['mean'], 'M': S[('antimuon,electron', 1)]['mean']}
pr(f"H851s (sonda z przelacznikiem: |srednie P / przyrost| <= 1e-10 n): {'TAK' if all(abs(v) <= 1e-10 for v in s1.values()) else 'NIE'}"
   f"  (wodor {s1['H']:+.3e}, mion {s1['M']:+.3e}; Ps {S[('electron,positron', 1)]['mean']:+.3e}, bez {S[('electron,positron', 0)]['mean']:+.3e})")
def parse(path):
    dec = []; cur = None; ev = []; ce = None
    for ln in open(path):
        if ln.startswith('DRIFT850'):
            t = ln.split(); d = {t[j]: t[j+1] for j in range(1, 31, 2)}
            cur = dict(nr=None, nf=None, D={k2: float(v) for k2, v in d.items() if k2 not in ('started',)})
            dec.append(cur)
        elif ln.startswith('LADDER E_rel'):
            d = kv(ln); ce = dict(n=float(d['n']), photon=float(d['photon'])); ev.append(ce)
        if ln.startswith('RULES'):
            d = kv(ln)
            if cur is not None and cur['nr'] is None: cur['nr'] = float(d['n']); cur['nf'] = int(d['n_final'])
            if ce is not None: ce['nr'] = float(d['n']); ce['nf'] = int(d['n_final'])
    nref = sum(1 for e in ev if e['photon'] == 0.0 and e.get('nr', e['n']) >= 1.5)
    got = any(e['photon'] > 0 and e.get('nf') == 1 for e in ev)
    return dec, nref, got
worst = 0.0; nphot = 0; R = {}
for s, lab in SYS.items():
    rows = []; refs = []; gots = []
    for i in SEEDS:
        dec, nref, got = parse(f'{O}/W/{s}_4f_{i}.err'); refs.append(nref); gots.append(got)
        for a, b in zip(dec, dec[1:]):
            if a['nf'] is None or b['nr'] is None or round(a['nr']) != 4 or a['nf'] != 3: continue
            dn = b['nr'] - a['nf']; d = b['D']; n = d['n_orb']; B = -d['E'] * n * n; f = n**3 / (2 * B)
            c = {x: d[x] * f for x in ('P', 'R', 'Dy', 'X')}
            cl = dn - sum(c.values()); worst = max(worst, abs(cl)); nphot += 1
            rows.append(dict(seed=i, dn=dn, c=c, cl=cl, incr=d['incr']))
    if not rows: pr(f'  {lab}: brak fotonow 4 -> 3'); continue
    R0 = rms([r['dn'] for r in rows]); per = [r['c']['P'] / r['incr'] for r in rows if r['incr'] > 0]; med = st.median(per)
    u = {x: rms([r['c'][x] for r in rows]) / R0 for x in ('P', 'R', 'Dy', 'X')}
    R[s] = dict(rms=R0, med=med, u=u)
    pr(f'  {lab}: fotony 4 -> 3 N {len(rows)}; rms dn {R0:.3e} (850 {RMS850[s]:.3e}, x{RMS850[s] / R0:.1f} mniej); max |dn| '
       f'{max(abs(r["dn"]) for r in rows):.3e}; mediana P / przyrost {med:+.2e}' + (f' (850 {MED850[s]:+.2e})' if s in MED850 else '')
       + f'; odrzucenia {refs}; do 1s {sum(gots)}/10')
    for r in rows:
        pr(f'    ziarno {r["seed"]}: dn {r["dn"]:+.3e} = P {r["c"]["P"]:+.3e} R {r["c"]["R"]:+.3e} Dy {r["c"]["Dy"]:+.3e} X {r["c"]["X"]:+.3e}'
           f' (domkn. {r["cl"]:+.1e}; przyrosty {r["incr"]:.0f})')
    pr('    udzialy u_c: ' + ', '.join(f'{x} {v:.3f}' for x, v in u.items()))
hd = all(s in R and abs(R[s]['med']) <= abs(MED850[s]) / 10 and R[s]['rms'] <= RMS850[s] / 5 for s in MED850)
pr(f"H851d (GLOWNA: |mediana P / przyrost| >= 10 x mniej i rms dn >= 5 x mniej niz 850 w H0, H1, M0, M1): {'TAK' if hd else 'NIE'}  ("
   + ', '.join(f"{s} P x{abs(MED850[s] / R[s]['med']) if R[s]['med'] else float('inf'):.0f}, rms x{RMS850[s] / R[s]['rms']:.1f}" for s in MED850 if s in R) + ')')
hp = all(s in R and 0.5 <= R[s]['rms'] / RMS850[s] <= 2 for s in ('pPs', 'oPs'))
pr(f"H851p (Ps rms dn w czynniku 2 od 850): {'TAK' if hp else 'NIE'}  (" + ', '.join(f"{s} {R[s]['rms'] / RMS850[s]:.2f}" for s in ('pPs', 'oPs') if s in R) + ')')
pr(f"H851z (domkniecie <= 1e-6 dla kazdego z {nphot} fotonow 4 -> 3): {'TAK' if worst <= 1e-6 else 'NIE'}  (max {worst:.2e})")

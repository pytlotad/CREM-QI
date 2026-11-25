# Audit 854: periapsis spin-orbit form (Thomas) sensitivity.  Usage: LC_ALL=C python3 ana854.py -> out/ana854.txt.  Written during the runs, before any result.
import os, re, math
SP = '/tmp/claude-1000/-home-teddy-Projekty-CREM-QI/4fa83d5c-8878-4cfd-898e-1a9b53a01ab2/scratchpad'
REG = os.environ.get('R854REG', SP + '/r854/reg'); pr = lambda *a: print(*a, flush=True)
rd = lambda t, x: open(f'{REG}/{t}/{x}', 'rb').read()
PM = re.compile(r'^(.*?)([-+]?\d[\d.]*(?:e[-+]?\d+)?)\s*\+/-\s*([\d.]+(?:e[-+]?\d+)?)\s*(\w*)')
def pm(t):
    out = []
    for ln in open(f'{REG}/{t}/out.txt'):
        m = PM.match(ln.strip())
        if m: out.append((m.group(1).strip(), float(m.group(2)), float(m.group(3)), m.group(4)))
    return out
def grab(t):
    s = open(f'{REG}/{t}/out.txt').read()
    a = re.search(r'lifetime at the final state ([\d.]+ \w+)', s); e = re.search(r'E\[T\] from the preparation \(cascade \+ annihilation\) ([\d.]+ \w+)', s)
    m = re.search(r'mean of completed runs\s+([\d.e+-]+) \+/- ([\d.e+-]+) (\w+)', s)
    return (float(m.group(1)), float(m.group(2)), m.group(3)) if m else None, a.group(1) if a else None, e.group(1) if e else None
k = all(rd('K_2', x) == rd('OFF_2', x) for x in ('out.txt', 'err.txt'))
on = {t: b'PAIR PERISO854 on' in rd(f'ON_{t}', 'err.txt') for t in ('1', '2', 'L2')}
off = {t: b'PERISO854' in rd(f'OFF_{t}', 'err.txt') for t in ('1', '2', 'L2')}
pr(f"H854k (K = OFF_2 bajt w bajt; PERISO854 w ON 2 i L2, brak w OFF): {'TAK' if k and on['2'] and on['L2'] and not any(off.values()) else 'NIE'}  (K {k}; ON {on}; OFF {off})")
for ln in rd('ON_2', 'err.txt').decode(errors='replace').split('\n'):
    if 'PERISO854' in ln: pr('  ' + ln.strip()); break
def compare(t, lab):
    o, n = pm(f'OFF_{t}'), pm(f'ON_{t}'); same = rd(f'OFF_{t}', 'out.txt') == rd(f'ON_{t}', 'out.txt')
    pr(f'  {lab} (zestaw {t}): out.txt ON = OFF bajt w bajt: {same}; linii "+/-": {len(o)} / {len(n)}')
    worst = 0.0
    for (lo, xo, so, u), (ln_, xn, sn, _) in zip(o, n):
        z = abs(xn - xo) / so if so > 0 else (0.0 if xn == xo else float('inf')); worst = max(worst, z)
        if xn != xo: pr(f'    {lo[:80]}: {xo:.6g} +/- {so:.3g} {u} -> {xn:.6g} +/- {sn:.3g}  ({(xn / xo - 1) * 100 if xo else float("nan"):+.3f} %, {z:.3f} Y_OFF)')
    go, gn = grab(f'OFF_{t}'), grab(f'ON_{t}')
    pr(f'    zapadanie (mean of completed runs) OFF {go[0]} ON {gn[0]}; lifetime at final state {go[1]} -> {gn[1]}; E[T] {go[2]} -> {gn[2]}')
    return worst, go, gn, same
wp, gp_o, gp_n, sp = compare('1', 'p-Ps (PHEN=1)')
pr(f"H854p (p-Ps: kazde X +/- Y ON = OFF do 0,1 Y; czasy bez zmiany): {'TAK' if wp <= 0.1 and gp_o[1:] == gp_n[1:] else 'NIE'}  (max {wp:.3f} Y)")
res = []
for t, lab in (('2', 'o-Ps (PHEN=2)'), ('L2', 'o-Ps n = 2 (PHEN=2, --level 2)')):
    w, go, gn, s = compare(t, lab)
    d = abs(gn[0][0] - go[0][0]) / go[0][1] if go[0] and gn[0] else float('nan'); res.append((t, d))
pr(f"H854m (GLOWNA: o-Ps |zapadanie ON - OFF| <= Y_OFF w 2 i L2): {'TAK' if all(d <= 1 for _, d in res) else 'NIE'}  (" + ', '.join(f'{t} {d:.3f} Y' for t, d in res) + ')')

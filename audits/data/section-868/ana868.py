# Audit 868 (prototype quasi-3D, E1): verdicts H868c, k, a, t, l, n -> out/ana868.txt.  Written before the results; definitions from ../section-773/ana773.py.
import os, math, numpy as np
from collections import Counter
D = os.path.dirname(os.path.abspath(__file__)); O = os.path.join(D, 'out', 'runs'); pr = lambda *a: print(*a, flush=True)
S = {}
for ln in open(os.path.join(D, '..', 'section-773', 'starts773.txt')):
    i, ini = ln.split(); S[i] = [(int(a), int(b), float(c)) for a, b, c, _ in (x.split(':') for x in ini.split(';'))]
def lattice(q3d):
    pts = []
    for m in range(-6, 7):
        Ls = ([0.5, -0.5] if m == 0 else [math.copysign(abs(m) + 0.5, m)]) if q3d else [float(m)]
        for L in Ls: pts.append((m, L))
    return pts
def nearest(j, l, q3d=True):
    best = (np.inf, None); nr = max(0, round(j - 0.5))
    for m, L in lattice(q3d):
        d = math.hypot(j - nr - 0.5, l - L)
        if d < best[0]: best = (d, (nr, m))
    return best
def ev(pref, i, per, GT):
    q = np.loadtxt(os.path.join(O, f"{pref}_{i}.txt"), comments='#', ndmin=2); assert len(q) == per, (pref, i, len(q))
    e = q[-20:, 2:4].mean(0); k = int(0.8 * per); b = q[k - 10:k + 10, 2:4].mean(0); dr = math.hypot(*(e - b))
    res = 500 * math.exp(-2 * GT * per) / q[-20:, 4].mean(); conv = dr < 0.005 and res < 1e-3
    d, t = nearest(*e); d2, t2 = nearest(*e, q3d=False); wm = q[-1, 5]; win = (int(q[-1, 6]), int(q[-1, 7])); one = wm >= 0.95 and d < 0.04
    fin = bool(np.all(np.isfinite(q))); rmin = float(q[:, 8].min()) if q.shape[1] > 8 else float('nan')
    return dict(j=e[0], l=e[1], dr=dr, res=res, conv=conv, d=d, tor=t, d2=d2, wm=wm, win=win, one=one, ok=conv and one, fin=fin, rmin=rmin)
def block(pref, ids, per, GT, nm):
    pr(f"=== {nm} ===");  pr("  id  J_g      L_g      odl.    torus     odl.2D  w_max   zwyc.     dryf    resztka   zb. sukces 1-mod  min r")
    R = {}
    for i in ids:
        r = ev(pref, i, per, GT); R[i] = r
        pr(f"  {i}  {r['j']:.5f}  {r['l']:+.5f}  {r['d']:.4f}  {str(r['tor']):9s} {r['d2']:.4f}  {r['wm']:.4f}  {str(r['win']):9s} {r['dr']:.4f}  {r['res']:.1e}  {'T' if r['conv'] else 'N'}   {'T' if r['ok'] else 'N'}      {'T' if r['one'] else 'N'}      {r['rmin']:.4f}")
    nc = sum(r['conv'] for r in R.values()); ns = sum(r['ok'] for r in R.values()); no = sum(r['one'] for r in R.values())
    pr(f"  zbiezne {nc}/{len(ids)}, uciete {len(ids) - nc}, sukces {ns}/{nc} = {100 * ns / max(nc, 1):.1f} %, jednomodowe na EBK (wszystkie) {no}/{len(ids)}")
    pr("  zwyciezcy-sukcesy: " + ', '.join(f"{k}: {v}" for k, v in Counter(r['win'] for r in R.values() if r['ok']).most_common()))
    return R, nc, ns, no
# H868c
ref = open(os.path.join(D, '..', 'section-773', 'out', 'runs', 'G_01.txt')).read().split('\n')
hdr = [x for x in ref if x.startswith('#')]; rec = [x for x in ref if x and not x.startswith('#')][:300]
same = open(os.path.join(O, 'C_01.txt')).read() == '\n'.join(hdr + rec) + '\n'
pr(f"H868c (C = naglowek + pierwsze 300 rekordow 773 G_01 bajt w bajt): {'TAK' if same else 'NIE'}")
# H868k
el = [x for x in open(os.path.join(O, 'E_01.txt')) if x.startswith('# norm') or x.startswith('# levels')]
for x in el: pr('  ' + x.strip())
nd = float(el[0].split('=')[1].split()[0]); ld = float(el[1].split('=')[1].split()[0])
pr(f"H868k (|norma - 1| < 1e-6, |E - Bohr| < 1e-12 przy beta 0): {'TAK' if nd < 1e-6 and ld < 1e-12 else 'NIE'}")
Q, ncq, nsq, noq = block('Q', [f"{i:02d}" for i in range(1, 41)], 12000, 0.001, '(Q) Q3D, beta 0,2, CAUST = CAUSTT = -pi/2, 12000 okresow, starty 773/01-40')
pr(f"H868a (GLOWNA; jednomodowe na EBK >= 24/40 i sukces >= 80 % zbieznych): {'TAK' if noq >= 24 and ncq > 0 and nsq >= 0.8 * ncq else 'NIE'}  (773 G: 24/40, 24/30)")
T, nct, nst, not_ = block('QT', [f"{i:02d}" for i in range(41, 61)], 6000, 0.001, '(QT) Q3D, beta 0,2, CAUSTT 0, 6000 okresow, starty 773/41-60')
pr(f"H868t (bez licznika theta: jednomodowe na EBK kwazi-3D <= 2/20): {'TAK' if not_ <= 2 else 'NIE'} ({not_}/20);  jednomodowe na siatce 2D (w_max >= 0,95, odl.2D < 0,04): "
   f"{sum(r['wm'] >= 0.95 and r['d2'] < 0.04 for r in T.values())}/20")
N, ncn, nsn, non = block('QN', [f"{i:02d}" for i in range(1, 41)], 12000, 0.001, '(QN) Q3D, beta -0,05, CAUST = CAUSTT = -pi/2, 12000 okresow, starty 773/01-40')
okl = all(r['fin'] and r['rmin'] > 0 for r in N.values())
pr(f"H868l (QN: skonczone i min r > 0 w 40/40): {'TAK' if okl else 'NIE'}  (min r min {min(r['rmin'] for r in N.values()):.4f})")
pr(f"H868n (QN: jednomodowe na EBK >= 24/40 i sukces >= 80 % zbieznych): {'TAK' if non >= 24 and ncn > 0 and nsn >= 0.8 * ncn else 'NIE'};  "
   f"l zwyciezcow-sukcesow: {dict(Counter(abs(r['win'][1]) for r in N.values() if r['ok']))}")

# Audit 869 (prototype quasi-3D, cause of the E1 failure): verdicts H869c, s, n, d, t, f, l -> out/ana869.txt.  Written before the results;
# definitions copied from ../section-868/ana868.py; 868 baselines on starts 773/01-20 from the prereg.
import os, math, tarfile, numpy as np
from collections import Counter
D = os.path.dirname(os.path.abspath(__file__)); O = os.path.join(D, 'out', 'runs'); pr = lambda *a: print(*a, flush=True)
def lattice():
    pts = []
    for m in range(-6, 7):
        for L in ([0.5, -0.5] if m == 0 else [math.copysign(abs(m) + 0.5, m)]): pts.append((m, L))
    return pts
def nearest(j, l):
    best = (np.inf, None); nr = max(0, round(j - 0.5))
    for m, L in lattice():
        d = math.hypot(j - nr - 0.5, l - L)
        if d < best[0]: best = (d, (nr, m))
    return best
def load(pref, i, per):
    q = np.loadtxt(os.path.join(O, f"{pref}_{i}.txt"), comments='#', ndmin=2); assert len(q) == per, (pref, i, len(q)); return q
def ev(pref, i, per=12000, GT=0.001):
    q = load(pref, i, per)
    e = q[-20:, 2:4].mean(0); k = int(0.8 * per); b = q[k - 10:k + 10, 2:4].mean(0); dr = math.hypot(*(e - b))
    res = 500 * math.exp(-2 * GT * per) / q[-20:, 4].mean(); conv = dr < 0.005 and res < 1e-3
    d, t = nearest(*e); wm = q[-1, 5]; win = (int(q[-1, 6]), int(q[-1, 7])); one = wm >= 0.95 and d < 0.04
    fin = bool(np.all(np.isfinite(q))); rmin = float(q[:, 8].min()); m0 = float((q[:1000, 7] == 0).mean())
    return dict(j=e[0], l=e[1], dr=dr, res=res, conv=conv, d=d, tor=t, wm=wm, win=win, one=one, ok=conv and one, fin=fin, rmin=rmin, m0=m0)
def block(pref, ids, nm):
    pr(f"=== {nm} ===");  pr("  id  J_g      L_g      odl.    torus     w_max   zwyc.     dryf    zb. sukces 1-mod  min r   m0(1000)")
    R = {}
    for i in ids:
        r = ev(pref, i); R[i] = r
        pr(f"  {i}  {r['j']:.5f}  {r['l']:+.5f}  {r['d']:.4f}  {str(r['tor']):9s} {r['wm']:.4f}  {str(r['win']):9s} {r['dr']:.4f}  {'T' if r['conv'] else 'N'}   {'T' if r['ok'] else 'N'}      {'T' if r['one'] else 'N'}      {r['rmin']:.4f}  {r['m0']:.3f}")
    nc = sum(r['conv'] for r in R.values()); ns = sum(r['ok'] for r in R.values()); no = sum(r['one'] for r in R.values())
    pr(f"  zbiezne {nc}/{len(ids)}, sukces {ns}/{nc} = {100 * ns / max(nc, 1):.1f} %, jednomodowe na EBK {no}/{len(ids)}")
    pr("  zwyciezcy-sukcesy: " + ', '.join(f"{k}: {v}" for k, v in Counter(r['win'] for r in R.values() if r['ok']).most_common()))
    pr("  zwyciezcy wszyscy: " + ', '.join(f"{k}: {v}" for k, v in Counter(r['win'] for r in R.values()).most_common()))
    return R, nc, ns, no
def klasa(no, nc, ns, ref):
    if no >= 12 and nc > 0 and ns >= 0.8 * nc: return 'PRZYCZYNA'
    if abs(no - ref) <= 3: return 'BEZ WPLYWU'
    return 'POSREDNIE'
ids20 = [f"{i:02d}" for i in range(1, 21)]; ids10 = ids20[:10]
# H869c
with tarfile.open(os.path.join(D, '..', 'section-868', 'runs', 'runs868.tar.xz')) as tf: ref = tf.extractfile('runs/Q_01.txt').read().decode().split('\n')
hdr = [x for x in ref if x.startswith('#')]; rec = [x for x in ref if x and not x.startswith('#')][:300]
same = open(os.path.join(O, 'C_01.txt')).read() == '\n'.join(hdr + rec) + '\n'
pr(f"H869c (C = naglowek + pierwsze 300 rekordow 868 Q_01 bajt w bajt): {'TAK' if same else 'NIE'}")
S, ncs, nss, nos = block('S', ids20, '(S) THSMOOTH, beta 0,2, starty 773/01-20  [868 Q: 6/20, sukces 4/12]')
pr(f"H869s (GLOWNA; jednomodowe >= 12/20 i sukces >= 80 % zbieznych): {'TAK' if nos >= 12 and ncs > 0 and nss >= 0.8 * ncs else 'NIE'}")
N, ncn, nsn, non = block('SN', ids20, '(SN) THSMOOTH, beta -0,05, starty 773/01-20  [868 QN: 5/20, sukces 4/8]')
pr(f"H869n (jednomodowe >= 12/20 i sukces >= 80 % zbieznych): {'TAK' if non >= 12 and ncn > 0 and nsn >= 0.8 * ncn else 'NIE'};  l zwyciezcow-sukcesow: {dict(Counter(abs(r['win'][1]) for r in N.values() if r['ok']))}")
Dd, ncd, nsd, nod = block('D', ids20, '(D) SIGDOM, CAUSTT -pi/2, beta -0,05, starty 773/01-20  [868 QN: 5/20]')
kd = klasa(nod, ncd, nsd, 5); pr(f"H869d (D klasa BEZ WPLYWU): {'TAK' if kd == 'BEZ WPLYWU' else 'NIE'} (klasa {kd})")
T, nct, nst, not_ = block('T', ids20, '(T) TR2D, CAUSTT -pi/2, beta 0,2, starty 773/01-20  [868 Q: 6/20]')
kt = klasa(not_, nct, nst, 6); pr(f"H869t (T klasa BEZ WPLYWU): {'TAK' if kt == 'BEZ WPLYWU' else 'NIE'} (klasa {kt})")
P3, _, _, no3 = block('P3', ids10, '(P3) CAUSTT -pi/3, beta 0,2, starty 773/01-10')
P23, _, _, no23 = block('P23', ids10, '(P23) CAUSTT -2pi/3, beta 0,2, starty 773/01-10')
pr(f"H869f (P3 i P23 jednomodowe <= 3/10 kazde): {'TAK' if no3 <= 3 and no23 <= 3 else 'NIE'} ({no3}/10, {no23}/10)")
Z = {}
for pref in ('Z0', 'ZS'):
    Z[pref], *_ = block(pref, ids10, f"({pref}) starty l = 0 (starts869), beta -0,05, {'licznik theta -pi/2' if pref == 'Z0' else 'THSMOOTH'}")
okl = all(all(r['fin'] and r['rmin'] > 0 for r in R.values()) and np.median([r['m0'] for r in R.values()]) >= 0.5 for R in Z.values())
for pref, R in Z.items(): pr(f"  {pref}: min r min {min(r['rmin'] for r in R.values()):.4f}, mediana m0(1000) {np.median([r['m0'] for r in R.values()]):.3f}, zwyciezcy z m = 0 na koncu {sum(r['win'][1] == 0 for r in R.values())}/10")
pr(f"H869l (Z0 i ZS: skonczone, min r > 0 i mediana m0(1000) >= 0,5): {'TAK' if okl else 'NIE'}")

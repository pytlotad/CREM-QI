# Audit 778: is the lead taken by the mode resonant with the CURRENT torus?  Verdicts H778a-e -> out/ana778.txt
import os, numpy as np
from collections import defaultdict
D = os.path.dirname(os.path.abspath(__file__)); A = os.path.join(D, '..'); O = os.path.join(D, 'out', 'runs'); Q = os.path.join(D, 'out', 'pred')
pr = lambda *a: print(*a, flush=True)
MODES = [(nr, m) for m in range(-6, 7) for nr in range(13)]       # order of guide772 (m outer, nr inner)
ORIG = {'774G': 'section-774/out/runs/G_{}.txt', '775B05': 'section-775/out/runs/B05_{}.txt', '775B80': 'section-775/out/runs/B80_{}.txt'}
def nearest(j, l):
    best = (np.inf, None)
    for m in range(-6, 7):
        nr = max(0, round(j - 0.5)); d = np.hypot(j - nr - 0.5, l - m)
        if d < best[0]: best = (d, (nr, m))
    return best
# H778e
c = open(os.path.join(O, 'C_ctrl.txt')).read(); c0 = open(os.path.join(A, 'section-775', 'out', 'runs', 'C_ctrl.txt')).read()
okC = c == c0
J = [l.split() for l in open(os.path.join(D, 'jobs778.txt'))]
same = 0
for a, i, b, kd, _ in J:
    r = open(os.path.join(O, f"{a}_{i}.txt")).read().splitlines(); o = open(os.path.join(A, ORIG[a].format(i))).read().splitlines()
    same += r == o[:len(r)] and len(r) == int(kd) + 22
pr(f"H778e (kontrola kodu): bez DUMP = 775 C_ctrl bajt w bajt: {'TAK' if okC else 'NIE'}; powtorki = poczatek oryginalu bajt w bajt: {same}/{len(J)} -> {'TAK' if okC and same == len(J) else 'NIE'}")
def pred(a, i, off):
    fn = os.path.join(Q, f"{a}_{i}_m{off}.txt")
    if not os.path.exists(fn): return None
    return [tuple(int(x) for x in ln.split()[:2]) for ln in open(fn) if not ln.startswith('#')]
R = []
for a, i, b, kd, _ in J:
    kd = int(kd); q = np.loadtxt(os.path.join(O, f"{a}_{i}.txt"), comments='#', ndmin=2); W = np.loadtxt(os.path.join(O, f"{a}_{i}.dump"), ndmin=2)[:, 1:]
    o = np.loadtxt(os.path.join(A, ORIG[a].format(i)), comments='#', ndmin=2); win = nearest(*o[-20:, 2:4].mean(0))[1]; iw = MODES.index(win)
    assert (int(q[kd - 1, 6]), int(q[kd - 1, 7])) == win
    kp = kd - 10
    def at(k):
        w = W[k - 1]; order = list(np.argsort(-w)); return dict(w=w[iw], rank=order.index(iw) + 1, top=[MODES[j] for j in order[:3]], lead=MODES[order[0]], wl=w[order[0]])
    r = dict(a=a, i=i, b=float(b), kd=kd, win=win, pre=at(kp), dec=at(kd), m50=at(kd - 50) if kd > 50 else None,
             near=nearest(*q[kp - 1, 2:4])[1], dom=(int(q[kp - 1, 6]), int(q[kp - 1, 7])), PD={off: pred(a, i, off) for off in (50, 10, 0)})
    prev = r['pre']['lead']; ip = MODES.index(prev)
    r['dwin'] = W[kd - 1, iw] - W[kp - 1, iw]; r['dprev'] = W[kd - 1, ip] - W[kp - 1, ip] if prev != win else 0.0; r['prev'] = prev
    r['P2'] = r['pre']['top'][1]
    R.append(r)
n = len(R)
pr(f"Biegi w mianowniku: {n} (11 <= k_d <= 2500; poza: patrz prereg/jobs)")
pr("=== szczegoly ===")
pr("  ramie  id  zwyc.     k_d   P_D(k-10) P_near    P_dom     P_2       ranga_PD  w_zw(k-10) ranga  w_zw(k_d) lider(k-10) w_lid  dw_zw   dw_lid")
for r in R:
    pd = r['PD'][10]; rk = pd.index(r['win']) + 1 if pd else -1
    pr(f"  {r['a']:6s} {r['i']}  {str(r['win']):9s} {r['kd']:5d}  {str(pd[0] if pd else None):9s} {str(r['near']):9s} {str(r['dom']):9s} {str(r['P2']):9s} {rk:5d}     "
       f"{r['pre']['w']:.3f}      {r['pre']['rank']:3d}    {r['dec']['w']:.3f}     {str(r['prev']):9s}  {r['pre']['wl']:.3f}  {r['dwin']:+.3f}  {r['dprev']:+.3f}")
h = lambda f: sum(f(r) for r in R)
hPD = h(lambda r: r['PD'][10] is not None and r['PD'][10][0] == r['win']); hN = h(lambda r: r['near'] == r['win']); hDo = h(lambda r: r['dom'] == r['win']); h2 = h(lambda r: r['P2'] == r['win'])
pr(f"Trafienia w k_d - 10: P_D {hPD}/{n} = {100 * hPD / n:.1f} %, P_near {hN}/{n} = {100 * hN / n:.1f} %, P_dom {hDo}/{n}, P_2 {h2}/{n} = {100 * h2 / n:.1f} %")
ok = hPD >= 0.6 * n and 100 * (hPD - hN) / n >= 15
pr(f"H778a (P_D(k_d-10) >= 60 % i >= P_near + 15 pp): {'TAK' if ok else 'NIE'}")
pr("H778b (raport):")
pr(f"  zwyciezca w top-3 wag w k_d - 10: {h(lambda r: r['pre']['rank'] <= 3)}/{n}; ranga mediana {np.median([r['pre']['rank'] for r in R]):.0f}; waga mediana {np.median([r['pre']['w'] for r in R]):.3f}")
M = [r for r in R if r['m50']]
if M: pr(f"  w k_d - 50 ({len(M)} biegow): top-3 {sum(r['m50']['rank'] <= 3 for r in M)}/{len(M)}, waga mediana {np.median([r['m50']['w'] for r in M]):.3f}")
pr(f"  w k_d: waga zwyciezcy mediana {np.median([r['dec']['w'] for r in R]):.3f}")
up = h(lambda r: r['dwin'] > 0 and r['dprev'] >= -abs(r['dwin'])); dn = h(lambda r: r['dprev'] < 0 and -r['dprev'] > abs(r['dwin']))
pr(f"  przejecie prowadzenia k_d-10 -> k_d: wzrost zwyciezcy dominuje {up}/{n}, spadek lidera dominuje {dn}/{n}; mediana dw_zw {np.median([r['dwin'] for r in R]):+.3f}, dw_lid {np.median([r['dprev'] for r in R]):+.3f}")
pr("H778c (raport):")
for off in (50, 10, 0):
    S = [r for r in R if r['PD'][off] is not None]
    pr(f"  P_D w k_d - {off}: {sum(r['PD'][off][0] == r['win'] for r in S)}/{len(S)}; top-3 {sum(r['win'] in r['PD'][off][:3] for r in S)}/{len(S)}; mediana rangi {np.median([r['PD'][off].index(r['win']) + 1 for r in S]):.0f}")
for a in ('774G', '775B05', '775B80'):
    S = [r for r in R if r['a'] == a]
    if S: pr(f"  {a}: n {len(S)}, P_D {sum(r['PD'][10][0] == r['win'] for r in S)}, P_near {sum(r['near'] == r['win'] for r in S)}, P_2 {sum(r['P2'] == r['win'] for r in S)}")

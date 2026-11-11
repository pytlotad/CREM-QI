# Audit 780: scale test -- verdicts H780a-c, e -> out/ana780.txt
import os, numpy as np
from collections import Counter
D = os.path.dirname(os.path.abspath(__file__)); A = os.path.join(D, '..'); O = os.path.join(D, 'out', 'runs')
pr = lambda *a: print(*a, flush=True)
def nearest(j, l):
    best = (np.inf, None)
    for m in range(-6, 7):
        nr = max(0, round(j - 0.5)); d = np.hypot(j - nr - 0.5, l - m)
        if d < best[0]: best = (d, (nr, m))
    return best
load = lambda fn: np.loadtxt(fn, comments='#', ndmin=2)
def state(q, k):                                    # k = 1-based period
    j, l = q[k - 20:k, 2:4].mean(0); d, tor = nearest(j, l); dom = (int(q[k - 1, 6]), int(q[k - 1, 7])); w = q[k - 1, 5]
    kd = k
    while kd > 1 and (int(q[kd - 2, 6]), int(q[kd - 2, 7])) == dom: kd -= 1
    return dict(j=j, l=l, d=d, tor=tor, dom=dom, w=w, one=w >= 0.95 and d < 0.04, Lh=abs(l - tor[1]) < 0.01, kd=kd)
# H780e
k0 = open(os.path.join(O, 'K_ctrl.txt')).read(); k1 = open(os.path.join(A, 'section-775', 'out', 'runs', 'C_ctrl.txt')).read()
pr(f"H780e (kontrola kodu: K = 775 C_ctrl bajt w bajt): {'TAK' if k0 == k1 else 'NIE'}")
# A
pr("=== A: zbieznosc kroku (start z czystego modu, GT 0,0003, 4000 okresow) ===")
pr("  beta  mod      e      d(2000)  d(4000)  d(8000)  skok(2000/4000/8000)  koniec(8000)")
HI = [(0.05, 1, 0), (0.05, 2, 0), (0.2, 2, 0)]; CT = [(0.05, 0, 0), (0.05, 0, 2), (0.2, 1, 0), (0.2, 0, 2)]
RA = {}
for b, nr, m in HI + CT:
    Lam = np.sqrt(m * m + 2 * b); e = np.sqrt(1 - (Lam / (nr + 0.5 + Lam)) ** 2); r = {}
    for s in (2000, 4000, 8000):
        q = load(os.path.join(O, f"A_b{b}_{nr}_{m}_s{s}.txt")); assert len(q) == 4000
        j, l = q[-20:, 2:4].mean(0); dom = (int(q[-1, 6]), int(q[-1, 7]))
        r[s] = dict(d=np.hypot(j - nr - 0.5, l - m), jump=dom != (nr, m), j=j, l=l, w=q[-1, 5])
    RA[(b, nr, m)] = r
    pr(f"  {b:<5} {str((nr, m)):8s} {e:.3f}  {r[2000]['d']:.5f}  {r[4000]['d']:.5f}  {r[8000]['d']:.5f}  {'/'.join('T' if r[s]['jump'] else 'N' for s in (2000, 4000, 8000))}"
       f"                 J {r[8000]['j']:.4f} L {r[8000]['l']:+.4f} w {r[8000]['w']:.3f}")
hi = sum(not any(RA[t][s]['jump'] for s in RA[t]) and abs(RA[t][2000]['d'] - RA[t][8000]['d']) > 0.005 for t in HI)
cts = [t for t in CT if not any(RA[t][s]['jump'] for s in RA[t])]
ct = sum(abs(RA[t][2000]['d'] - RA[t][8000]['d']) < 0.002 for t in cts)
okA = hi >= 2 and len(cts) >= 2 and ct == len(cts)
pr(f"H780a: duzy mimosrod |d(2000) - d(8000)| > 0,005 w {hi}/3; kontrole bez skoku {len(cts)}/4, z nich < 0,002: {ct}/{len(cts)} -> {'TAK' if okA else 'NIE'}")
# B, C
ids = [f"{k:02d}" for k in range(1, 21)]
def arm(pat, k):
    return {i: state(load(pat.format(i)), k) for i in ids}
ARMS = {'773G (GT 0,001)': (os.path.join(A, 'section-773', 'out', 'runs', 'G_{}.txt'), 3000, 1000),
        '774G (GT 0,0003)': (os.path.join(A, 'section-774', 'out', 'runs', 'G_{}.txt'), 10000, 3333.33),
        'B (GT 0,0001)': (os.path.join(O, 'B_{}.txt'), 30000, 10000),
        '775B05 (beta 0,05, GT 0,0003)': (os.path.join(A, 'section-775', 'out', 'runs', 'B05_{}.txt'), 10000, 3333.33),
        'C (beta 0,0125, GT 0,000075)': (os.path.join(O, 'C_{}.txt'), 40000, 13333.33)}
R = {a: arm(p, k) for a, (p, k, _) in ARMS.items()}
S = {a: sum(r['one'] for r in R[a].values()) for a in R}
pr("=== B, C: stan w 3 tau (starty 01-20) ===")
for a, (p, k, tau) in ARMS.items():
    V = list(R[a].values()); O1 = [r for r in V if r['one']]
    pr(f"  {a:30s}: jednomodowe na EBK {S[a]}/20; w_max >= 0,95 {sum(r['w'] >= 0.95 for r in V)}/20; mediana odl. jednomodowych {np.median([r['d'] for r in O1]) if O1 else float('nan'):.4f}; "
       f"trafienia L {sum(r['Lh'] for r in O1)}/{len(O1)}; mediana k_d/tau {np.median([r['kd'] / tau for r in V]):.3f}")
    pr(f"      zwyciezcy: " + ', '.join(f"{x}: {c}" for x, c in Counter(r['dom'] for r in V).most_common()))
s1, s3, s10 = S['773G (GT 0,001)'], S['774G (GT 0,0003)'], S['B (GT 0,0001)']
pr(f"H780b (niemalejace 0,001 -> 0,0003 -> 0,0001 i S(0,0001) >= S(0,001) + 3): {s1} -> {s3} -> {s10} -> {'TAK' if s1 <= s3 <= s10 and s10 >= s1 + 3 else 'NIE'}")
sc, sb = S['C (beta 0,0125, GT 0,000075)'], S['775B05 (beta 0,05, GT 0,0003)']
pr(f"H780c (|S_C - S_B05| <= 3): S_C {sc}, S_B05 {sb} -> {'TAK' if abs(sc - sb) <= 3 else 'NIE'}")
pr("Raport: pary startow (zwyciezca w 3 tau): 773G / 774G / B ; 775B05 / C   (* = jednomodowy na EBK)")
f = lambda r: f"{r['dom']}{'*' if r['one'] else ' '}"
for i in ids: pr(f"  {i}  {f(R['773G (GT 0,001)'][i])}  {f(R['774G (GT 0,0003)'][i])}  {f(R['B (GT 0,0001)'][i])}   ;  {f(R['775B05 (beta 0,05, GT 0,0003)'][i])}  {f(R['C (beta 0,0125, GT 0,000075)'][i])}")
same_b = sum(R['774G (GT 0,0003)'][i]['dom'] == R['B (GT 0,0001)'][i]['dom'] for i in ids); same_c = sum(R['775B05 (beta 0,05, GT 0,0003)'][i]['dom'] == R['C (beta 0,0125, GT 0,000075)'][i]['dom'] for i in ids)
pr(f"  ten sam zwyciezca 774G / B: {same_b}/20; 775B05 / C: {same_c}/20")
x = np.log10([0.001, 0.0003, 0.0001]); y = np.array([s1, s3, s10]); c = np.polyfit(x, y, 1)
md = [np.median([r['d'] for r in R[a].values() if r['one']] or [np.nan]) for a in ('773G (GT 0,001)', '774G (GT 0,0003)', 'B (GT 0,0001)')]
pr(f"Ekstrapolacja (opis): S(log10 GT) liniowo: nachylenie {c[0]:+.2f} na dekade; w GT 7,6e-7 -> {np.polyval(c, np.log10(7.6e-7)):.1f} (obciete do 20); mediana odl. jednomodowych {md[0]:.4f} / {md[1]:.4f} / {md[2]:.4f}")

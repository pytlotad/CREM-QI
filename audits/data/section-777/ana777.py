# Audit 777: when is the winning mode decided?  Verdicts H777a-d -> out/ana777.txt (offline, existing runs of 772-775)
import os, numpy as np
from collections import defaultdict
D = os.path.dirname(os.path.abspath(__file__)); A = os.path.join(D, '..')
pr = lambda *a: print(*a, flush=True)
def starts(fn):
    S = {}
    for ln in open(fn):
        i, ini = ln.split(); S[i] = [(int(a), int(b), float(c)) for a, b, c, _ in (x.split(':') for x in ini.split(';'))]
    return S
S2, S3 = starts(os.path.join(A, 'section-772', 'starts772.txt')), starts(os.path.join(A, 'section-773', 'starts773.txt'))
def nearest(j, l):
    best = (np.inf, None)
    for m in range(-6, 7):
        nr = max(0, round(j - 0.5)); d = np.hypot(j - nr - 0.5, l - m)
        if d < best[0]: best = (d, (nr, m))
    return best
ARMS = [('772 G', 'section-772/out/runs/G_{}.txt', S2, sorted(S2), 3000, 0.001),
        ('772 W', 'section-772/out/runs/W_{}.txt', S2, sorted(S2)[:12], 10000, 0.0003),
        ('773 G', 'section-773/out/runs/G_{}.txt', S3, [f"{k:02d}" for k in range(1, 41)], 12000, 0.001),
        ('773 W', 'section-773/out/runs/W_{}.txt', S2, ['05', '08', '14', '15', '36'], 20000, 0.0003),
        ('774 G', 'section-774/out/runs/G_{}.txt', S3, [f"{k:02d}" for k in range(1, 41)], 20000, 0.0003),
        ('775 B05', 'section-775/out/runs/B05_{}.txt', S3, [f"{k:02d}" for k in range(1, 21)], 20000, 0.0003),
        ('775 B80', 'section-775/out/runs/B80_{}.txt', S3, [f"{k:02d}" for k in range(1, 21)], 20000, 0.0003)]
GRID = [0.01, 0.03, 0.1, 0.3, 1.0]
rows = []
for arm, pat, S, ids, per, GT in ARMS:
    tau = 1 / GT
    for i in ids:
        q = np.loadtxt(os.path.join(A, pat.format(i)), comments='#', ndmin=2); assert len(q) == per, (arm, i, len(q))
        e = q[-20:, 2:4].mean(0); d, win = nearest(*e); wm = q[-1, 5]
        if not (wm >= 0.95 and d < 0.04): continue
        dom = [(int(a), int(b)) for a, b in q[:, 6:8]]
        J0 = sum(w * (a + 0.5) for a, b, w in S[i]) / sum(w for *_, w in S[i]); L0 = sum(w * b for a, b, w in S[i]) / sum(w for *_, w in S[i])
        dw = lambda j, l: np.hypot(j - win[0] - 0.5, l - win[1])
        r = dict(arm=arm, i=i, win=win, tau=tau, per=per, endok=dom[-1] == win, dw0=dw(J0, L0))
        if r['endok']:
            kd = per
            while kd > 1 and dom[kd - 2] == win: kd -= 1        # kd = 1-based first period of the final run of the winner
            r['kd'] = kd; r['sw'] = sum(dom[k] != dom[k - 1] for k in range(1, kd)); r['wkd'] = q[kd - 1, 5]; r['dwkd'] = dw(*q[kd - 1, 2:4])
            k95 = next((k for k in range(kd - 1, per) if q[k, 5] >= 0.95), None); r['t95'] = None if k95 is None else k95 + 1 - kd
        for g in GRID:
            k = max(1, int(round(g * tau)))
            if k <= per: r[('near', g)] = nearest(*q[k - 1, 2:4])[1] == win; r[('dom', g)] = dom[k - 1] == win
        rows.append(r)
n = len(rows); E = [r for r in rows if r['endok']]
pr(f"Konce jednomodowe na EBK: {n}; mod dominujacy na koncu = zwyciezca: {len(E)}; inne: " + (', '.join(f"{r['arm']} {r['i']}" for r in rows if not r['endok']) or 'brak'))
pr("=== szczegoly ===")
pr("  ramie    id  zwyc.     k_d     k_d/tau  zmian  w(k_d)  d_w(k_d) d_w(0)  k_d->w95")
for r in E: pr(f"  {r['arm']:8s} {r['i']}  {str(r['win']):9s} {r['kd']:6d}  {r['kd'] / r['tau']:7.3f}  {r['sw']:5d}  {r['wkd']:.3f}   {r['dwkd']:.3f}    {r['dw0']:.3f}   {r['t95']}")
pr("=== krzywe trafien (wszystkie ramiona; mianownik = biegi, w ktorych chwila miesci sie w biegu) ===")
for g in GRID:
    R = [r for r in rows if ('near', g) in r]
    pr(f"  k = {g:5.2f} tau: P_near {sum(r[('near', g)] for r in R)}/{len(R)} = {100 * sum(r[('near', g)] for r in R) / len(R):.1f} %, P_dom {sum(r[('dom', g)] for r in R)}/{len(R)} = {100 * sum(r[('dom', g)] for r in R) / len(R):.1f} %")
R = [r for r in rows if ('near', 0.1) in r]; hn = sum(r[('near', 0.1)] for r in R); hd = sum(r[('dom', 0.1)] for r in R)
pr(f"H777a (k = 0,1 tau: P_near lub P_dom >= 80 %): P_near {100 * hn / n:.1f} %, P_dom {100 * hd / n:.1f} % (na {n}) -> {'TAK' if max(hn, hd) >= 0.8 * n else 'NIE'}")
a = sum(r['dwkd'] < 0.15 for r in E); b = sum(r['dwkd'] < r['dw0'] for r in E)
pr(f"H777b (d_w(k_d) < 0,15 w >= 80 % i d_w(k_d) < d_w(0) w >= 80 %): {a}/{len(E)} = {100 * a / len(E):.1f} %, {b}/{len(E)} = {100 * b / len(E):.1f} % -> {'TAK' if a >= 0.8 * len(E) and b >= 0.8 * len(E) else 'NIE'}")
pr("H777c (raport) na ramie: n, mediana k_d (okresy), mediana k_d/tau [min-max], mediana zmian przed k_d, mediana w(k_d), mediana d_w(k_d), mediana k_d->w95")
for arm, *_ in ARMS:
    R = [r for r in E if r['arm'] == arm]
    if not R: continue
    kt = [r['kd'] / r['tau'] for r in R]; t9 = [r['t95'] for r in R if r['t95'] is not None]
    pr(f"  {arm:8s} n {len(R):2d}: {np.median([r['kd'] for r in R]):7.0f}  {np.median(kt):.3f} [{min(kt):.3f}-{max(kt):.3f}]  {np.median([r['sw'] for r in R]):4.0f}  {np.median([r['wkd'] for r in R]):.3f}  {np.median([r['dwkd'] for r in R]):.3f}  {np.median(t9) if t9 else float('nan'):.0f}")
kt = [r['kd'] / r['tau'] for r in E]
pr(f"  razem: mediana k_d/tau {np.median(kt):.3f}, kwartyle {np.percentile(kt, 25):.3f}-{np.percentile(kt, 75):.3f}; k_d/tau < 0,1 w {sum(x < 0.1 for x in kt)}/{len(E)}; "
   f"zmian przed k_d: mediana {np.median([r['sw'] for r in E]):.0f}, 0 zmian w {sum(r['sw'] == 0 for r in E)}; w(k_d) < 0,5 w {sum(r['wkd'] < 0.5 for r in E)}/{len(E)}")
P = {(r['arm'], r['i']): r for r in E}
pairs = [(P[('773 G', i)], P[('774 G', i)]) for i in [f"{k:02d}" for k in range(1, 41)] if ('773 G', i) in P and ('774 G', i) in P and P[('773 G', i)]['win'] == P[('774 G', i)]['win']]
if pairs:
    x = [a['kd'] / a['tau'] for a, _ in pairs]; y = [b['kd'] / b['tau'] for _, b in pairs]
    pr(f"H777d (raport): pary 773 G / 774 G z tym samym zwyciezca: {len(pairs)}; mediana k_d/tau {np.median(x):.3f} / {np.median(y):.3f}; mediana k_d {np.median([a['kd'] for a, _ in pairs]):.0f} / {np.median([b['kd'] for _, b in pairs]):.0f} okresow")

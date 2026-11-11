# Audit 781: classical radiation reaction taken from the field -- verdicts H781a-c, e -> out/ana781.txt
import os, numpy as np
from collections import Counter
D = os.path.dirname(os.path.abspath(__file__)); A = os.path.join(D, '..'); O = os.path.join(D, 'out', 'runs')
pr = lambda *a: print(*a, flush=True)
BETA = 0.2
def nearest(j, l):
    best = (np.inf, None)
    for m in range(-6, 7):
        nr = max(0, round(j - 0.5)); d = np.hypot(j - nr - 0.5, l - m)
        if d < best[0]: best = (d, (nr, m))
    return best
Nm = lambda t: t[0] + 0.5 + np.sqrt(t[1] ** 2 + 2 * BETA)
def end(fn, per):
    q = np.loadtxt(fn, comments='#', ndmin=2); assert len(q) == per, (fn, len(q))
    j, l = q[-20:, 2:4].mean(0); d, tor = nearest(j, l); w = q[-1, 5]; dom = (int(q[-1, 6]), int(q[-1, 7]))
    dfc = np.median(q[-1000:, 8]) if q.shape[1] > 8 else np.nan
    return dict(j=j, l=l, d=d, tor=tor, w=w, dom=dom, one=w >= 0.95 and d < 0.04, dfc=dfc, q=q)
c0 = open(os.path.join(A, 'section-775', 'out', 'runs', 'C_ctrl.txt')).read()
ok = all(open(os.path.join(O, f)).read() == c0 for f in ('K_norad.txt', 'K_rad0.txt'))
pr(f"H781e (kontrola kodu: bez RAD i RAD=0 = 775 C_ctrl bajt w bajt): {'TAK' if ok else 'NIE'}")
pr("=== P: czyste mody, 10000 okresow ===")
pr("  RAD  start     koniec: torus   odl.    w_max  dominujacy  N(start) N(koniec)  deficyt (mediana ost. 1000)  1-mod")
P = {}
for r in (1, 10):
    for md in ((0, 0), (0, 1), (0, 2), (1, 1), (1, -2), (2, 0)):
        e = end(os.path.join(O, f"P{r}_{md[0]}_{md[1]}.txt"), 10000); P[(r, md)] = e
        pr(f"  {r:3d}  {str(md):9s}         {str(e['tor']):9s} {e['d']:.4f}  {e['w']:.3f}  {str(e['dom']):9s}   {Nm(md):.3f}    {Nm(e['tor']):.3f}     {e['dfc']:.3e}                  {'T' if e['one'] else 'N'}")
okA = all(P[(r, (0, 0))]['one'] and P[(r, (0, 0))]['tor'] == (0, 0) for r in (1, 10))
pr(f"H781a (czysty (0,0) zostaje jednomodowo na (0,0) przy RAD 1 i 10): {'TAK' if okA else 'NIE'}")
EX = [(0, 1), (0, 2), (1, 1), (1, -2), (2, 0)]
nb = sum(P[(10, md)]['one'] and Nm(P[(10, md)]['tor']) < Nm(md) - 1e-9 for md in EX)
pr(f"H781b (RAD 10: z 5 wzbudzonych >= 4 jednomodowo na EBK o nizszej energii): {nb}/5 -> {'TAK' if nb >= 4 else 'NIE'}")
nb1 = sum(P[(1, md)]['one'] and Nm(P[(1, md)]['tor']) < Nm(md) - 1e-9 for md in EX)
pr(f"  raport RAD 1: {nb1}/5; zostaje na wlasnym torusie: RAD 1 {sum(P[(1, md)]['one'] and P[(1, md)]['tor'] == md for md in EX)}/5, RAD 10 {sum(P[(10, md)]['one'] and P[(10, md)]['tor'] == md for md in EX)}/5")
pr("  kaskada (mod dominujacy przy zmianach, okres):")
for r in (1, 10):
    for md in EX:
        q = P[(r, md)]['q']; seq = []; prev = None
        for k in range(len(q)):
            dm = (int(q[k, 6]), int(q[k, 7]))
            if dm != prev: seq.append(f"{dm}@{k + 1}"); prev = dm
        pr(f"    RAD {r} {md}: " + ' -> '.join(seq[:12]) + (' ...' if len(seq) > 12 else ''))
pr("=== R1, R10: starty 773/01-20, 20000 okresow, wobec 774 G (bez RAD) ===")
ids = [f"{k:02d}" for k in range(1, 21)]
G = {i: end(os.path.join(A, 'section-774', 'out', 'runs', f"G_{i}.txt"), 20000) for i in ids}
R = {r: {i: end(os.path.join(O, f"R{r}_{i}.txt"), 20000) for i in ids} for r in (1, 10)}
S0 = sum(G[i]['one'] for i in ids); S = {r: sum(R[r][i]['one'] for i in ids) for r in R}
pr("  id  774G          R1                    R10            (* = jednomodowy na EBK; odl.)")
for i in ids:
    f = lambda e: f"{str(e['tor'])}{'*' if e['one'] else ' '} {e['d']:.4f}"
    pr(f"  {i}  {f(G[i])}   {f(R[1][i])}  def {R[1][i]['dfc']:.2e}   {f(R[10][i])}  def {R[10][i]['dfc']:.2e}")
pr(f"  jednomodowe na EBK: 774G {S0}/20, R1 {S[1]}/20, R10 {S[10]}/20; ten sam torus co 774G: R1 {sum(R[1][i]['tor'] == G[i]['tor'] for i in ids)}/20, R10 {sum(R[10][i]['tor'] == G[i]['tor'] for i in ids)}/20")
for r in (1, 10):
    O1 = [R[r][i] for i in ids if R[r][i]['one']]
    pr(f"  RAD {r}: mediana odl. jednomodowych {np.median([e['d'] for e in O1]) if O1 else float('nan'):.4f}; mediana N konca {np.median([Nm(R[r][i]['tor']) for i in ids]):.3f} (774G {np.median([Nm(G[i]['tor']) for i in ids]):.3f}); "
       f"zwyciezcy: " + ', '.join(f"{x}: {c}" for x, c in Counter(R[r][i]['tor'] for i in ids).most_common()))
pr(f"H781c (R1 >= 774G - 3): {S[1]} >= {S0 - 3} -> {'TAK' if S[1] >= S0 - 3 else 'NIE'}")

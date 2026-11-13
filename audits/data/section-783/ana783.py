# Audit 783: radiative widths (Weisskopf-Wigner) instead of uniform Gamma, with / without the field-dipole term -- verdicts H783e, k, a, b, c, f -> out/ana783.txt
import os, numpy as np
from collections import Counter
D = os.path.dirname(os.path.abspath(__file__)); A = os.path.join(D, '..'); O = os.path.join(D, 'out', 'runs')
pr = lambda *a: print(*a, flush=True)
BETA = 0.2; GREF = 0.0003; N0 = 0.5 + np.sqrt(1 + 2 * BETA); Tr0 = 2 * np.pi * N0 ** 3
Nm = lambda t: t[0] + 0.5 + np.sqrt(t[1] ** 2 + 2 * BETA); En = lambda t: -0.5 / Nm(t) ** 2
LN = [(0,0),(0,1),(0,-1),(1,0),(0,2),(0,-2),(1,1),(1,-1),(2,0),(1,2),(1,-2),(2,1),(2,-1)]
def nearest(j, l):
    best = (np.inf, None)
    for m in range(-6, 7):
        nr = max(0, round(j - 0.5)); d = np.hypot(j - nr - 0.5, l - m)
        if d < best[0]: best = (d, (nr, m))
    return best
def W(fn):
    k, w = [], []
    for line in open(fn):
        if line.startswith('#W'): p = line.split(); k.append(int(p[1])); w.append([float(x) for x in p[2:]])
    return np.array(k), np.array(w)
def Gtab(fn):
    return {(int(l.split()[1]), int(l.split()[2])): float(l.split()[4]) for l in open(fn) if l.startswith('#G')}
def settled(q, L=200):
    dm = [(int(a), int(b)) for a, b in q[:, 6:8]]; seq = []; i = 0
    while i < len(dm):
        j = i
        while j < len(dm) and dm[j] == dm[i]: j += 1
        if j - i >= L and (not seq or seq[-1] != dm[i]): seq.append(dm[i])
        i = j
    return seq
def end(fn, per, start=None):
    q = np.loadtxt(fn, comments='#', ndmin=2); assert len(q) == per, (fn, len(q))
    j, l = q[-20:, 2:4].mean(0); d, tor = nearest(j, l); w = q[-1, 5]; dom = (int(q[-1, 6]), int(q[-1, 7]))
    ch = next((i + 1 for i in range(len(q)) if (int(q[i, 6]), int(q[i, 7])) != start), None) if start else None
    return dict(j=j, l=l, d=d, tor=tor, w=w, dom=dom, one=w >= 0.95 and d < 0.04, q=q, seq=settled(q), ch=ch, Pmax=q[:, 4].max())
# H783e
c0 = open(os.path.join(A, 'section-775', 'out', 'runs', 'C_ctrl.txt')).read()
okE = all(open(os.path.join(O, f)).read() == c0 for f in ('E_nogw.txt', 'E_gw0.txt'))
pr(f"H783e (kontrola kodu: bez GW i GW=0 = 775 C_ctrl bajt w bajt): {'TAK' if okE else 'NIE'}")
# H783k
k, w = W(os.path.join(O, 'K_01_00.txt')); G10 = Gtab(os.path.join(O, 'K_01_00.txt'))
c = np.polyfit(k * Tr0, np.log(w[:, 0] / w[:, 1]), 1); g01 = G10[(0, 1)] * GREF / Tr0
pr(f"H783k (nachylenie ln(w00/w01) = Gamma_(0,1) = {G10[(0, 1)]:.4f} Gamma_ref): stosunek {c[0] / g01:.5f} -> {'TAK' if abs(c[0] / g01 - 1) <= 0.01 else 'NIE'}")
# P
MODES = [(0, 1), (0, 2), (1, 1), (1, -2), (2, 0), (0, 0)]; EX = MODES[:5]
G30 = Gtab(os.path.join(O, 'PWD_0_1.txt'))
pr("Gamma_n / Gamma_ref przy s = 30: " + ', '.join(f"{m}: {G30[m]:.3f}" for m in LN))
P = {}
for arm in ('WD', 'W', 'D'):
    pr(f"=== P, ramie {arm} (s = 30, 30000 okresow) ===")
    pr("  start     koniec: torus  odl.    w_max  dominujacy  N(start) N(koniec) 1-mod  pierwsza zmiana  x = Gamma_n t   max P      stany osiadle")
    for md in MODES:
        e = end(os.path.join(O, f"P{arm}_{md[0]}_{md[1]}.txt"), 30000, md); P[(arm, md)] = e
        x = G30[md] * GREF * e['ch'] if e['ch'] else None
        pr(f"  {str(md):9s}        {str(e['tor']):9s} {e['d']:.4f}  {e['w']:.3f}  {str(e['dom']):9s}   {Nm(md):.3f}    {Nm(e['tor']):.3f}    {'T' if e['one'] else 'N'}     "
           f"{str(e['ch']):>8s}       {('%.2f' % x) if x is not None else '-':>8s}   {e['Pmax']:.2e}   {' -> '.join(map(str, e['seq']))}")
    nb = sum(P[(arm, md)]['one'] and Nm(P[(arm, md)]['tor']) < Nm(md) - 1e-9 for md in EX)
    pr(f"  {arm}: jednomodowo na nizszym EBK {nb}/5; (0,0) zostaje: {P[(arm, (0, 0))]['one'] and P[(arm, (0, 0))]['tor'] == (0, 0)}")
    for md in EX:
        _, ww = W(os.path.join(O, f"P{arm}_{md[0]}_{md[1]}.txt")); top = sorted(zip(ww[-1], LN), reverse=True)[:3]
        pr(f"    {str(md):8s} koncowe wagi: " + ', '.join(f"{m}: {v:.3f}" for v, m in top))
okA = P[('WD', (0, 0))]['one'] and P[('WD', (0, 0))]['tor'] == (0, 0)
pr(f"H783a (czysty (0,0) w WD jednomodowo na (0,0)): {'TAK' if okA else 'NIE'}")
nb = sum(P[('WD', md)]['one'] and Nm(P[('WD', md)]['tor']) < Nm(md) - 1e-9 for md in EX)
pr(f"H783b (WD: z 5 wzbudzonych >= 4 jednomodowo na EBK o nizszej energii): {nb}/5 -> {'TAK' if nb >= 4 else 'NIE'}")
xs = [(md, G30[md] * GREF * P[('WD', md)]['ch']) for md in EX if P[('WD', md)]['ch']]
okF = len(xs) >= 3 and all(1 <= x <= 30 for _, x in xs)
pr(f"H783f (WD: Gamma_n x czas pierwszej zmiany w [1, 30] dla wszystkich opuszczajacych start, >= 3): " + ', '.join(f"{m}: {x:.2f}" for m, x in xs) + f" -> {'TAK' if okF else 'NIE'}")
# R
pr("=== R (WD, s = 1): starty 773/01-20, 20000 okresow, wobec 774 G, 781 R1, 782 R1 ===")
ids = [f"{i:02d}" for i in range(1, 21)]
G = {i: end(os.path.join(A, 'section-774', 'out', 'runs', f"G_{i}.txt"), 20000) for i in ids}
R = {i: end(os.path.join(O, f"R_{i}.txt"), 20000) for i in ids}
L1 = {i: end(os.path.join(A, 'section-781', 'out', 'runs', f"R1_{i}.txt"), 20000) for i in ids}
D1 = {i: end(os.path.join(A, 'section-782', 'out', 'runs', f"R1_{i}.txt"), 20000) for i in ids}
f = lambda e: f"{str(e['tor'])}{'*' if e['one'] else ' '} {e['d']:.4f} w {e['w']:.3f}"
pr("  id  774G                       783 R                      781 R1    782 R1   (* = jednomodowy na EBK)")
for i in ids: pr(f"  {i}  {f(G[i])}   {f(R[i])}   {str(L1[i]['tor']):9s} {str(D1[i]['tor']):9s}")
S0 = sum(G[i]['one'] for i in ids); S = sum(R[i]['one'] for i in ids)
pr(f"  jednomodowe na EBK: 774G {S0}/20, 783 {S}/20 (781 {sum(L1[i]['one'] for i in ids)}/20, 782 {sum(D1[i]['one'] for i in ids)}/20); ten sam torus co 774G {sum(R[i]['tor'] == G[i]['tor'] for i in ids)}/20")
pr(f"  mediana N konca 783 {np.median([Nm(R[i]['tor']) for i in ids]):.3f} (774G {np.median([Nm(G[i]['tor']) for i in ids]):.3f}); max P {max(R[i]['Pmax'] for i in ids):.2e}; zwyciezcy: "
   + ', '.join(f"{x}: {c}" for x, c in Counter(R[i]['tor'] for i in ids).most_common()))
pr(f"H783c (R >= 17/20): {S} -> {'TAK' if S >= 17 else 'NIE'}")

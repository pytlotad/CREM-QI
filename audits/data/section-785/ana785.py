# Audit 785: exclusion thesis (torus selection by damping vs radiative decay), DIP scan 30 / 100 / 300 at uniform Gamma -- verdicts H785e, a, p, s, w -> out/ana785.txt
import os, numpy as np
from collections import Counter
D = os.path.dirname(os.path.abspath(__file__)); A = os.path.join(D, '..'); O = os.path.join(D, 'out', 'runs')
pr = lambda *a: print(*a, flush=True)
BETA = 0.2
Nm = lambda t: t[0] + 0.5 + np.sqrt(t[1] ** 2 + 2 * BETA)
def nearest(j, l):
    best = (np.inf, None)
    for m in range(-6, 7):
        nr = max(0, round(j - 0.5)); d = np.hypot(j - nr - 0.5, l - m)
        if d < best[0]: best = (d, (nr, m))
    return best
def settled(q, L=200):
    dm = [(int(a), int(b)) for a, b in q[:, 6:8]]; seq = []; i = 0
    while i < len(dm):
        j = i
        while j < len(dm) and dm[j] == dm[i]: j += 1
        if j - i >= L and (not seq or seq[-1] != dm[i]): seq.append(dm[i])
        i = j
    return seq
def end(q):
    j, l = q[-20:, 2:4].mean(0); d, tor = nearest(j, l); w = q[-1, 5]
    return dict(d=d, tor=tor, w=w, one=w >= 0.95 and d < 0.04)
def dwell(q, L=200):
    # excited stays: >= L consecutive records single-mode (w >= 0.95, record distance < 0.04) on the same EBK torus != (0,0); returns list of (torus, start period, length)
    out = []; cur = None; s0 = 0; n = 0
    for k in range(len(q)):
        d, tor = nearest(q[k, 2], q[k, 3]); ok = q[k, 5] >= 0.95 and d < 0.04 and tor != (0, 0)
        key = tor if ok else None
        if key == cur and key is not None: n += 1
        else:
            if cur is not None and n >= L: out.append((cur, s0, n))
            cur, s0, n = key, k + 1, 1
    if cur is not None and n >= L: out.append((cur, s0, n))
    return out
load = lambda fn: np.loadtxt(fn, comments='#', ndmin=2)
# H785e
c0 = open(os.path.join(A, 'section-775', 'out', 'runs', 'C_ctrl.txt')).read()
ref = []
for line in open(os.path.join(A, 'section-782', 'out', 'runs', 'P10_0_1.txt')):
    ref.append(line)
    if line.startswith('#W 300 '): break
e10 = [l for l in open(os.path.join(O, 'E10.txt')) if not l.startswith('# 783 GW')]
okE = open(os.path.join(O, 'E_nodip.txt')).read() == c0 and ''.join(e10) == ''.join(ref)
pr(f"H785e (kontrola kodu: E = 775 C_ctrl, E10 = 782 P10_0_1 do okresu 300 bez linii '# 783 GW', bajt w bajt): {'TAK' if okE else 'NIE'}")
DS = [30, 100, 300]; MD = [(0, 1), (0, 2), (1, 1), (1, -2), (2, 0), (0, 0)]; EX = MD[:5]
P = {}; casc = {}
for dp in DS:
    pr(f"=== P, DIP {dp} (15000 okresow) ===")
    pr("  start     koniec: torus  odl.    w_max  1-mod  nizszy EBK  stany osiadle (|dm|)")
    for m in MD:
        q = load(os.path.join(O, f"P{dp}_{m[0]}_{m[1]}.txt")); assert len(q) == 15000; e = end(q); sq = settled(q); P[(dp, m)] = e
        low = e['one'] and Nm(e['tor']) < Nm(m) - 1e-9; e['low'] = low
        pr(f"  {str(m):9s}         {str(e['tor']):9s} {e['d']:.4f}  {e['w']:.3f}    {'T' if e['one'] else 'N'}      {'T' if low else 'N'}       {' -> '.join(map(str, sq))} {[abs(b[1] - a[1]) for a, b in zip(sq, sq[1:])]}")
    casc[dp] = sum(P[(dp, m)]['low'] for m in EX); pr(f"  Kaskada P(DIP {dp}) = {casc[dp]}/5")
okA = all(P[(dp, (0, 0))]['one'] and P[(dp, (0, 0))]['tor'] == (0, 0) for dp in DS)
pr(f"H785a (czysty (0,0) zostaje przy DIP 30, 100, 300): {'TAK' if okA else 'NIE'}")
okP = casc[30] <= casc[100] <= casc[300] and casc[300] >= 4
pr(f"H785p (Kaskada P niemalejaca i P(300) >= 4): {casc[30]} -> {casc[100]} -> {casc[300]} -> {'TAK' if okP else 'NIE'}")
ids = [f"{i:02d}" for i in range(1, 21)]
def arm(pat, per):
    res = {}
    for i in ids:
        q = load(pat.format(i))[:per]; assert len(q) == per, (pat, i, len(q)); res[i] = dict(dw=dwell(q), **end(q))
    return res
REF = {'782 R1 (DIP 1)': arm(os.path.join(A, 'section-782', 'out', 'runs', 'R1_{}.txt'), 10000),
       '782 R10 (DIP 10)': arm(os.path.join(A, 'section-782', 'out', 'runs', 'R10_{}.txt'), 10000)}
R = {dp: arm(os.path.join(O, f"R{dp}_" + "{}.txt"), 10000) for dp in DS}
Sx = {}
pr("=== R: starty 773/01-20, 10000 okresow ===")
for name, r in list(REF.items()) + [(f"785 DIP {dp}", R[dp]) for dp in DS]:
    s = sum(bool(r[i]['dw']) for i in ids); first = [r[i]['dw'][0] for i in ids if r[i]['dw']]
    pr(f"  {name:17s}: S_exc {s:2d}/20; konce jednomodowe na EBK {sum(r[i]['one'] for i in ids)}/20, na (0,0) {sum(r[i]['one'] and r[i]['tor'] == (0, 0) for i in ids)}/20; "
       f"mediana startu 1. pobytu {np.median([f[1] for f in first]) if first else float('nan'):.0f}, dlugosci {np.median([f[2] for f in first]) if first else float('nan'):.0f} okresow; "
       f"tory 1. pobytu: " + ', '.join(f"{t}: {c}" for t, c in Counter(f[0] for f in first).most_common()))
    if name.startswith('785'): Sx[int(name.split()[-1])] = s
pr("  pierwsze pobyty wzbudzone (torus@start+dlugosc) i koniec:")
for i in ids:
    pr(f"  {i}  " + '   '.join(f"D{dp}: {(lambda d: f'{d[0]}@{d[1]}+{d[2]}' if d else '-')(R[dp][i]['dw'][0] if R[dp][i]['dw'] else None):22s} kon {R[dp][i]['tor']}{'*' if R[dp][i]['one'] else ' '}" for dp in DS))
okS = Sx[30] >= Sx[100] >= Sx[300] and Sx[300] <= 10
pr(f"H785s (S_exc nierosnacy i S_exc(300) <= 10): {Sx[30]} -> {Sx[100]} -> {Sx[300]} -> {'TAK' if okS else 'NIE'}")
win = [dp for dp in DS if casc[dp] >= 4 and Sx[dp] >= 17]
pr(f"H785w (teza wykluczenia: zaden DIP z Kaskada P >= 4 i S_exc >= 17): okno {win if win else 'brak'} -> {'TAK' if not win else 'NIE'}")

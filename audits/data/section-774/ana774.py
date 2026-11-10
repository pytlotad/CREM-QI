# Audit 774: Q scaling, GT 0.0003, starts of 773 -- verdicts H774a-d -> out/ana774.txt (based on ../section-773/ana773.py)
import os, numpy as np
from collections import Counter
D = os.path.dirname(os.path.abspath(__file__)); O = os.path.join(D, 'out', 'runs'); pr = lambda *a: print(*a, flush=True)
def starts(fn):
    S = {}
    for ln in open(os.path.join(D, fn)):
        i, ini = ln.split(); S[i] = [(int(a), int(b), float(c), float(d)) for a, b, c, d in (x.split(':') for x in ini.split(';'))]
    return S
ST = starts('../section-773/starts773.txt')
def nearest(j, l):
    best = (np.inf, None)
    for m in range(-6, 7):
        nr = max(0, round(j - 0.5)); d = np.hypot(j - nr - 0.5, l - m)
        if d < best[0]: best = (d, (nr, m))
    return best
def ev(pref, i, per, GT):
    q = np.loadtxt(os.path.join(O, f"{pref}_{i}.txt"), comments='#', ndmin=2); assert len(q) == per, (pref, i, len(q))
    e = q[-20:, 2:4].mean(0); k = int(0.8 * per); b = q[k - 10:k + 10, 2:4].mean(0); dr = np.hypot(*(e - b))
    res = 500 * np.exp(-2 * GT * per) / q[-20:, 4].mean(); conv = dr < 0.005 and res < 1e-3
    d, t = nearest(*e); wm = q[-1, 5]; win = (int(q[-1, 6]), int(q[-1, 7])); one = wm >= 0.95 and d < 0.04
    sec = win[0] == 1 and abs(win[1]) == 1 and 0.08 <= d <= 0.2 and wm < 0.95
    return dict(j=e[0], l=e[1], dr=dr, res=res, conv=conv, d=d, tor=t, wm=wm, win=win, one=one, ok=conv and one, sec=sec)
def block(pref, S, ids, per, GT, nm):
    pr(f"=== {nm} ===");  pr("  id  J_g      L_g      odl.    torus     w_max   zwyc.     dryf    resztka   zb. sukces 1-mod  start (nr,m,w)")
    R = {}
    for i in ids:
        r = ev(pref, i, per, GT); R[i] = r; s = ' '.join(f"({a},{b},{c:.2f})" for a, b, c, _ in S[i])
        pr(f"  {i}  {r['j']:.5f}  {r['l']:+.5f}  {r['d']:.4f}  {str(r['tor']):9s} {r['wm']:.4f}  {str(r['win']):9s} {r['dr']:.4f}  {r['res']:.1e}  {'T' if r['conv'] else 'N'}   {'T' if r['ok'] else 'N'}      {'T' if r['one'] else 'N'}      {s}")
    nc = sum(r['conv'] for r in R.values()); ns = sum(r['ok'] for r in R.values()); no = sum(r['one'] for r in R.values())
    pr(f"  zbiezne {nc}/{len(ids)}, uciete {len(ids) - nc}, sukces {ns}/{nc} = {100 * ns / max(nc, 1):.1f} %, jednomodowe na EBK (wszystkie) {no}/{len(ids)}")
    return R, nc, ns, no
G, ncg, nsg, _ = block('G', ST, [f"{i:02d}" for i in range(1, 41)], 20000, 0.0003, '(G) GT 0,0003, CAUST -pi/2, 20000 okresow, starty 01-40')
pr(f"H774a (sukces >= 36/40, uciete = porazka): {'TAK' if nsg >= 36 else 'NIE'} ({nsg}/40; zbiezne {ncg}/40)")
K, nck, nsk, nok = block('K', ST, [f"{i:02d}" for i in range(41, 61)], 10000, 0.0003, '(K) GT 0,0003, CAUST 0, 10000 okresow, starty 41-60')
pr(f"H774b (kontrola: jednomodowe na EBK <= 2/20): {'TAK' if nok <= 2 else 'NIE'} ({nok}/20);  w_max kontroli min {min(r['wm'] for r in K.values()):.3f}, mediana {np.median([r['wm'] for r in K.values()]):.3f}, max {max(r['wm'] for r in K.values()):.3f}")
O3 = os.path.join(D, '..', 'section-773', 'out', 'runs'); Ob = O
O = O3; G3 = {f"{i:02d}": ev('G', f"{i:02d}", 12000, 0.001) for i in range(1, 41)}; O = Ob
fs = sum(1 for i in G if not G3[i]['ok'] and G[i]['ok']); sf = sum(1 for i in G if G3[i]['ok'] and not G[i]['ok']); same = sum(1 for i in G if G3[i]['ok'] and G[i]['ok'] and G3[i]['win'] == G[i]['win'])
pr(f"H774c (raport): 773 porazka -> 774 sukces {fs}; 773 sukces -> 774 porazka {sf}; oba sukces {sum(1 for i in G if G3[i]['ok'] and G[i]['ok'])} (ten sam zwyciezca {same})")
pr("  porazki 774: " + '; '.join(f"{i} {G[i]['win']} odl. {G[i]['d']:.3f} w {G[i]['wm']:.3f} zb. {'T' if G[i]['conv'] else 'N'} (773: {G3[i]['win']} {G3[i]['d']:.3f})" for i in G if not G[i]['ok']))
pr("  zwyciezcy-sukcesy: " + ', '.join(f"{k}: {v}" for k, v in Counter(r['win'] for r in G.values() if r['ok']).most_common()))
pr("H774d (raport): przesuniecia (0,+-2) i (1,+-2): 774 " + ', '.join(f"{i} {G[i]['win']} {G[i]['d']:.4f}" for i in G if G[i]['ok'] and abs(G[i]['win'][1]) == 2 and G[i]['win'][0] <= 1) + ";  773 " + ', '.join(f"{i} {G3[i]['win']} {G3[i]['d']:.4f}" for i in G3 if G3[i]['ok'] and abs(G3[i]['win'][1]) == 2 and G3[i]['win'][0] <= 1))

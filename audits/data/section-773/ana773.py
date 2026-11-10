# Audit 773: decisive repeat with new starts, long runs -- verdicts H773a-d -> out/ana773.txt
import os, numpy as np
from collections import Counter
D = os.path.dirname(os.path.abspath(__file__)); O = os.path.join(D, 'out', 'runs'); pr = lambda *a: print(*a, flush=True)
def starts(fn):
    S = {}
    for ln in open(os.path.join(D, fn)):
        i, ini = ln.split(); S[i] = [(int(a), int(b), float(c), float(d)) for a, b, c, d in (x.split(':') for x in ini.split(';'))]
    return S
ST = starts('starts773.txt'); SW = starts('startsW.txt')
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
G, ncg, nsg, _ = block('G', ST, [f"{i:02d}" for i in range(1, 41)], 12000, 0.001, '(G) GT 0,001, CAUST -pi/2, 12000 okresow, starty 01-40')
pr(f"H773a (uciete <= 4/40 i sukces >= 80 % zbieznych): {'TAK' if (40 - ncg) <= 4 and ncg > 0 and nsg >= 0.8 * ncg else 'NIE'}")
K, nck, nsk, nok = block('K', ST, [f"{i:02d}" for i in range(41, 61)], 6000, 0.001, '(K) GT 0,001, CAUST 0, 6000 okresow, starty 41-60')
pr(f"H773b (kontrola: jednomodowe na EBK <= 2/20): {'TAK' if nok <= 2 else 'NIE'} ({nok}/20);  w_max kontroli min {min(r['wm'] for r in K.values()):.3f}, mediana {np.median([r['wm'] for r in K.values()]):.3f}, max {max(r['wm'] for r in K.values()):.3f}")
pr("H773c (raport): zwyciezcy-sukcesy (G): " + ', '.join(f"{k}: {v}" for k, v in Counter(r['win'] for r in G.values() if r['ok']).most_common()))
ins = sum(1 for i, r in G.items() if r['ok'] and r['win'] in [(a, b) for a, b, _, _ in ST[i]]); mx = sum(1 for i, r in G.items() if r['ok'] and r['win'] == max(ST[i], key=lambda x: x[2])[:2])
pr(f"  zwyciezca wsrod modow startu {ins}/{nsg}; = mod o najwiekszej wadze startowej {mx}/{nsg};  drugi atraktor (wszystkie konce): {sum(r['sec'] for r in G.values())}/40 (zbiezne {sum(r['sec'] and r['conv'] for r in G.values())})")
pr("  odleglosci sukcesow: " + ', '.join(f"{r['win']} {r['d']:.4f}" for r in G.values() if r['ok']))
W, _, _, _ = block('W', SW, sorted(SW), 20000, 0.0003, '(W) GT 0,0003, 20000 okresow, starty 772 drugiego atraktora')
pr("H773d (raport): " + ', '.join(f"{i}: odl. {r['d']:.4f}, w_max {r['wm']:.4f}, {r['win']}, zb. {'T' if r['conv'] else 'N'}" for i, r in W.items()))

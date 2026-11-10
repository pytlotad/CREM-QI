# Audit 772: confirmation of 771 POST HOC (random starts) -- verdicts H772a-e -> out/ana772.txt
import os, filecmp, numpy as np
D = os.path.dirname(os.path.abspath(__file__)); O = os.path.join(D, 'out', 'runs'); pr = lambda *a: print(*a, flush=True)
ST = {}
for ln in open(os.path.join(D, 'starts772.txt')):
    i, ini = ln.split(); ST[i] = [(int(a), int(b), float(c), float(d)) for a, b, c, d in (x.split(':') for x in ini.split(';'))]
same = filecmp.cmp(os.path.join(O, 'C_ctrl.txt'), os.path.join(D, '..', 'section-770', 'out', 'runs', 'G_m2_A_P500_f0.3.txt'), shallow=False)
pr(f"H772a (kontrola kodu = 770 G_m2_A_P500_f0.3 bajt w bajt): {'TAK' if same else 'NIE'}")
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
    d, t = nearest(*e); wm = q[-1, 5]; win = (int(q[-1, 6]), int(q[-1, 7]))
    return dict(j=e[0], l=e[1], dr=dr, res=res, conv=conv, d=d, tor=t, wm=wm, win=win, ok=conv and wm >= 0.95 and d < 0.04)
def block(pref, ids, per, GT, nm):
    pr(f"=== {nm} ===");  pr("  id  J_g      L_g      odl.    torus     w_max   zwyc.     dryf    resztka   zbiezny sukces  start (nr,m,w)")
    R = {}
    for i in ids:
        r = ev(pref, i, per, GT); R[i] = r
        s = ' '.join(f"({a},{b},{c:.2f})" for a, b, c, _ in ST[i])
        pr(f"  {i}  {r['j']:.5f}  {r['l']:+.5f}  {r['d']:.4f}  {str(r['tor']):9s} {r['wm']:.4f}  {str(r['win']):9s} {r['dr']:.4f}  {r['res']:.1e}  {'T' if r['conv'] else 'N'}       {'T' if r['ok'] else 'N'}       {s}")
    nc = sum(r['conv'] for r in R.values()); ns = sum(r['ok'] for r in R.values())
    pr(f"  zbiezne {nc}/{len(ids)}, uciete {len(ids) - nc}, sukces {ns}/{nc} = {100 * ns / max(nc, 1):.1f} %")
    return R, nc, ns
ids = sorted(ST)
G, ncg, nsg = block('G', ids, 3000, 0.001, '(G) GT 0,001, CAUST -pi/2, 40 startow')
pr(f"H772b (uciete <= 8/40 i sukces >= 80 % zbieznych): {'TAK' if (40 - ncg) <= 8 and nsg >= 0.8 * ncg and ncg > 0 else 'NIE'}")
K, nck, nsk = block('K', ids[:20], 3000, 0.001, '(K) GT 0,001, CAUST 0, starty 1-20')
pr(f"H772c (kontrola: sukces <= 10 % zbieznych): {'TAK' if nck > 0 and nsk <= 0.1 * nck else 'NIE'};  w_max kontroli: min {min(r['wm'] for r in K.values()):.3f}, mediana {np.median([r['wm'] for r in K.values()]):.3f}, max {max(r['wm'] for r in K.values()):.3f}")
pr("H772d (raport): zwyciezcy w (G) (sukcesy): ")
from collections import Counter
cw = Counter(r['win'] for r in G.values() if r['ok']); pr("  " + ', '.join(f"{k}: {v}" for k, v in cw.most_common()))
ins = sum(1 for i, r in G.items() if r['ok'] and r['win'] in [(a, b) for a, b, _, _ in ST[i]]); mx = sum(1 for i, r in G.items() if r['ok'] and r['win'] == max(ST[i], key=lambda x: x[2])[:2])
nsu = max(nsg, 1); pr(f"  zwyciezca wsrod modow startu: {ins}/{nsg};  zwyciezca = mod o najwiekszej wadze startowej: {mx}/{nsg}")
pr("  odleglosc od EBK (sukcesy): " + ', '.join(f"{r['win']} {r['d']:.4f}" for r in G.values() if r['ok']))
W, ncw, nsw = block('W', ids[:12], 10000, 0.0003, '(W) GT 0,0003, CAUST -pi/2, starty 1-12')
pr("H772e (raport): (0,2) w (W): " + (', '.join(f"{i} {r['d']:.4f}" for i, r in W.items() if r['win'] == (0, 2)) or 'brak') + ";  (0,2) w (G): " + (', '.join(f"{i} {r['d']:.4f}" for i, r in G.items() if r['win'] == (0, 2)) or 'brak'))

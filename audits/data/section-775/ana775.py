# Audit 775: second clock (precession) -- beta scan 0.05 / 0.2 (= 774 G_01..20) / 0.8; verdicts H775a-e -> out/ana775.txt
import os, numpy as np
from collections import Counter
D = os.path.dirname(os.path.abspath(__file__)); O = os.path.join(D, 'out', 'runs'); O4 = os.path.join(D, '..', 'section-774', 'out', 'runs')
pr = lambda *a: print(*a, flush=True)
ST = {}
for ln in open(os.path.join(D, '..', 'section-773', 'starts773.txt')):
    i, ini = ln.split()
    if int(i) <= 20: ST[i] = [(int(a), int(b), float(c)) for a, b, c, _ in (x.split(':') for x in ini.split(';'))]
ids = sorted(ST); PER, GT = 20000, 0.0003
# H775e: code control
c = open(os.path.join(O, 'C_ctrl.txt')).read().splitlines(); g = open(os.path.join(O4, 'G_01.txt')).read().splitlines()
pr(f"H775e (kontrola kodu: C = pierwsze {len(c)} wierszy 774 G_01 bajt w bajt): {'TAK' if len(c) == 302 and c == g[:len(c)] else 'NIE'}")
def nearest(j, l):
    best = (np.inf, None)
    for m in range(-6, 7):
        nr = max(0, round(j - 0.5)); d = np.hypot(j - nr - 0.5, l - m)
        if d < best[0]: best = (d, (nr, m))
    return best
def ev(fn, beta):
    q = np.loadtxt(fn, comments='#', ndmin=2); assert len(q) == PER, (fn, len(q))
    e = q[-20:, 2:4].mean(0); k = int(0.8 * PER); b = q[k - 10:k + 10, 2:4].mean(0); dr = np.hypot(*(e - b))
    res = 500 * np.exp(-2 * GT * PER) / q[-20:, 4].mean(); conv = dr < 0.005 and res < 1e-3
    d, (nr, m) = nearest(*e); wm = q[-1, 5]; one = wm >= 0.95 and d < 0.04
    N = e[0] + np.sqrt(e[1] ** 2 + 2 * beta); Ne = nr + 0.5 + np.sqrt(m * m + 2 * beta); Lam = np.sqrt(m * m + 2 * beta)
    return dict(j=e[0], l=e[1], dr=dr, conv=conv, d=d, tor=(nr, m), wm=wm, one=one, ok=conv and one,
                dN=N - Ne, dL=e[1] - m, Nh=abs(N - Ne) < 0.01, Lh=abs(e[1] - m) < 0.01, prec=1 - abs(m) / Lam)
ARMS = [(0.05, lambda i: os.path.join(O, f"B05_{i}.txt")), (0.2, lambda i: os.path.join(O4, f"G_{i}.txt")), (0.8, lambda i: os.path.join(O, f"B80_{i}.txt"))]
R = {}
for beta, fn in ARMS:
    pr(f"=== beta {beta} {'(= 774 G, nie liczone ponownie)' if beta == 0.2 else ''} ===")
    pr("  id  J_g      L_g      torus     odl.    dN       dL       N-tr L-tr w_max   zb. sukces 1-mod  precesja")
    R[beta] = {}
    for i in ids:
        r = ev(fn(i), beta); R[beta][i] = r
        pr(f"  {i}  {r['j']:.5f}  {r['l']:+.5f}  {str(r['tor']):9s} {r['d']:.4f}  {r['dN']:+.4f}  {r['dL']:+.4f}  {'T' if r['Nh'] else 'N'}    {'T' if r['Lh'] else 'N'}    {r['wm']:.4f}  {'T' if r['conv'] else 'N'}   {'T' if r['ok'] else 'N'}      {'T' if r['one'] else 'N'}      {r['prec']:.4f}")
    V = R[beta].values(); n = len(ids)
    pr(f"  N-trafienia {sum(r['Nh'] for r in V)}/{n}, L-trafienia {sum(r['Lh'] for r in V)}/{n}, torus < 0,01 {sum(r['d'] < 0.01 for r in V)}/{n}, "
       f"sukces {sum(r['ok'] for r in V)}/{n}, uciete {sum(not r['conv'] for r in V)}, jednomodowe na EBK {sum(r['one'] for r in V)}/{n}")
    pr("  zwyciezcy (torus odniesienia): " + ', '.join(f"{k}: {v}" for k, v in Counter(r['tor'] for r in V).most_common()))
L = {b: sum(r['Lh'] for r in R[b].values()) for b in R}; Nn = {b: sum(r['Nh'] for r in R[b].values()) for b in R}
pr(f"H775a (L niemalejace w beta i L(0,8) - L(0,05) >= 5): L = {L[0.05]}, {L[0.2]}, {L[0.8]} -> {'TAK' if L[0.05] <= L[0.2] <= L[0.8] and L[0.8] - L[0.05] >= 5 else 'NIE'}")
pr(f"H775b (N przy 0,05 i 0,8 w +-4 od N(0,2) = {Nn[0.2]}): N = {Nn[0.05]}, {Nn[0.2]}, {Nn[0.8]} -> {'TAK' if abs(Nn[0.05] - Nn[0.2]) <= 4 and abs(Nn[0.8] - Nn[0.2]) <= 4 else 'NIE'}")
pr("H775c (raport): parowanie startow (torus odniesienia przy beta 0,05 / 0,2 / 0,8):")
for i in ids: pr(f"  {i}  " + '  '.join(f"{R[b][i]['tor']}{'*' if R[b][i]['ok'] else ' '}" for b in (0.05, 0.2, 0.8)) + "   (* sukces)")
pr("H775d (raport): L-trafienia wzgledem szybkosci precesji zwyciezcy (m != 0, wszystkie beta):")
rows = sorted((round(r['prec'], 4), b, abs(r['tor'][1]), r['Lh'], abs(r['dL'])) for b in R for r in R[b].values() if r['tor'][1] != 0)
for p in sorted(set(x[0] for x in rows)):
    s = [x for x in rows if x[0] == p]
    pr(f"  precesja {p:.4f} ({', '.join(sorted(set(f'beta {x[1]} |m| {x[2]}' for x in s)))}): L-trafienia {sum(x[3] for x in s)}/{len(s)}, mediana |dL| {np.median([x[4] for x in s]):.4f}")
z = [(b, r) for b in R for r in R[b].values() if r['tor'][1] == 0]
pr(f"  m = 0: L-trafienia {sum(r['Lh'] for _, r in z)}/{len(z)} (" + ', '.join(f"beta {b}: |dL| {abs(r['dL']):.4f}" for b, r in z) + ")")

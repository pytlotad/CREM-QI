# Audit 784: field norm fixed (P = P_fix) with the field-dipole term, no width -- verdicts H784e, n, a, b, t, f, c -> out/ana784.txt
import os, numpy as np
from collections import Counter
D = os.path.dirname(os.path.abspath(__file__)); A = os.path.join(D, '..'); O = os.path.join(D, 'out', 'runs')
pr = lambda *a: print(*a, flush=True)
BETA = 0.2; GREF = 0.0003
Nm = lambda t: t[0] + 0.5 + np.sqrt(t[1] ** 2 + 2 * BETA); En = lambda t: -0.5 / Nm(t) ** 2
# Gamma_n / Gamma_ref at s = 30 (../section-783/out/post783.txt, #G); per period = x GREF
GN = {(0, 1): 30.0, (0, 2): 10.073, (1, 1): 4.949, (1, -2): 3.729, (2, 0): 1.625, (0, 0): 0.0}
GP = {k: v * GREF for k, v in GN.items()}
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
    rp = np.median(q[:(ch - 1 if ch and ch > 1 else len(q)), 10]) if q.shape[1] > 10 else np.nan
    return dict(j=j, l=l, d=d, tor=tor, w=w, dom=dom, one=w >= 0.95 and d < 0.04, q=q, seq=settled(q), ch=ch, rp=rp)
# H784e
c0 = open(os.path.join(A, 'section-775', 'out', 'runs', 'C_ctrl.txt')).read()
ref = []
for line in open(os.path.join(A, 'section-783', 'out', 'runs', 'PD_0_1.txt')):
    ref.append(line)
    if line.startswith('#W 300 '): break
okE = open(os.path.join(O, 'E_nonorm.txt')).read() == c0 and open(os.path.join(O, 'E_D.txt')).read() == ''.join(ref)
pr(f"H784e (kontrola kodu: E_nonorm = 775 C_ctrl, E_D = 783 PD_0_1 do okresu 300, bajt w bajt): {'TAK' if okE else 'NIE'}")
PF = ['1e4', '1e6', '1e8']; MD = [(0, 1), (0, 2), (1, 1), (1, -2), (2, 0), (0, 0)]; EX = MD[:5]
N = {(p, m): end(os.path.join(O, f"N{p}_{m[0]}_{m[1]}.txt"), 30000, m) for p in PF for m in MD}
dev = max(np.max(np.abs(N[(p, m)]['q'][:, 4] / float(p) - 1)) for p in PF for m in MD)
pr(f"H784n (|P / P_fix - 1| < 1e-9 we wszystkich rekordach N): max {dev:.2e} -> {'TAK' if dev < 1e-9 else 'NIE'}")
for p in PF:
    pr(f"=== N, P_fix {p} (30000 okresow) ===")
    pr("  start     koniec: torus  odl.    w_max  dominujacy  N(start) N(koniec) 1-mod  1. zmiana  R_p/okres   R_p/Gamma_n  zejscie  stany osiadle (|dm|)")
    for m in MD:
        e = N[(p, m)]; nxt = next((s for s in e['seq'] if s != m), None); e['down'] = nxt is not None and En(nxt) < En(m)
        rat = e['rp'] / GP[m] if GP[m] > 0 else np.inf; e['rat'] = rat
        dms = [abs(b[1] - a[1]) for a, b in zip(e['seq'], e['seq'][1:])]
        pr(f"  {str(m):9s}         {str(e['tor']):9s} {e['d']:.4f}  {e['w']:.3f}  {str(e['dom']):9s}   {Nm(m):.3f}    {Nm(e['tor']):.3f}    {'T' if e['one'] else 'N'}   {str(e['ch']):>8s}   "
           f"{e['rp']:.3e}   {rat:8.3f}     {'T' if e['down'] else 'N'}     {' -> '.join(map(str, e['seq']))} {dms}")
okA = all(N[(p, (0, 0))]['one'] and N[(p, (0, 0))]['tor'] == (0, 0) for p in PF)
pr(f"H784a (czysty (0,0) jednomodowo na (0,0) przy 1e4, 1e6, 1e8): {'TAK' if okA else 'NIE'}")
nb = sum(N[('1e8', m)]['one'] and Nm(N[('1e8', m)]['tor']) < Nm(m) - 1e-9 for m in EX)
pr(f"H784b (P_fix 1e8: >= 4/5 jednomodowo na nizszym EBK): {nb}/5 -> {'TAK' if nb >= 4 else 'NIE'}")
for p in PF[:2]: pr(f"  raport P_fix {p}: {sum(N[(p, m)]['one'] and Nm(N[(p, m)]['tor']) < Nm(m) - 1e-9 for m in EX)}/5")
hit = sum((N[(p, m)]['rat'] < 1) == N[(p, m)]['down'] for p in PF for m in EX)
pr(f"H784t (zejscie <=> R_p/Gamma_n < 1; 15 biegow): trafnosc {hit}/15 -> {'TAK' if hit >= 13 else 'NIE'}")
for p in PF: pr(f"  {p}: " + ', '.join(f"{m} {N[(p, m)]['rat']:.2f}/{'Z' if N[(p, m)]['down'] else '-'}" for m in EX))
k, w = W(os.path.join(O, 'N1e8_0_1.txt'))
i1 = next((i for i in range(len(k)) if w[i, 0] >= 0.01), None); i2 = next((i for i in range(len(k)) if w[i, 0] >= 0.5), None)
if i1 is not None and i2 is not None and i2 > i1 + 1:
    y = np.log(w[i1:i2 + 1, 0] / w[i1:i2 + 1, 1]); sl = np.polyfit(k[i1:i2 + 1], y, 1)[0]; r = sl / GP[(0, 1)]
    pr(f"H784f (P_fix 1e8, (0,1): nachylenie ln(w00/w01) miedzy w00 0,01 (okres {k[i1]}) a 0,5 (okres {k[i2]}) = {sl:.5f} / okres; / Gamma_(0,1) = {r:.3f}): {'TAK' if abs(r - 1) <= 0.15 else 'NIE'}")
else:
    pr(f"H784f: w00 nie przechodzi 0,01 -> 0,5 (i1 {i1}, i2 {i2}) -> NIE")
pr("=== R6, R8: starty 773/01-20, 15000 okresow, wobec 774 G (koniec 20000) ===")
ids = [f"{i:02d}" for i in range(1, 21)]
G = {i: end(os.path.join(A, 'section-774', 'out', 'runs', f"G_{i}.txt"), 20000) for i in ids}
R = {p: {i: end(os.path.join(O, f"R{p[2]}_{i}.txt"), 15000) for i in ids} for p in ('1e6', '1e8')}
f = lambda e: f"{str(e['tor'])}{'*' if e['one'] else ' '} {e['d']:.4f} w {e['w']:.3f}"
pr("  id  774G                      R6                        R8           (* = jednomodowy na EBK)")
for i in ids: pr(f"  {i}  {f(G[i])}   {f(R['1e6'][i])}   {f(R['1e8'][i])}")
okC = True
for p in ('1e6', '1e8'):
    S = sum(R[p][i]['one'] for i in ids); okC &= S >= 17
    pr(f"  P_fix {p}: jednomodowe na EBK {S}/20; na (0,0) {sum(R[p][i]['one'] and R[p][i]['tor'] == (0, 0) for i in ids)}/20; ten sam torus co 774G {sum(R[p][i]['tor'] == G[i]['tor'] for i in ids)}/20; "
       f"mediana N konca {np.median([Nm(R[p][i]['tor']) for i in ids]):.3f} (774G {np.median([Nm(G[i]['tor']) for i in ids]):.3f}); mediana R_p {np.median([np.median(R[p][i]['q'][:, 10]) for i in ids]):.3e}; "
       f"zwyciezcy: " + ', '.join(f"{x}: {c}" for x, c in Counter(R[p][i]['tor'] for i in ids).most_common()))
    pr(f"    H784c P_fix {p} (>= 17/20): {'TAK' if S >= 17 else 'NIE'}")
pr(f"H784c (R6 i R8 osobno >= 17/20): {'TAK' if okC else 'NIE (patrz osobno)'}")

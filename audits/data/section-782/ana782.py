# Audit 782: Abraham-Lorentz reaction on the field dipole (Crisp-Jaynes) in the guidance prototype -- verdicts H782e, k, a-d -> out/ana782.txt
import os, numpy as np
from collections import Counter
D = os.path.dirname(os.path.abspath(__file__)); A = os.path.join(D, '..'); O = os.path.join(D, 'out', 'runs')
pr = lambda *a: print(*a, flush=True)
BETA = 0.2; GT = 0.0003; N0 = 0.5 + np.sqrt(1 + 2 * BETA); Tr0 = 2 * np.pi * N0 ** 3; Gam = GT / Tr0
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
def end(fn, per):
    q = np.loadtxt(fn, comments='#', ndmin=2); assert len(q) == per, (fn, len(q))
    j, l = q[-20:, 2:4].mean(0); d, tor = nearest(j, l); w = q[-1, 5]; dom = (int(q[-1, 6]), int(q[-1, 7]))
    return dict(j=j, l=l, d=d, tor=tor, w=w, dom=dom, one=w >= 0.95 and d < 0.04, q=q,
                erad=np.median(q[-1000:, 8]) if q.shape[1] >= 10 else np.nan, w5=q[-1, 9] if q.shape[1] >= 10 else np.nan)
def settled(q, L=200):
    dm = [(int(a), int(b)) for a, b in q[:, 6:8]]; seq = []; i = 0
    while i < len(dm):
        j = i
        while j < len(dm) and dm[j] == dm[i]: j += 1
        if j - i >= L and (not seq or seq[-1] != dm[i]): seq.append(dm[i])
        i = j
    return seq
def steps(seq): return [(a, b) for a, b in zip(seq, seq[1:]) if En(b) < En(a) - 1e-12]
# H782e
c0 = open(os.path.join(A, 'section-775', 'out', 'runs', 'C_ctrl.txt')).read()
okE = all(open(os.path.join(O, f)).read() == c0 for f in ('E_nodip.txt', 'E_dip0.txt'))
pr(f"H782e (kontrola kodu: bez DIP i DIP=0 = 775 C_ctrl bajt w bajt): {'TAK' if okE else 'NIE'}")
# H782k
k, w = W(os.path.join(O, 'K_01_00.txt')); w00 = w[:, 0] / (w[:, 0] + w[:, 1]); t = k * Tr0
sel = (w00 > 0.02) & (w00 < 0.98); c = np.polyfit(t[sel], np.log(w00[sel] / (1 - w00[sel])), 1); GF = 10 * Gam
k2, w2 = W(os.path.join(O, 'K_20_00.txt')); dev = np.max(np.abs(w2[:, 0] - w2[0, 0]))
okK = abs(c[0] / GF - 1) <= 0.02 and dev < 1e-6
pr(f"H782k: (0,1)+(0,0) Gamma_fit/Gamma_F = {c[0] / GF:.4f} (punkty {sel.sum()}, s0 z dopasowania {1 / (1 + np.exp(-c[1])):.4f}); (2,0)+(0,0) max |w00 - w00(0)| = {dev:.2e} -> {'TAK' if okK else 'NIE'}")
q = np.loadtxt(os.path.join(O, 'K_01_00.txt'), comments='#', ndmin=2)
pr(f"  K_01_00: suma E_rad {q[:, 8].sum():.6e} wobec spadku <E> (E01 - E00) * (w00 koniec - w00 start) {(En((0,1)) - En((0,0))) * (w00[-1] - 0.01):.6e}")
# P
pr("=== P: czyste mody, 10000 okresow ===")
pr("  DIP  start     koniec: torus   odl.    w_max  dominujacy  N(start) N(koniec)  E_rad/okres (med. ost. 1000)  w(|m|>=5)  1-mod  stany osiadle")
P = {}
MODES = [(0, 0), (0, 1), (0, 2), (1, 1), (1, -2), (2, 0), (1, 0)]
for r in (1, 10):
    for md in MODES:
        e = end(os.path.join(O, f"P{r}_{md[0]}_{md[1]}.txt"), 10000); P[(r, md)] = e; e['seq'] = settled(e['q'])
        pr(f"  {r:3d}  {str(md):9s}         {str(e['tor']):9s} {e['d']:.4f}  {e['w']:.3f}  {str(e['dom']):9s}   {Nm(md):.3f}    {Nm(e['tor']):.3f}     {e['erad']:.3e}               {e['w5']:.2e}   {'T' if e['one'] else 'N'}    {' -> '.join(map(str, e['seq']))}")
okA = all(P[(r, (0, 0))]['one'] and P[(r, (0, 0))]['tor'] == (0, 0) for r in (1, 10))
pr(f"H782a (czysty (0,0) jednomodowo na (0,0) przy DIP 1 i 10): {'TAK' if okA else 'NIE'}")
EX = [(0, 1), (0, 2), (1, 1), (1, -2), (2, 0)]
nb = sum(P[(10, md)]['one'] and Nm(P[(10, md)]['tor']) < Nm(md) - 1e-9 for md in EX)
pr(f"H782b (DIP 10: z 5 wzbudzonych >= 4 jednomodowo na EBK o nizszej energii): {nb}/5 -> {'TAK' if nb >= 4 else 'NIE'}")
pr(f"  raport DIP 1: {sum(P[(1, md)]['one'] and Nm(P[(1, md)]['tor']) < Nm(md) - 1e-9 for md in EX)}/5; (1,0): DIP 1 -> {P[(1, (1, 0))]['tor']}, DIP 10 -> {P[(10, (1, 0))]['tor']}")
pr("  pierwsza zmiana modu dominujacego (okres) i wagi koncowe 13 najnizszych modow (DIP 10):")
for r in (1, 10):
    for md in MODES[1:]:
        q = P[(r, md)]['q']; ch = next((i + 1 for i in range(len(q)) if (int(q[i, 6]), int(q[i, 7])) != md), None)
        _, ww = W(os.path.join(O, f"P{r}_{md[0]}_{md[1]}.txt"))
        top = sorted(zip(ww[-1], LN), reverse=True)[:4]
        pr(f"    DIP {r} {str(md):8s} pierwsza zmiana @{ch}; koncowe wagi: " + ', '.join(f"{m}: {x:.3f}" for x, m in top))
# R
pr("=== R1, R10: starty 773/01-20, 20000 okresow, wobec 774 G (bez promieniowania) i 781 R (Larmor orbity) ===")
ids = [f"{i:02d}" for i in range(1, 21)]
G = {i: end(os.path.join(A, 'section-774', 'out', 'runs', f"G_{i}.txt"), 20000) for i in ids}
R = {r: {i: end(os.path.join(O, f"R{r}_{i}.txt"), 20000) for i in ids} for r in (1, 10)}
L81 = {r: {i: end(os.path.join(A, 'section-781', 'out', 'runs', f"R{r}_{i}.txt"), 20000) for i in ids} for r in (1, 10)}
S0 = sum(G[i]['one'] for i in ids); S = {r: sum(R[r][i]['one'] for i in ids) for r in R}
f = lambda e: f"{str(e['tor'])}{'*' if e['one'] else ' '} {e['d']:.4f}"
pr("  id  774G               782 R1             782 R10            781 R1        781 R10   (* = jednomodowy na EBK)")
for i in ids: pr(f"  {i}  {f(G[i])}   {f(R[1][i])}   {f(R[10][i])}   {str(L81[1][i]['tor']):9s} {str(L81[10][i]['tor']):9s}")
pr(f"  jednomodowe na EBK: 774G {S0}/20, R1 {S[1]}/20, R10 {S[10]}/20 (781: {sum(L81[1][i]['one'] for i in ids)}/20, {sum(L81[10][i]['one'] for i in ids)}/20); "
   f"ten sam torus co 774G: R1 {sum(R[1][i]['tor'] == G[i]['tor'] for i in ids)}/20, R10 {sum(R[10][i]['tor'] == G[i]['tor'] for i in ids)}/20")
for r in (1, 10):
    pr(f"  DIP {r}: mediana N konca {np.median([Nm(R[r][i]['tor']) for i in ids]):.3f} (774G {np.median([Nm(G[i]['tor']) for i in ids]):.3f}); max w(|m|>=5) {max(R[r][i]['w5'] for i in ids):.2e}; "
       "zwyciezcy: " + ', '.join(f"{x}: {c}" for x, c in Counter(R[r][i]['tor'] for i in ids).most_common()))
pr(f"H782c (R1 >= 774G - 3): {S[1]} >= {S0 - 3} -> {'TAK' if S[1] >= S0 - 3 else 'NIE'}")
# H782d
def frac(runs):
    st = [s for q in runs for s in steps(settled(q))]; n1 = sum(abs(a[1] - b[1]) == 1 for a, b in st)
    return n1, len(st), Counter(abs(a[1] - b[1]) for a, b in st)
n1, nt, cnt = frac([P[k]['q'] for k in P] + [R[r][i]['q'] for r in R for i in ids])
okD = nt > 0 and n1 / nt >= 0.8
pr(f"H782d (kroki kaskady w dol z |dm| = 1 we wszystkich P i R 782 >= 80 %): {n1}/{nt} = {n1 / max(nt, 1):.2f}; |dm|: {dict(sorted(cnt.items()))} -> {'TAK' if okD else 'NIE'}")
P81 = [np.loadtxt(os.path.join(A, 'section-781', 'out', 'runs', f"P{r}_{a}_{b}.txt"), comments='#', ndmin=2) for r in (1, 10) for a, b in EX + [(0, 0)]]
a1, at, ac = frac(P81 + [L81[r][i]['q'] for r in L81 for i in ids]); b1, bt, bc = frac([G[i]['q'] for i in ids])
pr(f"  raport: 781 (Larmor orbity) {a1}/{at} {dict(sorted(ac.items()))}; 774G (bez promieniowania) {b1}/{bt} {dict(sorted(bc.items()))}")
pr("  kroki 782: " + '; '.join(f"{a}->{b}" for q in [P[k]['q'] for k in P] + [R[r][i]['q'] for r in R for i in ids] for a, b in steps(settled(q))))

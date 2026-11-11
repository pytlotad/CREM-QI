# Audit 779: closed-form resonance condition -- verdicts H779a-d -> out/ana779.txt (offline, one-period integrals from overlap779.cpp)
import os, sys, numpy as np
D = os.path.dirname(os.path.abspath(__file__)); A = os.path.join(D, '..'); sys.path.insert(0, D)
from phase779 import Phi, wrap, fejer, gammaT, Nmode
pr = lambda *a: print(*a, flush=True)
GT, K = 0.0003, 300
def blocks(fn):
    B = {}; key = None
    for ln in open(os.path.join(D, 'out', fn)):
        p = ln.split()
        if p[0] == '##': key = p[1]; h = dict(beta=float(p[2]), caust=float(p[3]), J=float(p[4]), L=float(p[5]), N=float(p[6]), T=float(p[7])); B[key] = (h, [])
        else: B[key][1].append((int(p[0]), int(p[1]), float(p[4]), float(p[5])))
    out = {}
    for k, (h, rows) in B.items():
        a = np.array(rows); h['nr'] = a[:, 0].astype(int); h['m'] = a[:, 1].astype(int); h['I2'] = a[:, 2]; h['Pn'] = a[:, 3]
        h['Pc'] = Phi(h['J'], h['L'], h['nr'], h['m'], h['beta'], h['caust']); h['modes'] = list(zip(h['nr'], h['m']))
        h['F'] = h['I2'] * fejer(h['Pc'], K); h['E'] = h['I2'] / (4 * np.sin(h['Pc'] / 2) ** 2 + gammaT(GT, h['beta'], h['N']) ** 2); out[k] = h
    return out
def near3(J, L):
    T = sorted((np.hypot(J - nr - 0.5, L - m), (nr, m)) for nr in range(13) for m in range(-6, 7)); return [t[1] for t in T[:3]]
def ruleZ(h):
    c = near3(h['J'], h['L']); idx = [h['modes'].index(x) for x in c]; return c[int(np.argmin([abs(h['Pc'][i]) for i in idx]))]
am = lambda h, v: h['modes'][int(np.argmax(v))]
S = blocks('I_spectra.txt'); E = blocks('I_ebk.txt')
def pd(key):
    s, f = key.split(':'); fn = os.path.join(A, f'section-{s}', 'out', 'pred', f + '.txt')
    return [tuple(int(x) for x in ln.split()[:2]) for ln in open(fn) if not ln.startswith('#')]
# H779a
dmax = max(np.max(np.abs(wrap(h['Pc'] - h['Pn']))) for h in S.values())
rec = {k: am(h, h['F']) == pd(k)[0] for k, h in S.items()}; nrec = sum(rec.values()); n = len(S)
pr(f"H779a (i) max |Phi_zamk - Phi_num| (mod 2 pi) na {n} torusach: {dmax:.2e} (prog 1e-3) -> {'TAK' if dmax < 1e-3 else 'NIE'}")
for c, lab in ((-1.5707963267948966, 'CAUST -pi/2'), (0.0, 'CAUST 0')):
    R = [k for k, h in S.items() if abs(h['caust'] - c) < 1e-5]; pr(f"   rekonstrukcja {lab}: {sum(rec[k] for k in R)}/{len(R)}")
pr(f"H779a (ii) argmax |I|^2 F_300(Phi_zamk) = argmax P_D: {nrec}/{n} = {100 * nrec / n:.1f} % (prog 95 %) -> {'TAK' if nrec >= 0.95 * n else 'NIE'}")
pr(f"H779a -> {'TAK' if dmax < 1e-3 and nrec >= 0.95 * n else 'NIE'}")
bad = [k for k in S if not rec[k]]
for k in bad: h = S[k]; P = pd(k); pr(f"   niezgodne: {k}: rekonstr. {am(h, h['F'])}, P_D {P[0]}; ranga rekonstr. w P_D {P.index(am(h, h['F'])) + 1}")
# H779b
pr("H779b (punkty stale na torusach EBK n_r 0..8, |m| <= 5; argmax obwiedni E_n, GT 0,0003):")
okc = {}
for c, lab in (('-1.5707963267948966', 'CAUST -pi/2'), ('0', 'CAUST 0')):
    tot = 0; own = 0
    for b in ('0.05', '0.2', '0.8'):
        R = [h for k, h in E.items() if k.split(':')[1] == b and k.split(':')[2] == c]
        o = [am(h, h['E']) == (int(round(h['J'] - 0.5)), int(round(h['L']))) for h in R]; tot += len(R); own += sum(o)
        miss = [(int(round(h['J'] - 0.5)), int(round(h['L'])), am(h, h['E'])) for h, x in zip(R, o) if not x]
        oF = sum(am(h, h['F']) == (int(round(h['J'] - 0.5)), int(round(h['L']))) for h in R)
        pr(f"   {lab} beta {b}: wlasny mod = argmax E {sum(o)}/{len(R)}; = argmax F_300 {oF}/{len(R)}" + (f"; wyjatki (torus -> argmax): {miss[:12]}{' ...' if len(miss) > 12 else ''}" if miss and sum(o) > len(R) / 2 else ""))
        if c == '0': pr(f"      CAUST 0: Phi wlasnego modu na wlasnym torusie: {sorted(set(round(float(abs(h['Pc'][h['modes'].index((int(round(h['J'] - 0.5)), int(round(h['L']))))])), 6) for h in R))}")
    okc[c] = (own, tot)
o1, t1 = okc['-1.5707963267948966']; o0, t0 = okc['0']
pr(f"H779b: CAUST -pi/2 {o1}/{t1} (>= 95 %), CAUST 0 {o0}/{t0} (<= 5 %) -> {'TAK' if o1 >= 0.95 * t1 and o0 <= 0.05 * t0 else 'NIE'}")
# H779c / H779d(i): runs of 778 at k_d - 10
ORIG = {'774G': 'section-774/out/runs/G_{}.txt', '775B05': 'section-775/out/runs/B05_{}.txt', '775B80': 'section-775/out/runs/B80_{}.txt'}
def nearest(j, l): return near3(j, l)[0]
J8 = [l.split() for l in open(os.path.join(A, 'section-778', 'jobs778.txt'))]
hz = he = hp = hn = hf = 0; rows = []
for a, i, b, kd, _ in J8:
    o = np.loadtxt(os.path.join(A, ORIG[a].format(i)), comments='#', ndmin=2); win = nearest(*o[-20:, 2:4].mean(0))
    h = S[f"778:{a}_{i}_m10"]; z = ruleZ(h); e = am(h, h['E']); p = pd(f"778:{a}_{i}_m10")[0]; nn = nearest(h['J'], h['L']); f = am(h, h['F'])
    hz += z == win; he += e == win; hp += p == win; hn += nn == win; hf += f == win; rows.append((a, i, win, z, e, p, nn))
n8 = len(J8)
pr(f"H779c (Z: 3 najblizsze torusy EBK, min |Phi_zamk|) w k_d - 10: Z {hz}/{n8} = {100 * hz / n8:.1f} %, P_near {hn}/{n8} = {100 * hn / n8:.1f} % -> "
   f"{'TAK' if hz >= 0.6 * n8 and 100 * (hz - hn) / n8 >= 15 else 'NIE'}")
pr(f"H779d (i) obwiednia E (z calkami) {he}/{n8}; rekonstrukcja F_300 {hf}/{n8}; P_D (778) {hp}/{n8}")
for arm in ORIG:
    R = [r for r in rows if r[0] == arm]
    pr(f"   {arm} n {len(R)}: Z {sum(r[3] == r[2] for r in R)}, E {sum(r[4] == r[2] for r in R)}, P_D {sum(r[5] == r[2] for r in R)}, P_near {sum(r[6] == r[2] for r in R)}")
pr("   szczegoly (ramie id zwyc. Z E P_D P_near):")
for r in rows: pr(f"     {r[0]:6s} {r[1]} {str(r[2]):9s} {str(r[3]):9s} {str(r[4]):9s} {str(r[5]):9s} {str(r[6])}")
# H779d (ii), (iii)
for c, lab in ((-1.5707963267948966, 'CAUST -pi/2'), (0.0, 'CAUST 0')):
    R = [k for k, h in S.items() if abs(h['caust'] - c) < 1e-5]
    pr(f"H779d (ii) {lab}: Z = argmax P_D {sum(ruleZ(S[k]) == pd(k)[0] for k in R)}/{len(R)};  (iii) argmax F_300 != argmax E w {sum(am(S[k], S[k]['F']) != am(S[k], S[k]['E']) for k in R)}/{len(R)}")
# H779d (iv) ends of 774 / 775 that are not single-mode on EBK or have offset > 0.01
pr("H779d (iv) konce 774 G / 775 B05 / 775 B80 z odl. > 0,01 lub w_max < 0,95: Phi_zamk (CAUST -pi/2) trzech najblizszych modow na torusie konca")
for arm, pat, beta, ids in (('774G', 'section-774/out/runs/G_{}.txt', 0.2, range(1, 41)), ('775B05', 'section-775/out/runs/B05_{}.txt', 0.05, range(1, 21)), ('775B80', 'section-775/out/runs/B80_{}.txt', 0.8, range(1, 21))):
    for k in ids:
        q = np.loadtxt(os.path.join(A, pat.format(f"{k:02d}")), comments='#', ndmin=2); j, l = q[-20:, 2:4].mean(0); c3 = near3(j, l); d = np.hypot(j - c3[0][0] - 0.5, l - c3[0][1])
        if d > 0.01 or q[-1, 5] < 0.95:
            pr(f"   {arm} {k:02d}: J {j:.4f} L {l:+.4f} w_max {q[-1, 5]:.3f} odl. {d:.4f}; " + ', '.join(f"{x} {float(Phi(j, l, x[0], x[1], beta, -np.pi / 2)):+.4f}" for x in c3))
# H779d (v) linear off-resonant residual offset on EBK tori (CAUST -pi/2)
pr("H779d (v) liniowe wymuszenie nierezonansowe na torusie EBK (GT 0,0003): suma wag pozostalych modow i przesuniecie torusa")
for b in ('0.05', '0.2', '0.8'):
    R = [h for k, h in E.items() if k.split(':')[1] == b and k.split(':')[2] == '-1.5707963267948966']
    best = []
    for h in R:
        own = h['modes'].index((int(round(h['J'] - 0.5)), int(round(h['L'])))); w = h['E'] / h['E'][own]; w[own] = 0
        dJ = np.sum(w * (h['nr'] + 0.5 - h['J'])) / (1 + w.sum()); dL = np.sum(w * (h['m'] - h['L'])) / (1 + w.sum()); best.append((np.hypot(dJ, dL), w.sum(), h['modes'][own]))
    s = sorted(best, reverse=True); z2 = [x for x in best if x[2][0] == 0 and abs(x[2][1]) == 2]
    pr(f"   beta {b}: max przesuniecie {s[0][0]:.2e} przy {s[0][2]} (suma wag {s[0][1]:.2e}); mediana {np.median([x[0] for x in best]):.2e}; (0,+-2): " + ', '.join(f"{x[0]:.2e}" for x in z2))

# Audit 874: prospective test of the w / r_peri rule -- verdicts H874c, q, n, o, g, l -> out/ana874.txt.  Written after prereg874.txt and the launch, before any result.
import os, glob, math, numpy as np
D = os.path.dirname(os.path.abspath(__file__)); O = os.path.join(D, 'out', 'runs'); pr = lambda *a: print(*a, flush=True)
def load(p):
    hdr = {}; rows = []; stop = False
    for ln in open(p):
        if ln.startswith('# 871'):
            t = ln.split()[2:]; hdr = {t[i]: float(t[i + 1]) for i in range(0, len(t) - 1, 2)}
        elif ln.startswith('# STOP'): stop = True
        elif not ln.startswith('#') and ln.strip(): rows.append([float(x) for x in ln.split()])
    return hdr, (np.array(rows) if rows else np.zeros((0, 13))), stop
def summ(p):
    name = os.path.basename(p)[:-4]; hdr, q, stop = load(p); n = len(q)
    N0, e, w = hdr['N0'], hdr['E0'], hdr['W']; mu = hdr['M1'] * hdr['M2'] / (hdr['M1'] + hdr['M2']); a = N0 ** 2 / mu
    T = 2 * math.pi * N0 ** 3 / mu; rp = a * (1 - e)
    r = dict(name=name, sp='H' if hdr['M1'] > 2 else 'Ps', N0=N0, e=e, w=w, dx=hdr['DX'], kappa=hdr['KAPPA'], gam=hdr['GAM'], q=w / rp, GT=hdr['GAM'] * T,
             n=n, stop=stop, T=T)
    if n < 10: r.update(ok=False, cls='SPADA' if stop else '?', lock=False, orb=(q[-1, 1] / T if n else 0)); return r
    k0 = int(0.8 * n); L0 = q[0, 4]; Lk = q[k0:, 4].mean(); J = q[k0:, 2].mean()
    drift = abs(q[int(0.9 * n):, 2].mean() - q[int(0.7 * n):k0, 2].mean()); fin = bool(np.all(np.isfinite(q))) and not stop
    rat = Lk / L0; cls = 'SPADA' if (stop or rat < 0.2) else ('TRZYMA' if rat > 0.8 else 'POSREDNIE')
    half = np.where(q[:, 4] < 0.5 * L0)[0]; oh = q[half[0], 1] / T if len(half) else float('nan')
    lock = fin and drift < 0.002 and abs(J - round(N0)) < 0.005
    r.update(ok=True, L0=L0, Lk=Lk, rat=rat, rmin=q[:, 9].min(), J=J, drift=drift, lock=lock, cls=cls, ohalf=oh, orb=q[-1, 1] / T,
             A0=q[0, 3], Ak=q[-1, 3], E0=q[0, 6], Ek=q[-1, 6])
    return r
def pred(r): return 'TRZYMA' if r['q'] <= 1.5 else ('SPADA' if r['q'] >= 3 else None)
def row(r):
    s = f"  {r['name']:9s} {r['sp']:2s} n {round(r['N0'])} e {r['e']:.2f} q {r['q']:.2f} w {r['w']:.3f} dx {r['dx']:.4f} Gamma {r['gam']:.2e} GammaT {r['GT']:.4f} obiegi {r['orb']:6.0f}"
    if not r['ok']: return s + f"  cykli {r['n']} -- {r['cls']}"
    return s + (f"  L0 {r['L0']:.4f} L_k/L0 {r['rat']:.4f} rmin {r['rmin']:.4f} J {r['J']:.6f} dryf {r['drift']:.1e} LOCK {'T' if r['lock'] else 'N'}"
                f"  obieg(L0/2) {r['ohalf']:.0f}  {r['cls']}{'  STOP' if r['stop'] else ''}  przewidywanie {pred(r)}")
R = {os.path.basename(p)[:-4]: summ(p) for p in sorted(glob.glob(os.path.join(O, '*.txt')))}
pr("=== wszystkie biegi 874 ==="); [pr(row(r)) for r in R.values()]
k = R['K']; dl, da, de = abs(k['Lk'] / k['L0'] - 1), abs(k['Ak'] / k['A0'] - 1), abs(k['Ek'] / k['E0'] - 1)
pr(f"H874c (kappa 0, e 0,97: L {dl:.1e}, A {da:.1e}, E {de:.1e} < 1e-6): {'TAK' if (dl < 1e-6 and da < 1e-6 and de < 1e-6 and not k['stop']) else 'NIE'}")
main = [R[x] for x in R if x.split('_')[0] in ('A', 'P', 'N2', 'N3')]
for sp in ('H', 'Ps'):
    sub = [r for r in main if r['sp'] == sp]; pr(f"  {sp}: zgodne z regula {sum(r['cls'] == pred(r) for r in sub)}/{len(sub)}; niezgodne {[(r['name'], r['cls']) for r in sub if r['cls'] != pred(r)]}")
g = sum(r['cls'] == pred(r) for r in main)
pr(f"H874q (GLOWNA; zgodne >= 11/12 -> {g}/{len(main)}): {'TAK' if (len(main) == 12 and g >= 11) else 'NIE'}")
nn = [r for r in main if round(r['N0']) >= 2]; gn = sum(r['cls'] == 'SPADA' for r in nn)
pr(f"H874n (n >= 2 przy q 4: SPADA 4/4 -> {gn}/{len(nn)} {[(r['name'], r['cls']) for r in nn]}): {'TAK' if (len(nn) == 4 and gn == 4) else 'NIE'}")
n1 = [r for r in main if round(r['N0']) == 1]; g1 = sum(r['cls'] == pred(r) for r in n1)
pr(f"H874o (n = 1 zgodne 8/8 -> {g1}/{len(n1)}): {'TAK' if (len(n1) == 8 and g1 == 8) else 'NIE'}")
G1, G2 = R['G1'], R['G2']
pr(f"H874g (G1 TRZYMA -> {G1['cls']}, G2 SPADA -> {G2['cls']}): {'TAK' if (G1['cls'] == 'TRZYMA' and G2['cls'] == 'SPADA') else 'NIE'}")
tr = [r for r in R.values() if r['cls'] == 'TRZYMA' and r['name'] not in ('K', 'G2')]
pr(f"H874l (TRZYMA -> LOCK: {sum(r['lock'] for r in tr)}/{len(tr)}; bez LOCK {[r['name'] for r in tr if not r['lock']]}): {'TAK' if all(r['lock'] for r in tr) else 'NIE'}")
pr(f"  G2 (bez reguly): LOCK {'T' if G2.get('lock') else 'N'}")

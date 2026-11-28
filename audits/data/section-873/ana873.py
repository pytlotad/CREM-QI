# Audit 873: n = 1 l ~ 0 fall -- verdicts H873c, g, w, x, s -> out/ana873.txt.  Written after prereg873.txt and the launch, before any result.
import os, glob, math, tarfile, numpy as np
D = os.path.dirname(os.path.abspath(__file__)); O = os.path.join(D, 'out', 'runs'); pr = lambda *a: print(*a, flush=True)
S72 = '/tmp/claude-1000/-home-teddy-Projekty-CREM-QI/4fa83d5c-8878-4cfd-898e-1a9b53a01ab2/scratchpad/r872'
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
    b = a * math.sqrt(1 - e * e); pc = a * (1 - e * e)
    r = dict(name=name, N0=N0, e=e, w=w, dx=hdr['DX'], kappa=hdr['KAPPA'], bw=b / w, pw=pc / w, n=n, stop=stop, sp='H' if hdr['M1'] > 2 else 'Ps')
    if n < 10: r.update(ok=False, cls='SPADA' if stop else '?', lock=False); return r
    k0 = int(0.8 * n); L0 = q[0, 4]; Lk = q[k0:, 4].mean(); J = q[k0:, 2].mean()
    drift = abs(q[int(0.9 * n):, 2].mean() - q[int(0.7 * n):k0, 2].mean()); fin = bool(np.all(np.isfinite(q))) and not stop
    rat = Lk / L0; cls = 'SPADA' if (stop or rat < 0.2) else ('TRZYMA' if rat > 0.8 else 'POSREDNIE')
    half = np.where(q[:, 4] < 0.5 * L0)[0]; th = q[half[0], 1] if len(half) else float('nan')
    lock = fin and drift < 0.002 and abs(J - round(N0)) < 0.005
    r.update(ok=True, L0=L0, Lk=Lk, rat=rat, rmin=q[:, 9].min(), J=J, drift=drift, lock=lock, cls=cls, thalf=th,
             A0=q[0, 3], Ak=q[-1, 3], E0=q[0, 6], Ek=q[-1, 6], Lend=q[-1, 4])
    return r
def row(r):
    if not r['ok']: return f"  {r['name']:12s} {r['sp']:2s} cykli {r['n']} -- {r['cls']}"
    return (f"  {r['name']:12s} {r['sp']:2s} N0 {r['N0']:.2f} e {r['e']:.2f} w {r['w']:.3f} dx {r['dx']:.3f}  L0 {r['L0']:.4f} L_k {r['Lk']:.5f} "
            f"L_k/L0 {r['rat']:.4f}  rmin {r['rmin']:.4f}  J {r['J']:.6f} dryf {r['drift']:.1e} LOCK {'T' if r['lock'] else 'N'}  b/w {r['bw']:.2f} p/w {r['pw']:.2f}"
            f"  t(L0/2) {r['thalf']:.0f}  {r['cls']}{'  STOP' if r['stop'] else ''}")
R = {os.path.basename(p)[:-4]: summ(p) for p in sorted(glob.glob(os.path.join(O, '*.txt'))) if not os.path.basename(p).startswith('XL_')}  # XL_ = POST HOC (UZUPELNIENIE 1), outside the prereg set
pr("=== wszystkie biegi 873 ===")
for r in R.values(): pr(row(r))
# H873c
k = R['K']; dl, da, de = abs(k['Lk'] / k['L0'] - 1), abs(k['Ak'] / k['A0'] - 1), abs(k['Ek'] / k['E0'] - 1)
okc = dl < 1e-6 and da < 1e-6 and de < 1e-6 and not k['stop']
pr(f"H873c (kappa 0: L {dl:.1e}, A {da:.1e}, E {de:.1e} < 1e-6): {'TAK' if okc else 'NIE'}")
Dd = [R[x] for x in R if x.startswith('D_')]; okg = sum(r['cls'] == 'SPADA' for r in Dd) == 2
pr(f"H873g (dx 0,025: SPADA 2/2 -> {[r['cls'] for r in Dd]}): {'TAK' if okg else 'NIE'}")
W5 = [R[x] for x in R if x.startswith('W05_')]; okw = True
for sp in ('H', 'Ps'):
    sub = [r for r in W5 if r['sp'] == sp]; g = sum(r['cls'] == 'TRZYMA' and r['lock'] for r in sub)
    pr(f"  W05 {sp}: TRZYMA i LOCK {g}/{len(sub)}"); okw &= g == 2
pr(f"H873w (GLOWNA; w 0,05: TRZYMA i LOCK 4/4, wodor 2/2 i pozytonium 2/2): {'TAK' if okw else 'NIE'}")
Xs = [R[x] for x in R if x.startswith('X_')]; okx = sum(r['cls'] == 'SPADA' for r in Xs) == 2
pr(f"H873x (n = 2 przy duzym w: SPADA 2/2 -> {[r['cls'] for r in Xs]}): {'TAK' if okx else 'NIE'}")
# H873s
os.makedirs(S72, exist_ok=True)
with tarfile.open(os.path.join(D, '..', 'section-872', 'runs', 'runs872.tar.xz')) as t:
    mem = [m for m in t.getmembers() if os.path.basename(m.name).startswith('GZ') and m.name.endswith('.txt')]
    t.extractall(S72, members=mem, filter='data')
R72 = [summ(p) for p in sorted(glob.glob(os.path.join(S72, 'runs', 'GZ*.txt')))]
pr("=== 872 GZ (w 0,25, dx 0,05) wg tych samych definicji ===")
for r in R72: pr(row(r))
pool = [r for r in list(R.values()) + R72 if r['kappa'] != 0 and r['cls'] in ('TRZYMA', 'SPADA')]
oks = False
for var, thr in (('bw', 3.0), ('pw', 1.2)):
    good = sum((r[var] >= thr) == (r['cls'] == 'TRZYMA') for r in pool); bad = [r['name'] for r in pool if (r[var] >= thr) != (r['cls'] == 'TRZYMA')]
    pr(f"  {var} prog {thr}: poprawnie {good}/{len(pool)} ({100 * good / len(pool):.0f} %); bledne {bad}"); oks |= good >= 0.9 * len(pool)
pr(f"  POSREDNIE pominiete: {[r['name'] for r in list(R.values()) + R72 if r['cls'] == 'POSREDNIE']}")
pr(f"H873s (>= 90 % jedna ze zmiennych): {'TAK' if oks else 'NIE'}")

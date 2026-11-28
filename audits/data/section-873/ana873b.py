# Audit 873 POST HOC (written after out/ana873.txt): (1) orbit-count censoring of arm X -> XL reruns with TEND x5;
# (2) early per-orbit L loss rate; (3) candidate post-hoc variable w / r_peri (r_peri = a (1 - e)).  -> out/ana873b.txt
import os, glob, math, importlib.util, contextlib, io
D = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location('a873', os.path.join(D, 'ana873.py')); A = importlib.util.module_from_spec(spec)
with contextlib.redirect_stdout(io.StringIO()): spec.loader.exec_module(A)   # reuses load/summ and R, R72 exactly as scored
pr = A.pr
def extra(p, r):
    hdr, q, stop = A.load(p); mu = hdr['M1'] * hdr['M2'] / (hdr['M1'] + hdr['M2']); a = hdr['N0'] ** 2 / mu
    T = 2 * math.pi * math.sqrt(mu * a ** 3); n = len(q); m = max(int(0.25 * n), 2)
    import numpy as np
    s = np.polyfit(q[:m, 1], q[:m, 4], 1)[0] if n > 3 else float('nan')
    r.update(T=T, orb=q[-1, 1] / T if n else 0, rate=s * T / q[0, 4] if n > 3 else float('nan'),
             orbhalf=r.get('thalf', float('nan')) / T, wr=hdr['W'] / (a * (1 - hdr['E0'])))
    return r
files = {os.path.basename(p)[:-4]: p for p in glob.glob(os.path.join(A.O, '*.txt'))}
files.update({os.path.basename(p)[:-4]: p for p in glob.glob(os.path.join(A.S72, 'runs', 'GZ*.txt'))})
XL = {k: A.summ(files[k]) for k in sorted(files) if k.startswith('XL_')}
allr = {**A.R, **{r['name']: r for r in A.R72}, **XL}
pr("=== POST HOC: obiegi, tempo utraty L na obieg (pierwsza 1/4 biegu), w / r_peri ===")
for k in sorted(allr):
    r = allr[k]
    if r.get('kappa', 1) == 0 or not r['ok']: continue
    extra(files[k], r)
    pr(f"  {k:12s} {r['sp']:2s} w {r['w']:.2f}  obiegi {r['orb']:6.0f}  t(L0/2) w obiegach {r['orbhalf']:6.0f}  dL/L0 na obieg {r['rate']:+.2e}  "
       f"L_k/L0 {r['rat']:.4f}  LOCK {'T' if r['lock'] else 'N'}  w/r_peri {r['wr']:.2f}  b/w {r['bw']:.2f} p/w {r['pw']:.2f}  {r['cls']}")
pr("=== POST HOC: X (TEND 6e4) wobec XL (TEND 3e5 / 6e5) ===")
for k in XL: r = XL[k]; pr(A.row(r))
pool = [r for r in list(A.R.values()) + A.R72 if r['kappa'] != 0 and r['cls'] in ('TRZYMA', 'SPADA') and not r['name'].startswith('X_')] + \
       [r for r in XL.values() if r['cls'] in ('TRZYMA', 'SPADA')]
for k, r in XL.items(): extra(files[k], r)
pr(f"=== POST HOC: klasyfikacja puli z X zastapionym przez XL ({len(pool)} biegow TRZYMA/SPADA) ===")
for var, thr in (('bw', 3.0), ('pw', 1.2), ('wr', 2.2)):
    if var == 'wr': ok = [(r[var] < thr) == (r['cls'] == 'TRZYMA') for r in pool]
    else: ok = [(r[var] >= thr) == (r['cls'] == 'TRZYMA') for r in pool]
    pr(f"  {var} prog {thr}{' (POST HOC, z tych danych)' if var == 'wr' else ' (prereg)'}: poprawnie {sum(ok)}/{len(pool)}; bledne {[r['name'] for r, o in zip(pool, ok) if not o]}")

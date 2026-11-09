# Audit 754: analysis -> out/ana754.txt  (surowe przebiegi biegu 2: tar -xJf out/runs754.tar.xz -C out)
import os, numpy as np
from scipy import stats
D = os.path.dirname(os.path.abspath(__file__)); O = os.path.join(D, 'out', 'runs'); pr = lambda *a: print(*a, flush=True)
au_t = 2.4188843265857e-17; alpha = 7.2973525643e-3; R1 = 128 * np.sqrt(6) / 243 / np.sqrt(3); GF = 4 / 3 * alpha**3 * 0.375**3 * R1**2
tauH = 1 / (GF / (1 + 1 / 1836.15267343)) * au_t * 1e9; tauPs = 1 / (GF * 0.5) * au_t * 1e9
import glob
def L(f):   # bieg 2: tekst z bohm754.cpp, czesci _p0.._p3 jednego zespolu (lub jeden plik)
    fs = sorted(glob.glob(os.path.join(O, f[:-4] + '_p*.txt'))) or [os.path.join(O, f[:-4] + '.txt')]
    ws, te, cr, vmax = [], {}, {}, 0.0
    for x in fs:
        a = np.loadtxt(x, comments='#', ndmin=2); t, Pe = a[:, 0], a[:, 1]; ws.append(a[:, 2:])
        for line in open(x):
            p = line.split()
            if line.startswith('# TE'): te[int(p[2])] = float(p[3])
            elif line.startswith('# CR'): cr[int(p[2])] = int(p[3])
            elif line.startswith('# steps'): vmax = max(vmax, float(p[4]))
    k = sorted(te); return dict(t=t, Pe=Pe, w=np.hstack(ws), te=np.array([te[i] for i in k]), cross=np.array([cr[i] for i in k]), vmax=vmax)
pr(f"skale: tau wodoru (Fermi) {tauH:.5f} ns (pomiar 1,600(4) ns, Bickel i Goodman 1966); Ps (p-Ps i o-Ps, linia E1 bez spinu) {tauPs:.4f} ns")
pr("  linie: wodor 10,19872 eV / 121,5684 nm / 2,46604e15 Hz; Ps 5,10213 eV / 243,0045 nm / 1,23369e15 Hz")
g = L('g_N800.npz'); pr(f"H754b: start |g, proznia>: max|v| = {float(g['vmax']):.1e} (steps line), min w = {g['w'].min():.6f} -> {'TAK' if g['vmax'] == 0 and g['w'].min() == 1 else 'NIE'}")
def H(f, name):
    d = L(f); t, w, Pe, te, cr = d['t'], d['w'].astype(float), d['Pe'], d['te'], d['cross']; M = w.shape[1]
    mw = w.mean(1); se = w.std(1) / np.sqrt(M); dev = np.abs(mw - (1 - Pe)); dexp = np.abs(Pe - np.exp(-t))
    j = np.argmin(abs(t - 1)); sharp = np.mean((w[j] < 0.1) | (w[j] > 0.9))
    pr(f"=== {name} ({f}; M {M}) ===")
    pr(f"  max |<w> - (1 - |c_e|^2)| = {dev.max():.4f} (przy t {t[dev.argmax()]:.2f}, blad sredniej {se[dev.argmax()]:.4f}; max dev/blad {np.max(dev[1:]/se[1:]):.2f}); max ||c_e|^2 - e^-t| = {dexp.max():.2e}")
    pr(f"  udzial biegow z w < 0,1 lub w > 0,9 przy t = tau: {sharp:.3f}")
    return d, dev.max(), dexp.max(), sharp
d, dv, de, sh = H('e_N800.npz', 'glowny, N 800, dt 0,25')
pr(f"H754a: {'TAK' if dv < 0.03 and de < 0.01 else 'NIE'};  H754c: {'TAK' if sh >= 0.8 else 'NIE'}")
t, w, te, cr = d['t'], d['w'].astype(float), d['te'], d['cross']; M = w.shape[1]
for tq in [0.25, 0.5, 1, 2, 3, 5]:
    j = np.argmin(abs(t - tq)); h = np.histogram(w[j], bins=[0, .1, .3, .5, .7, .9, 1.0001])[0] / M
    pr(f"  t = {tq} tau: <w> {w[j].mean():.4f} (1 - e^-t {1-np.exp(-tq):.4f}); kwantyle 10/50/90 % " + ' / '.join(f"{x:.3f}" for x in np.percentile(w[j], [10, 50, 90]))
       + "; histogram [0;0,1;0,3;0,5;0,7;0,9;1] " + ' '.join(f"{x:.3f}" for x in h))
u = te[te >= 0]; nc = int((te < 0).sum()); T = t[-1]; mle = (u.sum() + nc * T) / len(u); cv = u.std() / u.mean()
cdf = lambda x: (1 - np.exp(-x / mle)) / (1 - np.exp(-T / mle)); ks = stats.kstest(u, cdf)
ks1 = stats.kstest(u, lambda x: (1 - np.exp(-x)) / (1 - np.exp(-T)))
pr(f"=== H754d: t_e (pierwsze w >= 1/2) ===")
pr(f"  przejsc {len(u)}/{M}, cenzurowanych {nc}; srednia MLE {mle:.4f} tau; CV (nieenzur.) {cv:.3f}; min {u.min():.4f} tau")
pr(f"  kwantyle 10/25/50/75/90 %: " + ' / '.join(f"{q:.3f}" for q in np.percentile(u, [10, 25, 50, 75, 90])) + f" (wykladnik tau: " + ' / '.join(f"{-np.log(1-p):.3f}" for p in [.1, .25, .5, .75, .9]) + ")")
pr(f"  KS wobec wykladnika (srednia MLE): D {ks.statistic:.4f}, p {ks.pvalue:.2e};  wobec wykladnika ze srednia tau: D {ks1.statistic:.4f}, p {ks1.pvalue:.2e}")
pr(f"H754d: {'TAK' if ks.pvalue > 0.01 and abs(mle - 1) < 0.1 else 'NIE'}")
pr(f"  wodor: mediana t_e {np.median(u)*tauH:.3f} ns; Ps (p-Ps i o-Ps): {np.median(u)*tauPs:.3f} ns")
pr(f"=== H754e (raport) ===")
pr(f"  przejscia w przez 1/2 na bieg: " + ', '.join(f"{k}: {np.mean(cr == k):.3f}" for k in range(0, 8)) + f"; srednio {cr.mean():.2f}")
w1 = w[np.argmin(abs(t - 1))]; pr(f"  korelacja w(tau) z w(2 tau) w biegach: {np.corrcoef(w1, w[np.argmin(abs(t - 2))])[0,1]:.3f}")
for f, name in [('e_N800_f025.npz', 'kontrola krok f 0,25'), ('e_N400.npz', 'N 400'), ('e_N1600.npz', 'N 1600')]: H(f, name)
dd = L('e_N800_f025.npz'); m = dd['w'].shape[1]; dw = np.abs(dd['w'] - w[:, :m]).max(0)
pr(f"  f 0,25 wobec 0,5 (te same {m} biegow): max |w roznica| {dw.max():.3e}; biegow z roznica > 0,01: {(dw > 0.01).sum()}; t_e: max roznica {np.abs(dd['te'] - te[:m]).max():.4f} tau, mediana {np.median(np.abs(dd['te'] - te[:m])):.1e}")

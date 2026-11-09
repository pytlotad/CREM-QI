# Audit 755: analysis -> out/ana755.txt  (surowe przebiegi: tar -xJf out/runs755.tar.xz -C out, jesli spakowane)
import os, numpy as np
from scipy import stats
D = os.path.dirname(os.path.abspath(__file__)); O = os.path.join(D, 'out', 'runs'); pr = lambda *a: print(*a, flush=True)
au_t = 2.4188843265857e-17; alpha = 7.2973525643e-3; R1 = 128 * np.sqrt(6) / 243 / np.sqrt(3); GF = 4 / 3 * alpha**3 * 0.375**3 * R1**2
tauH = 1 / (GF / (1 + 1 / 1836.15267343)) * au_t * 1e9; tauPs = 1 / (GF * 0.5) * au_t * 1e9; kap, Dl, tau = 0.15, 2.0, 100.0; delay = Dl / kap / tau
pr(f"skale: tau wodoru (Fermi) {tauH:.5f} ns (pomiar 1,600(4) ns, Bickel i Goodman 1966); Ps (p-Ps i o-Ps, linia E1 bez spinu) {tauPs:.4f} ns")
pr("  linie: wodor 10,19872 eV / 121,5684 nm / 2,46604e15 Hz; Ps 5,10213 eV / 243,0045 nm / 1,23369e15 Hz")
d = np.load(os.path.join(O, 'main_h005.npz')); c = np.load(os.path.join(O, 'ctrl_h0025.npz'))
t, X, Pe = d['t'], d['X'], d['Pe']; T = t[-1]
m = (t >= 0.5) & (t <= 4); Geff = -np.polyfit(t[m], np.log(Pe[m]), 1)[0]; tauE = 1 / Geff
pr(f"P_e: " + ' / '.join(f"{Pe[np.argmin(abs(t-x))]:.4f}" for x in [0.1, 0.5, 1, 2, 3]) + " przy 0,1/0,5/1/2/3 tau (e^-t " + ' / '.join(f"{np.exp(-x):.4f}" for x in [0.1, 0.5, 1, 2, 3])
   + f"); Gamma_eff/Gamma = {Geff:.4f} (dopasowanie [0,5; 4] tau); norma na koncu {float(d['norm_end']):.8f}")
def ksx(xs, rho):
    cdf = np.cumsum(rho); cdf /= cdf[-1]; return stats.kstest(xs, lambda z: np.interp(z, X, cdf))
okA = True
for name, key in [('P', 'XP'), ('F', 'XF')]:
    for tt, rk in [(1, 'rho1'), (2, 'rho2')]:
        r = ksx(d[key][np.argmin(abs(t - tt))].astype(float), d[rk]); okA &= r.pvalue > 0.01
        pr(f"  ekwiwariancja {name}, t = {tt} tau: KS x wobec rho: D {r.statistic:.4f}, p {r.pvalue:.3f}")
WF = d['WF'].astype(float); dev = np.abs(WF.mean(1) - (1 - Pe)); se = WF.std(1) / np.sqrt(WF.shape[1])
pr(f"  F: max |<w> - (1 - P_e)| = {dev.max():.4f} (przy {t[dev.argmax()]:.2f} tau; max dev/blad {np.max(dev[1:]/np.maximum(se[1:],1e-12)):.2f})")
okA &= dev.max() < 0.03
pr(f"H755a: {'TAK' if okA else 'NIE'}")
def expo(tr, name):
    te = tr - delay; reg = tr >= 0; u = te[reg]; nc = int((~reg).sum()); Tc = T - delay
    mle = (u.sum() + nc * Tc) / len(u); cdf = lambda z: (1 - np.exp(-z / mle)) / (1 - np.exp(-Tc / mle))
    ks = stats.kstest(u, cdf); ksg = stats.kstest(u, lambda z: (1 - np.exp(-z / tauE)) / (1 - np.exp(-Tc / tauE)))
    pr(f"=== {name}: t_e = t_r - {delay:.4f} tau; rejestracji {reg.sum()}/{len(tr)}, cenzurowanych {nc} ===")
    pr(f"  srednia MLE {mle:.4f} tau, 1/Gamma_eff {tauE:.4f} tau (stosunek {mle/tauE:.4f}); CV {u.std()/u.mean():.3f}; min t_e {u.min():+.4f} tau")
    pr(f"  kwantyle 10/25/50/75/90 %: " + ' / '.join(f"{q:.3f}" for q in np.percentile(u, [10, 25, 50, 75, 90])) + " (wykladnik 1/Gamma_eff: "
       + ' / '.join(f"{-tauE*np.log(1-p):.3f}" for p in [.1, .25, .5, .75, .9]) + ")")
    pr(f"  KS wobec wykladnika (MLE): D {ks.statistic:.4f}, p {ks.pvalue:.3f};  wobec 1/Gamma_eff: D {ksg.statistic:.4f}, p {ksg.pvalue:.3f}")
    pr(f"  wodor: mediana t_e {np.median(u)*tauH:.3f} ns (1/Gamma_eff {tauE*tauH:.4f} ns); Ps (p-Ps i o-Ps): {np.median(u)*tauPs:.3f} ns ({tauE*tauPs:.4f} ns)")
    return ks.pvalue > 0.01 and abs(mle / tauE - 1) < 0.05
okB = expo(d['trP'], 'P (tylko wskaznik)'); pr(f"H755b: {'TAK' if okB else 'NIE'}")
okC = expo(d['trF'], 'F (pole + wskaznik)'); pr(f"H755c: {'TAK' if okC else 'NIE'}")
j = np.argmin(abs(t - 1)); sh = np.mean((WF[j] < 0.1) | (WF[j] > 0.9))
pr(f"H755d: udzial w < 0,1 lub > 0,9 przy tau = {sh:.3f} -> {'TAK' if sh >= 0.8 else 'NIE'}")
for tt in [0.25, 0.5, 1, 2]:
    k = np.argmin(abs(t - tt)); pr(f"  t = {tt} tau: udzial ostrych {np.mean((WF[k] < 0.1) | (WF[k] > 0.9)):.3f}; histogram w [0;0,1;0,5;0,9;1] " + ' '.join(f"{x:.3f}" for x in np.histogram(WF[k], [0, .1, .5, .9, 1.0001])[0] / WF.shape[1]))
pr("=== H755e (raport) ===")
for name, tr, x0 in [('P', d['trP'], d['xP0']), ('F', d['trF'], d['xF0'])]:
    r = tr >= 0; rho = stats.spearmanr(x0[r], tr[r]).statistic; pr(f"  {name}: korelacja rang t_r z x(0) = {rho:+.4f} (R^2 rang {rho**2:.4f})")
n = len(d['trF']); a, b = d['trF'], d['trP'][:n]; r = (a >= 0) & (b >= 0)
pr(f"  F wobec P bieg w bieg (te same x(0), {r.sum()} biegow): mediana |t_r^F - t_r^P| {np.median(abs(a[r]-b[r])):.4f} tau, 90 % {np.percentile(abs(a[r]-b[r]), 90):.4f} tau")
m2 = len(c['trF']); a2 = c['trF']; r2 = (a2 >= 0) & (a[:m2] >= 0); dd = abs(a2[r2] - a[:m2][r2])
pr(f"  kontrola h 0,025 (F, {m2} biegow): |t_r roznica| mediana {np.median(dd):.2e}, max {dd.max():.4f} tau, >0,01 tau: {(dd > 0.01).sum()}; P: max {np.abs(c['trP'] - d['trP'][:len(c['trP'])]).max():.2e} tau")

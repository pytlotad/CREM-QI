# Audit 755 POST HOC (after ana755): (1) tail of t_e (memorylessness) for t_e > 0.5 tau; (2) absolute threshold x >= x_th: Bohm fraction vs QM
import os, numpy as np
from scipy import stats
D = os.path.dirname(os.path.abspath(__file__)); O = os.path.join(D, 'out', 'runs'); pr = lambda *a: print(*a, flush=True)
d = np.load(os.path.join(O, 'main_h005.npz')); t, X, Pe = d['t'], d['X'], d['Pe']; delay = 2 / 0.15 / 100
m = (t >= 0.5) & (t <= 4); Ge = -np.polyfit(t[m], np.log(Pe[m]), 1)[0]
pr(f"Gamma_eff/Gamma (P_e, [0,5; 4] tau) = {Ge:.4f}; bez detektora (754, N 800, ta sama metoda) nachylenie 1,016")
for name, key in [('P', 'trP'), ('F', 'trF')]:
    te = d[key] - delay; reg = d[key] >= 0
    for t0 in [0.5, 1.0]:
        u = te[reg & (te > t0)] - t0; nc = int((~reg).sum()); Tc = t[-1] - delay - t0
        mle = (u.sum() + nc * Tc) / len(u); ks = stats.kstest(u, lambda z: (1 - np.exp(-z / mle)) / (1 - np.exp(-Tc / mle)))
        pr(f"  {name}: ogon t_e > {t0} tau ({len(u)} biegow): srednia MLE nadwyzki {mle:.4f} tau (1/Gamma_eff {1/Ge:.4f}, stosunek {mle*Ge:.4f}); KS p {ks.pvalue:.3f}")
# absolute threshold: fraction of trajectories with x(t) >= x_th vs QM integral of rho over x >= x_th (P and F); grid of saved times
for xth in [3.0, 4.0]:
    pr(f"  prog bezwzgledny x_th = {xth}: P(x >= x_th) Bohm P / Bohm F / QM (int rho) przy t = 1 i 2 tau:")
    for tt, rk in [(1, 'rho1'), (2, 'rho2')]:
        k = np.argmin(abs(t - tt)); rho = d[rk]; q = rho[X >= xth].sum() / rho.sum()
        pr(f"    t {tt}: {np.mean(d['XP'][k] >= xth):.4f} / {np.mean(d['XF'][k] >= xth):.4f} / {q:.4f}")
    # registration time with absolute threshold (P: monotone motion -> first crossing); compare its distribution to exponential with detector delay fitted
    XP = d['XP'].astype(float); first = np.array([t[np.argmax(XP[:, i] >= xth)] if (XP[:, i] >= xth).any() else -1 for i in range(XP.shape[1])])
    u = first[first >= 0]; pr(f"    P: rejestracja x >= {xth}: {len(u)}/{len(first)}; kwantyle 10/50/90 % {np.percentile(u,10):.3f} / {np.percentile(u,50):.3f} / {np.percentile(u,90):.3f} tau")

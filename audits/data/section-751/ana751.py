# Audit 751: analysis -> out/ana751.txt
import os, glob, numpy as np
D = os.path.dirname(os.path.abspath(__file__)); O = os.path.join(D, 'out', 'runs'); pr = lambda *a: print(*a, flush=True)
alpha = 7.2973525643e-3; R1 = 128 * np.sqrt(6) / 243 / np.sqrt(3); w21 = 0.375
GF = 4 / 3 * alpha**3 * w21**3 * R1**2; au_t = 2.4188843265857e-17; Ha = 27.211386245988; hc = 1239.84198
mH = 1 / (1 + 1 / 1836.15267343); mPs = 0.5
ld = lambda f: np.loadtxt(f, comments='#')
def logistic(d):
    t, P = d[:, 0], d[:, 1]; m = (P > 0.05) & (P < 0.95); y = np.log((1 - P[m]) / P[m]); a, b = np.polyfit(t[m], y, 1); return a, b
def freq(d, P0, P1):
    t, P = d[:, 0], d[:, 1]; m = (P < P0) & (P > P1); ph = np.unwrap(np.angle(d[m, 4] + 1j * d[m, 5])); return abs(np.polyfit(t[m], ph, 1)[0])
pr(f"Regula Fermiego (mu = 1): Gamma = {GF:.6e} a.u.; tau wodoru (mu {mH:.6f}) = {1/(GF*mH)*au_t*1e9:.5f} ns; Ps = {1/(GF*mPs)*au_t*1e9:.4f} ns")
pr("=== H751a, H751b: s0 = 1e-2 ===")
Gs, fs, Ks = [], [], []
for K in [1e5, 1e4, 1e3]:
    d = ld(os.path.join(O, f'K{K:.0e}'.replace('e+0', 'e') + '_s1e-2.txt')); a, b = logistic(d)
    f_early, f_mid, f_late = freq(d, 1.0, 0.9), freq(d, 0.6, 0.4), freq(d, 0.1, 0.02)
    Gs.append(a / K); fs.append(f_mid); Ks.append(K)
    pr(f"  K {K:.0e}: Gamma_fit/K = {a/K:.6e} (/Fermi - 1 = {a/K/GF-1:+.2e}); s0 z wyrazu wolnego {1/(1+np.exp(-b)):.4e}; "
       f"czestosc (P 1-0,9 / 0,6-0,4 / 0,1-0,02): {f_early:.8f} / {f_mid:.8f} / {f_late:.8f} a.u. (/omega_21 - 1: {f_early/w21-1:+.2e} / {f_mid/w21-1:+.2e} / {f_late/w21-1:+.2e})")
c = np.polyfit(Ks, np.array(fs) / w21 - 1, 1)
pr(f"  przesuniecie czestosci (srodek) liniowo w K: nachylenie {c[0]:.3e} / K, wyraz wolny {c[1]:+.2e} -> przy K = 1: {c[0]+c[1]:+.2e}")
okA = abs(freq(ld(os.path.join(O, 'K1e3_s1e-2.txt')), 0.6, 0.4) / w21 - 1) < 1e-4 and abs(c[0] + c[1]) < 1e-8
okB = all(abs(g / GF - 1) < 0.01 for g in Gs)
for name, mu in [('wodor', mH), ('Ps (p-Ps i o-Ps)', mPs)]:
    E = w21 * mu * Ha; tau = 1 / (np.mean(Gs) * mu) * au_t * 1e9
    pr(f"  {name}: linia {E:.5f} eV = {hc/E:.4f} nm = {E/4.135667696e-15:.5e} Hz; tau z dopasowania {tau:.5f} ns")
pr("  pomiar wodoru 2p: 1,600(4) ns (Bickel i Goodman 1966)")
pr(f"H751a: {'TAK' if okA else 'NIE'};  H751b: {'TAK' if okB else 'NIE'}")
d0 = ld(os.path.join(O, 'K1e4_s0.txt')); dev = np.abs(d0[:, 1] - 1).max()
pr(f"H751c: s0 = 0, 20 tau: max |P_2p - 1| = {dev:.2e}, max P_1s = {d0[:,2].max():.2e} -> {'TAK' if dev < 1e-12 else 'NIE'}")
s0 = float(open(os.path.join(D, 'out', 's0.txt')).read()); ds = ld(os.path.join(O, 'K1e4_sorbit.txt'))
th = np.log((1 - s0) / s0) / (GF * 1e4); i = np.argmin(abs(ds[:, 1] - 0.5))
pr(f"H751d: s0 z orbity = {s0:.4f}; t_1/2 przewidywane {th*GF*1e4:.4f} tau, z biegu {ds[i,0]*GF*1e4:.4f} tau (P_2p(0) = {ds[0,1]:.4f})")
d = ld(os.path.join(O, 'K1e4_s1e-2.txt')); t = d[:, 0] * GF * 1e4; P = d[:, 1]
best = min(((np.abs(P - np.where(t > t0, np.exp(-(t - t0)), 1)).max(), t0) for t0 in np.linspace(0, 6, 601)))
pr(f"H751e: s0 1e-2, najlepszy wykladnik z opoznieniem t0 = {best[1]:.2f} tau: max |P - exp| = {best[0]:.3f} -> {'TAK (nie wykladnik)' if best[0] > 0.05 else 'NIE'}")
fs_ = sorted(glob.glob(os.path.join(O, 'ens', 's*.txt'))); arr = [ld(f) for f in fs_]; n = min(len(a) for a in arr)
tt = arr[0][:n, 0] * GF * 1e5; Pm = np.mean([a[:n, 1] for a in arr], axis=0)
for tq in [0.5, 1, 2, 3, 5]:
    j = np.argmin(abs(tt - tq)); pr(f"  zespol s0 ~ U(0,1), 100 biegow: t = {tq} tau: <P_2p> = {Pm[j]:.4f}, e^(-t) = {np.exp(-tq):.4f}, analitycznie {((1+np.expm1(tq))/np.expm1(tq)**2*np.log(np.exp(tq)) - 1/np.expm1(tq)):.4f}")

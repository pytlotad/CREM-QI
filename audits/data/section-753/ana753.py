# Audit 753: analysis -> out/ana753.txt  (surowe przebiegi: tar -xJf out/runs753.tar.xz -C out)
import os, numpy as np
from scipy import stats
D = os.path.dirname(os.path.abspath(__file__)); O = os.path.join(D, 'out', 'runs'); pr = lambda *a: print(*a, flush=True)
alpha = 7.2973525643e-3; R1 = 128 * np.sqrt(6) / 243 / np.sqrt(3); w21 = 0.375; GF = 4 / 3 * alpha**3 * w21**3 * R1**2
au_t = 2.4188843265857e-17; mH = 1 / (1 + 1 / 1836.15267343); mPs = 0.5; Ha = 27.211386245988; hc = 1239.84198
def load(f):
    te = {}; rows = []
    for line in open(os.path.join(O, f)):
        if line.startswith('# TE'): _, _, i, t = line.split(); te[int(i)] = float(t)
        elif not line.startswith('#'): rows.append([float(x) for x in line.split()])
    a = np.array(rows); runs = int(a[:, 0].max()) + 1; a = a.reshape(runs, -1, 5)
    return a, np.array([te[i] for i in range(runs)])
pr(f"Fermi (mu = 1) Gamma = {GF:.6e} a.u.; tau wodoru {1/(GF*mH)*au_t*1e9:.5f} ns, Ps {1/(GF*mPs)*au_t*1e9:.5f} ns (p-Ps i o-Ps, ta sama linia E1)")
for name, mu in [('wodor', mH), ('Ps (p-Ps i o-Ps)', mPs)]:
    E = w21 * mu * Ha; pr(f"  linia {name}: {E:.5f} eV = {hc/E:.4f} nm = {E/4.135667696e-15:.5e} Hz")
c, _ = load('C_2p_nozpf.txt'); pr(f"kontrola C (bez proznii, czysty 2p): max P_1s = {c[0,:,3].max():.2e} -> {'OK' if c[0,:,3].max() < 1e-10 else 'ZLE'}")
d, _ = load('D_s1e-2_nozpf.txt'); t, P = d[0, :, 1], d[0, :, 2]; m = (P > 0.05) & (P < 0.95)
a = np.polyfit(t[m], np.log((1 - P[m]) / P[m]), 1)[0]
pr(f"kontrola D (bez proznii, s0 1e-2): Gamma_fit/Gamma_Fermi - 1 = {a-1:+.2e} (751 K 1e5: -1,5e-5)")
A, te = load('A_2p_zpf.txt'); B, _ = load('B_1s_zpf.txt'); t = A[0, :, 1]; M = len(A)
PA = A[:, :, 2].mean(0); sA = A[:, :, 2].std(0) / np.sqrt(M); PB = B[:, :, 2].mean(0) + B[:, :, 4].mean(0)
pr(f"=== zespol A (start czysty 2p(m=+1), proznia; {M} biegow) ===")
mk = t <= 4 + 1e-9; dev = np.abs(PA[mk] - np.exp(-t[mk])); j = np.argmax(dev)
s = t <= 0.2 + 1e-9; slope = np.polyfit(t[s], PA[s], 1)[0]
pr(f"  max |<P_2p> - e^-t| na [0; 4 tau] = {dev.max():.4f} przy t = {t[mk][j]:.2f} (blad sredniej {sA[mk][j]:.4f}); nachylenie na [0; 0,2 tau] = {slope:+.4f} Gamma")
okA = dev.max() < 0.04 and abs(slope + 1) < 0.1; pr(f"H753a: {'TAK' if okA else 'NIE'}")
for tq in [0.1, 0.25, 0.5, 1, 2, 3, 4, 5]:
    i = np.argmin(abs(t - tq))
    pr(f"  t = {tq} tau: A <P_2p> {PA[i]:.4f} +- {sA[i]:.4f}, <P_1s> {A[:,i,3].mean():.4f}, <P_2s> {A[:,i,4].mean():.4f} (e^-t {np.exp(-tq):.4f});"
       f"  B <P_2p> {B[:,i,2].mean():.4f}, <P_1s> {B[:,i,3].mean():.4f}, <P_2s> {B[:,i,4].mean():.4f}")
pr(f"=== zespol B (start 1s, proznia) ===")
pr(f"  max <P_2p + P_2s> na [0; 5 tau] = {PB.max():.4f}; nachylenie na [0; 0,2 tau] = {np.polyfit(t[s], PB[s], 1)[0]:+.4f} Gamma (przewidywanie 3/2)")
pr(f"H753b: {'TAK' if PB.max() < 0.02 else 'NIE'}")
w = t >= 4 - 1e-9
pr(f"H753d: poziom ustalony (srednia [4; 5] tau): A <P_2p> {PA[w].mean():.4f}, <P_2s> {A[:,w,4].mean():.4f}; B <P_2p> {B[:,w,2].mean():.4f}, <P_2s> {B[:,w,4].mean():.4f}")
pr(f"  rozklad P_2p w pojedynczych biegach przy 5 tau: A kwantyle 10/50/90 % " + ' / '.join(f"{q:.3f}" for q in np.percentile(A[:, -1, 2], [10, 50, 90]))
   + "; B " + ' / '.join(f"{q:.3f}" for q in np.percentile(B[:, -1, 2], [10, 50, 90])))
pr(f"=== H753c: chwile emisji t_e (pierwsze P_2p < 1/2), zespol A ===")
u = te[te >= 0]; nc = int((te < 0).sum()); T = t[-1]
mle = (u.sum() + nc * T) / len(u)              # MLE sredniej wykladnika z cenzura w T
cv = u.std() / u.mean()
# KS wobec wykladnika uciętego do [0, T] (warunkowo na t_e < T)
cdf = lambda x: (1 - np.exp(-x / mle)) / (1 - np.exp(-T / mle))
ks = stats.kstest(u, cdf)
pr(f"  przejsc {len(u)}/{M}, cenzurowanych {nc}; srednia (nieenzur.) {u.mean():.3f} tau, MLE sredniej z cenzura {mle:.3f} tau; CV {cv:.3f}")
pr(f"  kwantyle t_e 10/25/50/75/90 %: " + ' / '.join(f"{q:.3f}" for q in np.percentile(u, [10, 25, 50, 75, 90])) + f"; min {u.min():.3f}")
pr(f"  KS wobec wykladnika (srednia MLE, ucietego w {T:.0f} tau): D = {ks.statistic:.4f}, p = {ks.pvalue:.2e}")
pr(f"H753c: {'TAK' if ks.pvalue > 0.01 and 0.8 <= cv <= 1.2 else 'NIE'}")
for name, mu in [('wodor', mH), ('Ps (p-Ps i o-Ps)', mPs)]:
    pr(f"  {name}: 1 tau = {1/(GF*mu)*au_t*1e9:.4f} ns; mediana t_e {np.median(u)/(GF*mu)*au_t*1e9:.3f} ns")

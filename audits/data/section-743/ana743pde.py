# Audit 743 (F): PDE results -- H743f controls, H743c (PDE vs M), H743d (response scan), H743e (back-reaction plateaus).
import sys, os, glob, filecmp, numpy as np
D = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, D)
from modes743 import torus, modes2d, delta, wrap, response
O = os.path.join(D, 'out', 'pde'); pr = lambda *a: print(*a, flush=True)
beta = 0.2
ld = lambda f: np.loadtxt(os.path.join(O, f), comments='#', ndmin=2)

pr("=== H743f: kontrole ===")
for a, b in [('ctrl446_eps0.txt', 'ctrl_eps0.txt'), ('ctrl446_eps3e-6_chi0.txt', 'eps3e-6_chi0.txt')]:
    pr(f"  field743 (sciezka 446) {a} wobec section-446/out/{b}: {'bajt w bajt' if filecmp.cmp(os.path.join(O, a), os.path.join(D, '..', 'section-446', 'out', b), shallow=False) else 'ROZNE'}")
c = ld('ctrl_eps0_Jr0.6_L1.1.txt')
pr(f"  eps = 0, 300 okresow, start (0,6; 1,1): max |J_r - J_r(1)| = {np.ptp(c[:, 2]):.2e}, max |L - L(1)| = {np.ptp(c[:, 3]):.2e} (kryterium 1e-5)")

pr("=== H743c: zrodlo przez jeden okres (wiersz 1 = okres emisji; Z_k = wiersze 2..4) ===")
for Jr, L, b in [(0.5, 1.0, 0.2), (0.6, 1.1, 0.2), (0.4, 0.9, 0.2), (0.5, 1.0, 0.0)]:
    d = ld(f'onep_Jr{Jr}_L{L}_b{b}.txt'); Z = d[:, 4] + 1j * d[:, 5]
    aF = [np.angle(Z[2] / Z[1]), np.angle(Z[3] / Z[2])]
    Zm = response(torus(Jr, L, b), modes2d(b, 12, 6), 0.0, onep=True)
    aM = [np.angle(Zm[1] / Zm[0]), np.angle(Zm[2] / Zm[1])]
    pr(f"  ({Jr}; {L}) beta {b}: F arg = {aF[0]:+.4f} {aF[1]:+.4f}; M arg = {aM[0]:+.4f} {aM[1]:+.4f}; |F - M| = "
       f"{abs(wrap(aF[0] - aM[0])):.3f} {abs(wrap(aF[1] - aM[1])):.3f} rad; |arg F| - pi = {abs(aF[0]) - np.pi:+.3f} {abs(aF[1]) - np.pi:+.3f}")

pr("=== H743d: skan odpowiedzi PDE (L = 1, beta 0,2, Gamma = 0,1/T_r(0,5; 1), srednia ostatnich 10 okresow) ===")
G = 0.1 / torus(0.5, 1.0, beta)['Tr']
Js = [0.20, 0.275, 0.35, 0.425, 0.50, 0.575, 0.65, 0.725, 0.80]
PF = np.array([ld(f'scan_Jr{J:.3f}'.rstrip('0') + '.txt' if False else f"scan_Jr{('%g' % J) if len(('%g' % J)) > 3 else '%.2f' % J}.txt")[-10:, 6].mean() for J in Js])
M2 = modes2d(beta, 12, 6)
PM = np.array([response(torus(J, 1.0, beta), M2, G)[0] for J in Js])
fine = np.round(np.arange(0.2, 0.8 + 1e-9, 0.005), 4); PMf = np.array([response(torus(J, 1.0, beta), M2, G)[0] for J in fine])
nF, nM = PF / PF.max(), PM / PM.max()
for J, a, b in zip(Js, nF, nM): pr(f"  J_r {J:.3f}: P_F/max = {a:.4f}  P_M/max = {b:.4f}")
pr(f"  argmax F {Js[nF.argmax()]}, M {Js[nM.argmax()]}; argmin F {Js[nF.argmin()]}, M {Js[nM.argmin()]}; korelacja log P: {np.corrcoef(np.log(PF), np.log(PM))[0, 1]:.3f}")
mx = [fine[i] for i in range(1, len(fine) - 1) if PMf[i] > PMf[i - 1] and PMf[i] >= PMf[i + 1]]
pr(f"  M (siatka 0,005) maksima: {mx}; P_M(0,5)/max_M = {response(torus(0.5, 1.0, beta), M2, G)[0] / PMf.max():.3f}")

pr("=== H743e: sila zwrotna (eps 3e-5, A 1, Gamma T_r0 = 0,1, 300 okresow) ===")
M2b = modes2d(beta, 12, 6)
for f in sorted(glob.glob(os.path.join(O, 'br_*.txt'))):
    d = np.loadtxt(f, comments='#', ndmin=2); n = len(d) // 10 * 10
    J = d[:n, 2].reshape(-1, 10).mean(1); L = d[:n, 3].reshape(-1, 10).mean(1)
    rJ, rL = np.ptp(J[-20:]), np.ptp(L[-20:])
    lock = rJ < 1e-3 and rL < 1e-3 and abs(J[-1] - 0.5) < 1e-3 and abs(L[-1] - 1) < 1e-3
    rate = np.hypot(np.diff(J), np.diff(L)); i = rate.argmin()
    Jm, Lm = 0.5 * (J[i] + J[i + 1]), 0.5 * (L[i] + L[i + 1]); T = torus(Jm, Lm, beta)
    near = min(M2b, key=lambda md: abs(wrap(delta(md, T))) if abs(md['m']) <= 4 else 9)
    pr(f"  {os.path.basename(f)}: okresy {len(d)}; koniec J_r {J[-1]:.4f} L {L[-1]:.4f}; zakres ost. 200: J {rJ:.2e} L {rL:.2e}; "
       f"zatrzasniecie: {'TAK' if lock else 'NIE'}; najwolniej przy blokach {i}-{i + 1} (J_r {Jm:.4f}, L {Lm:.4f}, tempo {rate[i]:.2e}/10 okr., "
       f"srednio {rate.mean():.2e}); tam min |delta_n| = {abs(wrap(delta(near, T))):.3f} rad (n_r {near['nr']}, m {near['m']})")

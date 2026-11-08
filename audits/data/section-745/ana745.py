# Audit 745: analysis of out/runs -> out/ana745.txt
import os, glob, numpy as np
D = os.path.dirname(os.path.abspath(__file__)); O = os.path.join(D, 'out', 'runs'); pr = lambda *a: print(*a, flush=True)
a = np.loadtxt(os.path.join(O, 'ctrl744_br_Jr0.5_L1.0_chi0.txt'))[:, :7]
b = np.loadtxt(os.path.join(D, '..', 'section-744', 'out', 'runs', 'br_Jr0.5_L1.0_chi0.txt'))
pr(f"kontrola: lock745 bez FORCE = lock744 (kolumny 1-7): max |roznica| = {np.abs(a - b).max():.1e}")
pr("=== sila eps grad|psi_mem|^2 ===")
anyLock = {}
for f in sorted(glob.glob(os.path.join(O, 'I_*.txt'))):
    d = np.loadtxt(f, comments='#', ndmin=2); n = len(d) // 10 * 10; name = os.path.basename(f)[:-4]
    J = d[:n, 2].reshape(-1, 10).mean(1); L = d[:n, 3].reshape(-1, 10).mean(1)
    rJ, rL = np.ptp(J[-20:]), np.ptp(L[-20:]); lock = rJ < 1e-3 and rL < 1e-3 and abs(J[-1] - .5) < 1e-3 and abs(L[-1] - 1) < 1e-3
    j0, l0 = d[0, 2], d[0, 3]; jm, lm = d[-100:, 2].mean(), d[-100:, 3].mean()
    d0 = np.hypot(j0 - .5, l0 - 1); dm = np.hypot(jm - .5, lm - 1)
    Eo = d[:, 7] + 0  # E_tot
    # E_orb from J_r, L: E = -1/(2 N^2), N = J_r + sqrt(L^2 + 0.4)
    Eorb = -0.5 / (d[:, 2] + np.sqrt(d[:, 3]**2 + 0.4))**2
    dEt, dEo = np.ptp(Eo[-200:]), np.ptp(Eorb[-200:])
    pr(f"  {name}: okresy {len(d)}; J_r " + ' '.join(f"{x:.4f}" for x in J[::5]) + " | L " + ' '.join(f"{x:.4f}" for x in L[::5]))
    pr(f"      zakres ost. 200: J_r {rJ:.2e}, L {rL:.2e}; zatrzasniecie {'TAK' if lock else 'NIE'}; srednia ost. 100 ({jm:.4f}; {lm:.4f}), "
       f"odl. od EBK {d0:.4f} -> {dm:.4f}; zakres E_tot {dEt:.2e} wobec E_orb {dEo:.2e} (stosunek {dEt/max(dEo,1e-30):.3f}); <|psi_mem|> {d[-100:,6].mean():.3f}")

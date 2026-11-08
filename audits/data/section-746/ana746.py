# Audit 746: drift-rate map around EBK -> out/ana746.txt
import os, glob, filecmp, numpy as np
D = os.path.dirname(os.path.abspath(__file__)); O = os.path.join(D, 'out', 'runs'); pr = lambda *a: print(*a, flush=True)
def meas(f):
    d = np.loadtxt(f, comments='#', ndmin=2); n = len(d) // 10 * 10
    J = d[:n, 2].reshape(-1, 10).mean(1); L = d[:n, 3].reshape(-1, 10).mean(1)
    R = np.hypot(J - J[0], L - L[0]).max()
    a = np.array([d[:20, 2].mean(), d[:20, 3].mean()]); b = np.array([d[-20:, 2].mean(), d[-20:, 3].mean()])
    return R, a, b, len(d)
c = filecmp.cmp(os.path.join(O, 'A_eps1e-3_Jr0.50_L1.00.txt'), os.path.join(D, '..', 'section-745', 'out', 'runs', 'I_k0_Jr0.5_L1.0_eps1e-3.txt'), shallow=False)
pr(f"kontrola: start EBK eps +1e-3 = bieg 745 (I_k0_Jr0.5_L1.0_eps1e-3): {'bajt w bajt' if c else 'ROZNE (naglowek/format?)'}")
Js = [f"{x:.2f}" for x in np.arange(0.2, 0.81, 0.1)]; Ls = [f"{x:.2f}" for x in np.arange(0.7, 1.31, 0.1)]
for E in ['1e-3', '-1e-3']:
    pr(f"=== siatka A, eps {E}: R (najwieksza odleglosc srednich 10-okresowych od pierwszej) ===")
    Rm = np.zeros((len(Js), len(Ls))); toward = {}
    for i, J in enumerate(Js):
        for k, L in enumerate(Ls):
            R, a, b, n = meas(os.path.join(O, f'A_eps{E}_Jr{J}_L{L}.txt')); Rm[i, k] = R
            e = np.array([0.5, 1.0]) - a; dn = np.linalg.norm(e)
            if 0.05 < dn < 0.15: toward[(J, L)] = (np.dot(b - a, e / dn), np.linalg.norm(b - a))
            if n < 300: pr(f"  uwaga: {J} {L} tylko {n} okresow")
    pr("  J_r \\ L " + ' '.join(f"{L:>7}" for L in Ls))
    for i, J in enumerate(Js): pr(f"  {J}    " + ' '.join(f"{x:7.4f}" for x in Rm[i]))
    iE, kE = Js.index('0.50'), Ls.index('1.00'); RE = Rm[iE, kE]; med = np.median(Rm); rank = int((Rm < RE).sum()) + 1
    imin = np.unravel_index(Rm.argmin(), Rm.shape)
    pr(f"  R(EBK) = {RE:.4f}; miejsce {rank}/49; minimum {Rm.min():.4f} przy ({Js[imin[0]]}; {Ls[imin[1]]}); mediana {med:.4f}; R(EBK)/mediana = {RE/med:.3f}")
    pos = sum(1 for v in toward.values() if v[0] > 0)
    pr(f"  sasiedzi EBK ({len(toward)}): skladowa D w strone EBK > 0: {pos}/{len(toward)}; " + '; '.join(f"({k[0]};{k[1]}) {v[0]:+.4f}/{v[1]:.4f}" for k, v in sorted(toward.items())))
    if E == '1e-3':
        pr(f"  H746a: {'TAK' if rank == 1 and RE < 0.2 * med else 'NIE'};  H746c: {'TAK' if pos >= 2 * len(toward) / 3 else 'NIE'}")
    else:
        pr(f"  H746d (R(EBK) nie najmniejsze przy eps < 0): {'TAK' if rank > 1 else 'NIE'}")
pr("=== siatka B (EBK (0,5; 2,0)), eps +1e-3 ===")
RB = {}
for f in sorted(glob.glob(os.path.join(O, 'B_*.txt'))):
    k = os.path.basename(f)[len('B_eps1e-3_'):-4]; RB[k] = meas(f)[0]; pr(f"  {k}: R = {RB[k]:.4f}")
e = RB.pop('Jr0.50_L2.00'); pr(f"  H746b: {'TAK' if e < min(RB.values()) else 'NIE'} (R(EBK) {e:.4f}, min sasiadow {min(RB.values()):.4f})")

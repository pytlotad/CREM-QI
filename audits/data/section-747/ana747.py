# Audit 747: rescaled-force drift map around EBK (0.5; 2) -> out/ana747.txt
import os, numpy as np
D = os.path.dirname(os.path.abspath(__file__)); O = os.path.join(D, 'out', 'runs'); pr = lambda *a: print(*a, flush=True)
EBK = np.array([0.5, 2.0])
def meas(f):
    d = np.loadtxt(f, comments='#', ndmin=2); n = len(d) // 10 * 10
    J = d[:n, 2].reshape(-1, 10).mean(1); L = d[:n, 3].reshape(-1, 10).mean(1)
    a = np.array([d[:20, 2].mean(), d[:20, 3].mean()]); b = np.array([d[-20:, 2].mean(), d[-20:, 3].mean()])
    return np.hypot(J - J[0], L - L[0]).max(), a, b, len(d)
def grid(prefix, Js, Ls, title):
    pr(f"=== {title} ===")
    Rm = np.zeros((len(Js), len(Ls))); tow = {}
    for i, J in enumerate(Js):
        for k, L in enumerate(Ls):
            R, a, b, n = meas(os.path.join(O, f'{prefix}Jr{J}_L{L}.txt')); Rm[i, k] = R
            e = EBK - a; dn = np.linalg.norm(e)
            if 0.05 < dn < 0.15: tow[(J, L)] = (np.dot(b - a, e / dn), np.linalg.norm(b - a))
            if n < 300: pr(f"  uwaga: ({J}; {L}) tylko {n} okresow")
    pr("  J_r \\ L " + ' '.join(f"{L:>7}" for L in Ls))
    for i, J in enumerate(Js): pr(f"  {J}    " + ' '.join(f"{x:7.4f}" for x in Rm[i]))
    RE = Rm[Js.index('0.50'), Ls.index('2.00')]; rank = int((Rm < RE).sum()) + 1; med = np.median(Rm)
    im = np.unravel_index(Rm.argmin(), Rm.shape)
    pr(f"  R(EBK) = {RE:.4f}; miejsce {rank}/{Rm.size}; minimum {Rm.min():.4f} przy ({Js[im[0]]}; {Ls[im[1]]}); mediana {med:.4f}; R(EBK)/mediana {RE/med:.3f}")
    pos = sum(1 for v in tow.values() if v[0] > 0)
    pr(f"  sasiedzi ({len(tow)}): w strone EBK {pos}/{len(tow)}; " + '; '.join(f"({k[0]};{k[1]}) {v[0]:+.4f}/{v[1]:.4f}" for k, v in sorted(tow.items())))
    return rank, RE, med, pos, len(tow)
Js = [f"{x:.2f}" for x in np.arange(0.2, 0.81, 0.1)]; Ls = [f"{x:.2f}" for x in np.arange(1.7, 2.31, 0.1)]
r, RE, med, pos, nt = grid('C_', Js, Ls, 'siatka C, eps_B = 3,945e-4')
pr(f"H747a: {'TAK' if r == 1 and RE < 0.2 * med else 'NIE'};  H747b: {'TAK' if pos >= 6 else 'NIE'} ({pos}/{nt})")
ok = True
for E in ['1.972e-4', '7.889e-4']:
    r, *_ = grid(f'D_eps{E}_', ['0.40', '0.50', '0.60'], ['1.90', '2.00', '2.10'], f'siatka D, eps {E}'); ok &= r == 1
pr(f"H747c: {'TAK' if ok else 'NIE'}")

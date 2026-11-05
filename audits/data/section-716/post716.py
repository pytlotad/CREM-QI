# Audit 716, POST HOC (after reading ana716): Ps on both sums (lrl = retarded, ins = instantaneous), the p-Ps spin pair's |s1 + s2| growth, and a rough
# Ps check of Delta psi = kappa (1 + 3 S_par/L) - 1 and apsidal 1 - 2 kappa S_par/L (units pi alpha^2/L^2), kappa = 3/4 (Breit-Pauli L.S of Ps, g = 2),
# S_par = (s1 + s2).L^ / 2 in hbar; spin-spin torques not included.
# usage: python3 post716.py <716 runs dir>
import sys, os, math
import numpy as np
from ana716 import parse, psi_track, ALPHA
R = sys.argv[1]; out = []
for D, lab in (('lrl', 'retarded'), ('ins', 'instantaneous')):
    for ph, ch in (('1', 'p-Ps'), ('2', 'o-Ps')):
        for s in ('0', '1'):
            d = parse(f'{R}/{D}/ps_A{ph}_{s}.err'); L = d['L']; u = math.pi*ALPHA**2/L**2
            a = psi_track(d, 1); b = psi_track(d, 2); aps = (np.mean(d['k']) - 2*math.pi)/u
            Sp = 0.5*(a[1] + b[1]); k = 0.75
            S = [np.linalg.norm(s1 + s2) for _, _, s1, s2 in d['lrl']]
            grow = np.polyfit(np.arange(len(S)), S, 1)[0]
            out.append(f'{lab:13s} {ch} seed {s}: Delta psi 1 / 2 {a[0]/u:+.4f} / {b[0]/u:+.4f} (mean {(a[0] + b[0])/(2*u):+.4f}) vs kappa(1 + 3 S_par/L) - 1'
                       f' {k*(1 + 3*Sp/L) - 1:+.4f}; apsidal {aps:+.4f} vs 1 - 2 kappa S_par/L {1 - 2*k*Sp/L:+.4f}; |s1 + s2| (unit spins) {S[0]:.2e} -> {S[-1]:.2e}'
                       f' ({grow:+.2e} per period)')
text = '\n'.join(out); print(text)
os.makedirs('out', exist_ok=True); open('out/post716.txt', 'w').write(text + '\n')

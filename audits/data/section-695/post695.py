# Audit 695 POST HOC (PHEN 2 only; the PHEN 1 trace projects a product state, see prereg695 uzupelnienie 2):
# time-mean of D (the offset a photon at a random instant would see) and the change of D over the photon's lag (13.3 / 84.3 ps).
# usage: python3 post695.py <695 runs dir with tr/>
import sys, glob, os
import numpy as np
R = sys.argv[1]


def parse(path):
    t, p = [], []; pf = None
    for line in open(path):
        s = line.split()
        if line.startswith('PAIR QM1S695 form'): pf = np.array([float(x) for x in s[17:29]])
        elif line.startswith('PAIR QM1S695 t'): t.append(float(s[3])); p.append([float(x) for x in s[5:17]])
    return np.array(t), np.array(p), pf


out = []
for f in sorted(glob.glob(f'{R}/tr/*_A2_*.err')):
    t, p, pf = parse(f); D = p - pf
    g = np.arange(t[0], t[-1], 1e-12); Dg = np.array([np.interp(g, t, D[:, k]) for k in range(12)]).T
    mean = Dg.mean(0); amp = Dg.std(0)
    lag = {L: np.sqrt(((Dg[L:] - Dg[:-L])**2).sum(1).mean()) for L in (13, 84)}
    out.append(f'{os.path.basename(f)[:-4]:26s} time-mean D max {np.abs(mean).max():.4f} (state {np.argmax(np.abs(mean))}) | oscillation std max {amp.max():.4f}'
               f' | rms |D(t + lag) - D(t)| 13 ps {lag[13]:.4f} 84 ps {lag[84]:.4f}')
text = '\n'.join(out); print(text); open('out/post695.txt', 'w').write(text + '\n')

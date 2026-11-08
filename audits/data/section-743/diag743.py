# Audit 743 POST HOC: |delta_n| of the nearest modes along the chi = pi/2 (start EBK) stop, blocks 12-17 -> out/diag743_przystanek.txt
import numpy as np, sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from modes743 import *
d = np.loadtxt(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'out/pde/br_Jr0.5_L1.0_chi1.5707963.txt'))
J = d[:, 2].reshape(-1, 10).mean(1); L = d[:, 3].reshape(-1, 10).mean(1); M = modes2d(0.2, 12, 6)
for b in [12, 13, 14, 15, 16, 17]:
    T = torus(J[b], L[b], 0.2); s = sorted(M, key=lambda m: abs(wrap(delta(m, T))))[:2]
    print(b, '%.4f %.4f' % (J[b], L[b]), [(m['nr'], m['m'], round(abs(wrap(delta(m, T))), 3)) for m in s])

# Audit 818, prediction before runs: nu_pred(scale) by the independent Plummer-softened integral (pred815) on the 815 starts (S025 8192).
import os, sys, math, tarfile
import numpy as np
D = os.path.dirname(os.path.abspath(__file__)); sys.path[:0] = [os.path.join(D, '..', 'section-815'), os.path.join(D, '..', 'section-814')]
import pred815 as P15
from pred814 import parse_start
x, w = np.polynomial.legendre.leggauss(5); u = 0.75 + 0.25*x[::-1]
NODES = dict(zip([f'u5{c}' for c in 'abcde'], zip(1/u, 0.25*w[::-1])))
nu_of = lambda Om: sum(w*L**2*Om[k] for k, (L, w) in NODES.items())/(2*math.pi)
G = {s: {} for s in (0.0, 0.01, 0.05)}
with tarfile.open(os.path.join(D, '..', 'section-815', 'runs', 'runs815.tar.xz')) as t:
    for k in NODES:
        m = [q for q in t.getnames() if q.endswith(f'S025/F1_ops_{k}_s8192_5.err')][0]
        tmp = os.path.join(D, 'out', f'_{k}.err'); open(tmp, 'wb').write(t.extractfile(m).read()); st = parse_start(tmp); os.remove(tmp)
        for s in G: G[s][k] = P15.omega(st, s)[0]
        print(k, ' '.join(f'{s}: {G[s][k]:.10e}' for s in G), flush=True)
n = {s: nu_of(G[s]) for s in G}
print(' '.join(f'nu_pred({s}) = {n[s]/1e6:.4f} MHz;' for s in G))
print(f'oslona 0,01 r*: {(n[0.01] - n[0.0])/1e6:+.4f} MHz; 0,05 r*: {(n[0.05] - n[0.0])/1e6:+.4f} MHz')

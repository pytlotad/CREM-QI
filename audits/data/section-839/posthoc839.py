# Audit 839 POST HOC: P0..P4 with the lab dipole set equal to the proper dipole (P839_NOBOOST=1) -- written after ana839.
import math, numpy as np, os
D = os.path.dirname(os.path.abspath(__file__))
def gauss(N):
    x, w = np.polynomial.legendre.leggauss(N); u = 0.75 + 0.25*x[::-1]; return 1/u, 0.25*w[::-1]
NODES = dict(zip([f'u5{c}' for c in 'abcde'], zip(*gauss(5))))
nu_of = lambda Om: sum(w*L**2*Om[k] for k, (L, w) in NODES.items())/(2*math.pi)
def om(f):
    return {int(l.split()[2]): -float(l.split()[4]) for l in open(f) if l.startswith('OMR')}
V = {}
for sk in ('ops', 'H0', 'H1'):
    for var, d in (('base', 'runs'), ('noboost', 'posthoc')):
        O = {nd: om(f'{D}/out/{d}/{sk}_{nd}_{var}.txt') for nd in NODES}
        V[sk, var] = {s: nu_of({nd: O[nd][s] for nd in NODES})/1e6 for s in (0, 1, 3, 7, 15)}
        print(sk, var, ' '.join(f'P{i} {V[sk, var][s]:+.4f}' for i, s in enumerate((0, 1, 3, 7, 15))), 'MHz')
for var in ('base', 'noboost'):
    print('wodor', var, ' '.join(f'P{i} parz {(V["H0", var][s] + V["H1", var][s])/2:+.4f} nieparz {(V["H1", var][s] - V["H0", var][s])/2:+.4f}'
                                 for i, s in enumerate((0, 1, 3, 7, 15))))

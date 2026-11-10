# Audit 772: 40 random starts (seed 772) -> starts772.txt  (id INIT)
import numpy as np, os
rng = np.random.default_rng(772); M = [(nr, m) for nr in range(3) for m in range(-3, 4)]
with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'starts772.txt'), 'w') as fo:
    for i in range(40):
        k = [2, 3, 4][i % 3]; sel = rng.choice(len(M), size=k, replace=False); w = rng.dirichlet(np.ones(k)); ph = rng.uniform(0, 2 * np.pi, k)
        fo.write(f"{i + 1:02d} " + ';'.join(f"{M[j][0]}:{M[j][1]}:{x:.6f}:{p:.6f}" for j, x, p in zip(sel, w, ph)) + "\n")

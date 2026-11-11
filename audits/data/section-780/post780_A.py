# Audit 780, POST HOC (nie zmienia werdyktow): przebieg czasowy ramienia A -- odleglosc od torusa startowego i w_max co 500 okresow -> out/post780_A.txt
import os, numpy as np
D = os.path.dirname(os.path.abspath(__file__))
for t in ("0.05_2_0", "0.2_1_0", "0.2_0_2", "0.05_0_2", "0.2_2_0"):
    b, nr, m = t.split('_'); nr = int(nr); m = int(m)
    for s in (2000, 4000, 8000):
        q = np.loadtxt(os.path.join(D, 'out', 'runs', f"A_b{t}_s{s}.txt"), comments='#')
        d = np.hypot(q[:, 2] - nr - 0.5, q[:, 3] - m)
        print(f"beta {b} {(nr, m)} STEPS {s}: d/w co 500: " + ' '.join(f"{d[k - 1]:.3f}/{q[k - 1, 5]:.2f}" for k in (1, 500, 1000, 1500, 2000, 2500, 3000, 3500, 4000)))

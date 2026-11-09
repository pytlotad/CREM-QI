# Audit 763: table rho3(n, l) = <r^-3> (pair units) from the EBK wave at order hbar^2 (../section-737, k = 1, "in"), n = 2..8, l = 1..n-1.
import os, sys, io, contextlib, numpy as np
D = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(D, '..', 'section-737'))
with contextlib.redirect_stdout(io.StringIO()): import ana737 as w
rows = []
for n in range(2, 9):
    r = w.grid(n)
    for l in range(1, n):
        u, ph = w.ebk(n, l, r, 1, (n, 1), False); v = np.trapezoid(u * u / r**3, r); q = 1 / (n**3 * l * (l + 0.5) * (l + 1))
        rows.append((n, l, v, v / q - 1)); print(f"{n} {l} {v:.17e} {v/q-1:+.3e}", flush=True)
    w._cache.clear()
with open(os.path.join(D, 'out', 'rho3_table.txt'), 'w') as f:
    f.write("# n l rho3 (pair units, EBK wave order hbar^2, 737 k=1 in) rho3/QM-1\n")
    for n, l, v, d in rows: f.write(f"{n} {l} {v:.17e} {d:+.3e}\n")

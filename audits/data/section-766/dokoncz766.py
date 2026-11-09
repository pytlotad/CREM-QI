# Audit 766: finish state 8d in a fresh process (runs 1-3 of ana766.py were OOM-killed, see prereg766.txt addendum); same ana766.run, same order of cases.
import os, sys, time, resource
D0 = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, D0)
import ana766 as a
pr = lambda s: print(s, flush=True); t0 = time.time(); run = a.run; n, l = 8, 2; T = lambda: f"[{time.time()-t0:.0f} s, maxrss {resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1e6:.2f} GB]"
k1_1 = run(n, l, 1); K3_1 = run(n, l, 1, K3=True); k2_1 = run(n, l, 2)
pr(f"  8d fine 1: k1 {k1_1:+.4e}  |K3-k1| {abs(K3_1-k1_1):.1e}  k2 {k2_1:+.4e} {T()}")
k1_2 = run(n, l, 1, fine=2); K3_2 = run(n, l, 1, K3=True, fine=2); k2_2 = run(n, l, 2, fine=2)
pr(f"  8d fine 2: k1 {k1_2:+.4e}  |K3-k1| {abs(K3_2-k1_2):.1e}  k2 {k2_2:+.4e}  (k1 zmiana {k1_2-k1_1:+.1e}, k2 zmiana {k2_2-k2_1:+.1e}) {T()}")
a6 = run(n, l, 2, rmin=1e-6); a4 = run(n, l, 2, rmin=1e-4)
pr(f"  8d rmin 1e-6: k2 {a6:+.4e} ({a6-k2_1:+.1e});  rmin 1e-4: k2 {a4:+.4e} ({a4-k2_1:+.1e}) {T()}")
g1 = run(n, l, 1, glob=True); g2 = run(n, l, 2, glob=True)
pr(f"  8d gl: k1 {g1:+.4e} ({g1-k1_1:+.1e})  k2 {g2:+.4e} (gl - in {g2-k2_1:+.1e}) {T()}")

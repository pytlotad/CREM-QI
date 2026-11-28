# Audit 875 POST HOC (after out/ana875.txt): H875r used s0 at the start n, but in collapsing runs J falls and the rate grows as 1/J^2;
# here the local rate over the first 10 % of cycles (min 5) is compared with -4 pi / (3 c^3 J^2) at the mean J of that window,
# and the lock offset J_k - n of the surviving runs with the 872 memory offset -0.024 (Gamma T)^2 / n.  -> out/ana875b.txt
import os, glob, math, importlib.util, contextlib, io, numpy as np
D = os.path.dirname(os.path.abspath(__file__))
sp_ = importlib.util.spec_from_file_location('a875', os.path.join(D, 'ana875.py')); A = importlib.util.module_from_spec(sp_)
with contextlib.redirect_stdout(io.StringIO()): sp_.loader.exec_module(A)
pr = A.pr
pr("=== POST HOC: tempo lokalne (pierwsze 10 % cykli) wobec -4 pi / (3 c^3 J^2) ===")
for k in sorted(A.R):
    r = A.R[k]
    if not k.startswith('K_'): continue
    hdr, q, stop, fall = A.load(os.path.join(A.O, k + '.txt')); m = max(len(q) // 10, 5)
    s = np.polyfit(q[:m, 0], q[:m, 2], 1)[0]; Jm = q[:m, 2].mean(); sl = -4 * math.pi / (3 * r['c'] ** 3 * Jm ** 2)
    pr(f"  {k:11s} {r['sp']:2s} n {r['n']} c {r['c']:8.3f}: J sr. {Jm:.4f}, s {s:+.4e}, s_lok {sl:+.4e}, s / s_lok {s / sl:.4f}")
pr("=== POST HOC: przesuniecie blokady J_k - n wobec pamieci 872 (-0,024 (Gamma T)^2 / n) ===")
for k in sorted(A.R):
    r = A.R[k]
    if r['cls'] != 'PRZETRWA': continue
    mem = -0.024 * (1e-3 * r['T']) ** 2 / r['n']
    pr(f"  {k:11s} {r['sp']:2s} {r['ch']:9s} n {r['n']} c {r['c']:8.3f}: J_k - n {r['J'] - r['n']:+.3e}, pamiec {mem:+.3e}, reszta {r['J'] - r['n'] - mem:+.3e}, s0 {r['s0']:+.2e}")

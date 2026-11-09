# Audit 767: one state (n l) at fine = 2 in its own process -> out/stany/n_l.txt (json).  Uses ana766.run unchanged.
import os, sys, json, time, resource, io, contextlib
D0 = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(D0, '..', 'section-766'))
with contextlib.redirect_stdout(io.StringIO()): import ana766 as a
n, l = int(sys.argv[1]), int(sys.argv[2]); t0 = time.time()
k1 = a.run(n, l, 1, fine=2); K3 = a.run(n, l, 1, N=80, W=0.06, K3=True, fine=2)
k2 = {f"{W}_{N}": a.run(n, l, 2, N=N, W=W, fine=2) for W, N in ((0.06, 80), (0.06, 160), (0.10, 80), (0.10, 160))}
res = dict(n=n, l=l, k1=k1, K3=K3, k2=k2, seconds=time.time() - t0, maxrss_GB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6)
p = os.path.join(D0, 'out', 'stany', f'{n}_{l}.txt'); json.dump(res, open(p + '.tmp', 'w')); os.replace(p + '.tmp', p)
print(f"{n} {l} k1 {k1:+.4e} |K3-k1| {abs(K3-k1):.1e} k2 {' '.join(f'{v:+.2e}' for v in k2.values())} [{time.time()-t0:.0f} s, {res['maxrss_GB']:.2f} GB]", flush=True)

# Audit 787 POST HOC (after scoring): t1 identity P vs D (same first threshold), photon-2 overlap with photon-1 window, M at photon 1 by quarter, H tau bias P vs D -> out/post787.txt
import os, numpy as np
D = os.path.dirname(os.path.abspath(__file__)); O = os.path.join(D, 'out', 'runs')
def ph(fn):
    r = []
    for l in open(fn):
        if l.startswith('PHASE787 photon'):
            p = l.split(); d = {p[i]: p[i + 1] for i in range(2, len(p) - 1, 2)}; r.append({k: float(v) for k, v in d.items()})
    return r
def period(fn):
    for l in open(fn):
        if l.startswith('PERIOD'): return float(l.split()[1])
for s in ('H', 'pPs', 'oPs'):
    same = sum(ph(f"{O}/P/{s}_{i}.err")[0]['t'] == ph(f"{O}/D/{s}_{i}.err")[0]['t'] for i in range(30))
    T3 = period(f"{O}/P/{s}_0.out"); T2 = T3 * 8 / 27
    P = [ph(f"{O}/P/{s}_{i}.err") for i in range(100)]; dt = np.array([p[1]['t'] - p[0]['t'] for p in P]); M = np.array([p[0]['M'] for p in P])
    q = np.histogram(M, bins=4, range=(0, 2 * np.pi))[0]
    dtP30 = dt[:30].mean(); dtD = np.mean([(lambda p: p[1]['t'] - p[0]['t'])(ph(f"{O}/D/{s}_{i}.err")) for i in range(30)])
    print(f"{s}: t1 P = D (ziarna 0-29) {same}/30; okres n = 3 {T3:.4e} s, n = 2 {T2:.4e} s; Dt < 1 obieg n = 2: {np.mean(dt < T2):.2f}; "
          f"srednie Dt/T2 P {dt.mean() / T2:.3f}; M fotonu 1 w cwiartkach {q.tolist()}; srednie Dt ziarna 0-29: P {dtP30 * 8e5 * 1e9:.4f} ns, D {dtD * 8e5 * 1e9:.4f} ns")

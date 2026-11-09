# Audit 765 (offline): preregistered replication of the 764 window method (k = 2) for n = 2..8, l = 1..n-1 -> out/ana765.txt, out/rho3_k2_table.txt
import os, sys, time, io, contextlib, numpy as np
D0 = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(D0, '..', 'section-764'))
with contextlib.redirect_stdout(io.StringIO()): import posthoc764 as ph
a = ph.a; w = ph.w
out = []; pr = lambda s: (print(s, flush=True), out.append(s))
k1 = {}
for L in open(os.path.join(D0, '..', 'section-763', 'out', 'rho3_table.txt')):
    if L.startswith('#'): continue
    n, l, v, _ = L.split(); k1[(int(n), int(l))] = float(v)
COMB = [(0.06, 80), (0.06, 160), (0.10, 80), (0.10, 160)]
def fsplit(r3, l, g=2 * (1 + a.AE)): return 0.5 * a.ALPHA**2 * r3 * (l + 0.5) * a.MU**3 * ((g - 1) + g / a.MP_ME) * a.HZ / 1e6
SOM = {(2, 1): 10968.704747, (3, 1): 3249.986592, (3, 2): 1083.328864, (4, 1): 1371.088093, (4, 2): 457.029364, (4, 3): 228.514682}
if __name__ == '__main__':
    t0 = time.time(); okA = okB = okC = True; res = {}
    pr('Audit 765: window method (764 POST HOC) preregistered for n = 2..8, l = 1..n-1; k = 2 mean of (W, N) = ' + str(COMB))
    for n in range(2, 9):
        r = w.grid(n)
        for l in range(1, n):
            q = 1 / (n**3 * l * (l + 0.5) * (l + 1)); e1 = k1[(n, l)] / q - 1
            K3 = ph.run(n, l, r, 1, 80, 0.06, K3=True)[0]; dK = abs(K3 - e1)
            vals = []; rv = []
            for W, N in COMB:
                R = ph.run(n, l, r, 2, N, W); vals.append(R[0]); rv.append(R[4])
            m = float(np.mean(vals)); sp = max(vals) - min(vals); rho = float(np.mean(rv))
            okA &= sp < 2e-6; okB &= dK < 1e-8; okC &= abs(m) < 5e-6 and abs(m) < abs(e1); res[(n, l)] = (rho, sp * q, m, q)
            pr(f"  {n}{'spdfghik'[l]}: k1 {e1:+.3e}  K3-k1 {dK:.1e}  k2 {' '.join(f'{v:+.2e}' for v in vals)}  srednia {m:+.3e} rozrzut {sp:.1e}"
               f"  {'ok' if (sp < 2e-6 and dK < 1e-8 and abs(m) < 5e-6 and abs(m) < abs(e1)) else 'NIE'}  [{time.time()-t0:.0f} s]")
    pr(f"H765a (rozrzut < 2e-6, 28/28): {'TAK' if okA else 'NIE'};  H765b (|K3 - k1| < 1e-8): {'TAK' if okB else 'NIE'};  "
       f"H765c (|k2| < 5e-6 i < |k1|): {'TAK' if okC else 'NIE'}")
    pr('H765d: rozszczepienia wodoru E(j=l+1/2) - E(j=l-1/2) [MHz]: k2 (rozrzut) / Sommerfeld silnika / <r^-3> QM')
    for (n, l), s in SOM.items():
        rho, spr, m, q = res[(n, l)]; f = fsplit(rho, l); fq = fsplit(q, l)
        pr(f"  {n}{'spdf'[l]}: {f:.3f} (+-{fsplit(spr, l):.3f}) / {s:.3f} ({f-s:+.3f}) / {fq:.3f} ({f-fq:+.3f})")
    f = fsplit(res[(2, 1)][0], 1)
    pr(f"  2P wobec pomiaru: {a.un(f)}; - pomiar {f-a.MEAS:+.3f} MHz ({f/a.MEAS-1:+.2e})")
    pr(f"  o-Ps 2^3P stala LS: {a.xi(res[(2,1)][0]):.3f} MHz wobec {a.XI_PS} (teoria QED) ({a.xi(res[(2,1)][0])/a.XI_PS-1:+.2e}); p-Ps 2^1P_1: bez LS")
    with open(os.path.join(D0, 'out', 'rho3_k2_table.txt'), 'w') as fo:
        fo.write("# n l rho3 (pair units, EBK wave order hbar^4, 764 window method, mean of 4 (W,N)) spread rho3/QM-1\n")
        for (n, l), (rho, spr, m, q) in res.items(): fo.write(f"{n} {l} {rho:.17e} {spr:.2e} {m:+.3e}\n")
    pr(f"runtime {time.time()-t0:.0f} s")
    open(os.path.join(D0, 'out', 'ana765.txt'), 'w').write('\n'.join(out) + '\n')

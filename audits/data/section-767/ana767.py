# Audit 767: verdicts and the candidate table from out/stany/ -> out/ana767.txt, out/rho3_k2_fine2_table.txt
import os, sys, json, io, contextlib, numpy as np
D0 = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(D0, '..', 'section-765'))
with contextlib.redirect_stdout(io.StringIO()): import ana765 as b
a = b.a; out = []; pr = lambda s: (print(s, flush=True), out.append(s))
f1 = {}
for L in open(os.path.join(D0, '..', 'section-765', 'out', 'rho3_k2_table.txt')):
    if L[0] != '#': p = L.split(); f1[(int(p[0]), int(p[1]))] = float(p[4])
okA = okB = okC = True; tab = {}; ND = 0; maxd = 0.0
pr('Audit 767: table k = 2 at fine = 2; n = 2..8, l = 1..n-1')
for n in range(2, 9):
    for l in range(1, n):
        R = json.load(open(os.path.join(D0, 'out', 'stany', f'{n}_{l}.txt'))); v = list(R['k2'].values()); m = float(np.mean(v)); sp = max(v) - min(v); dK = abs(R['K3'] - R['k1'])
        q = 1 / (n**3 * l * (l + 0.5) * (l + 1)); d = m - f1[(n, l)]; u = max(sp, abs(d)); A = dK < 1e-8; B = sp < 2e-6; C = abs(m) < 5e-6 and abs(m) < abs(R['k1'])
        okA &= A; okB &= B; okC &= C; ND += not (A and B and C); maxd = max(maxd, abs(d)); tab[(n, l)] = (q * (1 + m), u, m, q)
        pr(f"  {n}{'spdfghik'[l]}: k1 {R['k1']:+.3e}  K3 {dK:.1e}  k2 srednia {m:+.3e}  rozrzut {sp:.1e}  fine2-fine1 {d:+.1e}  {'ok' if A and B and C else 'NIE'}  [{R['seconds']:.0f} s, {R['maxrss_GB']:.2f} GB]")
pr(f"H767a (K3 < 1e-8): {'TAK' if okA else 'NIE'};  H767b (rozrzut < 2e-6): {'TAK' if okB else 'NIE'};  H767c (|k2| < 5e-6 i < |k1|): {'TAK' if okC else 'NIE'};  stanow z niespelnionym: {ND}")
pr(f"H767d (raport): max |fine 2 - fine 1| = {maxd:.1e}")
pr('H767e: rozszczepienia wodoru [MHz]: k2 (niepewnosc) / Sommerfeld silnika / <r^-3> QM')
for (n, l), s in b.SOM.items():
    rho, u, m, q = tab[(n, l)]; f = b.fsplit(rho, l); fq = b.fsplit(q, l)
    pr(f"  {n}{'spdf'[l]}: {f:.3f} (+-{b.fsplit(u * q, l):.3f}) / {s:.3f} ({f-s:+.3f}) / {fq:.3f} ({f-fq:+.3f})")
rho = tab[(2, 1)][0]; f = b.fsplit(rho, 1)
pr(f"  2P wobec pomiaru: {a.un(f)}; - pomiar {f-a.MEAS:+.3f} MHz ({f/a.MEAS-1:+.2e})")
pr(f"  o-Ps 2^3P stala LS: {a.xi(rho):.3f} MHz wobec {a.XI_PS} (teoria QED) ({a.xi(rho)/a.XI_PS-1:+.2e}); p-Ps 2^1P_1: bez LS")
with open(os.path.join(D0, 'out', 'rho3_k2_fine2_table.txt'), 'w') as fo:
    fo.write("# n l rho3 (pair units, EBK wave order hbar^4, window method, fine = 2, mean of 4 (W,N)) uncertainty(relative) rho3/QM-1\n")
    for (n, l), (rh, u, m, q) in tab.items(): fo.write(f"{n} {l} {rh:.17e} {u:.2e} {m:+.3e}\n")
open(os.path.join(D0, 'out', 'ana767.txt'), 'w').write('\n'.join(out) + '\n')

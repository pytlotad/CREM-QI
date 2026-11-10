# Audit 770: kinematic condition (first-order guidance) -- verdicts H770a-e -> out/ana770.txt
import os, numpy as np
D = os.path.dirname(os.path.abspath(__file__)); O = os.path.join(D, 'out', 'runs'); pr = lambda *a: print(*a, flush=True)
FS = [f"{x:.1f}" for x in np.arange(0, 1.01, 0.1)]; TOR = {'A': ((0.5, 1.0), (1.5, 1.0)), 'B': ((0.5, 1.0), (0.5, 2.0))}
def load(a, p, P0, f):
    q = np.loadtxt(os.path.join(O, f"G_{a}_{p}_P{P0}_f{f}.txt"), comments='#', ndmin=2); assert len(q) == 300, (a, p, P0, f, len(q)); return q
def nearest(j, l):
    best = (np.inf, None)
    for m in range(-6, 7):
        nr = max(0, round(j - 0.5)); d = np.hypot(j - nr - 0.5, l - m)
        if d < best[0]: best = (d, (nr + 0.5, float(m)))
    return best
def end(a, p, P0, f):
    q = load(a, p, P0, f); j, l = q[-20:, 2].mean(), q[-20:, 3].mean(); d, t = nearest(j, l)
    dd = np.array([nearest(x, y)[0] for x, y in q[:, 2:4]]); bad = np.where(dd >= 0.02)[0]
    tl = (bad[-1] + 2 if len(bad) else 1) if d < 0.02 else None
    return dict(j=j, l=l, d=d, tor=t, tl=tl, P=q[-20:, 4].mean(), wm=q[-1, 5], mm=(int(q[-1, 6]), int(q[-1, 7])))
R = {(a, p, P0, f): end(a, p, P0, f) for a in ['m2', 'z', 'a0'] for p in 'AB' for P0 in ['50', '500'] for f in FS}
for a, nm in [('m2', 'CAUST -pi/2, A 1 (glowne)'), ('z', 'CAUST 0, A 1'), ('a0', 'CAUST -pi/2, A 0')]:
    pr(f"=== ramie {nm} ===")
    for p in 'AB':
        for P0 in ['50', '500']:
            pr(f"  para {p}, P0 {P0}:   f   J_g      L_g      odl.   torus        t_lad  <P>        w_max (mod)")
            for f in FS:
                r = R[a, p, P0, f]
                pr(f"    {f}  {r['j']:.5f}  {r['l']:+.5f}  {r['d']:.4f}  {str(r['tor']):12s} {str(r['tl']):5s}  {r['P']:.3e}  {r['wm']:.4f} {r['mm']}")
okA = all(R['m2', p, P0, f]['d'] < 0.02 and R['m2', p, P0, f]['tor'] == TOR[p][0 if f == '0.0' else 1] for p in 'AB' for P0 in ['50', '500'] for f in ['0.0', '1.0'])
pr(f"H770a (czysty mod -> wlasny torus < 0,02): {'TAK' if okA else 'NIE'};  odleglosci: " + ', '.join(f"{p}{P0} f{f} {R['m2', p, P0, f]['d']:.4f}" for p in 'AB' for P0 in ['50', '500'] for f in ['0.0', '1.0']))
MID = FS[1:-1]; cnt = lambda a, p, P0, ok=lambda r, p: True: sum(R[a, p, P0, f]['d'] < 0.02 and ok(R[a, p, P0, f], p) for f in MID)
okB = all(cnt('m2', p, P0, lambda r, p: r['tor'] in TOR[p]) >= 8 for p in 'AB' for P0 in ['50', '500'])
pr(f"H770b (>= 8/9 ladowan na a lub b): {'TAK' if okB else 'NIE'};  " + ', '.join(f"{p}{P0} {cnt('m2', p, P0, lambda r, p: r['tor'] in TOR[p])}/9" for p in 'AB' for P0 in ['50', '500']))
c0 = {(p, P0): cnt('a0', p, P0) for p in 'AB' for P0 in ['50', '500']}; cz = {(p, P0): cnt('z', p, P0) for p in 'AB' for P0 in ['50', '500']}
okC = all(v == 0 for v in c0.values()) and all(v <= 2 for v in cz.values())
pr(f"H770c (A 0: 0/9; CAUST 0: <= 2/9): {'TAK' if okC else 'NIE'};  A 0: " + ', '.join(f"{k[0]}{k[1]} {v}/9" for k, v in c0.items()) + ";  CAUST 0: " + ', '.join(f"{k[0]}{k[1]} {v}/9" for k, v in cz.items()))
pr("H770d (raport): zwyciezca wzdluz f (a / b / - brak ladowania / inny torus), ramie glowne")
for p in 'AB':
    for P0 in ['50', '500']:
        s = ''.join('-' if R['m2', p, P0, f]['d'] >= 0.02 else ('a' if R['m2', p, P0, f]['tor'] == TOR[p][0] else ('b' if R['m2', p, P0, f]['tor'] == TOR[p][1] else 'x')) for f in FS)
        ch = sum(1 for i in range(len(s) - 1) if s[i] in 'ab' and s[i + 1] in 'ab' and s[i] != s[i + 1])
        tl = [R['m2', p, P0, f]['tl'] for f in MID if R['m2', p, P0, f]['tl']]
        pr(f"  {p} P0 {P0}: f 0,0..1,0 -> {s}  zmiany zwyciezcy {ch}  czas ladowania (okresy) {('mediana ' + str(int(np.median(tl)))) if tl else 'brak'}")
pr("H770e (raport): czysty mod, ramie glowne: przesuniecie (J_g - J_EBK, L_g - L_EBK), <P>, w_max")
for p in 'AB':
    for P0 in ['50', '500']:
        for f, k in [('0.0', 0), ('1.0', 1)]:
            r = R['m2', p, P0, f]; t = TOR[p][k]; pr(f"  {p} P0 {P0} f {f}: dJ {r['j'] - t[0]:+.5f}  dL {r['l'] - t[1]:+.5f}  <P> {r['P']:.3e}  w_max {r['wm']:.4f} {r['mm']}")

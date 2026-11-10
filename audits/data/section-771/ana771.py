# Audit 771: kinematic condition at narrower resonance (GT 0.01, 0.001) -- verdicts H771a-e -> out/ana771.txt
import os, numpy as np
D = os.path.dirname(os.path.abspath(__file__)); O = os.path.join(D, 'out', 'runs'); pr = lambda *a: print(*a, flush=True)
FS = ['0.0', '0.1', '0.3', '0.5', '0.7', '0.9', '1.0']; MID = FS[1:-1]; TOR = {'A': ((0.5, 1.0), (1.5, 1.0)), 'B': ((0.5, 1.0), (0.5, 2.0))}
GS = {'0.01': (1000, '50000'), '0.001': (3000, '5000000')}
def nearest(j, l):
    best = (np.inf, None)
    for m in range(-6, 7):
        nr = max(0, round(j - 0.5)); d = np.hypot(j - nr - 0.5, l - m)
        if d < best[0]: best = (d, (nr + 0.5, float(m)))
    return best
def load(g, a, p, P0, f):
    q = np.loadtxt(os.path.join(O, f"G_g{g}_{a}_{p}_P{P0}_f{f}.txt"), comments='#', ndmin=2); assert len(q) == GS[g][0], (g, a, p, P0, f, len(q)); return q
def end(g, a, p, P0, f):
    q = load(g, a, p, P0, f); j, l = q[-20:, 2].mean(), q[-20:, 3].mean(); d, t = nearest(j, l)
    return dict(q=q, j=j, l=l, d=d, tor=t, P=q[-20:, 4].mean(), wm=q[-1, 5], mm=(int(q[-1, 6]), int(q[-1, 7])))
for g, (per, PB) in GS.items():
    P0s = ['500', PB]; R = {(a, p, P0, f): end(g, a, p, P0, f) for a in ['m2', 'z'] for p in 'AB' for P0 in P0s for f in FS}
    pr(f"################ GT {g} ({per} okresow) ################")
    for a, nm in [('m2', 'CAUST -pi/2 (glowne)'), ('z', 'CAUST 0 (kontrola)')]:
        pr(f"=== ramie {nm} ===")
        for p in 'AB':
            for P0 in P0s:
                pr(f"  para {p}, P0 {P0}:   f   J_g      L_g      odl.   torus        <P>        w_max (mod)")
                for f in FS:
                    r = R[a, p, P0, f]; pr(f"    {f}  {r['j']:.5f}  {r['l']:+.5f}  {r['d']:.4f}  {str(r['tor']):12s} {r['P']:.3e}  {r['wm']:.4f} {r['mm']}")
    okA = all(R['m2', p, P0, f]['d'] < 0.02 and R['m2', p, P0, f]['tor'] == TOR[p][0 if f == '0.0' else 1] for p in 'AB' for P0 in P0s for f in ['0.0', '1.0'])
    pr(f"H771a GT {g} (czysty mod -> wlasny torus < 0,02): {'TAK' if okA else 'NIE'};  " + ', '.join(f"{p}{P0} f{f} {R['m2', p, P0, f]['d']:.4f}" for p in 'AB' for P0 in P0s for f in ['0.0', '1.0']))
    cnt = lambda a, p, P0, tor: sum(R[a, p, P0, f]['d'] < 0.02 and (not tor or R[a, p, P0, f]['tor'] in TOR[p]) for f in MID)
    okB = all(cnt('m2', p, P0, True) >= 4 for p in 'AB' for P0 in P0s)
    pr(f"H771b GT {g} (>= 4/5 ladowan na a lub b): {'TAK' if okB else 'NIE'};  " + ', '.join(f"{p}{P0} {cnt('m2', p, P0, True)}/5" for p in 'AB' for P0 in P0s))
    okC = all(cnt('z', p, P0, False) <= 1 for p in 'AB' for P0 in P0s)
    pr(f"H771c GT {g} (CAUST 0 <= 1/5): {'TAK' if okC else 'NIE'};  " + ', '.join(f"{p}{P0} {cnt('z', p, P0, False)}/5" for p in 'AB' for P0 in P0s))
    pr(f"H771d GT {g} (raport): czysty mod (0,1), ramie glowne: (okres: w_max, odleglosc od (0,5; 1))")
    for P0 in P0s:
        q = R['m2', 'A', P0, '0.0']['q']
        pr(f"  P0 {P0}: " + ';  '.join(f"{k}: {q[k - 1, 5]:.4f}, {np.hypot(q[k - 1, 2] - 0.5, q[k - 1, 3] - 1):.4f}" for k in [1, 10, 100, per]))
    for p in 'AB':
        for P0 in P0s:
            s = ''.join('-' if R['m2', p, P0, f]['d'] >= 0.02 else ('a' if R['m2', p, P0, f]['tor'] == TOR[p][0] else ('b' if R['m2', p, P0, f]['tor'] == TOR[p][1] else 'x')) for f in FS)
            pr(f"  zwyciezca {p} P0 {P0}: f {' '.join(FS)} -> {s}")
    pr(f"H771e GT {g} (raport): odleglosc koncowa czystego modu (ramie glowne): " + ', '.join(f"{p}{P0} f{f} {R['m2', p, P0, f]['d']:.4f}" for p in 'AB' for P0 in P0s for f in ['0.0', '1.0']) + "  (770, GT 0,1: 0,36-0,41)")

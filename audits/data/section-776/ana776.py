# Audit 776: does the drive spectrum on the frozen start torus (P_D) predict the winning mode?  Verdicts H776a-c -> out/ana776.txt
import os, numpy as np
from collections import Counter, defaultdict
D = os.path.dirname(os.path.abspath(__file__)); A = os.path.join(D, '..'); P = os.path.join(D, 'out', 'pred')
pr = lambda *a: print(*a, flush=True)
def starts(fn):
    S = {}
    for ln in open(fn):
        i, ini = ln.split(); S[i] = [(int(a), int(b), float(c)) for a, b, c, _ in (x.split(':') for x in ini.split(';'))]
    return S
S2, S3 = starts(os.path.join(A, 'section-772', 'starts772.txt')), starts(os.path.join(A, 'section-773', 'starts773.txt'))
def nearest(j, l):
    best = (np.inf, None)
    for m in range(-6, 7):
        nr = max(0, round(j - 0.5)); d = np.hypot(j - nr - 0.5, l - m)
        if d < best[0]: best = (d, (nr, m))
    return best
def pred(key, c, i):
    fn = os.path.join(P, f"{key}_{c}_{i}.txt"); h = open(fn).readline().split()
    J0, L0 = float(h[h.index('J0') + 1]), float(h[h.index('L0') + 1])
    q = [tuple(int(x) for x in ln.split()[:2]) for ln in open(fn) if not ln.startswith('#')]
    return q, (J0, L0)
#        arm      run file pattern                                   starts ids                         per    pred key
ARMS = [('772 G', 'section-772/out/runs/G_{}.txt', S2, sorted(S2), 3000, 'S772_b02'),
        ('772 W', 'section-772/out/runs/W_{}.txt', S2, sorted(S2)[:12], 10000, 'S772_b02'),
        ('773 G', 'section-773/out/runs/G_{}.txt', S3, [f"{k:02d}" for k in range(1, 41)], 12000, 'S773_b02'),
        ('773 W', 'section-773/out/runs/W_{}.txt', S2, ['05', '08', '14', '15', '36'], 20000, 'S772_b02'),
        ('774 G', 'section-774/out/runs/G_{}.txt', S3, [f"{k:02d}" for k in range(1, 41)], 20000, 'S773_b02'),
        ('775 B05', 'section-775/out/runs/B05_{}.txt', S3, [f"{k:02d}" for k in range(1, 21)], 20000, 'S773_b005'),
        ('775 B80', 'section-775/out/runs/B80_{}.txt', S3, [f"{k:02d}" for k in range(1, 21)], 20000, 'S773_b08')]
rows = []
for arm, pat, S, ids, per, key in ARMS:
    for i in ids:
        q = np.loadtxt(os.path.join(A, pat.format(i)), comments='#', ndmin=2); assert len(q) == per, (arm, i, len(q))
        e = q[-20:, 2:4].mean(0); d, tor = nearest(*e); wm = q[-1, 5]
        pd, (J0, L0) = pred(key, 'c1', i); p0, _ = pred(key, 'c0', i)
        rows.append(dict(arm=arm, i=i, key=key, one=wm >= 0.95 and d < 0.04, win=tor, wm=wm, d=d, PD=pd, PD0=p0,
                         B1=nearest(J0, L0)[1], B2=max(S[i], key=lambda x: x[2])[:2], start=[x[:2] for x in S[i]], J0=J0, L0=L0))
OK = [r for r in rows if r['one']]
B3 = {arm: Counter(r['win'] for r in OK if r['arm'] == arm).most_common(1)[0][0] for arm, *_ in ARMS}
for r in OK: r['B3'] = B3[r['arm']]
pr(f"Biegi: {len(rows)}; konce jednomodowe na EBK (mianownik): {len(OK)}; pozostale {len(rows) - len(OK)}")
pr("=== szczegoly (konce jednomodowe) ===")
pr("  ramie    id  zwyc.     P_D       P_D0      B1        B2        ranga  top5 P_D")
for r in OK:
    rk = r['PD'].index(r['win']) + 1
    pr(f"  {r['arm']:8s} {r['i']}  {str(r['win']):9s} {str(r['PD'][0]):9s} {str(r['PD0'][0]):9s} {str(r['B1']):9s} {str(r['B2']):9s} {rk:5d}  {' '.join(str(x) for x in r['PD'][:5])}")
hit = lambda k: sum((r[k][0] if k in ('PD', 'PD0') else r[k]) == r['win'] for r in OK)
n = len(OK); h = {k: hit(k) for k in ('PD', 'PD0', 'B1', 'B2', 'B3')}
pr("=== trafienia (wszystkie ramiona, na biegach) ===")
for k in ('PD', 'PD0', 'B1', 'B2', 'B3'): pr(f"  {k:4s} {h[k]}/{n} = {100 * h[k] / n:.1f} %")
best = max(h['B1'], h['B2'], h['B3'])
ok = h['PD'] >= 0.5 * n and 100 * (h['PD'] - best) / n >= 15
pr(f"H776a (P_D >= 50 % i >= najlepsza z B1-B3 + 15 pp): P_D {100 * h['PD'] / n:.1f} %, najlepsza linia bazowa {100 * best / n:.1f} % -> {'TAK' if ok else 'NIE'}")
pr(f"H776b (raport): P_D0 (CAUST 0) {h['PD0']}/{n} wobec P_D {h['PD']}/{n}")
pr("H776c (raport): na ramie -- P_D / P_D0 / B1 / B2 / B3 / top-3 P_D / mediana rangi")
for arm, *_ in ARMS:
    R = [r for r in OK if r['arm'] == arm]; m = len(R)
    if not m: pr(f"  {arm}: brak koncow jednomodowych"); continue
    f = lambda k: sum((r[k][0] if k in ('PD', 'PD0') else r[k]) == r['win'] for r in R)
    t3 = sum(r['win'] in r['PD'][:3] for r in R); rk = np.median([r['PD'].index(r['win']) + 1 for r in R])
    pr(f"  {arm:8s} n {m:2d}: {f('PD')} / {f('PD0')} / {f('B1')} / {f('B2')} / {f('B3')} / {t3} / {rk:.0f}")
t3 = sum(r['win'] in r['PD'][:3] for r in OK); pr(f"  razem top-3 P_D {t3}/{n}; mediana rangi {np.median([r['PD'].index(r['win']) + 1 for r in OK]):.0f}")
U = defaultdict(list)
for r in OK: U[(r['key'], r['i'])].append(r)
cons = {k: v for k, v in U.items() if len(set(r['win'] for r in v)) == 1}
pr(f"  unikalne (start, beta): {len(U)}; z jednym zwyciezca we wszystkich ramionach {len(cons)}; P_D trafia w {sum(v[0]['PD'][0] == v[0]['win'] for v in cons.values())}/{len(cons)}; "
   f"zwyciezca rozny miedzy ramionami (GT) w {len(U) - len(cons)}")
ch = [i for i in [f"{k:02d}" for k in range(1, 21)]]
W = {(r['arm'], r['i']): r for r in OK}
nb = sum(1 for i in ch if all((a, i) in W for a in ('775 B05', '774 G', '775 B80')))
same_pred = 0; both = 0; diff = 0
for i in ch:
    rr = [W.get((a, i)) for a in ('775 B05', '774 G', '775 B80')]
    if None in rr: continue
    if len(set(r['win'] for r in rr)) > 1:
        diff += 1; both += all(r['PD'][0] == r['win'] for r in rr)
pr(f"  775: starty z koncem jednomodowym przy trzech beta {nb}; zwyciezca zmienia sie z beta w {diff}; P_D trafia wszystkie trzy w {both}/{diff}")
pr("  pozostale konce (nie jednomodowe na EBK): " + ', '.join(f"{r['arm']} {r['i']} w {r['wm']:.2f} odl. {r['d']:.3f} P_D {r['PD'][0]}" for r in rows if not r['one']))

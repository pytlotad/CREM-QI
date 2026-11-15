# Audit 794: why 3p misses (793) -- relativistically matched torus start (CREM_PAIR_TORQUE_RELSTART), fixed-time window vs fit estimator,
# and hydrogen F = 0 (PHEN=1) next to F = 1 (PHEN=2).  All runs CREM_NO_DIPOLE_FORCE, 512 steps.  usage: python3 ana794.py
# 793 G references unpacked from ../section-793/runs/runs793.tar.xz into $R793 (default scratchpad copy).
import sys, os, math
import numpy as np
D = os.path.dirname(os.path.abspath(__file__)); O = os.path.join(D, 'out', 'runs')
sys.path.insert(0, os.path.join(D, '..', 'section-699')); import ana699
R793 = os.environ.get('R793', '/tmp/claude-1000/-home-teddy-Projekty-CREM-QI/4fa83d5c-8878-4cfd-898e-1a9b53a01ab2/scratchpad/r793x/runs')
HEV, C = 4.135667696e-15, 299792458.0; FS = 10969.045e6
pr = lambda *a: print(*a, flush=True)
units = lambda f: f'{f/1e6:.2f} MHz = {f*HEV*1e6:.4f} ueV = {C/f*100:.4f} cm'
TORI = [(3, 1), (3, 1.5), (3, 2), (3, 2.5), (3, 3), (2, 1), (2, 1.5), (2, 2)]
lab = lambda L: '%g' % L
ecc = lambda n, L: math.sqrt(max(0.0, 1 - L*L/(n*n)))
TYPES = ('proton', 'antiproton'); SEEDS = '456'
def omega(path):
    m = ana699.parse(path)['inc']['mechB0']
    return np.linalg.norm(m['s21'] - m['s20'])/(m['dt']*np.linalg.norm(np.cross(m['L0'], m['s20'])))/(2*math.pi), m
def trace(path):
    rows = []
    for ln in open(path):
        if not ln.startswith('PAIR TORQUE792'): continue
        w = ln.split(); g = lambda k: float(w[w.index(k) + 1]); v = lambda k: np.array([float(x) for x in w[w.index(k) + 1:w.index(k) + 4]])
        rows.append(dict(t=g('t'), E=g('E'), eps=g('eps'), L=v('L'), s2=v('s2')))
    return rows
def relstart(path):
    for ln in open(path):
        if ln.startswith('PAIR TORQUE794'):
            w = ln.split(); return {w[i]: float(w[i + 1]) for i in range(3, len(w) - 1, 2)}
def fit(path):   # least-squares slope of the accumulated cone angle over the 200 one-orbit pieces
    om, m = omega(path); rows = trace(path); s_prev, L_prev = m['s20'], m['L0']; t = [0.0]; a = [0.0]
    for r in rows:
        a.append(a[-1] + np.linalg.norm(r['s2'] - s_prev)/np.linalg.norm(np.cross(L_prev, s_prev))); t.append(r['t']); s_prev, L_prev = r['s2'], r['L']
    return np.polyfit(np.array(t), np.array(a), 1)[0]/(2*math.pi)
def rule(W, tag, okD, okM, out):
    for n, l, st in ((2, 1, '2p'), (3, 1, '3p'), (3, 2, '3d')):
        F = lambda L: L*W[(n, L)]
        od = F(l) - F(l + 1); q = od/W[(n, l + 0.5)]; qt = (l + 0.5)**2/(l*(l + 1)); okD[0] &= abs(q/qt - 1) <= 0.01
        rf = FS*16/(n**3*l*(l + 1)); out.setdefault(st, []).append(q/qt)
        pr(f'  {tag} {st}: Omega_rozn/2pi {units(od)} | {q/qt:.5f} x regula | wobec skalowania pomiaru 2P {rf/1e6:.2f} MHz -> {od/rf:.5f} ({(od - rf)/1e6:+.2f} MHz)')
        if st == '2p': okM[0] &= abs(od/FS - 1) <= 0.01; out.setdefault('2P', []).append(od)
def kn(W): return {k: W[k]*k[0]**3*k[1]**2 for k in W}
# --- H794k
okK = True
for n, L in TORI:
    a = ''.join(l for l in open(f'{O}/A/proton_n{n}L{lab(L)}_5.err') if not l.startswith('PAIR TORQUE792'))
    okK &= a == open(f'{R793}/G/proton_n{n}L{lab(L)}_5.err').read()
pr(f"H794k (A bez linii TORQUE792 = 793 G bajt w bajt, 8 plikow): {'TAK' if okK else 'NIE'}")
# --- H794r
pr('=== start na torusie (B wodor ziarno 5; A = stary start) ===')
okR = True
for n, L in TORI:
    p = f'{O}/B/proton_n{n}L{lab(L)}_5.err'; rs = relstart(p); tb = trace(p); ta = trace(f'{O}/A/proton_n{n}L{lab(L)}_5.err')
    d1 = tb[0]['E']/rs['ES'] - 1; okR &= abs(d1) <= 3e-4
    pr(f'  n {n} L {L:3} e {ecc(n, L):.3f}: stary start E0/E_S - 1 {rs["E0"]/rs["ES"] - 1:+.3e} | r1/r0 - 1 {rs["r1"]/rs["r0"] - 1:+.2e} bracket {int(rs["bracket"])}'
       f' | B: E(1)/E_S - 1 {d1:+.2e}, E(200)/E_S - 1 {tb[-1]["E"]/rs["ES"] - 1:+.2e} | A: E(1)/E_S - 1 {ta[0]["E"]/rs["ES"] - 1:+.2e}, E(200)/E_S - 1 {ta[-1]["E"]/rs["ES"] - 1:+.2e}')
pr(f"H794r (|E(1)/E_S - 1| <= 3e-4, 8 torusow): {'TAK' if okR else 'NIE'}")
# --- B, C rule
RES = {}
for cfg, F in (('B', 'F = 1'), ('C', 'F = 0')):
    okP = [True]; okD = [True]; okM = [True]; out = {}; ALL = {}
    pr(f'=== {cfg} ({F}, nowy start) ===')
    for t in TYPES:
        for s in SEEDS:
            W = {k: omega(f'{O}/{cfg}/{t}_n{k[0]}L{lab(k[1])}_{s}.err')[0] for k in TORI}
            for k in TORI: ALL.setdefault(k, []).append(W[k])
            K = kn(W); sp = max(K.values())/min(K.values()) - 1; okP[0] &= sp <= 0.01
            pr(f'  {t} ziarno {s}: rozrzut Omega n^3 L^2 {sp*100:.3f} % | ' + ', '.join(f'n{k[0]}L{lab(k[1])} {K[k]/1e6:.3f}' for k in TORI) + ' MHz')
            rule(W, f'{t} {s}', okD, okM, out)
    M = {k: np.mean(v) for k, v in ALL.items()}; KM = kn(M); ref = KM[(3, 3)]
    pr(f'  srednie 6 biegow: Omega n^3 L^2 wobec n3L3 -- ' + ', '.join(f'n{k[0]}L{lab(k[1])} (e {ecc(*k):.3f}) {KM[k]/ref - 1:+.4%}' for k in TORI))
    rule(M, 'srednie', [True], [True], {})
    RES[cfg] = dict(okP=okP[0], okD=okD[0], okM=okM[0], out=out, M=M, ALL=ALL, KM=KM)
    pr(f"  {cfg}: H794p rozrzut <= 1 %: {'TAK' if okP[0] else 'NIE'}; regula 2p {min(out['2p']):.5f}-{max(out['2p']):.5f}, 3p {min(out['3p']):.5f}-{max(out['3p']):.5f},"
       f" 3d {min(out['3d']):.5f}-{max(out['3d']):.5f}; 2P {min(out['2P'])/1e6:.2f}-{max(out['2P'])/1e6:.2f} MHz ({min(out['2P'])/FS:.4f}-{max(out['2P'])/FS:.4f} pomiaru)")
KB = RES['B']['KM']; okE = abs(KB[(3, 1)]/KB[(3, 3)] - 1) <= 3e-3
pr(f"H794e (B: n3L1 / n3L3 Omega n^3 L^2 w 0,3 %: {KB[(3, 1)]/KB[(3, 3)] - 1:+.4%}): {'TAK' if okE else 'NIE'}")
pr(f"H794d (B, F = 1: 2p, 3p, 3d w 1 %): {'TAK' if RES['B']['okD'] else 'NIE'}")
pr(f"H794m (B, F = 1: 2P w 1 % pomiaru): {'TAK' if RES['B']['okM'] else 'NIE'}")
pr(f"H794p (B i C: rozrzut <= 1 %): {'TAK' if RES['B']['okP'] and RES['C']['okP'] else 'NIE'}")
# --- H794w
pr('=== estymator dopasowania (wodor ziarno 5) wobec okna 699 ===')
okW = True
for cfg in ('A', 'B'):
    Wf = {}; Wo = {}
    for k in TORI:
        p = f'{O}/{cfg}/proton_n{k[0]}L{lab(k[1])}_5.err'; Wf[k] = fit(p); Wo[k] = omega(p)[0]
    Kf = kn(Wf); Ko = kn(Wo); ref_f = Kf[(3, 3)]; ref_o = Ko[(3, 3)]
    rest = [Kf[k] for k in TORI if k != (3, 1)]; sp = max(rest)/min(rest) - 1; okW &= sp <= 1e-3
    pr(f'  {cfg}: rozrzut dopasowania bez n3L1 {sp*100:.3f} %')
    for k in TORI: pr(f'    n{k[0]}L{lab(k[1])} (e {ecc(*k):.3f}): okno {Ko[k]/ref_o - 1:+.4%} | dopasowanie {Kf[k]/ref_f - 1:+.4%} | dopasowanie / okno - 1 {Wf[k]/Wo[k] - 1:+.2e}')
    rule(Wf, f'{cfg} dopasowanie', [True], [True], {})
pr(f"H794w (dopasowanie: rozrzut bez n3L1 <= 0,1 % w A i B): {'TAK' if okW else 'NIE'}")
# --- H794s, singlet rule
okS = True; mx = 0.0
for t in TYPES:
    for s in SEEDS:
        for k in TORI:
            q = RES['C']['ALL'][k][TYPES.index(t)*3 + SEEDS.index(s)]/RES['B']['ALL'][k][TYPES.index(t)*3 + SEEDS.index(s)] - 1
            mx = max(mx, abs(q)); okS &= abs(q) <= 3e-3
pr(f"H794s (|Omega(F = 0)/Omega(F = 1) - 1| <= 0,3 %, max {mx*100:.3f} %): {'TAK' if okS else 'NIE'}")
pr(f"H794d0 (C, F = 0: 2p, 3p, 3d w 1 %): {'TAK' if RES['C']['okD'] else 'NIE'}")
pr(f"H794m0 (C, F = 0: 2P w 1 % pomiaru): {'TAK' if RES['C']['okM'] else 'NIE'}")
# --- report
pr('=== RAPORT ===')
for cfg, F in (('B', 'F = 1'), ('C', 'F = 0')):
    o = RES[cfg]['out']['2P']
    pr(f'  {F}: 2P srednio {units(np.mean(o))} -> {np.mean(o)/FS:.5f} pomiaru; wodor {np.mean(o[:3])/1e6:.2f}, antywodor {np.mean(o[3:])/1e6:.2f} MHz')
c2 = 0.25*np.mean(RES['C']['out']['2P']) + 0.75*np.mean(RES['B']['out']['2P'])
pr(f'  centroid 1/4 F = 0 + 3/4 F = 1 (wyprowadzone): 2P {units(c2)} -> {c2/FS:.5f} pomiaru')
R21 = [x*9/8 for x in RES['B']['out']['2p'] + RES['C']['out']['2p']]
pr(f'  o-Ps: R (2p) {min(R21):.5f}-{max(R21):.5f}; 8/9 x R x 0,9998..0,9999 = {8/9*min(R21)*0.9998:.5f}-{8/9*max(R21)*0.9999:.5f} multipletu (Hagena i in. 1993); p-Ps bez czesci spinowej')
pr('Pomiar: 2P3/2 - 2P1/2 wodoru 10 969,045(15) MHz (Hagley i Pipkin, PRL 72, 1172 (1994); Lundeen i Pipkin, PRL 46, 232 (1981)).')

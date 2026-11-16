# Audit 799: spin precession per radial orbit (CREM_PAIR_TORQUE_PASS, apoapsis samples of the mechB0 run) instead of per unit time.  usage: python3 ana799.py
# 798 references unpacked from ../section-798/runs/runs798.tar.xz into $R798.
import sys, os, math
import numpy as np
D = os.path.dirname(os.path.abspath(__file__)); O = os.path.join(D, 'out', 'runs')
sys.path.insert(0, os.path.join(D, '..', 'section-699')); import ana699
R798 = os.environ.get('R798', '/tmp/claude-1000/-home-teddy-Projekty-CREM-QI/4fa83d5c-8878-4cfd-898e-1a9b53a01ab2/scratchpad/r798x/runs')
HEV, C = 4.135667696e-15, 299792458.0; FS = 10969.045e6
pr = lambda *a: print(*a, flush=True)
units = lambda f: f'{f/1e6:.2f} MHz = {f*HEV*1e6:.4f} ueV = {C/f*100:.4f} cm'
TORI = [(3, 1), (3, 1.5), (3, 2), (3, 2.5), (3, 3), (2, 1), (2, 1.5), (2, 2)]
TYPES = ('proton', 'antiproton'); SEEDS = '456'
lab = lambda L: '%g' % L
fn = lambda d, t, k, s: f'{d}/{t}_n{k[0]}L{lab(k[1])}_{s}.err'
def tomega(path):
    m = ana699.parse(path)['inc']['mechB0']
    return np.linalg.norm(m['s21'] - m['s20'])/(m['dt']*np.linalg.norm(np.cross(m['L0'], m['s20'])))/(2*math.pi)
CACHE = {}
def orb(path):   # per-orbit precession: Omega_obr, delta, K periapses, T_r/T_n, J apoapses
    if path in CACHE: return CACHE[path]
    L = []; S = []; per = None; peri = None
    for ln in open(path):
        if not ln.startswith('PAIR TORQUE799 mechB0'): continue
        w = ln.split()
        if w[3] == 'apo':
            i = w.index('L'); L.append(np.array([float(x) for x in w[i + 1:i + 4]])); i = w.index('s2'); S.append(np.array([float(x) for x in w[i + 1:i + 4]]))
        elif w[3] == 'sum':
            per = float(w[w.index('period') + 1]); k = int(w[w.index('peri') + 1]); peri = [float(x) for x in w[w.index('peri') + 2:w.index('peri') + 2 + k]]
    d = [np.linalg.norm(S[j] - S[j - 1])/np.linalg.norm(np.cross(L[j - 1], S[j - 1])) for j in range(1, len(S))]
    delta = float(np.mean(d)); tr = (peri[-1] - peri[0])/(len(peri) - 1)
    CACHE[path] = r = dict(om=delta/(2*math.pi*per), delta=delta, K=len(peri), trn=tr/per, J=len(S)); return r
def rule(W, tag, out):
    ok_d = ok_m = True
    for n, l, st in ((2, 1, '2p'), (3, 1, '3p'), (3, 2, '3d')):
        od = l*W[(n, l)] - (l + 1)*W[(n, l + 1)]; q = od/W[(n, l + 0.5)]/((l + 0.5)**2/(l*(l + 1))); ok_d &= abs(q - 1) <= 0.01
        out.setdefault(st, []).append(q)
        if st == '2p': ok_m &= abs(od/FS - 1) <= 0.01; out.setdefault('2P', []).append(od)
        pr(f'    {tag} {st}: {units(od)} | {q:.5f} x regula' + (f' | {od/FS:.5f} pomiaru' if st == '2p' else ''))
    return ok_d, ok_m
# H799k
okK = True; nK = 0
for sub in ('A', 'A0', 'AZ', 'AK'):
    for f in sorted(os.listdir(f'{O}/{sub}')):
        a = ''.join(l for l in open(f'{O}/{sub}/{f}') if not l.startswith('PAIR TORQUE799')); okK &= a == open(f'{R798}/{sub}/{f}').read(); nK += 1
pr(f"H799k ({nK} plikow bez linii TORQUE799 = 798 bajt w bajt): {'TAK' if okK else 'NIE'}")
# H799n, H799q
okN = okQ = True; trd = []
pr('=== AK (stary start) wobec A: skok czasowy x 200 | K_AK - K_A | T_r/T_n - 1 (AK, A) | Omega_obr AK/A - 1 x 200 ===')
for k in ((3, 1), (2, 1)):
    for s in SEEDS:
        for t in TYPES:
            a = orb(fn(f'{O}/A', t, k, s)); b = orb(fn(f'{O}/AK', t, k, s))
            j = (tomega(fn(f'{O}/AK', t, k, s))/tomega(fn(f'{O}/A', t, k, s)) - 1)*200; dk = b['K'] - a['K']
            if abs(j) >= 0.5: okN &= dk != 0 and np.sign(dk) == np.sign(j)
            if k == (3, 1): okN &= abs(b['trn']/a['trn'] - 1) > 1e-3
            q = (b['om']/a['om'] - 1)*200; okQ &= abs(q) <= 0.3
            pr(f'  n{k[0]}L1 {t} {s}: skok {j:+.2f} | dK {dk:+d} (K {b["K"]} / {a["K"]}) | T_r/T_n - 1 {b["trn"] - 1:+.3e} / {a["trn"] - 1:+.3e} | {q:+.3f}')
pr(f"H799n (liczba peryaps wyjasnia skoki; n3L1 T_r rozne > 1e-3): {'TAK' if okN else 'NIE'}")
pr(f"H799q (|Omega_obr AK/A - 1| x 200 <= 0,3, 12 par): {'TAK' if okQ else 'NIE'}")
okZ = True
for t in TYPES:
    q = orb(fn(f'{O}/AZ', t, (3, 1), 5))['om']/orb(fn(f'{O}/A', t, (3, 1), 5))['om'] - 1; okZ &= abs(q) <= 1e-3
    pr(f'  AZ {t}: Omega_obr(2048)/Omega_obr(512) - 1 = {q:+.3e}')
pr(f"H799z (<= 1e-3): {'TAK' if okZ else 'NIE'}")
# rule
okP = okD = okM = okG = True; out = {}; ALL = {}; W = {}
pr('=== A (F = 1), Omega_obr ===')
for s in SEEDS:
    for t in TYPES:
        W[(t, s)] = {k: orb(fn(f'{O}/A', t, k, s))['om'] for k in TORI}; Wt = {k: tomega(fn(f'{O}/A', t, k, s)) for k in TORI}
        for k in TORI: ALL.setdefault(k, []).append(W[(t, s)][k])
        K = {k: W[(t, s)][k]*k[0]**3*k[1]**2 for k in TORI}; sp = max(K.values())/min(K.values()) - 1; okP &= sp <= 0.01
        pr(f'  {t} ziarno {s}: rozrzut Omega_obr n^3 L^2 {sp*100:.3f} % | Omega_obr / czasowe - 1 x 200: ' + ', '.join(f'n{k[0]}L{lab(k[1])} {(W[(t, s)][k]/Wt[k] - 1)*200:+.2f}' for k in TORI))
        pr('    T_r/T_n - 1: ' + ', '.join(f'n{k[0]}L{lab(k[1])} {orb(fn(f"{O}/A", t, k, s))["trn"] - 1:+.2e}' for k in TORI))
        d, m = rule(W[(t, s)], f'{t} {s}', out); okD &= d; okM &= m
    g = W[('proton', s)][(3, 1)]/W[('antiproton', s)][(3, 1)] - 1; okG &= abs(g) <= 0.005
    g8 = tomega(fn(f'{O}/A', 'proton', (3, 1), s))/tomega(fn(f'{O}/A', 'antiproton', (3, 1), s)) - 1
    pr(f'  luka n3L1 ziarno {s}: proton / antyproton - 1 = {g*100:+.3f} % (798 czasowo: {g8*100:+.3f} %)')
pr(f"H799p (rozrzut <= 1 %): {'TAK' if okP else 'NIE'}")
pr(f"H799d (2p {min(out['2p']):.5f}-{max(out['2p']):.5f}, 3p {min(out['3p']):.5f}-{max(out['3p']):.5f}, 3d {min(out['3d']):.5f}-{max(out['3d']):.5f} w 1 %): {'TAK' if okD else 'NIE'}")
pr(f"H799m (2P {min(out['2P'])/1e6:.2f}-{max(out['2P'])/1e6:.2f} MHz = {min(out['2P'])/FS:.4f}-{max(out['2P'])/FS:.4f} pomiaru): {'TAK' if okM else 'NIE'}")
pr(f"H799g (luka n3L1 <= 0,5 %): {'TAK' if okG else 'NIE'}")
pr('  rozrzut 6 biegow na torusie (max/min - 1, x 200): Omega_obr | czasowe (798):')
for k in TORI:
    a = ALL[k]; b = [tomega(fn(f'{O}/A', t, k, s)) for s in SEEDS for t in TYPES]
    pr(f'    n{k[0]}L{lab(k[1])}: {(max(a)/min(a) - 1)*200:.2f} | {(max(b)/min(b) - 1)*200:.2f}')
M = {k: np.mean(v) for k, v in ALL.items()}
pr('  srednie 6 biegow:'); rule(M, 'srednie', {})
okS = True; out0 = {}; okD0 = okM0 = True
pr('=== A0 (F = 0, PHEN=1, ziarno 5), Omega_obr ===')
for t in TYPES:
    W0 = {k: orb(fn(f'{O}/A0', t, k, 5))['om'] for k in TORI}
    qs = {k: W0[k]/W[(t, '5')][k] - 1 for k in TORI}; okS &= all(abs(q) <= 1e-3 for q in qs.values())
    pr(f'  {t}: F0/F1 - 1: ' + ', '.join(f'n{k[0]}L{lab(k[1])} {q:+.2e}' for k, q in qs.items()))
    pr(f'    dK (F0 - F1): ' + ', '.join(f'n{k[0]}L{lab(k[1])} {orb(fn(f"{O}/A0", t, k, 5))["K"] - orb(fn(f"{O}/A", t, k, 5))["K"]:+d}' for k in TORI))
    d, m = rule(W0, f'{t} F = 0', out0); okD0 &= d; okM0 &= m
pr(f"H799s (|F0/F1 - 1| <= 1e-3): {'TAK' if okS else 'NIE'}")
pr(f"H799s0 (F = 0: regula i 2P w 1 %): {'TAK' if okD0 and okM0 else 'NIE'}")
R21 = [x*9/8 for x in out['2p']]
pr(f'  o-Ps: 8/9 x R x 0,9998..0,9999 = {8/9*min(R21)*0.9998:.5f}-{8/9*max(R21)*0.9999:.5f} multipletu (Hagena i in., PRL 71, 2887 (1993)); p-Ps bez czesci spinowej')
pr('Pomiar: 2P3/2 - 2P1/2 wodoru 10 969,045(15) MHz (Hagley i Pipkin, PRL 72, 1172 (1994); Lundeen i Pipkin, PRL 46, 232 (1981)).')

# Audit 801: raw secular averages (Kepler K / relativistic R nodes) at the 699 start vs the per-orbit mechanics of 799/800.  usage: python3 ana801.py
# 799 runs unpacked into $R799 (../section-799/runs/runs799.tar.xz), 800 runs into $R800 (../section-800/runs/runs800.tar.xz).
import os, math
import numpy as np
D = os.path.dirname(os.path.abspath(__file__)); O = os.path.join(D, 'out', 'runs')
SP = '/tmp/claude-1000/-home-teddy-Projekty-CREM-QI/4fa83d5c-8878-4cfd-898e-1a9b53a01ab2/scratchpad'
R799 = os.environ.get('R799', SP + '/r799x/runs'); R800 = os.environ.get('R800', SP + '/r800x/runs')
HEV, C, ALPHA = 4.135667696e-15, 299792458.0, 7.2973525693e-3; FS = 10969.045e6
pr = lambda *a: print(*a, flush=True)
units = lambda f: f'{f/1e6:.2f} MHz = {f*HEV*1e6:.4f} ueV = {C/f*100:.4f} cm'
TORI = [(3, 1), (3, 1.5), (3, 2), (3, 2.5), (3, 3), (2, 1), (2, 1.5), (2, 2)]
TYPES = ('proton', 'antiproton'); SEEDS = {1: (4, 5, 6, 7, 8, 9, 10, 11, 12), 0: (4, 5, 6)}
lab = lambda L: '%g' % L
nm = lambda t, k, s: f'{t}_n{k[0]}L{lab(k[1])}_{s}.err'
def ref(F, t, k, s):
    if F == 1: return f'{R799}/A/{nm(t, k, s)}' if s in (4, 5, 6) else f'{R800}/B/{nm(t, k, s)}'
    return f'{R799}/A0/{nm(t, k, s)}' if s == 5 else f'{R800}/B0/{nm(t, k, s)}'
new = lambda F, t, k, s: f'{O}/S/F{F}_{nm(t, k, s)}'
def om(p):   # Omega_obr as in ana799 / ana800
    L = []; S = []; per = None
    for ln in open(p):
        if not ln.startswith('PAIR TORQUE799 mechB0'): continue
        w = ln.split()
        if w[3] == 'apo':
            i = w.index('L'); L.append(np.array([float(x) for x in w[i + 1:i + 4]])); i = w.index('s2'); S.append(np.array([float(x) for x in w[i + 1:i + 4]]))
        elif w[3] == 'sum': per = float(w[w.index('period') + 1])
    d = [np.linalg.norm(S[j] - S[j - 1])/np.linalg.norm(np.cross(L[j - 1], S[j - 1])) for j in range(1, len(S))]
    return float(np.mean(d))/(2*math.pi*per)
def v3(w, tag): i = w.index(tag); return np.array([float(x) for x in w[i + 1:i + 4]])
def secline(p):
    for ln in open(p):
        if ln.startswith('PAIR QM1S699 secB0'): return ln.split()
def sec801(p):
    out = {}
    for ln in open(p):
        if ln.startswith('PAIR SEC801'): w = ln.split(); out[w[2]] = (int(w[4]), v3(w, 'second'), int(w[w.index('nodes') + 1]))
    return out
def rate(omega, Lh, s): s = s/np.linalg.norm(s); Lh = Lh/np.linalg.norm(Lh); return np.linalg.norm(np.cross(omega, s))/np.linalg.norm(np.cross(Lh, s))/(2*math.pi)
# --- H801k
K = os.path.join(O, 'K', 'proton_n3L1_5.err')
kl = [ln for ln in open(K) if not ln.startswith('PAIR SEC801')]
okK = kl == list(open(f'{R799}/A/proton_n3L1_5.err'))
pr(f"H801k (K bez SEC801 = 799 A/proton_n3L1_5.err bajt w bajt): {'TAK' if okK else 'NIE'} ({len(kl)} linii)")
# --- per run
okN = sec801(K) == sec801(new(1, 'proton', (3, 1), 5)) if False else None
a, b = sec801(K), sec801(new(1, 'proton', (3, 1), 5))
okN = all(a[x][0] == b[x][0] and a[x][2] == b[x][2] and np.array_equal(a[x][1], b[x][1]) for x in 'KR')
okV = True; W = {}
for F, seeds in SEEDS.items():
    for s in seeds:
        for t in TYPES:
            for k in TORI:
                wn, wr = secline(new(F, t, k, s)), secline(ref(F, t, k, s))
                okN &= all(np.array_equal(v3(wn, x), v3(wr, x)) for x in ('L0', 's10', 's20'))
                e = sec801(new(F, t, k, s)); okV &= e['K'][0] == 1 and e['R'][0] == 1
                L0, s20, s21, dt = v3(wn, 'L0'), v3(wn, 's20'), v3(wn, 's21'), float(wn[wn.index('dt') + 1])
                full = np.linalg.norm(s21/np.linalg.norm(s21) - s20/np.linalg.norm(s20))/(dt*np.linalg.norm(np.cross(L0/np.linalg.norm(L0), s20/np.linalg.norm(s20))))/(2*math.pi)
                W[(F, t, k, s)] = dict(K=rate(e['K'][1], L0, s20), R=rate(e['R'][1], L0, s20), obr=om(ref(F, t, k, s)), nodes=e['K'][2], full=full)
pr(f"H801n (start L0, s10, s20 = 799/800 we wszystkich biegach; SEC801 K = S): {'TAK' if okN else 'NIE'};  wszystkie linie SEC801 valid: {'TAK' if okV else 'NIE'}")
pm = lambda F, k, s, f: 0.5*(f(W[(F, 'proton', k, s)]) + f(W[(F, 'antiproton', k, s)]))
def vp2(k):
    n, L = k; e = math.sqrt(max(0.0, 1 - L*L/(n*n))); return (ALPHA/n)**2*(1 + e)/(1 - e)
okC = okR = okX = okB = True; B2 = []
for F, seeds in SEEDS.items():
    pr(f'=== F = {F} (PHEN={2 if F else 1}) ===')
    pr('  torus  v_p^2/c^2 | R/K - 1 (zakres biegow) | p(K)/p(obr) - 1 (zakres ziaren) | p(R)/p(obr) - 1 | K/obr - 1 pojedyncze biegi | pelna secB0/obr | wezly')
    for k in TORI:
        rk = [W[(F, t, k, s)]['R']/W[(F, t, k, s)]['K'] - 1 for s in seeds for t in TYPES]
        ck = [pm(F, k, s, lambda x: x['K'])/pm(F, k, s, lambda x: x['obr']) - 1 for s in seeds]
        xk = [pm(F, k, s, lambda x: x['R'])/pm(F, k, s, lambda x: x['obr']) - 1 for s in seeds]
        sk = [W[(F, t, k, s)]['K']/W[(F, t, k, s)]['obr'] - 1 for s in seeds for t in TYPES]
        fk = [W[(F, t, k, s)]['full']/W[(F, t, k, s)]['obr'] for s in seeds for t in TYPES]
        okR &= all(abs(x) <= 1e-4 for x in rk); okC &= all(abs(x) <= 1e-4 for x in ck); okX &= all(abs(x) <= 1e-4 for x in xk)
        pr(f'  n{k[0]}L{lab(k[1]):<4} {vp2(k):.2e} | {min(rk):+.2e} .. {max(rk):+.2e} | {min(ck):+.2e} .. {max(ck):+.2e} (sr. {np.mean(ck):+.2e}) | '
           f'{min(xk):+.2e} .. {max(xk):+.2e} (sr. {np.mean(xk):+.2e}) | {min(sk):+.2e} .. {max(sk):+.2e} | {min(fk):.5f} .. {max(fk):.5f} | {W[(F, "proton", k, seeds[0])]["nodes"]}')
    for n, l, st in ((2, 1, '2p'), (3, 1, '3p'), (3, 2, '3d')):
        q = []
        for s in seeds:
            od = pm(F, (n, l), s, lambda x: l*x['obr']) - pm(F, (n, l + 1), s, lambda x: (l + 1)*x['obr'])
            born = pm(F, (n, l + 0.5), s, lambda x: x['K'])*(l + 0.5)**2/(l*(l + 1)); q.append(born/od - 1)
            if st == '2p': B2.append(born)
        okB &= all(abs(x) <= 3e-4 for x in q)
        pr(f'  H801b {st}: Born x p(Omega_K(l + 1/2)) / p(Omega_rozn) - 1 = {min(q):+.2e} .. {max(q):+.2e} (sr. {np.mean(q):+.2e})')
pr(f"H801c (p(K)/p(obr) w 1e-4): {'TAK' if okC else 'NIE'}")
pr(f"H801r (R/K w 1e-4, kazdy bieg): {'TAK' if okR else 'NIE'}")
pr(f"H801x (p(R)/p(obr) w 1e-4): {'TAK' if okX else 'NIE'}")
pr(f"H801b (Born sciezki sekularnej = regula mechaniki w 3e-4, 2p 3p 3d, F = 1 i F = 0): {'TAK' if okB else 'NIE'}")
m = np.mean(B2); pr(f'RAPORT 2P ze sciezki sekularnej (9/8 x p(Omega_K) na n2L1.5, {len(B2)} par): {units(m)} = {m/FS:.5f} pomiaru ({(m - FS)/1e6:+.2f} MHz); zakres {min(B2)/1e6:.2f}-{max(B2)/1e6:.2f} MHz')
pr(f'  o-Ps przez czynnik z wodoru (2P modelu / pomiar x 0,9998..0,9999): {min(B2)/FS*0.9998:.5f}-{max(B2)/FS*0.9999:.5f} multipletu (Hagena i in., PRL 71, 2887 (1993)); p-Ps bez czesci spinowej')
pr('Pomiar: 2P3/2 - 2P1/2 wodoru 10 969,045(15) MHz (Hagley i Pipkin, PRL 72, 1172 (1994); Lundeen i Pipkin, PRL 46, 232 (1981)).')

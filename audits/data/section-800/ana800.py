# Audit 800: difference rule on the per-orbit measure (799) with 9 seeds F = 1 and 3 seeds F = 0; classification of Born's 9/8 import.  usage: python3 ana800.py
# 799 A / A0 unpacked from ../section-799/runs/runs799.tar.xz into $R799.
import sys, os, math
import numpy as np
D = os.path.dirname(os.path.abspath(__file__)); O = os.path.join(D, 'out', 'runs')
R799 = os.environ.get('R799', '/tmp/claude-1000/-home-teddy-Projekty-CREM-QI/4fa83d5c-8878-4cfd-898e-1a9b53a01ab2/scratchpad/r799x/runs')
HEV, C = 4.135667696e-15, 299792458.0; FS = 10969.045e6
pr = lambda *a: print(*a, flush=True)
units = lambda f: f'{f/1e6:.2f} MHz = {f*HEV*1e6:.4f} ueV = {C/f*100:.4f} cm'
TORI = [(3, 1), (3, 1.5), (3, 2), (3, 2.5), (3, 3), (2, 1), (2, 1.5), (2, 2)]
TYPES = ('proton', 'antiproton')
lab = lambda L: '%g' % L
def path(F, t, k, s):
    if F == 1: return f'{R799}/A/{t}_n{k[0]}L{lab(k[1])}_{s}.err' if s in (4, 5, 6) else f'{O}/B/{t}_n{k[0]}L{lab(k[1])}_{s}.err'
    return f'{R799}/A0/{t}_n{k[0]}L{lab(k[1])}_{s}.err' if s == 5 else f'{O}/B0/{t}_n{k[0]}L{lab(k[1])}_{s}.err'
def om(p):   # Omega_obr as in ana799
    L = []; S = []; per = None
    for ln in open(p):
        if not ln.startswith('PAIR TORQUE799 mechB0'): continue
        w = ln.split()
        if w[3] == 'apo':
            i = w.index('L'); L.append(np.array([float(x) for x in w[i + 1:i + 4]])); i = w.index('s2'); S.append(np.array([float(x) for x in w[i + 1:i + 4]]))
        elif w[3] == 'sum': per = float(w[w.index('period') + 1])
    d = [np.linalg.norm(S[j] - S[j - 1])/np.linalg.norm(np.cross(L[j - 1], S[j - 1])) for j in range(1, len(S))]
    return float(np.mean(d))/(2*math.pi*per)
RULES = ((2, 1, '2p'), (3, 1, '3p'), (3, 2, '3d'))
def analyse(F, seeds, tag):
    okP = okD = okM = okG = True; R = {st: [] for _, _, st in RULES}; P2 = []; W = {}
    pr(f'=== {tag} ===')
    for s in seeds:
        for t in TYPES:
            w = W[(t, s)] = {k: om(path(F, t, k, s)) for k in TORI}
            K = {k: w[k]*k[0]**3*k[1]**2 for k in TORI}; sp = max(K.values())/min(K.values()) - 1; okP &= sp <= 0.01
            cells = []
            for n, l, st in RULES:
                od = l*w[(n, l)] - (l + 1)*w[(n, l + 1)]; q = od/w[(n, l + 0.5)]/((l + 0.5)**2/(l*(l + 1))); okD &= abs(q - 1) <= 0.01; R[st].append(q)
                cells.append(f'{st} {q:.5f}')
                if st == '2p': okM &= abs(od/FS - 1) <= 0.01; P2.append(od)
            pr(f'  {t} ziarno {s}: rozrzut {sp*100:.3f} % | ' + ', '.join(cells) + f' | 2P {P2[-1]/1e6:.2f} MHz = {P2[-1]/FS:.5f} pomiaru')
        if F == 1:
            g = W[('proton', s)][(3, 1)]/W[('antiproton', s)][(3, 1)] - 1; okG &= abs(g) <= 0.005
            pr(f'    luka n3L1 ziarno {s}: {g*100:+.3f} %')
    A = {}
    for st, v in R.items():
        m = np.mean(v); se = np.std(v, ddof=1)/math.sqrt(len(v)); A[st] = (m, se)
        pr(f'  srednia {st} ({len(v)} ukladow): {m:.5f} +- {se:.5f} x regula (odchylka {(m - 1)/se:+.2f} bledu standardowego)')
    m = np.mean(P2); se = np.std(P2, ddof=1)/math.sqrt(len(P2))
    pr(f'  2P srednio: {units(m)} +- {se/1e6:.2f} MHz -> {m/FS:.5f} pomiaru ({(m - FS)/1e6:+.2f} MHz, {(m - FS)/se:+.2f} bledu standardowego modelu)')
    return dict(okP=okP, okD=okD, okM=okM, okG=okG, A=A, W=W, P2=P2, R2p=R['2p'])
F1 = analyse(1, (4, 5, 6, 7, 8, 9, 10, 11, 12), 'F = 1 (PHEN=2), ziarna 4-12, wodor i antywodor')
F0 = analyse(0, (4, 5, 6), 'F = 0 (PHEN=1), ziarna 4-6, wodor i antywodor')
okA = all(abs(m - 1) <= 2e-3 and abs(m - 1) <= 3*se for m, se in F1['A'].values())
okA0 = all(abs(m - 1) <= 3e-3 for m, _ in F0['A'].values())
pr(f"H800p (rozrzut <= 1 %, F = 1 i F = 0): {'TAK' if F1['okP'] and F0['okP'] else 'NIE'}")
pr(f"H800d (F = 1, 18 ukladow, 2p 3p 3d w 1 %): {'TAK' if F1['okD'] else 'NIE'}")
pr(f"H800d0 (F = 0, 6 ukladow): {'TAK' if F0['okD'] else 'NIE'}")
pr(f"H800a (F = 1 srednie w 2e-3 i 3 bledach standardowych): {'TAK' if okA else 'NIE'}")
pr(f"H800a0 (F = 0 srednie w 3e-3): {'TAK' if okA0 else 'NIE'}")
pr(f"H800m (2P kazdego ukladu F = 1 w 1 %): {'TAK' if F1['okM'] else 'NIE'}")
pr(f"H800g (luka n3L1 <= 0,5 %, 9 ziaren): {'TAK' if F1['okG'] else 'NIE'}")
okS = okL = True
pr('=== F = 0 / F = 1 (ziarna 4-6) ===')
for s in (4, 5, 6):
    for t in TYPES:
        q = {k: F0['W'][(t, s)][k]/F1['W'][(t, s)][k] - 1 for k in TORI}; okS &= all(abs(x) <= 1e-3 for x in q.values())
        r12 = q[(3, 1)]/q[(3, 2)]; r23 = q[(2, 1)]/q[(3, 1)] - 1; okL &= 6 <= r12 <= 10 and abs(r23) <= 0.1
        pr(f'  {t} {s}: ' + ', '.join(f'n{k[0]}L{lab(k[1])} {x:+.2e}' for k, x in q.items()) + f' | n3L1/n3L2 {r12:.2f}, n2L1/n3L1 - 1 {r23:+.3f}')
pr(f"H800s (|F0/F1 - 1| <= 1e-3; prawo 1/L^3: n3L1/n3L2 w [6, 10], n2L1 = n3L1 do 10 %): {'TAK' if okS and okL else 'NIE'} (granica {'TAK' if okS else 'NIE'}, prawo {'TAK' if okL else 'NIE'})")
cls = F1['okD'] and F0['okD'] and okA and okA0
pr(f"KLASYFIKACJA (H800d, d0, a, a0): {'czesc (A) importu Borna 9/8 = wynik dynamiki silnika (wodor F = 0 i F = 1); importem zostaje regula roznicowa (B)' if cls else 'status bez zmian'}")
R21 = [x*9/8 for x in F1['R2p']]
pr(f'  o-Ps: 8/9 x R x 0,9998..0,9999 = {8/9*min(R21)*0.9998:.5f}-{8/9*max(R21)*0.9999:.5f} multipletu (Hagena i in., PRL 71, 2887 (1993)); p-Ps bez czesci spinowej')
pr('Pomiar: 2P3/2 - 2P1/2 wodoru 10 969,045(15) MHz (Hagley i Pipkin, PRL 72, 1172 (1994); Lundeen i Pipkin, PRL 46, 232 (1981)).')

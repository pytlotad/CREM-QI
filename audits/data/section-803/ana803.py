# Audit 803: source of the +6.5e-5 excess of the mechanical difference rule -- the Sommerfeld alpha^2 term carried through the rule.  usage: python3 ana803.py
# n = 2, 3: mechanics from 799 / 800 ($R799, $R800), secular K / S from 802 ($R802); n = 4: new runs out/runs/M (mechanics and K / S in the same file).
import os, math
import numpy as np
D = os.path.dirname(os.path.abspath(__file__)); O = os.path.join(D, 'out', 'runs')
SP = '/tmp/claude-1000/-home-teddy-Projekty-CREM-QI/4fa83d5c-8878-4cfd-898e-1a9b53a01ab2/scratchpad'
R799 = os.environ.get('R799', SP + '/r799x/runs'); R800 = os.environ.get('R800', SP + '/r800x/runs'); R802 = os.environ.get('R802', SP + '/r802x/runs')
HEV, C, ALPHA = 4.135667696e-15, 299792458.0, 7.2973525693e-3; FS = 10969.045e6
pr = lambda *a: print(*a, flush=True)
units = lambda f: f'{f/1e6:.2f} MHz = {f*HEV*1e6:.4f} ueV = {C/f*100:.4f} cm'
TYPES = ('proton', 'antiproton'); lab = lambda L: '%g' % L
nm = lambda t, k, s: f'{t}_n{k[0]}L{lab(k[1])}_{s}.err'
SEEDS = {(1, 'old'): (4, 5, 6, 7, 8, 9, 10, 11, 12), (0, 'old'): (4, 5, 6), (1, 'new'): (4, 5, 6), (0, 'new'): (5,)}
RULES = {'old': ((2, 1, '2p'), (3, 1, '3p'), (3, 2, '3d')), 'new': ((4, 1, '4p'), (4, 2, '4d'), (4, 3, '4f'))}
def mech_path(F, t, k, s):
    if k[0] == 4: return f'{O}/M/F{F}_{nm(t, k, s)}'
    if F == 1: return f'{R799}/A/{nm(t, k, s)}' if s in (4, 5, 6) else f'{R800}/B/{nm(t, k, s)}'
    return f'{R799}/A0/{nm(t, k, s)}' if s == 5 else f'{R800}/B0/{nm(t, k, s)}'
sec_path = lambda F, t, k, s: mech_path(F, t, k, s) if k[0] == 4 else f'{R802}/S/F{F}_{nm(t, k, s)}'
def v3(w, tag): i = w.index(tag); return np.array([float(x) for x in w[i + 1:i + 4]])
def om(p):   # Omega_obr as in ana799-802
    L = []; S = []; per = None
    for ln in open(p):
        if not ln.startswith('PAIR TORQUE799 mechB0'): continue
        w = ln.split()
        if w[3] == 'apo': L.append(v3(w, 'L')); S.append(v3(w, 's2'))
        elif w[3] == 'sum': per = float(w[w.index('period') + 1])
    d = [np.linalg.norm(S[j] - S[j - 1])/np.linalg.norm(np.cross(L[j - 1], S[j - 1])) for j in range(1, len(S))]
    return float(np.mean(d))/(2*math.pi*per), len(S)
def sec(p):   # secular K / S rates about L at the start (as ana802)
    e = {}; wn = None
    for ln in open(p):
        if ln.startswith('PAIR SEC801'): w = ln.split(); e[w[2]] = (int(w[4]), v3(w, 'second'))
        elif ln.startswith('PAIR QM1S699 secB0') and wn is None: wn = ln.split()
    L0, s20 = v3(wn, 'L0'), v3(wn, 's20'); s = s20/np.linalg.norm(s20); Lh = L0/np.linalg.norm(L0)
    r = lambda o: np.linalg.norm(np.cross(o, s))/np.linalg.norm(np.cross(Lh, s))/(2*math.pi)
    return {x: r(e[x][1]) for x in 'KS'}, all(e[x][0] == 1 for x in 'KS')
Ef = lambda l: 1.5*ALPHA**2*((3*l*l + 3*l + 1)/(l*l*(l + 1)**2) - 1/(l + 0.5)**2)
W = {}; okV = True; okP = okD = True
for (F, g), seeds in SEEDS.items():
    for s in seeds:
        for t in TYPES:
            tori = sorted({(n, l + d) for n, l, _ in RULES[g] for d in (0, 0.5, 1)})
            for k in tori:
                (o, napo) = om(mech_path(F, t, k, s)); se, v = sec(sec_path(F, t, k, s)); okV &= v and napo >= 150
                W[(F, t, k, s)] = dict(obr=o, K=se['K'], S=se['S'])
            if g == 'new':
                Kn = {k: W[(F, t, k, s)]['obr']*k[0]**3*k[1]**2 for k in tori}; okP &= max(Kn.values())/min(Kn.values()) - 1 <= 0.01
                for n, l, st in RULES[g]:
                    w = lambda k: W[(F, t, k, s)]['obr']
                    q = (l*w((n, l)) - (l + 1)*w((n, l + 1)))/(w((n, l + 0.5))*(l + 0.5)**2/(l*(l + 1))); okD &= abs(q - 1) <= 0.01
pm = lambda F, k, s, x: 0.5*(W[(F, 'proton', k, s)][x] + W[(F, 'antiproton', k, s)][x])
def E(F, s, n, l, x):
    od = l*pm(F, (n, l), s, x) - (l + 1)*pm(F, (n, l + 1), s, x)
    return od/(pm(F, (n, l + 0.5), s, x)*(l + 0.5)**2/(l*(l + 1))) - 1, od
pr(f"linie SEC801 K / S valid i >= 150 apocentrow we wszystkich plikach: {'TAK' if okV else 'NIE'}")
okF = okM = okN = okM0 = okN0 = True; R = {}
for F in (1, 0):
    pr(f'=== F = {F} (PHEN={2 if F else 1}); srednie par (wodor + antywodor)/2 po ziarnach ===')
    pr('  regula | E_f        | E_S        E_S/E_f | E_K        | E_mech (sr. +- blad std; zakres ziaren)         | E_mech - E_S | E_mech - E_f')
    for g in ('old', 'new'):
        seeds = SEEDS[(F, g)]
        for n, l, st in RULES[g]:
            eS = np.mean([E(F, s, n, l, 'S')[0] for s in seeds]); eK = np.mean([E(F, s, n, l, 'K')[0] for s in seeds])
            em = [E(F, s, n, l, 'obr')[0] for s in seeds]; m = np.mean(em); se = np.std(em, ddof=1)/math.sqrt(len(em)) if len(em) > 1 else float('nan')
            R[(F, st)] = (m, eS)
            if F == 1: okF &= abs(eS/Ef(l) - 1) <= 0.1
            if F == 1 and g == 'old': okM &= abs(m - eS) <= 3e-5
            if F == 1 and g == 'new': okN &= abs(m - Ef(l)) <= 3e-5
            if F == 0 and g == 'old': okM0 &= abs(m - eS) <= 5e-5
            if F == 0 and g == 'new': okN0 &= abs(m - Ef(l)) <= 5e-5
            pr(f'  {st}     | {Ef(l):+.3e} | {eS:+.3e} {eS/Ef(l):.3f} | {eK:+.3e} | {m:+.3e} +- {se:.1e} ({min(em):+.2e} .. {max(em):+.2e}) | {m - eS:+.2e}    | {m - Ef(l):+.2e}')
    d42 = R[(F, '4p')][0] - R[(F, '2p')][0]; pr(f'  E_mech(4p) - E_mech(2p) = {d42:+.2e}')
    if F == 1: okN &= abs(d42) <= 3e-5
okM2 = True
for F in (1, 0):
    P = [E(F, s, 2, 1, 'obr')[1] for s in SEEDS[(F, 'old')]]; eS = R[(F, '2p')][1]; m = np.mean(P); c = m/(1 + eS)
    pr(f'RAPORT 2P F = {F}: regula roznicowa mechaniki {units(m)} = {m/FS:.6f} pomiaru; bez czlonu Sommerfelda reguly (/(1 + E_S(2p))) {units(c)} = {c/FS:.6f} pomiaru')
    if F == 1: R2 = (m, c)
pr(f"H803f (|E_S/E_f - 1| <= 0,1, 6 regul, F = 1): {'TAK' if okF else 'NIE'}")
pr(f"H803m (|E_mech - E_S| <= 3e-5, 2p 3p 3d, F = 1): {'TAK' if okM else 'NIE'}")
pr(f"H803n (slepo n = 4: |E_mech - E_f| <= 3e-5 dla 4p 4d 4f i |E(4p) - E(2p)| <= 3e-5, F = 1): {'TAK' if okN else 'NIE'}")
pr(f"H803p (n = 4, F = 1 i F = 0: rozrzut Omega n^3 L^2 <= 1 % i regula 4p 4d 4f w 1 % w kazdym biegu): {'TAK' if okP and okD else 'NIE'} (rozrzut {'TAK' if okP else 'NIE'}, regula {'TAK' if okD else 'NIE'})")
pr(f"H803s0 (F = 0: |E_mech - E_S| <= 5e-5 dla 2p 3p 3d; |E_mech - E_f| <= 5e-5 dla 4p 4d 4f): {'TAK' if okM0 and okN0 else 'NIE'}")
pr(f'  o-Ps przez czynnik z wodoru (2P / pomiar x 0,9998..0,9999): regula mechaniki {R2[0]/FS*0.9998:.5f}-{R2[0]/FS*0.9999:.5f}, bez czlonu Sommerfelda {R2[1]/FS*0.9998:.5f}-{R2[1]/FS*0.9999:.5f} multipletu (Hagena i in., PRL 71, 2887 (1993)); p-Ps bez czesci spinowej')
pr('Pomiar: 2P3/2 - 2P1/2 wodoru 10 969,045(15) MHz (Hagley i Pipkin, PRL 72, 1172 (1994); Lundeen i Pipkin, PRL 46, 232 (1981)).')

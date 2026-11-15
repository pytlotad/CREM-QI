# Audit 793: difference rule on tori that hold their energy (CREM_NO_DIPOLE_FORCE: no moment forces on the orbit, full precession) and the energy
# ledger of the drifting channel (PAIR TORQUE792 lines with E = conservedParticleEnergy, Ec = conservativeParticleEnergy).  usage: python3 ana793.py
# 792 trace references unpacked from ../section-792/runs/runs792.tar.xz into $R792 (default scratchpad copy).
import sys, os, math, re
import numpy as np
D = os.path.dirname(os.path.abspath(__file__)); O = os.path.join(D, 'out', 'runs')
sys.path.insert(0, os.path.join(D, '..', 'section-699')); import ana699
R792 = os.environ.get('R792', '/tmp/claude-1000/-home-teddy-Projekty-CREM-QI/4fa83d5c-8878-4cfd-898e-1a9b53a01ab2/scratchpad/r792x/runs')
HEV, C = 4.135667696e-15, 299792458.0; FS = 10969.045e6
pr = lambda *a: print(*a, flush=True)
units = lambda f: f'{f/1e6:.2f} MHz = {f*HEV*1e6:.4f} ueV = {C/f*100:.4f} cm'
TORI = [(3, 1), (3, 1.5), (3, 2), (3, 2.5), (3, 3), (2, 1), (2, 1.5), (2, 2)]
lab = lambda L: '%g' % L
def omega(path):
    m = ana699.parse(path)['inc']['mechB0']
    return np.linalg.norm(m['s21'] - m['s20'])/(m['dt']*np.linalg.norm(np.cross(m['L0'], m['s20'])))/(2*math.pi), m['Ln0']
def trace(path):
    out = []
    for ln in open(path):
        if not ln.startswith('PAIR TORQUE792'): continue
        w = ln.split(); g = lambda k: float(w[w.index(k) + 1])
        out.append(dict(eps=g('eps'), E=g('E'), Ec=g('Ec'), minsep=g('minsep')))
    return out
strip = lambda path: re.sub(r' E \S+ Ec \S+\n', '\n', open(path).read())
# --- H793k
okK = all(strip(f'{O}/E/{t}_E512.err') == open(f'{R792}/T/{t}_n3L1_5.err').read() for t in ('proton', 'antiproton'))
pr(f"H793k (E512 bez pol E, Ec = 792 T bajt w bajt): {'TAK' if okK else 'NIE'}")
# --- energy
pr('=== bilans energii (n = 3 L = 1, ziarno 5; wzgledne zmiany obieg 1 -> 200) ===')
TRc = {}
for t in ('proton', 'antiproton'):
    for w in ('E512', 'E2048', 'E0'):
        r = trace(f'{O}/E/{t}_{w}.err'); TRc[(t, w)] = r
        q = lambda k: r[-1][k]/r[0][k] - 1
        pr(f'  {t} {w}: eps {q("eps"):+.3e} | E {q("E"):+.3e} | Ec {q("Ec"):+.3e} | E(1) {r[0]["E"]/1.602176634e-19:.6f} eV, E(200) {r[-1]["E"]/1.602176634e-19:.6f} eV'
           f' | Omega {omega(f"{O}/E/{t}_{w}.err")[0]/1e6:.2f} MHz')
r = TRc[('proton', 'E512')]; okE = abs(r[-1]['E']/r[0]['E'] - 1) <= 1e-4
pr(f"H793e (proton E512: |E(200)/E(1) - 1| <= 1e-4): {'TAK' if okE else 'NIE'}")
okN = all(abs(TRc[(t, 'E0')][-1]['eps']/TRc[(t, 'E0')][0]['eps'] - 1) <= 1e-3 for t in ('proton', 'antiproton'))
pr(f"H793n (E0 obu kanalow: |eps(200)/eps(1) - 1| <= 1e-3): {'TAK' if okN else 'NIE'}")
# --- Z
okZ = True
for n, L in ((3, 1), (2, 1)):
    a = omega(f'{O}/G/proton_n{n}L{L}_4.err')[0]; b = omega(f'{O}/Z/proton_n{n}L{L}_4.err')[0]; okZ &= abs(a/b - 1) <= 1e-3
    pr(f'  n {n} L {L}: 512 {a/1e6:.3f} MHz, 2048 {b/1e6:.3f} MHz -> {a/b - 1:+.2e}')
pr(f"H793z (zbieznosc 512 / 2048 w 1e-3, bez sil od momentow): {'TAK' if okZ else 'NIE'}")
# --- rule
okP = okD = okM = True; R21 = []; ALL = {}
for t in ('proton', 'antiproton'):
    for s in '456':
        pr(f'=== {t} ziarno {s} (CREM_NO_DIPOLE_FORCE) ===')
        W = {}
        for n, L in TORI:
            om, ln = omega(f'{O}/G/{t}_n{n}L{lab(L)}_{s}.err'); W[(n, L)] = (om, ln); ALL.setdefault((n, L), []).append(om)
            pr(f'  n {n} L {L:3}: Ln0 {ln:.6f} | Omega/2pi {units(om)} | Omega n^3 L^2 {om*n**3*ln**2/1e6:.3f} MHz')
        k = [om*n**3*ln**2 for (n, L), (om, ln) in W.items()]; sp = max(k)/min(k) - 1; okP &= sp <= 0.01
        pr(f'  H793p rozrzut Omega n^3 L^2: {sp*100:.3f} %')
        for n, l, st in ((2, 1, '2p'), (3, 1, '3p'), (3, 2, '3d')):
            F = lambda L: W[(n, L)][1]*W[(n, L)][0]
            od = F(l) - F(l + 1); oc = W[(n, l + 0.5)][0]; q = od/oc; qt = (l + 0.5)**2/(l*(l + 1)); okD &= abs(q/qt - 1) <= 0.01
            rf = FS*16/(n**3*l*(l + 1))
            pr(f'  {st}: Omega_rozn/2pi {units(od)} | stosunek {q:.5f} wobec {qt:.5f} -> {q/qt:.5f} | wobec skalowania pomiaru {rf/1e6:.2f} MHz -> {od/rf:.5f} ({(od - rf)/1e6:+.2f} MHz)')
            if st == '2p': okM &= abs(od/FS - 1) <= 0.01; R21.append(q)
pr(f"H793p (Omega n^3 L^2 stale w 1 %): {'TAK' if okP else 'NIE'}")
pr(f"H793d (regula roznicowa w 1 %, 2p 3p 3d): {'TAK' if okD else 'NIE'}")
pr(f"H793m (wodor 2P = 10 969,045(15) MHz w 1 %): {'TAK' if okM else 'NIE'}")
pr('=== RAPORT ===')
for (n, L), v in ALL.items():
    v = np.array(v); pr(f'  n {n} L {L:3}: rozrzut 6 biegow {(v.max()/v.min() - 1)*100:.3f} % | srednia {v.mean()/1e6:.2f} MHz | Omega n^3 L^2 {v.mean()*n**3*L**2/1e6:.3f} MHz')
for s in '456':
    a = omega(f'{O}/G/proton_n3L1_{s}.err')[0]; b = omega(f'{O}/G/antiproton_n3L1_{s}.err')[0]
    pr(f'  luka kanalow n3L1 ziarno {s}: {a/b - 1:+.4f}')
pr(f'  o-Ps: R (2p) {min(R21):.5f}-{max(R21):.5f}; 8/9 x R x 0,9998..0,9999 = {8/9*min(R21)*0.9998:.5f}-{8/9*max(R21)*0.9999:.5f} multipletu (Hagena i in. 1993); p-Ps bez czesci spinowej')
pr('Pomiar: 2P3/2 - 2P1/2 wodoru 10 969,045(15) MHz (Hagley i Pipkin, PRL 72, 1172 (1994); Lundeen i Pipkin, PRL 46, 232 (1981)).')

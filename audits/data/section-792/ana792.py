# Audit 792: where the seed / channel dependence of the lepton spin precession rate on the L = 1 torus (791) comes from.
# K control, T chunk traces (PAIR TORQUE792, k = 1), S sector "charge", D no dipole-dipole force.  usage: python3 ana792.py
# 791 reference runs (G) from ../section-791/runs/runs791.tar.xz unpacked into $R791 (default scratchpad copy).
import sys, os, math, re
import numpy as np
D = os.path.dirname(os.path.abspath(__file__)); O = os.path.join(D, 'out', 'runs')
sys.path.insert(0, os.path.join(D, '..', 'section-699')); import ana699
R791 = os.environ.get('R791', '/tmp/claude-1000/-home-teddy-Projekty-CREM-QI/4fa83d5c-8878-4cfd-898e-1a9b53a01ab2/scratchpad/r791x/runs')
HEV, C = 4.135667696e-15, 299792458.0
pr = lambda *a: print(*a, flush=True)
units = lambda f: f'{f/1e6:.2f} MHz = {f*HEV*1e6:.4f} ueV = {C/f*100:.4f} cm'
def omega(path):
    m = ana699.parse(path)['inc']['mechB0']
    return np.linalg.norm(m['s21'] - m['s20'])/(m['dt']*np.linalg.norm(np.cross(m['L0'], m['s20'])))/(2*math.pi), m
def chunks(path):
    rows = []
    for ln in open(path):
        if not ln.startswith('PAIR TORQUE792'): continue
        w = ln.split(); g = lambda k: float(w[w.index(k) + 1]); v = lambda k: np.array([float(x) for x in w[w.index(k) + 1:w.index(k) + 4]])
        rows.append(dict(c=int(g('c')), t=g('t'), minsep=g('minsep'), Ln=g('Ln'), eps=g('eps'), e=g('e'), L=v('L'), s1=v('s1'), s2=v('s2'), A=v('A')))
    return rows
strip = lambda path: ''.join(l for l in open(path) if not l.startswith('PAIR TORQUE792'))
# --- H792k
okK = open(f'{O}/K/proton_n3L1_5.err').read() == open(f'{R791}/G/proton_n3L1_5.err').read()
okK &= all(strip(f'{O}/T/{t}_n3L1_5.err') == open(f'{R791}/G/{t}_n3L1_5.err').read() for t in ('proton', 'antiproton'))
pr(f"H792k (K = 791 bajt w bajt; T bez linii TORQUE792 = 791): {'TAK' if okK else 'NIE'}")
# --- H792c and trace report
okC = True; TR = {}
for t in ('proton', 'antiproton'):
    p = f'{O}/T/{t}_n3L1_5.err'; om, m = omega(p); rows = chunks(p)
    s_prev, L_prev, t_prev = m['s20'], m['L0'], 0.0; rate = []; ang = 0.0
    for r in rows:
        d = np.linalg.norm(r['s2'] - s_prev)/np.linalg.norm(np.cross(L_prev, s_prev)); ang += d
        rate.append(d/(r['t'] - t_prev)/(2*math.pi)); s_prev, L_prev, t_prev = r['s2'], r['L'], r['t']
    rate = np.array(rate); oms = ang/rows[-1]['t']/(2*math.pi); okC &= abs(oms/om - 1) <= 1e-3; TR[t] = (rows, rate, om)
    pr(f'  {t}: mechB0 {units(om)} | suma 200 kawalkow {oms/1e6:.2f} MHz -> {oms/om - 1:+.2e} | koniec: |s2 - s21| {np.linalg.norm(rows[-1]["s2"] - m["s21"]):.1e}, |L - L1| {np.linalg.norm(rows[-1]["L"] - m["L1"]):.1e}')
pr(f"H792c (suma kawalkow = mechB0 w 1e-3): {'TAK' if okC else 'NIE'}")
pr('=== RAPORT: slad T (k = 1, wodor / antywodor ziarno 5, n = 3 L = 1) ===')
for t, (rows, rate, om) in TR.items():
    ep = np.array([r['eps'] for r in rows]); ee = np.array([r['e'] for r in rows]); ln = np.array([r['Ln'] for r in rows]); ms = np.array([r['minsep'] for r in rows])
    As = np.array([abs(np.dot(r['A'], r['s2'])) for r in rows]); A0 = rows[0]['A']; tilt = [math.degrees(math.acos(max(-1, min(1, np.dot(r['A'], A0))))) for r in rows]
    pr(f'  {t}: tempo na obieg {rate.mean()/1e6:.2f} MHz (sd {rate.std(ddof=1)/1e6:.2f}, min {rate.min()/1e6:.2f}, max {rate.max()/1e6:.2f}); '
       f'eps {ep[0]:.6e} -> {ep[-1]:.6e} ({ep[-1]/ep[0] - 1:+.2e}); e {ee[0]:.6f} -> {ee[-1]:.6f}; Ln {ln[0]:.6f} -> {ln[-1]:.6f}; minsep {ms.min():.4e}-{ms.max():.4e} m; '
       f'obrot A^ {tilt[-1]:.3f} st.; |A^.s2^| {As[0]:.4f} -> {As[-1]:.4f}')
    for q, name in ((ee, 'e'), (ms, 'minsep'), (As, '|A.s|')):
        if np.std(q) > 0: pr(f'    korelacja tempa na obieg z {name}: {np.corrcoef(rate, q)[0, 1]:+.3f}')
    cum = np.cumsum(rate)/np.arange(1, len(rate) + 1)
    pr('    srednia narastajaca (obiegi 1, 2, 5, 10, 20, 50, 100, 150, 200): ' + ', '.join(f'{cum[i - 1]/1e6:.1f}' for i in (1, 2, 5, 10, 20, 50, 100, 150, 200)) + ' MHz')
rp, ra = TR['proton'][1], TR['antiproton'][1]
d = np.cumsum(rp - ra)/np.arange(1, len(rp) + 1)
pr('  roznica srednich narastajacych proton - antyproton (te same obiegi): ' + ', '.join(f'{d[i - 1]/1e6:+.1f}' for i in (1, 2, 5, 10, 20, 50, 100, 150, 200)) + ' MHz')
pr(f'  korelacja tempa na obieg proton / antyproton: {np.corrcoef(rp, ra)[0, 1]:+.3f}')
# --- H792s
pr('=== H792s: sektor "charge" (bez pola dipola partnera w precesji) ===')
S = {}
for t in ('proton', 'antiproton'):
    for s in '456':
        a = omega(f'{O}/S/{t}_n3L1_{s}.err')[0]; b = omega(f'{R791}/G/{t}_n3L1_{s}.err')[0]; S[(t, s)] = a
        pr(f'  {t} ziarno {s}: charge {units(a)} | 791 {b/1e6:.2f} MHz | roznica {(a - b)/1e6:+.2f} MHz ({a/b - 1:+.4f})')
v = np.array(list(S.values())); sp = v.max()/v.min() - 1
pr(f'  rozrzut 6 biegow: {sp*100:.3f} % (791: 1,92 %); Omega n^3 L^2 srednio {v.mean()*27/1e6:.3f} MHz (791 n3L1 173 712,7, n3L1,5 175 054,1)')
pr(f"H792s (rozrzut w sektorze charge <= 0,3 %): {'TAK' if sp <= 3e-3 else 'NIE'}")
for t in ('proton', 'antiproton'):
    a = omega(f'{O}/S/{t}_n3L1.5_5.err')[0]; b = omega(f'{R791}/G/{t}_n3L1.5_5.err')[0]
    pr(f'  raport L = 1,5 {t} ziarno 5: charge {a/1e6:.2f} MHz | 791 {b/1e6:.2f} MHz ({a/b - 1:+.4f}); Omega n^3 L^2 {a*27*2.25/1e6:.3f} MHz')
# --- H792f
pr('=== H792f: bez sily dipol-dipol ===')
a = omega(f'{O}/D/proton_n3L1_5.err')[0]; b = omega(f'{O}/D/antiproton_n3L1_5.err')[0]
pr(f'  proton {units(a)}, antyproton {units(b)} -> {a/b - 1:+.4f} (791: {omega(f"{R791}/G/proton_n3L1_5.err")[0]/omega(f"{R791}/G/antiproton_n3L1_5.err")[0] - 1:+.4f})')
pr(f"H792f (|p / pbar - 1| <= 0,3 %): {'TAK' if abs(a/b - 1) <= 3e-3 else 'NIE'}")
pr('Diagnostyka silnika -- pomiar odniesienia nie dotyczy; wynik dotyczy czesci spinowej o-Ps przez czynnik z wodoru (791); p-Ps bez czesci spinowej.')

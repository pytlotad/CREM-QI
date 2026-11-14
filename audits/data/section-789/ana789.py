# Audit 789: (A) canonical pair of the retarded sum -- candidate K4 = J + Jq + JTh + JT (Hamilton ranking, Ps o - p against the 2^3P_J multiplet,
# retarded K4 = instantaneous K3 on Ps n = 3 and hydrogen); (B) lepton spin rate of the mechanics over the mean field of H_2P = 4/3, L^ = 8/9.
# usage: python3 ana789.py   (runs in out/runs/{P,K,X3,XH,B})
import sys, os, math
import numpy as np
D = os.path.dirname(os.path.abspath(__file__)); O = os.path.join(D, 'out', 'runs')
for sec in ('section-711', 'section-710', 'section-709', 'section-702', 'section-699'): sys.path.insert(0, os.path.join(D, '..', sec))
import ana709, ana702
from ana709 import means, CAND, EV
from ana702 import A_MHZ, B_MHZ, HEV
import ana699
CAND['K4 J+Jq+JTh+JT'] = (('J', 'Jq', 'JTh', 'JT'), ('Lcan', 'Lq', 'LTh', 'LT'))
K3, K4 = 'K3 J+Jq+JTh', 'K4 J+Jq+JTh+JT'
pr = lambda *a: print(*a, flush=True)
C = {}
def P(p):
    if p not in C: C[p] = ana702.parse(p)
    return C[p]
def residual(base, cand, de):
    m, p = means(P(f'{base}_e-{de}.err'), cand), means(P(f'{base}_e{de}.err'), cand)
    m0, p0 = means(P(f'{base}_e-{de}.err'), 'K0 J'), means(P(f'{base}_e{de}.err'), 'K0 J')
    T = 0.5*(m['T'] + p['T']); dphi = 0.5*(m['dphi'] + p['dphi'])
    return (p['Jr'] - m['Jr'] - T*(p['E'] - m['E']) + dphi*(p['L'] - m['L']))/abs(dphi*(p0['L'] - m0['L']))
rich = lambda base, cand: (4*residual(base, cand, '1e-3') - residual(base, cand, '2e-3'))/3
def spinE(fa, fb, cand):   # energy of n_clock(a) - n_clock(b) [eV]
    a, b = P(fa), P(fb); ma, mb = means(a, cand), means(b, cand)
    n = math.sqrt(a['binding']/abs(ma['E'])); return -(2*a['binding']/n**3)*(ma['n'] - mb['n'])/EV
def multiplet(fo, fp):
    o, p = P(fo), P(fp); L, s1, s2 = o['L'], o['s1'], o['s2']
    LS = 1.5*float(L @ (0.5*s1 + 0.5*s2)); S12 = 0.5*float(s1 @ s2) - 1.5*float(s1 @ L)*float(s2 @ L)
    return 8/9*(A_MHZ*LS + B_MHZ*S12)*1e6*HEV - np.mean([k['Udd'] for k in p['k']])/EV, LS, S12
MHZ = lambda e: e/HEV/1e6
pr('=== CZESC A ===')
okK = True; REF = {('ins', '0'): 0.9993, ('ins', '1'): 0.9994, ('ret', '0'): 1.0001, ('ret', '1'): 1.0002}
for sm, cand in (('ins', K3), ('ret', K4)):
    for s in ('0', '1'):
        fo, fp = f'{O}/K/{sm}/ps_A2_{s}_e0.err', f'{O}/K/{sm}/ps_A1_{s}_e0.err'
        e = spinE(fo, fp, cand); pred = multiplet(fo, fp)[0]; okK &= abs(e/pred - REF[(sm, s)]) <= 1e-3
        pr(f'  K {sm} {cand} ziarno {s}: o - p {MHZ(e):+.2f} MHz wobec {MHZ(pred):+.2f} MHz -> {e/pred:.4f} (711: {REF[(sm, s)]})')
pr(f"H789k (kontrola HEAD wobec 711 w 1e-3): {'TAK' if okK else 'NIE'}")
okA = True
for s in ('2', '3'):
    sp = {}
    for cand in CAND:
        ro, rp = rich(f'{O}/P/ps_A2_{s}', cand), rich(f'{O}/P/ps_A1_{s}', cand); sp[cand] = ro - rp
        pr(f'  P ziarno {s} {cand:15s}: R/|dphi dL| o-Ps {ro:+.3e}, p-Ps {rp:+.3e} | czesc spinowa (o - p) {sp[cand]:+.3e}')
    order = sorted(CAND, key=lambda c: abs(sp[c])); ratio = abs(sp[order[1]])/max(abs(sp[order[0]]), 1e-300)
    okA &= order[0] == K4 and ratio >= 5; pr(f'   -> najmniejszy {order[0]}, nastepny {order[1]} ({ratio:.1f}x)')
pr(f"H789a (Hamilton: K4 wygrywa >= 5x, Ps ziarna 2 i 3, suma opozniona): {'TAK' if okA else 'NIE'}")
okB = True
for s in ('2', '3'):
    fo, fp = f'{O}/P/ps_A2_{s}_e0.err', f'{O}/P/ps_A1_{s}_e0.err'
    e = spinE(fo, fp, K4); pred, LS, S12 = multiplet(fo, fp); okB &= abs(e/pred - 1) <= 0.01
    a, b = P(fo), P(fp); n = lambda d: means(d, K4)['n']
    pr(f'  P ziarno {s}: L.S {LS:+.4f} S12 {S12:+.4f} | n_clock - n_E o-Ps {n(a):+.4e}, p-Ps {n(b):+.4e} | o - p {MHZ(e):+.2f} MHz = {e*1e6:+.4f} ueV'
       f' = {abs(299792458/(MHZ(e)*1e6))*100:.3f} cm wobec {MHZ(pred):+.2f} MHz -> {e/pred:.4f}')
pr(f"H789b (pomiar: o - p z K4 = 8/9 multipletu - U_dd(p) w 1 %): {'TAK' if okB else 'NIE'}")
okC = True
for lab, sub, pairs in (('Ps n = 3', 'X3', [(f'ps_A2_{s}', f'ps_A1_{s}') for s in ('0', '1')]),
                        ('wodor n = 2', 'XH', [(f'proton_{s}', f'antiproton_{s}') for s in ('2', '3')])):
    for a, b in pairs:
        er = spinE(f'{O}/{sub}/ret/{a}_e0.err', f'{O}/{sub}/ret/{b}_e0.err', K4)
        ei = spinE(f'{O}/{sub}/ins/{a}_e0.err', f'{O}/{sub}/ins/{b}_e0.err', K3)
        e0 = spinE(f'{O}/{sub}/ret/{a}_e0.err', f'{O}/{sub}/ret/{b}_e0.err', 'K0 J')
        okC &= abs(er/ei - 1) <= 0.02
        pr(f'  {lab} {a} - {b}: opozniona K4 {MHZ(er):+.2f} MHz ({er*1e6:+.4f} ueV), chwilowa K3 {MHZ(ei):+.2f} MHz -> {er/ei:.4f} (opozniona K0 {er and e0/ei:.4f})')
pr(f"H789c (opozniona K4 = chwilowa K3 w 2 %, Ps n = 3 i wodor): {'TAK' if okC else 'NIE'}")
pr('=== CZESC B ===')
def b0(f):
    d = ana699.parse(f); g = d['inc']['sec']; X0 = [g['L0'], g['s10'], g['s20']]; m = d['inc']['mechB0']
    mech = [np.linalg.norm(np.asarray(m[v]) - X0[i]) for i, v in enumerate(('L1', 's11', 's21'))]
    hz = ana699.hartree(d['H']['fs'] + d['H']['hfs'], X0, g['dt']); hb = [np.linalg.norm(np.asarray(hz[i]) - X0[i]) for i in range(3)]
    om = mech[2]/(g['dt']*np.linalg.norm(np.cross(m['L0'], m['s20'])))/(2*math.pi)
    return mech[0]/hb[0], mech[2]/hb[2], om, m['Ln0']
REF699 = {'proton': (0.892, 1.337), 'antiproton': (0.890, 1.336)}; FS = 10969.045e6
okQ = okS = okW = True
for t in ('proton', 'antiproton'):
    for s in ('0', '2', '3'):
        rl, rs, om, ln = b0(f'{O}/B/{t}_{s}.err'); pw = 8/9*FS*(1.5/ln)**2
        if s == '0': okQ &= abs(rl/REF699[t][0] - 1) <= 0.005 and abs(rs/REF699[t][1] - 1) <= 0.005
        else:
            okS &= abs(rl/(8/9) - 1) <= 0.01 and abs(rs/(4/3) - 1) <= 0.01; okW &= abs(om/pw - 1) <= 0.02
        pr(f'  {t:10s} ziarno {s}: mech/Hartree L^ {rl:.4f} (8/9 = 0,8889), lepton {rs:.4f} (4/3 = 1,3333) | Omega_lepton/2pi {om/1e6:.2f} MHz'
           f' = {om*HEV*1e6:.3f} ueV wobec (8/9) x 10 969,045 MHz x (1,5/Ln0)^2 = {pw/1e6:.2f} MHz (Ln0 {ln:.5f}) -> {om/pw:.4f}')
pr(f"H789q (kontrola ziarno 0 wobec 699 w 0,5 %): {'TAK' if okQ else 'NIE'}")
pr(f"H789s (mech/Hartree = 8/9 i 4/3 w 1 %, ziarna 2, 3): {'TAK' if okS else 'NIE'}")
pr(f"H789w (Omega_lepton = (8/9) x przedzial 2P w 2 %): {'TAK' if okW else 'NIE'}")
pr('Pomiar: 2P3/2 - 2P1/2 wodoru 10 969,045(15) MHz (Hagley i Pipkin 1994; Lundeen i Pipkin 1981); Ps 2^3P_J: Hagena i in., PRL 71, 2887 (1993).')

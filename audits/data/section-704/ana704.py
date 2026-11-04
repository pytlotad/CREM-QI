# Audit 704: the classical spin energies of the mechanical path against Breit-Pauli -- spin-orbit coefficient g/(2 m^2) (own motional dipole only)
# against BP (g-1)/(2 m^2) + g/(2 m1 m2), classical <r^-3> at L = 1.5 hbar against the quantum 2P value.
# usage: python3 ana704.py <704 runs dir with h/> <702 runs dir with s8192/>
import sys, os, math
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'section-702'))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'section-701'))
from ana702 import parse, A_MHZ, B_MHZ, HEV
from ana701 import derived, HBAR, H as HPLANCK
C = 299792458.0; ME = 9.1093837015e-31; MP = 1.67262192369e-27; ALPHA = 7.2973525693e-3; KC = ALPHA*HBAR*C
G = 2.00231930436256; EV = 1.602176634e-19
DE_H = 10969.045e6*HEV   # eV, 2P3/2 - 2P1/2 (Hagley & Pipkin 1994 + Lundeen & Pipkin 1981)


def r3_classical(mu, n, k):
    a = n*n*HBAR**2/(mu*KC); return n**3/(a**3*k**3)


def main():
    RH, R2 = sys.argv[1:3]; out = []
    # H704a
    g8 = G/(G - 1)*8/9
    out.append(f'H704a: g/(g-1) x 8/9 = {g8:.4f}; 701 ratios 1.775 / 1.774: {abs(1.775/g8 - 1) < 0.01 and abs(1.774/g8 - 1) < 0.01}')
    # H704b hydrogen
    mu = ME*MP/(ME + MP); okB = True; diffs = {}
    for s in ('0', '1'):
        for tag in ('proton', 'antiproton'):
            d = parse(f'{RH}/h/{tag}_B0.01_A2_{s}.err'); rows = derived(d)
            Uqd = np.mean([q['Uqd'] for q in d['k']]); Udd = np.mean([q['Udd'] for q in d['k']])
            n = np.mean([r['nE'] for r in rows]); k = np.mean([r['X'] for r in rows])
            lep = d['s1'] if tag == 'antiproton' and False else None
            # the lepton is the lighter particle: second for proton,electron / antiproton,positron in the probe's PAIR order
            sl = d['s2']; LS = (k*HBAR)*(0.5*HBAR)*float(d['L'] @ sl)
            pred = G/(2*ME*ME)*KC/(C*C)*r3_classical(mu, n, k)*LS
            okB &= abs(Uqd/pred - 1) <= 0.02
            ncl = np.mean([r['ncl'] - r['nE'] for r in rows])
            diffs[(s, tag)] = (ncl, Uqd + Udd, n, d['binding'])
            out.append(f'hydrogen-like {tag:10s} seed {s}: L.s {float(d["L"] @ sl):+.4f} <U_qd> {Uqd/EV*1e6:+.4f} ueV  own-motional prediction {pred/EV*1e6:+.4f} ueV'
                       f' (ratio {Uqd/pred:.4f}) | <U_dd> {Udd/EV*1e6:+.4f} ueV | BP spin-orbit (g-1)/(2m^2) with QM <r^-3> {pred*(G-1)/G*9/8/EV*1e6:+.4f} ueV')
    okB2 = True
    for s in ('0', '1'):
        (nh, uh, n, B), (na, ua, _, _) = diffs[(s, 'proton')], diffs[(s, 'antiproton')]
        lhs = nh - na; rhs = -(n**3/(2*B))*(uh - ua); okB2 &= abs(lhs/rhs - 1) <= 0.02
        out.append(f'seed {s}: (H - Hbar) clock {lhs:+.4e} vs -(n^3/2B) (H - Hbar) <U_qd + U_dd> {rhs:+.4e} (ratio {lhs/rhs:.4f}; default step -- clock error cancels in H - Hbar)')
    out.append(f'H704b (U_qd = own-motional g/(2m^2) prediction within 2 %; clock H - Hbar = U difference within 2 %): {okB} / {okB2}')
    # H704c Ps
    okC = True; okD = True
    pred_so = G/(2*G - 1)*8/9
    for s in ('0', '1', '2', '3'):
        o = parse(f'{R2}/s8192/ps_A2_{s}.err'); p = parse(f'{R2}/s8192/ps_A1_{s}.err')
        Uo = np.mean([q['Uqd'] for q in o['k']]); Do = np.mean([q['Udd'] for q in o['k']]); Dp = np.mean([q['Udd'] for q in p['k']])
        Up = np.mean([q['Uqd'] for q in p['k']])
        L, s1, s2 = o['L'], o['s1'], o['s2']
        LS = 1.5*float(L @ (0.5*s1 + 0.5*s2)); S12 = 0.5*float(s1 @ s2) - 1.5*float(s1 @ L)*float(s2 @ L)
        Ea = A_MHZ*1e6*HEV*LS*EV; Eb = B_MHZ*1e6*HEV*S12*EV
        rso = (Uo - Up)/Ea if Ea else float('nan'); rt = (Do - Dp)/Eb if Eb else float('nan')
        okC &= abs(rso/pred_so - 1) <= 0.05
        tot_model = (Uo - Up) + (Do - Dp); tot_meas = Ea + Eb
        out.append(f'Ps seed {s}: L.S {LS:+.4f} S12 {S12:+.4f} | spin-orbit model {(Uo - Up)/EV*1e6:+.4f} ueV / measured a L.S {Ea/EV*1e6:+.4f} ueV = {rso:.4f} (prediction {pred_so:.4f})'
                   f' | tensor model {(Do - Dp)/EV*1e6:+.4f} ueV / measured b S12 {Eb/EV*1e6:+.4f} ueV = {rt:.4f} | total {tot_model/tot_meas:.4f}')
    out.append(f'H704c (Ps spin-orbit part = g/(2g-1) x 8/9 = {pred_so:.4f} within 5 %, every seed): {okC}')
    text = '\n'.join(out); print(text)
    os.makedirs('out', exist_ok=True); open('out/ana704.txt', 'w').write(text + '\n')


if __name__ == '__main__':
    main()

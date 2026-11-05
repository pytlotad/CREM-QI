# Audit 707 POST HOC: the same clock with the energy function of the instantaneous sum (charge-dipole weight 0): E -> <E> - <U_qd> per period
# (both are period means over the same samples), for the 707 runs (L_T form) and the 705 runs (kappa = 1, default retarded sum).
# usage: python3 post707.py <707 runs> <705 runs>
import sys, os, math
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'section-706'))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'section-702'))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'section-701'))
from ana706 import DE_H, EV
from ana702 import parse, A_MHZ, B_MHZ, HEV
from ana701 import HBAR, H


def clk(path, extra):
    d = parse(path); n = []; nE = []
    for v in d['k']:
        Jr = v['J'] - v['Lcan']*v['dphi']; ncl = Jr/H + v['Lcan']/HBAR
        e = sum(v.get(k, 0.0) for k in extra)/H
        n.append(ncl + e - math.sqrt(d['binding']/abs(v['E'] - v['Uqd']))); nE.append(math.sqrt(d['binding']/abs(v['E'] - v['Uqd'])))
    return np.mean(n), np.mean(nE), d


def main():
    R7, R5 = sys.argv[1:3]; out = []
    for tag, R, extra in (('707 L_T form, J + Jq + JTh', R7, ('Jq', 'JTh')), ('705 kappa = 1, J only', R5, ()), ('705 kappa = 1, J + Jh + JT', R5, ('Jh', 'JT'))):
        for s in ('0', '1'):
            h, n, d = clk(f'{R}/h/proton_B0.01_A2_{s}.err', extra); a, _, _ = clk(f'{R}/h/antiproton_B0.01_A2_{s}.err', extra)
            cos = float(d['L'] @ d['s2']); r = -(2*d['binding']/n**3)*(h - a)/EV/(DE_H*cos)
            out.append(f'{tag:28s} hydrogen seed {s}: energy weight 0 -> (H - Hbar) / (Delta E cos) {r:.4f} (BP x 8/9 = {8/9:.4f})')
    for s in ('0', '1', '2', '3'):
        o, n, d = clk(f'{R7}/ps/ps_A2_{s}.err', ('Jq', 'JTh')); p, _, dp = clk(f'{R7}/ps/ps_A1_{s}.err', ('Jq', 'JTh'))
        L, s1, s2 = d['L'], d['s1'], d['s2']
        LS = 1.5*float(L @ (0.5*s1 + 0.5*s2)); S12 = 0.5*float(s1 @ s2) - 1.5*float(s1 @ L)*float(s2 @ L)
        Udd = np.mean([k['Udd'] for k in dp['k']])
        pred = 8/9*(A_MHZ*LS + B_MHZ*S12)*1e6*HEV - Udd/EV
        e = -(2*d['binding']/n**3)*(o - p)/EV
        out.append(f'707 L_T form Ps seed {s}: energy weight 0 -> o - p {e*1e6:+.4f} ueV vs 8/9 (a L.S + b S12) - U_dd(p) {pred*1e6:+.4f} ueV (ratio {e/pred:.4f})')
    text = '\n'.join(out); print(text)
    os.makedirs('out', exist_ok=True); open('out/post707.txt', 'w').write(text + '\n')


if __name__ == '__main__':
    main()

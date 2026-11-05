# Audit 706: the clock on the instantaneous sum (production mechanical path): canonical action = J (kinetic + Darwin) + int L_qmu + int L_T,
# energy with charge-dipole weight 0; H - Hbar and o-Ps - p-Ps against Breit-Pauli x 8/9 (measured intervals).
# usage: python3 ana706.py <706 runs dir with h/ ps/ ch/ cps/> <705 runs dir>
import sys, os, glob, math, filecmp
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'section-702'))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'section-701'))
from ana702 import parse, A_MHZ, B_MHZ, HEV
from ana701 import derived, H
EV = 1.602176634e-19; DE_H = 10969.045e6*HEV   # eV (Hagley & Pipkin 1994 + Lundeen & Pipkin 1981)


def clock(path):
    d = parse(path); rows = derived(d)
    base = np.mean([r['ncl'] - r['nE'] for r in rows]); q = np.mean([k['Jq'] for k in d['k']])/H; t = np.mean([k['JTh'] for k in d['k']])/H
    return dict(base=base, q=q, T=t, n=np.mean([r['nE'] for r in rows]), B=d['binding'], d=d, Udd=np.mean([k['Udd'] for k in d['k']]))


def main():
    R, R5 = sys.argv[1:3]; out = []
    ctl = [(f, f'{R5}/h/{os.path.basename(f)}') for f in sorted(glob.glob(f'{R}/ch/*.err'))] + [(f, f'{R5}/ps/{os.path.basename(f)}') for f in sorted(glob.glob(f'{R}/cps/*.err'))]
    same = sum(filecmp.cmp(a, b, shallow=False) for a, b in ctl)
    out.append(f'H706a (files): clock without CREM_PAIR_CLOCK_INSTANT = 705 byte for byte {same}/{len(ctl)}: {same == len(ctl) == 12}')
    okB = okC = True
    for s in ('0', '1'):
        h = clock(f'{R}/h/proton_B0.01_A2_{s}.err'); a = clock(f'{R}/h/antiproton_B0.01_A2_{s}.err')
        cos = float(h['d']['L'] @ h['d']['s2']); toE = lambda dn: -(2*h['B']/h['n']**3)*dn/EV
        b = h['base'] - a['base']; dq = h['q'] - a['q']; dT = h['T'] - a['T']
        r_all = toE(b + dq + dT)/(DE_H*cos); okB &= abs(r_all/(8/9) - 1) <= 0.03
        out.append(f'hydrogen seed {s}: cos {cos:+.4f} | H - Hbar in n: kinetic+Darwin {b:+.4e}, L_qmu {dq:+.4e}, L_T {dT:+.4e} | energy / (Delta E cos):'
                   f' canonical {r_all:.4f} (BP x 8/9 = {8/9:.4f}); without L_T {toE(b + dq)/(DE_H*cos):.4f}; kinetic only {toE(b)/(DE_H*cos):.4f}')
    for s in ('0', '1', '2', '3'):
        o = clock(f'{R}/ps/ps_A2_{s}.err'); p = clock(f'{R}/ps/ps_A1_{s}.err')
        L, s1, s2 = o['d']['L'], o['d']['s1'], o['d']['s2']
        LS = 1.5*float(L @ (0.5*s1 + 0.5*s2)); S12 = 0.5*float(s1 @ s2) - 1.5*float(s1 @ L)*float(s2 @ L)
        pred = 8/9*(A_MHZ*LS + B_MHZ*S12)*1e6*HEV - p['Udd']/EV
        toE = lambda dn: -(2*o['B']/o['n']**3)*dn/EV
        e = toE((o['base'] + o['q'] + o['T']) - (p['base'] + p['q'] + p['T']))
        okC &= abs(e/pred - 1) <= 0.05
        out.append(f'Ps seed {s}: L.S {LS:+.4f} S12 {S12:+.4f} | o - p canonical {e*1e6:+.4f} ueV vs 8/9 (a L.S + b S12) - U_dd(p) {pred*1e6:+.4f} ueV (ratio {e/pred:.4f})'
                   f' | parts: kinetic {toE(o["base"] - p["base"])*1e6:+.4f} L_qmu {toE(o["q"] - p["q"])*1e6:+.4f} L_T {toE(o["T"] - p["T"])*1e6:+.4f} ueV')
    out.append(f'H706b (hydrogen canonical = BP x 8/9 within 3 %, both seeds): {okB}')
    out.append(f'H706c (Ps o - p = 8/9 (a L.S + b S12) - U_dd(p) within 5 %, every seed): {okC}')
    text = '\n'.join(out); print(text)
    os.makedirs('out', exist_ok=True); open('out/ana706.txt', 'w').write(text + '\n')


if __name__ == '__main__':
    main()

# Audit 707: the clock on the retarded sum in its 569 form (L_T, CREM_LAB_SPIN_MOMENTA=1) -- canonical action J + int L_qmu + int L_T against
# Breit-Pauli x 8/9 (measured intervals), as 706 on the instantaneous sum; torque of the dipole sector against the precession (p567); periapsis
# precession H - Hbar for the three forms (707 L_T, 706 instantaneous, 705 kappa = 1).
# usage: python3 ana707.py <707 runs dir with h/ ps/ ch/ cps/ tq/> <706 runs dir> <705 runs dir>
import sys, os, glob, math, filecmp, re
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'section-706'))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'section-702'))
from ana706 import clock, DE_H, EV
from ana702 import parse, A_MHZ, B_MHZ, HEV


def prec(path):
    return np.mean([k['dphi'] - 2*math.pi for k in parse(path)['k']])


def torque(path):
    rows = []
    for line in open(path):
        m = re.search(r'PRODUCTION .*rel\s+(\S+) \(/\|dS\|\)', line)
        if m: rows.append(float(m.group(1)))
    return rows


def main():
    R, R6, R5 = sys.argv[1:4]; out = []
    ctl = [(f, f'{R6}/ch/{os.path.basename(f)}') for f in sorted(glob.glob(f'{R}/ch/*.err'))] + [(f, f'{R6}/cps/{os.path.basename(f)}') for f in sorted(glob.glob(f'{R}/cps/*.err'))]
    same = sum(filecmp.cmp(a, b, shallow=False) for a, b in ctl)
    out.append(f'H707a (files): clock without CREM_LAB_SPIN_MOMENTA = 706 controls byte for byte {same}/{len(ctl)}: {same == len(ctl) == 12}')
    okB = okC = True
    for s in ('0', '1'):
        h = clock(f'{R}/h/proton_B0.01_A2_{s}.err'); a = clock(f'{R}/h/antiproton_B0.01_A2_{s}.err')
        cos = float(h['d']['L'] @ h['d']['s2']); toE = lambda dn: -(2*h['B']/h['n']**3)*dn/EV
        b = h['base'] - a['base']; dq = h['q'] - a['q']; dT = h['T'] - a['T']
        r_all = toE(b + dq + dT)/(DE_H*cos); okB &= abs(r_all/(8/9) - 1) <= 0.03
        out.append(f'hydrogen seed {s}: cos {cos:+.4f} | H - Hbar in n: kinetic+Darwin {b:+.4e}, L_qmu {dq:+.4e}, L_T {dT:+.4e} | energy / (Delta E cos):'
                   f' J + Jq + JTh {r_all:.4f} (BP x 8/9 = {8/9:.4f}); J + Jq {toE(b + dq)/(DE_H*cos):.4f}; J only {toE(b)/(DE_H*cos):.4f}')
        pr = {}
        for tag, D in (('707 L_T', f'{R}/h'), ('706 instantaneous', f'{R6}/h'), ('705 kappa = 1', f'{R5}/h')):
            pr[tag] = prec(f'{D}/proton_B0.01_A2_{s}.err') - prec(f'{D}/antiproton_B0.01_A2_{s}.err')
        out.append(f'   periapsis precession per period H - Hbar: ' + ', '.join(f'{k} {v:+.4e} rad' for k, v in pr.items())
                   + f' | L_T / instantaneous {pr["707 L_T"]/pr["706 instantaneous"]:.4f}, kappa 1 / instantaneous {pr["705 kappa = 1"]/pr["706 instantaneous"]:.4f}')
    for s in ('0', '1', '2', '3'):
        o = clock(f'{R}/ps/ps_A2_{s}.err'); p = clock(f'{R}/ps/ps_A1_{s}.err')
        L, s1, s2 = o['d']['L'], o['d']['s1'], o['d']['s2']
        LS = 1.5*float(L @ (0.5*s1 + 0.5*s2)); S12 = 0.5*float(s1 @ s2) - 1.5*float(s1 @ L)*float(s2 @ L)
        pred = 8/9*(A_MHZ*LS + B_MHZ*S12)*1e6*HEV - p['Udd']/EV
        toE = lambda dn: -(2*o['B']/o['n']**3)*dn/EV
        e = toE((o['base'] + o['q'] + o['T']) - (p['base'] + p['q'] + p['T']))
        okC &= abs(e/pred - 1) <= 0.05
        out.append(f'Ps seed {s}: L.S {LS:+.4f} S12 {S12:+.4f} | o - p J + Jq + JTh {e*1e6:+.4f} ueV vs 8/9 (a L.S + b S12) - U_dd(p) {pred*1e6:+.4f} ueV (ratio {e/pred:.4f})'
                   f' | parts: kinetic {toE(o["base"] - p["base"])*1e6:+.4f} L_qmu {toE(o["q"] - p["q"])*1e6:+.4f} L_T {toE(o["T"] - p["T"])*1e6:+.4f} ueV'
                   f' | precession o - p: L_T {prec(f"{R}/ps/ps_A2_{s}.err") - prec(f"{R}/ps/ps_A1_{s}.err"):+.4e}'
                   f' instantaneous {prec(f"{R6}/ps/ps_A2_{s}.err") - prec(f"{R6}/ps/ps_A1_{s}.err"):+.4e} rad')
    okD = True; worst = {}
    for m in ('def', 'lab'):
        for f in sorted(glob.glob(f'{R}/tq/{m}_*.txt')):
            r = torque(f); w = max(r) if r else float('nan'); worst[(m, os.path.basename(f)[4:-4])] = w
            if m == 'lab': okD &= w <= 0.02
    for (m, k), w in worst.items():
        out.append(f'torque {m} {k:28s} max |<tau> + <dS/dt>|/|dS| (PRODUCTION) {w:.5f}')
    out.append(f'H707a (files): see first line; regression and validation -- addendum')
    out.append(f'H707b (hydrogen J + Jq + JTh = BP x 8/9 within 3 %, both seeds): {okB}')
    out.append(f'H707c (Ps o - p = 8/9 (a L.S + b S12) - U_dd(p) within 5 %, every seed): {okC}')
    out.append(f'H707d (L_T form: torque = precession within 2 %, every PRODUCTION row): {okD}')
    text = '\n'.join(out); print(text)
    os.makedirs('out', exist_ok=True); open('out/ana707.txt', 'w').write(text + '\n')


if __name__ == '__main__':
    main()

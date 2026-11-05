# Audit 705: the clock action with the spin momenta of the retarded sum (hidden m x E/c^2, Thomas q E x S/(m c^2)) -- H - Hbar and o-Ps - p-Ps against
# Breit-Pauli x 8/9 (measured intervals).
# usage: python3 ana705.py <705 runs dir with h/ ps/>
import sys, os, math
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'section-702'))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'section-701'))
from ana702 import parse, A_MHZ, B_MHZ, HEV
from ana701 import derived, HBAR, H
EV = 1.602176634e-19; G = 2.00231930436256; DE_H = 10969.045e6*HEV   # eV (Hagley & Pipkin 1994 + Lundeen & Pipkin 1981)


def clock(path):
    d = parse(path); rows = derived(d); out = []
    for r, q in zip(rows, d['k']):
        corr_h = q['Jh']/H + q['Lh']*(1 - q['dphi']/(2*math.pi))/HBAR
        corr_T = q['JT']/H + q['LT']*(1 - q['dphi']/(2*math.pi))/HBAR
        out.append((r['ncl'] - r['nE'], corr_h, corr_T))
    a = np.array(out); n = np.mean([r['nE'] for r in rows])
    return dict(base=a[:, 0].mean(), h=a[:, 1].mean(), T=a[:, 2].mean(), n=n, B=d['binding'], d=d,
                Udd=np.mean([q['Udd'] for q in d['k']]))


def main():
    R = sys.argv[1]; out = []; okB = True; okC = True
    for s in ('0', '1'):
        h = clock(f'{R}/h/proton_B0.01_A2_{s}.err'); a = clock(f'{R}/h/antiproton_B0.01_A2_{s}.err')
        cos = float(h['d']['L'] @ h['d']['s2'])
        dE = lambda dn: -(2*h['B']/h['n']**3)*dn/EV        # eV per unit of n difference (opposite sign of the action shift)
        base = h['base'] - a['base']; dh = h['h'] - a['h']; dT = h['T'] - a['T']
        ref = DE_H*abs(cos)
        r0 = dE(base)/(DE_H*cos) ; r1 = dE(base + dh + dT)/(DE_H*cos)
        okB &= abs(r1/((G - 1)/G*G/(G - 1)*8/9) - 1) <= 0.02
        out.append(f'hydrogen seed {s}: cos {cos:+.4f} | H - Hbar: n_clock - n_E {base:+.4e}, hidden {dh:+.4e}, Thomas {dT:+.4e} | energy / (Delta E cos):'
                   f' without spin momenta {r0:.4f} (704: 1.7757), with {r1:.4f} (BP x 8/9 = {8/9:.4f}); hidden alone {dE(base + dh)/(DE_H*cos):.4f}, Thomas alone {dE(base + dT)/(DE_H*cos):.4f}')
    for s in ('0', '1', '2', '3'):
        o = clock(f'{R}/ps/ps_A2_{s}.err'); p = clock(f'{R}/ps/ps_A1_{s}.err')
        L, s1, s2 = o['d']['L'], o['d']['s1'], o['d']['s2']
        LS = 1.5*float(L @ (0.5*s1 + 0.5*s2)); S12 = 0.5*float(s1 @ s2) - 1.5*float(s1 @ L)*float(s2 @ L)
        Emeas = (A_MHZ*LS + B_MHZ*S12)*1e6*HEV                              # eV
        pred = 8/9*Emeas - p['Udd']/EV                                     # eV; the p-Ps classical U_dd is in the clock's p-Ps value
        dEf = lambda dn: -(2*o['B']/o['n']**3)*dn/EV
        e0 = dEf(o['base'] - p['base']); e1 = dEf((o['base'] + o['h'] + o['T']) - (p['base'] + p['h'] + p['T']))
        okC &= abs(e1/pred - 1) <= 0.03
        out.append(f'Ps seed {s}: L.S {LS:+.4f} S12 {S12:+.4f} | o - p energy without spin momenta {e0*1e6:+.4f} ueV, with {e1*1e6:+.4f} ueV'
                   f' | prediction 8/9 (a L.S + b S12) - U_dd(p) {pred*1e6:+.4f} ueV (ratio {e1/pred:.4f}) | hidden {dEf(o["h"] - p["h"])*1e6:+.4f} Thomas {dEf(o["T"] - p["T"])*1e6:+.4f} ueV')
    out.append(f'H705b (hydrogen: with spin momenta = BP x 8/9 = 0.8889 within 2 %, both seeds): {okB}')
    out.append(f'H705c (Ps: o - p with spin momenta = 8/9 (a L.S + b S12) - U_dd(p) within 3 %, every seed): {okC}')
    text = '\n'.join(out); print(text)
    os.makedirs('out', exist_ok=True); open('out/ana705.txt', 'w').write(text + '\n')


if __name__ == '__main__':
    main()

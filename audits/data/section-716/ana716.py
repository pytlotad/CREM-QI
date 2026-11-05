# Audit 716: angle psi of each spin's in-plane component measured from the periapsis (Runge-Lenz direction) per radial period on the mechanical path
# (CREM_PAIR_CLOCK_LRL=1), against the secular prediction Delta psi = (pi alpha^2/L^2)[(g - 2) + 3 (g - 1) S_par/L] (hydrogen, lighter particle).
# usage: python3 ana716.py <716 runs dir with lrl/ ctl/>
import sys, os, math, re
import numpy as np

ALPHA = 7.2973525693e-3; G_E = 2.00231930436256


def parse(path):
    d = {'k': [], 'lrl': []}
    for line in open(path):
        if line.startswith('PAIR CLOCK701 start'):
            w = line.split(); d['L'] = float(w[w.index('Xlock') + 1])
        elif line.startswith('PAIR CLOCK701 k'):
            w = line.split(); d['k'].append(float(w[w.index('dphi') + 1]))
        elif line.startswith('PAIR LRL716'):
            w = line.split(); v = lambda t: np.array([float(x) for x in w[w.index(t) + 1:w.index(t) + 4]])
            d['lrl'].append((v('peri'), v('L'), v('s1'), v('s2')))
    return d


def psi_track(d, which):
    ps, cs = [], []
    for p, L, s1, s2 in d['lrl']:
        s = s1 if which == 1 else s2
        pp = p - L*(p @ L); pp /= np.linalg.norm(pp); q = np.cross(L, pp)
        ps.append(math.atan2(s @ q, s @ pp)); cs.append(s @ L)
    ps = np.unwrap(np.array(ps)); k = np.arange(len(ps))
    slope = np.polyfit(k, ps, 1)[0]
    return slope, float(np.mean(cs)), float(np.std(cs)), ps


def main():
    R = sys.argv[1]; out = []
    same = []
    for name in ('proton_2', 'ps_A1_0'):
        a = ''.join(l for l in open(f'{R}/lrl/{name}.err') if not l.startswith('PAIR LRL716'))
        same.append(a == open(f'{R}/ctl/{name}.err').read())
    out.append(f'H716a (without the switch = with it minus the LRL716 lines, byte for byte): hydrogen {same[0]}, Ps {same[1]}')
    okB = okC = True
    for t in ('proton', 'antiproton'):
        for s in ('2', '3'):
            d = parse(f'{R}/lrl/{t}_{s}.err'); L = d['L']; unit = math.pi*ALPHA**2/L**2
            dpsi, cos, cstd, _ = psi_track(d, 2); Spar = 0.5*cos
            pred = unit*((G_E - 2) + 3*(G_E - 1)*Spar/L)
            aps = np.mean(d['k']) - 2*math.pi
            okB &= abs(dpsi/pred - 1) <= 0.15
            if abs(cos) > 0.2: okC &= abs(dpsi)/unit > 0.1
            out.append(f'{t:10s} seed {s}: L/hbar {L:.6f}, periods {len(d["lrl"])}, cos(L, s) {cos:+.5f} (sd {cstd:.1e}); Delta psi {dpsi:+.4e} rad/period'
                       f' vs predicted {pred:+.4e} (ratio {dpsi/pred:.4f}); /(pi alpha^2/L^2) {dpsi/unit:+.4f}; apsidal advance {aps:+.4e} ({aps/unit:+.4f} units)')
    out.append(f'H716b (hydrogen Delta psi = prediction within 15 %): {okB}')
    out.append(f'H716c (|Delta psi| > 0.1 pi alpha^2/L^2 wherever |cos| > 0.2: no lock to the Runge-Lenz vector): {okC}')
    for ph, ch in (('1', 'p-Ps'), ('2', 'o-Ps')):
        for s in ('0', '1'):
            d = parse(f'{R}/lrl/ps_A{ph}_{s}.err'); L = d['L']; unit = math.pi*ALPHA**2/L**2
            r1 = psi_track(d, 1); r2 = psi_track(d, 2); aps = np.mean(d['k']) - 2*math.pi
            out.append(f'{ch} seed {s}: L/hbar {L:.6f}; spin 1 cos {r1[1]:+.5f} Delta psi {r1[0]:+.4e} ({r1[0]/unit:+.4f}); spin 2 cos {r2[1]:+.5f}'
                       f' Delta psi {r2[0]:+.4e} ({r2[0]/unit:+.4f}); apsidal advance {aps:+.4e} ({aps/unit:+.4f}) [units pi alpha^2/L^2 = {unit:.4e}]')
    for lab, Lq, sq in (('rule Q', 1.5, 1.0), ('QM', math.sqrt(2), math.sqrt(0.75))):
        out.append(f'lock direction S_par = 0: |L + s| = sqrt(L^2 + s^2) = {math.hypot(Lq, sq):.4f} ({lab} lengths) -- not an integer k = j + 1/2')
    text = '\n'.join(out); print(text)
    os.makedirs('out', exist_ok=True); open('out/ana716.txt', 'w').write(text + '\n')


if __name__ == '__main__':
    main()

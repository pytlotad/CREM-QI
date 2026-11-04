# Audit 700: the 699 start increments with the fine-structure midpoint (528) and the Born inverse-cube rule switched off -- secular against
# mechanical at zero field, L^ and the lepton spin; c = <dD_mech, dD_sec>/<dD_sec, dD_sec>.
# usage: python3 ana700.py <700 runs dir with both/ fs/>
import sys, glob, os
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'section-699'))
from ana699 import parse, pops, unit


def main():
    R = sys.argv[1]; out = []; ok = {'both': True, 'fs': True}; okC = True
    for v, target in (('both', 1.0), ('fs', 9/8)):
        for f in sorted(glob.glob(f'{R}/{v}/*.err')):
            d = parse(f); g = d['inc']['sec']; X0 = [g['L0'], g['s10'], g['s20']]
            r = [np.linalg.norm(np.asarray(d['inc']['secB0'][a]) - g[b])/np.linalg.norm(np.asarray(d['inc']['mechB0'][a]) - g[b])
                 for a, b in (('L1', 'L0'), ('s21', 's20'), ('s11', 's10'))]
            ok[v] &= all(abs(x/target - 1) <= 0.02 for x in r[:2])
            w, V = np.linalg.eigh(d['H']['all']); p0 = pops(V, X0)
            dD = {k: pops(V, [d['inc'][k]['L1'], d['inc'][k]['s11'], d['inc'][k]['s21']]) - p0 for k in ('sec', 'mech')}
            c = float(dD['mech'] @ dD['sec']/(dD['sec'] @ dD['sec']))
            zs = np.linalg.norm(np.asarray(d['inc']['sec']['s21']) - np.asarray(d['inc']['secB0']['s21']))
            zm = np.linalg.norm(np.asarray(d['inc']['mech']['s21']) - np.asarray(d['inc']['mechB0']['s21']))
            if v == 'both': okC &= abs(c - 1) <= 0.10
            out.append(f'{v:4s} {os.path.basename(f)[:-4]:26s} L.s {g["L0"] @ g["s20"]:+.3f} secB0/mechB0 L^ {r[0]:.4f} lepton {r[1]:.4f} heavy {r[2]:.4f}'
                       f' (target {target:.4f}) | c {c:.3f} | Zeeman lepton sec/mech {zs/zm:.4f}')
    out.append(f'H700a (both off: secB0/mechB0 = 1 within 2 %, L^ and lepton, every run): {ok["both"]}')
    out.append(f'H700b (528 off, Born on: = 9/8 within 2 %, every run): {ok["fs"]}')
    out.append(f'H700c (both off: c within 10 % of 1, every run): {okC}')
    text = '\n'.join(out); print(text)
    os.makedirs('out', exist_ok=True); open('out/ana700.txt', 'w').write(text + '\n')


if __name__ == '__main__':
    main()

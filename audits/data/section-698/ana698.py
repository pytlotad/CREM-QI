# Audit 698: the formation state of the 2P pair on the mechanical path (CREM_PAIR_MECH_TRACE, 2.5 ps, no reaction) against the secular trace
# of 695 (CREM_PAIR_QM1S_TRACE) from the same state -- D_k(t) = p_k(t) - p_k(form) of the 12 formation 2P eigenstates, PHEN 2.
# usage: python3 ana698.py <698 runs dir with mech/ ctl/>
import sys, glob, math, os
import numpy as np
R = sys.argv[1] if __name__ == "__main__" else None
PAIRS = (('antiproton', 'antihydrogen'), ('proton', 'hydrogen'))
FIELDS = ('1.0329', '0.1', '0.01')


def nums(t, key, n):
    i = t.index(key); return [float(x) for x in t[i+1:i+1+n]]


def parse(path):
    d = {'sec': [], 'mech': []}
    for line in open(path):
        t = line.split()
        if line.startswith('PAIR QM1S695 form'): d['form'] = np.array(nums(t, 'p', 12))
        elif line.startswith('PAIR QM1S695 t'): d['sec'].append([float(t[3])] + nums(t, 'p', 12))
        elif line.startswith('PAIR QM1S698 t'):
            d['mech'].append([float(t[3])] + nums(t, 'p', 12) + [float(t[t.index(k)+1]) for k in ('L', 'm1', 'm2')])
        elif line.startswith('PAIR QM1S698 end'): d['end698'] = line.strip()
        elif line.startswith('END'): d['end'] = (t[t.index('stop')+1], t[t.index('photons')+1])
    d['sec'] = np.array(d['sec']); d['mech'] = np.array(d['mech'])
    return d


def at(tr, t):   # D at time t by linear interpolation, 12 states
    return np.array([np.interp(t, tr[:, 0], tr[:, 1+k]) for k in range(12)])



def main():
    out = []
    fs = sorted(glob.glob(f'{R}/mech/*.err'))
    same = sum(''.join(l for l in open(f) if not l.startswith('PAIR QM1S698')) == open(f'{R}/ctl/{os.path.basename(f)}').read() for f in fs)
    out.append(f'H698a (files): with CREM_PAIR_MECH_TRACE without QM1S698 = without it {same}/{len(fs)}: {same == len(fs) == 12}')
    okB = True; okF = okG = True; okE = True
    for tag, name in PAIRS:
        for B in FIELDS:
            for s in ('0', '1'):
                d = parse(f'{R}/mech/{tag}_B{B}_A2_{s}.err'); f0 = d['form']; M = d['mech']; S = d['sec']
                m0 = M[0, 1:13] - f0; L = M[:, 13]; m1 = M[:, 14]; m2 = M[:, 15]
                b = float(np.max(np.abs(m0))) <= 1e-3 and abs(L[0] - 1.5) <= 1e-3
                okB &= b
                tmax = min(M[-1, 0], S[-1, 0])
                grid = [x for x in (1e-12, 2e-12) if x <= tmax + 1e-15]
                Dm = np.array([at(M, x) - f0 for x in grid]); Ds = np.array([at(S, x) - f0 for x in grid])
                rd = math.sqrt(((Dm - Ds)**2).mean()); rs = math.sqrt((Ds**2).mean())
                g2m = float(np.max(np.abs(Dm[-1]))); g2s = float(np.max(np.abs(Ds[-1])))
                mxm = float(np.max(np.abs(M[:, 1:13] - f0)))
                if B == '1.0329':
                    okF &= rd < 0.3*rs; okG &= g2m > 0.3*g2s; okE &= mxm > 0.02
                out.append(f'{name:12s} B {B:>6} T seed {s} | start |D_mech| {np.max(np.abs(m0)):.1e} L0 {L[0]:.6f} | L {L.min():.5f}-{L.max():.5f}'
                           f' |m1| {m1.min():.6f}-{m1.max():.6f} |m2| {m2.min():.6f}-{m2.max():.6f} | points mech {len(M)} sec {len(S)} to {tmax*1e12:.2f} ps'
                           f' | t = ' + ','.join(f'{x*1e12:.0f}' for x in grid) + f' ps: rms|D_mech - D_sec| {rd:.4f} rms|D_sec| {rs:.4f} (ratio {rd/rs if rs > 0 else float("nan"):.2f})'
                           f' | max|D(2 ps)| mech {g2m:.4f} sec {g2s:.4f} | max|D_mech| window {mxm:.4f} | {d.get("end698", "")[17:]} end {d.get("end")}')
                out.append('     D_mech(2 ps) ' + ' '.join(f'{x:+.4f}' for x in Dm[-1]))
                out.append('     D_sec (2 ps) ' + ' '.join(f'{x:+.4f}' for x in Ds[-1]))
    out.append(f'H698b (start = formation to 1e-3, L0 = 1.5 hbar): {okB}')
    out.append(f'H698f (1.0329 T: rms |D_mech - D_sec| < 0.3 rms |D_sec| at 1, 2 ps, every run): {okF}')
    out.append(f'H698g (1.0329 T: max |D_mech(2 ps)| > 0.3 max |D_sec(2 ps)|, every run): {okG}')
    out.append(f'H698e report (1.0329 T: max |D_mech| > 0.02 in the window, every run): {okE}')
    text = '\n'.join(out); print(text)
    os.makedirs('out', exist_ok=True); open('out/ana698.txt', 'w').write(text + '\n')


if __name__ == '__main__':
    main()

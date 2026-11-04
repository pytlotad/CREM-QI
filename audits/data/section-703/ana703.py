# Audit 703: the nutation lock CREM_WAVE_LOCK_J=1 (wave law on X = |(l + 1/2) L^ + s^|, rungs k = l, l + 1) against the default, same seeds.
# usage: python3 ana703.py <703 runs dir with on/ off/>
import sys, glob, os, math
import numpy as np
PAIRS = (('antiproton', 'antihydrogen'), ('proton', 'hydrogen'))
FIELDS = ('0.01', '0.1', '1.0329')
EV = 1.602176634e-19; HEV = 4.135667696e-15


def vec(t, key, n):
    i = t.index(key); return [float(x) for x in t[i+1:i+1+n]]


def parse(path):
    d = {'tr': []}
    for line in open(path):
        t = line.split()
        if line.startswith('PAIR JLOCK703 photon'): d['X'] = float(t[4]); d['cos'] = float(t[6])
        elif line.startswith('PAIR j projection jcl'): d['jcl'] = t[t.index('jcl')+1]
        elif line.startswith('PAIR QM1S695 form'): d['form'] = np.array([float(x) for x in t[17:29]])
        elif line.startswith('PAIR QM1S695 t'): d['tr'].append([float(x) for x in t[5:17]])
        elif line.startswith('PAIR QM1S684'): d['v3'] = vec(t, 'v3', 4)
        elif line.startswith('PAIR QM1S path'): d['qm'] = vec(t, 'qm', 4); d['sc'] = vec(t, 'sc', 4)
        elif line.startswith('CREM_LEDGER t'): d['Eg'] = float(t[t.index('Eg')+1])
        elif line.startswith('PAIR eigen photon'): d['jump'] = float(t[t.index('jump')+1])
        elif line.startswith('END'): d['end'] = (t[t.index('stop')+1], t[t.index('photons')+1])
    d['tr'] = np.array(d['tr'])
    return d


def main():
    R = sys.argv[1]; out = []; okB = True; same = tot = 0; okD = True
    for tag, name in PAIRS:
        for B in FIELDS:
            for s in ('0', '1', '2', '3'):
                key = f'{tag}_B{B}_A2_{s}.err'; a = parse(f'{R}/on/{key}'); b = parse(f'{R}/off/{key}')
                off = abs(a['X'] - round(a['X'])) if 'X' in a else float('nan')
                if B != '1.0329': okB &= off < 1e-3
                tot += 1; same += a.get('jcl') == b.get('jcl')
                Da = np.abs(a['tr'] - a['form']).max() if len(a['tr']) else float('nan')
                Db = np.abs(b['tr'] - b['form']).max() if len(b['tr']) else float('nan')
                if B == '0.01': okD &= Da < 0.3*Db
                Aa = np.array(a['v3']) - np.array(a['qm']); Ab = np.array(b['v3']) - np.array(b['qm'])
                dEg = (a['Eg'] - b['Eg'])/EV
                out.append(f'{name:12s} B {B:>6} T seed {s} | X at photon {a.get("X", float("nan")):.6f} (|X - k| {off:.1e}) cos {a.get("cos", float("nan")):+.4f}'
                           f' | jcl on {a.get("jcl")} off {b.get("jcl")} | max|D| on {Da:.4f} off {Db:.4f} | (A) on ' + ' '.join(f'{x:+.4f}' for x in Aa)
                           + ' off ' + ' '.join(f'{x:+.4f}' for x in Ab) + f' | photon on - off {dEg*1e6:+.3f} ueV ({dEg/HEV/1e9:+.4f} GHz) | end {a.get("end")} {b.get("end")}')
    out.append(f'H703b (0.01 and 0.1 T: |X - round(X)| < 1e-3 at the photon, every run): {okB}')
    out.append(f'H703c (jcl on = off in >= 95 % of runs): {same}/{tot}: {same >= 0.95*tot}')
    out.append(f'H703d (0.01 T: max|D| on < 0.3 x off, every run): {okD}')
    text = '\n'.join(out); print(text)
    os.makedirs('out', exist_ok=True); open('out/ana703.txt', 'w').write(text + '\n')


if __name__ == '__main__':
    main()

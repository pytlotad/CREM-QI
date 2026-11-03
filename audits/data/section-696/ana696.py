# Audit 696: the 695 traces (PHEN 2) with the wave lock off (CREM_WAVE_LOCK=0) -- does the lock damp or shift the 2P population nutation?
# usage: python3 ana696.py <696 runs dir with tr/> <695 runs dir with tr/>
import sys, glob, os
import numpy as np
R, R5 = sys.argv[1:3]


def parse(path):
    d = dict(t=[], p=[])
    for line in open(path):
        s = line.split()
        if line.startswith('PAIR QM1S695 form'): d['form'] = line; d['pf'] = np.array([float(x) for x in s[17:29]])
        elif line.startswith('PAIR QM1S695 t'): d['t'].append(float(s[3])); d['p'].append([float(x) for x in s[5:17]])
        elif line.startswith('PAIR QM1S695 photon'): d['ts'] = float(s[4])
        elif line.startswith('CREM_LEDGER t'): d['tau'] = float(s[2])
        elif line.startswith('END'): d['end'] = (s[s.index('stop')+1], s[s.index('photons')+1])
    d['t'] = np.array(d['t']); d['p'] = np.array(d['p']); return d


def grid(d, g): return np.array([np.interp(g, d['t'], d['p'][:, k] - d['pf'][k]) for k in range(12)]).T


out = []; okA = okB = okC = True
for f in sorted(glob.glob(f'{R}/tr/*.err')):
    b = os.path.basename(f); a = parse(f); z = parse(f'{R5}/tr/{b}')
    okA &= a.get('end') == ('4', '1') and a['form'] == z['form']
    m = {}
    for k, d in (('off', a), ('on', z)):
        g = np.arange(d['t'][0], d['t'][-1], 1e-12); Dg = grid(d, g); m[k] = (np.abs(Dg.mean(0)).max(), Dg.std(0).max(), np.abs(Dg).max())
    T = min(a['t'][-1], z['t'][-1]); g = np.arange(max(a['t'][0], z['t'][0]), T, 1e-12)
    dd = np.sqrt(((grid(a, g) - grid(z, g))**2).sum(1)); early = dd[g < g[0] + 50e-12]
    close = lambda x, y: abs(x - y) <= max(0.2*abs(y), 0.005)
    okB &= close(m['off'][0], m['on'][0]) and close(m['off'][1], m['on'][1]); okC &= np.sqrt((early**2).mean()) < 0.01
    out.append(f'{b[:-4]:26s} mean|D| off {m["off"][0]:.4f} on {m["on"][0]:.4f} | std off {m["off"][1]:.4f} on {m["on"][1]:.4f} | max|D| off {m["off"][2]:.4f} on {m["on"][2]:.4f}'
               f' | rms|D_off - D_on| first 50 ps {np.sqrt((early**2).mean()):.2e} all {np.sqrt((dd**2).mean()):.2e} ({T*1e9:.2f} ns)'
               f' | photon off {a["tau"]*1e12:.1f} on {z["tau"]*1e12:.1f} ps, spins off {a["ts"]*1e12:.1f} on {z["ts"]*1e12:.1f} ps')
out.append(f'H696a (12/12 stop 4 after 1 photon, same formation line as 695): {okA}')
out.append(f'H696b (time-mean and std of D without the lock within 20 % or 0.005 of 695, every run): {okB}')
out.append(f'H696c (rms |D_off - D_on| over the first 50 ps < 0.01, every run): {okC}')
text = '\n'.join(out); print(text)
os.makedirs('out', exist_ok=True); open('out/ana696.txt', 'w').write(text + '\n')

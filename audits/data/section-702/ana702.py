# Audit 702: the clock phase of 701 for p-Ps and o-Ps (n = 2, L = 1.5 hbar, no field) -- clock = lock; n_clock - n_E against the semiclassical Breit
# spin-independent shift (kinetic + Darwin orbit-orbit, two equal masses) and the spin part o-Ps - p-Ps against the measured 2^3P_J multiplet.
# Recipe (checked on Sommerfeld in 701): n_clock - n_E = -(n^3/(2B)) <H_pert>, <> the time average over the Kepler orbit of actions (n, k).
# usage: python3 ana702.py <702 runs dir with s512/ s2048/ s8192/ ctl/>
import sys, glob, os, math
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'section-701'))
from ana701 import derived, HBAR
C = 299792458.0; ME = 9.1093837015e-31; ALPHA = 7.2973525693e-3; ALPHAP = ALPHA*HBAR*C   # CODATA 2018; alpha' = e^2/(4 pi eps0)
HEV = 4.135667696e-15
A_MHZ, B_MHZ = 2742.9, 457.4   # from Hagena et al., PRL 71, 2887 (1993) (prereg702)


def parse(path):
    d = {'k': []}
    for line in open(path):
        t = line.split()
        if line.startswith('PAIR CLOCK701 start'):
            d['nlock'] = float(t[4]); d['Xlock'] = float(t[6]); d['binding'] = float(t[8])
            for key in ('L', 's1', 's2', 'peri'):
                i = t.index(key); d[key] = np.array([float(x) for x in t[i+1:i+4]])
        elif line.startswith('PAIR CLOCK701 k'):
            d['k'].append({t[i]: float(t[i+1]) for i in range(2, len(t) - 1, 2)})
    return d


def kepler_avg(f, n, k, mu, B):
    """time average of f(r, p2, pr2) over the Kepler orbit of energy -B/n^2 and L = k hbar (reduced mass mu)"""
    a_mu = HBAR**2/(mu*ALPHAP); a = n*n*a_mu; e = math.sqrt(max(0.0, 1 - (k/n)**2)); E = -B/(n*n); L = k*HBAR
    u = (np.arange(20000) + 0.5)*2*math.pi/20000; r = a*(1 - e*np.cos(u)); w = r/r.sum()
    p2 = 2*mu*(E + ALPHAP/r); pr2 = p2 - L*L/(r*r)
    return float((f(r, p2, pr2)*w).sum())


def main():
    R = sys.argv[1]; out = []
    fs = sorted(glob.glob(f'{R}/s512/*.err'))
    cf = sorted(glob.glob(f'{R}/s128/*.err'))
    same = sum(''.join(l for l in open(f) if not l.startswith('PAIR CLOCK701')) == open(f'{R}/ctl/{os.path.basename(f)}').read() for f in cf)
    out.append(f'H702a (files): default step with CREM_PAIR_CLOCK without CLOCK701 = without it {same}/{len(cf)}: {same == len(cf) == 8}')
    # self-test of the recipe on the test particle: Sommerfeld
    mu_t = ME; B_t = mu_t*ALPHAP**2/(2*HBAR**2)
    som = -(8/(2*B_t))*kepler_avg(lambda r, p2, pr2: -p2*p2/(8*mu_t**3*C*C), 2, 1.5, mu_t, B_t)
    out.append(f'recipe self-test (test particle, n 2, k 1.5): {som:.6e} vs alpha^2 (1/2k - 3/8n) {ALPHA**2*(1/3 - 3/16):.6e}')
    res = {}
    for f in fs:
        key = os.path.basename(f); P = key.split('_A')[1][0]; s = key[:-4].split('_')[-1]
        v = {}; dd = {}
        for st in ('512', '2048', '8192'):
            d = parse(f'{R}/s{st}/{key}'); rows = derived(d)
            v[st] = np.mean([r['ncl'] - r['nE'] for r in rows])
            if st == '8192':
                rows8 = rows; d8 = d
        r1, r2 = v['512'] - v['2048'], v['2048'] - v['8192']
        p = math.log(r1/r2)/math.log(4.0) if r1*r2 > 0 else float('nan')
        x0 = v['8192'] - r2/(4.0**p - 1.0) if math.isfinite(p) else float('nan')
        n = np.mean([r['nE'] for r in rows8]); k = np.mean([r['X'] for r in rows8]); B = d8['binding']; mu = ME/2
        Hkin = kepler_avg(lambda r, p2, pr2: -p2*p2/(4*ME**3*C*C), n, k, mu, B)
        Hoo = kepler_avg(lambda r, p2, pr2: -(ALPHAP/(2*ME*ME*C*C*r))*(p2 + pr2), n, k, mu, B)
        pred = -(n**3/(2*B))*(Hkin + Hoo)
        Udd = np.mean([q['Udd'] for q in d8['k']]); Uqd = np.mean([q['Uqd'] for q in d8['k']])
        dn_dd = -(n**3/(2*B))*(Udd + Uqd)
        cl = max(max(abs(r['ncl'] - d8['nlock']) for r in rows8), max(abs(r['X'] - d8['Xlock']) for r in rows8))
        L, s1, s2 = d8['L'], d8['s1'], d8['s2']
        if P == '2':
            LS = 1.5*float(L @ (0.5*s1 + 0.5*s2))                     # (L.S)/hbar^2: |L| = 1.5 hbar, each |s| = hbar/2
            S12 = 0.5*float(s1 @ s2) - 1.5*float(s1 @ L)*float(s2 @ L)
            Espin = (A_MHZ*LS + B_MHZ*S12)*1e6*HEV
        else:
            LS = S12 = 0.0; Espin = 0.0                                # QM singlet: no a, b terms
        res[(P, s)] = dict(x0=x0, p=p, pred=pred, dn_dd=dn_dd, cl=cl, Espin=Espin, LS=LS, S12=S12, n=n, B=B, Udd=Udd, Uqd=Uqd)
        out.append(f'{"p-Ps" if P == "1" else "o-Ps"} seed {s}: n_clock - n_E 512 {v["512"]:+.4e} 2048 {v["2048"]:+.4e} 8192 {v["8192"]:+.4e} | p {p:.3f} x0 {x0:+.4e}'
                   f' | max |clock - lock| {cl:.1e} | Breit kin+oo {pred:+.4e} | classical spin terms <Udd> {Udd/1.602176634e-19*1e6:+.3f} ueV <Uqd> {Uqd/1.602176634e-19*1e6:+.3f} ueV'
                   f' -> dn {dn_dd:+.4e} | L.S {LS:+.3f} S12 {S12:+.3f}')
    okB = all(r['cl'] < 1e-3 for r in res.values()); okC = all(abs(r['p'] - 2) < 0.2 for r in res.values())
    okD = all(abs((r['x0'] - r['dn_dd'])/r['pred'] - 1) <= 0.2 for (P, s), r in res.items() if P == '1')
    out.append('p-Ps: (x0 - dn_dd) / Breit kin+oo: ' + ' '.join(f'seed {s} {(r["x0"] - r["dn_dd"])/r["pred"]:.3f}' for (P, s), r in sorted(res.items()) if P == '1'))
    okE = True
    for s in ('0', '1', '2', '3'):
        o, pp = res[('2', s)], res[('1', s)]
        dE = -(2*o['B']/o['n']**3)*(o['x0'] - pp['x0'])/1.602176634e-19       # eV
        dEm = o['Espin'] - pp['Espin']
        dEm_model = -(2*o['B']/o['n']**3)*(o['dn_dd'] - pp['dn_dd'])/1.602176634e-19
        okE &= abs(dE/dEm - 1) <= 0.3 if dEm != 0 else False
        out.append(f'seed {s}: o-Ps - p-Ps from the clock {dE*1e6:+.3f} ueV ({dE/HEV/1e9:+.3f} GHz) | from the classical spin energies <Udd + Uqd> {dEm_model*1e6:+.3f} ueV'
                   f' | measured a, b at the classical vectors {dEm*1e6:+.3f} ueV ({dEm/HEV/1e9:+.3f} GHz) | ratio clock / measured {dE/dEm if dEm else float("nan"):.3f}')
    out.append(f'H702b (|clock - lock| < 1e-3, p-Ps and o-Ps): {okB}')
    out.append(f'H702c (order p within 0.2 of 2): {okC}')
    out.append(f'H702d (p-Ps: x0 - dn_dd = Breit kin+oo within 20 %, every seed): {okD}')
    out.append(f'H702e (o-Ps - p-Ps = measured a, b within 30 %, every seed): {okE}')
    text = '\n'.join(out); print(text)
    os.makedirs('out', exist_ok=True); open('out/ana702.txt', 'w').write(text + '\n')


if __name__ == '__main__':
    main()

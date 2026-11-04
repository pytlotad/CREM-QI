# Audit 701: the clock's action per radial period on the mechanical path against the level n and L/hbar of the wave lock.
# J_r = J - L_kan dphi, n_clock = J_r/h + L_kan/hbar, n_E = sqrt(binding/|E|); Sommerfeld (test particle, Coulomb):
# n_clock - n_E = alpha^2 (1/(2k) - 3/(8n)).
# usage: python3 ana701.py <701 runs dir with s128/ s512/ s2048/ ctl/>
import sys, glob, os, math
import numpy as np
HBAR = 1.054571817e-34; H = 2*math.pi*HBAR; ALPHA = 7.2973525693e-3   # CODATA 2018
PAIRS = (('antiproton', 'antihydrogen'), ('proton', 'hydrogen'))


def parse(path):
    d = {'k': []}
    for line in open(path):
        t = line.split()
        if line.startswith('PAIR CLOCK701 start'):
            d['nlock'] = float(t[4]); d['Xlock'] = float(t[6]); d['binding'] = float(t[8])
        elif line.startswith('PAIR CLOCK701 k'):
            v = {t[i]: float(t[i+1]) for i in range(2, len(t) - 1, 2)}; d['k'].append(v)
    return d


def derived(d):
    rows = []
    for v in d['k']:
        Jr = v['J'] - v['Lcan']*v['dphi']; ncl = Jr/H + v['Lcan']/HBAR; nE = math.sqrt(d['binding']/abs(v['E']))
        full = (v['J'] - v['intE'])/HBAR
        rows.append(dict(ncl=ncl, nE=nE, X=v['Lcan']/HBAR, nr=Jr/H, prec=v['dphi'] - 2*math.pi, full=math.fmod(full, 2*math.pi)))
    return rows


def main():
    R = sys.argv[1]; out = []
    fs = sorted(glob.glob(f'{R}/s128/*.err'))
    same = sum(''.join(l for l in open(f) if not l.startswith('PAIR CLOCK701')) == open(f'{R}/ctl/{os.path.basename(f)}').read() for f in fs)
    out.append(f'H701a (files): default step with CREM_PAIR_CLOCK without CLOCK701 = without it {same}/{len(fs)}: {same == len(fs) == 4}')
    okB = okC = True; res = {}
    for s in ('0', '1'):
        for tag, name in PAIRS:
            key = f'{tag}_B0.01_A2_{s}.err'; m = {}
            for st in ('128', '512', '2048'):
                d = parse(f'{R}/s{st}/{key}'); rows = derived(d)
                ncl = np.array([r['ncl'] for r in rows]); nE = np.array([r['nE'] for r in rows]); X = np.array([r['X'] for r in rows])
                m[st] = dict(ncl=ncl.mean(), sd=ncl.std(ddof=1), dif=(ncl - nE).mean(), difsd=(ncl - nE).std(ddof=1), nE=nE.mean(), X=X.mean(),
                             nr=np.mean([r['nr'] for r in rows]), prec=np.mean([r['prec'] for r in rows]), n=len(rows),
                             nlock=d['nlock'], Xlock=d['Xlock'], full=[r['full'] for r in rows[:3]])
                okB &= max(abs(ncl - d['nlock']).max(), abs(X - d['Xlock']).max()) < 1e-3
                out.append(f'{name:12s} seed {s} steps {st:>4} periods {len(rows)} | n_clock {m[st]["ncl"]:.9f} (sd {m[st]["sd"]:.1e}) n_E {m[st]["nE"]:.9f} nlock {d["nlock"]:.9f}'
                           f' | n_clock - n_E {m[st]["dif"]:+.3e} (sd {m[st]["difsd"]:.1e}) | L/hbar {m[st]["X"]:.9f} Xlock {d["Xlock"]:.9f} | J_r/h {m[st]["nr"]:.9f}'
                           f' | precession/period {m[st]["prec"]:+.3e} rad | full phase mod 2pi ' + ' '.join(f'{x:+.3f}' for x in m[st]['full']))
            c512, c2048 = m['512'], m['2048']
            okC &= abs(c512['dif'] - c2048['dif']) < 1e-6 and c2048['difsd'] < 1e-6
            rich = c2048['dif'] + (c2048['dif'] - c512['dif'])/15.0      # h^2: (16 x - y)/15
            k, n = c2048['X'], c2048['nE']; pred = ALPHA**2*(1/(2*k) - 3/(8*n))
            res[(s, name)] = (rich, pred)
            out.append(f'   Richardson (512, 2048) n_clock - n_E {rich:+.4e}; Sommerfeld alpha^2 (1/(2k) - 3/(8n)) {pred:+.4e} (k {k:.6f}, n {n:.6f})')
    okD = True
    for s in ('0', '1'):
        (ra, pa), (rh, ph) = res[(s, 'antihydrogen')], res[(s, 'hydrogen')]
        mean = 0.5*(ra + rh); pred = 0.5*(pa + ph); okD &= abs(mean/pred - 1) <= 0.2
        diff = rh - ra
        dEdn = 2*13.598/8                                              # eV per unit n at n = 2 (hydrogen, reduced mass)
        e = diff*dEdn
        out.append(f'seed {s}: mean (H, Hbar) n_clock - n_E {mean:+.4e} vs Sommerfeld {pred:+.4e} (ratio {mean/pred:.3f}) | H - Hbar {diff:+.4e}'
                   f' -> {e*1e6:+.2f} ueV = {e/4.135667696e-15/1e9:+.3f} GHz (dE/dn = 2 Ry/n^3 = {dEdn:.4f} eV)')
    out.append(f'H701b (|n_clock - nlock|, |L/hbar - Xlock| < 1e-3 every period): {okB}')
    out.append(f'H701c (512 vs 2048 difference < 1e-6, scatter at 2048 < 1e-6): {okC}')
    out.append(f'H701d (mean of H and Hbar = Sommerfeld within 20 %, both seeds): {okD if okC else "not decidable (H701c failed)"}')
    text = '\n'.join(out); print(text)
    os.makedirs('out', exist_ok=True); open('out/ana701.txt', 'w').write(text + '\n')


if __name__ == '__main__':
    main()

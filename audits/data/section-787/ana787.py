# Audit 787: deterministic emission threshold from the orbital phase on the mechanical path -- verdicts H787e, u, x, t, m, d -> out/ana787.txt
import os, glob, numpy as np
from scipy import stats
D = os.path.dirname(os.path.abspath(__file__)); O = os.path.join(D, 'out', 'runs')
pr = lambda *a: print(*a, flush=True)
SCALE = 8e5; MP = 1.67262192595e-27; ME = 9.1093837139e-31
TAU = {'H': 1.59572e-9, 'pPs': 1.59572e-9 * 2 * MP / (MP + ME), 'oPs': 1.59572e-9 * 2 * MP / (MP + ME)}   # engine rules (738), Ps by mu
def photons(fn):
    ph = []
    for l in open(fn):
        if l.startswith('PHASE787 photon'):
            p = l.split(); d = {p[i]: p[i + 1] for i in range(2, len(p) - 1, 2)}
            ph.append(dict(t=float(d['t']), U=float(d['U']) if 'U' in d else np.nan))
    return ph
# H787e
okE = all(open(os.path.join(O, 'F', f'head_{i}.{x}')).read() == open(os.path.join(O, 'F', f'new_{i}.{x}')).read() for i in range(10) for x in ('out', 'err'))
pr(f"H787e (bez przelacznikow: sonda HEAD = sonda 787, wodor, 10 ziaren, stdout i stderr bajt w bajt): {'TAK' if okE else 'NIE'}")
# H787d
okD = all(open(os.path.join(O, 'R', f'H_{i}.err')).read() == open(os.path.join(O, 'P', f'H_{i}.err')).read() for i in range(10))
pr(f"H787d (powtorzenie P, wodor, 10 ziaren: bit w bit): {'TAK' if okD else 'NIE'}")
res = {}
for v, n in (('P', 100), ('D', 30)):   # poprawka 2 przed ocena
    for s in ('H', 'pPs', 'oPs'):
        rows = []; cens = 0
        for i in range(n):
            ph = photons(os.path.join(O, v, f'{s}_{i}.err'))
            if len(ph) < 2: cens += 1; continue
            rows.append((ph[0]['t'], ph[1]['t'] - ph[0]['t'], ph[0]['U']))
        res[(v, s)] = (np.array(rows), cens)
okU = okX = okT = okM = True
for s in ('H', 'pPs', 'oPs'):
    pr(f"=== {s} ===")
    for v in ('P', 'D'):
        a, cens = res[(v, s)]; N = len(a); t1, dt, U = a[:, 0], a[:, 1], a[:, 2]
        tau = dt.mean(); x = dt / tau
        pU = stats.kstest(U, 'uniform').pvalue; pX = stats.kstest(x, 'expon').pvalue
        fs = np.mean(x < 0.1); f0 = 1 - np.exp(-0.1); sf = np.sqrt(f0 * (1 - f0) / N)
        tphys = tau * SCALE; st = tphys / np.sqrt(N); dev = (tphys - TAU[s]) / st
        rho = stats.spearmanr(t1, dt).correlation
        pr(f"  {v}: N {N} (cenzura {cens}); KS U p {pU:.3f}; KS Dt/tau ~ Exp(1) p {pX:.3f}; Dt < 0,1 tau {fs:.3f} (oczek. {f0:.3f} +- {sf:.3f}); "
           f"tau fiz. {tphys * 1e9:.4f} +- {st * 1e9:.4f} ns (regula {TAU[s] * 1e9:.4f}, {dev:+.2f} sigma); Spearman(t1, Dt) {rho:+.3f} (prog {3 / np.sqrt(N):.3f})")
        if v == 'P':
            okU &= pU > 0.01; okX &= pX > 0.01 and abs(fs - f0) <= 3 * sf; okT &= abs(dev) <= 3; okM &= abs(rho) < 3 / np.sqrt(N)
pr(f"H787u (U w chwili fotonu 3d -> 2p jednostajne, KS p > 0,01, kazdy uklad): {'TAK' if okU else 'NIE'}")
pr(f"H787x (Dt/tau ~ Exp(1), KS p > 0,01 i krotkie czasy w 3 sigma): {'TAK' if okX else 'NIE'}")
pr(f"H787t (tau x skala = regula w 3 sigma): {'TAK' if okT else 'NIE'}")
pr(f"H787m (|Spearman(t1, Dt)| < 3/sqrt(N)): {'TAK' if okM else 'NIE'}")
pr("Pomiar: wodor 2p 1,600(4) ns (Bickel i Goodman 1966); Ps 2P ~3,19 ns (teoria; p-Ps i o-Ps osobno).")

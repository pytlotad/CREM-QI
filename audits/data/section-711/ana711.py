# Audit 711: independent confirmation of the 710 ranking -- Hamilton residual on the velocity-scaling family (CREM_PAIR_CLOCK_VSCALE), spin part
# hydrogen - antihydrogen (seeds 2-3) and o-Ps - p-Ps (seeds 0-1); a candidate wins a seed when its |spin part| is the smallest and >= 5x below the
# next.  At eps = 0: Ps o - p with the winner against 8/9 (a L.S + b S12) - U_dd(p) (Hagena et al. 1993), as 706/707.
# usage: python3 ana711.py <711 runs dir with ret/ ins/>
import sys, os, math
import numpy as np
here = os.path.dirname(os.path.abspath(__file__))
for sec in ('section-710', 'section-709', 'section-702'): sys.path.insert(0, os.path.join(here, '..', sec))
from ana709 import means, CAND, HBAR, H, EV
from ana702 import parse, A_MHZ, B_MHZ, HEV

SYS = {'hydrogen': (('2', '3'), lambda s: (f'proton_{s}', f'antiproton_{s}'), 'H - Hbar'),
       'Ps': (('0', '1'), lambda s: (f'ps_A2_{s}', f'ps_A1_{s}'), 'o - p')}
CACHE = {}


def P(path):
    if path not in CACHE: CACHE[path] = parse(path)
    return CACHE[path]


def residual(D, name, cand, de):
    m, p = means(P(f'{D}/{name}_e-{de}.err'), cand), means(P(f'{D}/{name}_e{de}.err'), cand)
    m0, p0 = means(P(f'{D}/{name}_e-{de}.err'), 'K0 J'), means(P(f'{D}/{name}_e{de}.err'), 'K0 J')
    T = 0.5*(m['T'] + p['T']); dphi = 0.5*(m['dphi'] + p['dphi'])
    return (p['Jr'] - m['Jr'] - T*(p['E'] - m['E']) + dphi*(p['L'] - m['L']))/abs(dphi*(p0['L'] - m0['L']))


def rich(D, name, cand):
    r1 = residual(D, name, cand, '1e-3'); r2 = residual(D, name, cand, '2e-3')
    return (4*r1 - r2)/3


def main():
    R = sys.argv[1]; out = []; win = {}; spread = {}
    for sumtag, D in (('instantaneous (weight 0)', 'ins'), ('retarded (weight 1)', 'ret')):
        for sysname, (seeds, names, lab) in SYS.items():
            for s in seeds:
                a, b = names(s); sp = {}
                for cand in CAND:
                    ra, rb = rich(f'{R}/{D}', a, cand), rich(f'{R}/{D}', b, cand); sp[cand] = ra - rb
                    out.append(f'{sumtag:24s} {sysname:8s} seed {s} {cand:12s}: R/|dphi dL| {a} {ra:+.3e}, {b} {rb:+.3e} | spin part ({lab}) {sp[cand]:+.3e}')
                order = sorted(CAND, key=lambda c: abs(sp[c])); ratio = abs(sp[order[1]])/max(abs(sp[order[0]]), 1e-300)
                w = order[0] if ratio >= 5 else None; win[(D, sysname, s)] = w
                spread[(D, sysname, s)] = {c: sp[c] - sp['K0 J'] for c in CAND}
                out.append(f'   -> smallest {order[0]}, next {order[1]} ({ratio:.1f}x): winner {w}')
    okA = all(win[('ins', y, s)] == 'K3 J+Jq+JTh' for y, (ss, _, _) in SYS.items() for s in ss)
    okB = all(win[('ret', 'hydrogen', s)] == 'K2 J+JTh' for s in SYS['hydrogen'][0])
    okC = all(win[('ret', 'Ps', s)] == 'K2 J+JTh' for s in SYS['Ps'][0])
    okD = True
    for D in ('ins', 'ret'):
        for s in SYS['hydrogen'][0]:
            d = spread[(D, 'hydrogen', s)]   # 710: K2 - K0 = +1.0e-5 (JTh), K3 - K2 = -2.0e-5 (Jq) on both sums
            jq = d['K3 J+Jq+JTh'] - d['K2 J+JTh']; jth = d['K2 J+JTh'] - d['K0 J']
            okD &= 0.5 <= abs(jq)/2.0e-5 <= 2 and 0.5 <= abs(jth)/1.0e-5 <= 2
            out.append(f'{D} hydrogen seed {s}: candidate differences Jq term (K3 - K2) {jq:+.3e}, JTh term (K2 - K0) {jth:+.3e} (710: -2e-5, +1e-5)')
    out.append(f'H711a (instantaneous: K3 wins in both hydrogen and both Ps seeds): {okA}')
    out.append(f'H711b (retarded hydrogen seeds 2-3: K2 wins): {okB}')
    out.append(f'H711c (retarded Ps seeds 0-1: K2 wins): {okC}')
    out.append(f'H711d (hydrogen candidate differences as 710 within x2): {okD}')
    okE = True
    for D, cand in (('ret', 'K2 J+JTh'), ('ins', 'K3 J+Jq+JTh')):
        for s in SYS['Ps'][0]:
            o, p = P(f'{R}/{D}/ps_A2_{s}_e0.err'), P(f'{R}/{D}/ps_A1_{s}_e0.err')
            mo, mp = means(o, cand), means(p, cand)
            L, s1, s2 = o['L'], o['s1'], o['s2']
            LS = 1.5*float(L @ (0.5*s1 + 0.5*s2)); S12 = 0.5*float(s1 @ s2) - 1.5*float(s1 @ L)*float(s2 @ L)
            Udd_p = np.mean([k['Udd'] for k in p['k']])/EV
            pred = 8/9*(A_MHZ*LS + B_MHZ*S12)*1e6*HEV - Udd_p
            n = math.sqrt(o['binding']/abs(mo['E'])); toE = lambda dn: -(2*o['binding']/n**3)*dn/EV
            e = toE(mo['n'] - mp['n']); okE &= abs(e/pred - 1) <= 0.05
            eo, ep = toE(mo['n']), toE(mp['n'])
            out.append(f'{D} Ps seed {s} {cand}: L.S {LS:+.4f} S12 {S12:+.4f} | o-Ps {eo*1e6:+.4f} ueV, p-Ps {ep*1e6:+.4f} ueV (n_clock - n_E as energy) | o - p'
                       f' {e*1e6:+.4f} ueV = {e/HEV/1e6:+.2f} MHz vs 8/9 (a L.S + b S12) - U_dd(p) {pred*1e6:+.4f} ueV (ratio {e/pred:.4f})')
    out.append(f'H711e (Ps o - p with the sum\'s winner = 8/9 multiplet - U_dd(p) within 5 %, both sums, both seeds): {okE}')
    text = '\n'.join(out); print(text)
    os.makedirs('out', exist_ok=True); open('out/ana711.txt', 'w').write(text + '\n')


if __name__ == '__main__':
    main()

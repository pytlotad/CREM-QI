# Audit 709: which action is the adiabatic invariant of each sum, paired with its conserved energy (instantaneous: charge-dipole weight 0,
# retarded: weight 1, audit 708).  Hamilton's relation on a family of orbits: dJ_r = T dE - dphi dL, residual R = Delta J_r - <T> Delta E + <dphi> Delta L,
# family by scaling the clock's start (CREM_PAIR_CLOCK_SCALE = eps, central differences +-1e-3, +-2e-3, Richardson); spin part = hydrogen - antihydrogen.
# usage: python3 ana709.py <709 runs dir with ret/ ins/ ctl/>
import sys, os, re, math
import numpy as np
HBAR = 1.054571817e-34; H = 2*math.pi*HBAR; EV = 1.602176634e-19; HEV = 4.135667696e-15
DE_H = 10969.045e6*HEV   # eV, 2P3/2 - 2P1/2 (Hagley & Pipkin 1994 + Lundeen & Pipkin 1981)
CAND = {'K0 J': (('J',), ('Lcan',)), 'K1 J+Jh+JT': (('J', 'Jh', 'JT'), ('Lcan', 'Lh', 'LT')),
        'K2 J+JTh': (('J', 'JTh'), ('Lcan', 'LTh')), 'K3 J+Jq+JTh': (('J', 'Jq', 'JTh'), ('Lcan', 'Lq', 'LTh'))}


def parse(path):
    d = {'k': []}
    for line in open(path):
        t = line.split()
        if line.startswith('PAIR CLOCK701 start'):
            d['binding'] = float(t[8]); d['L'] = np.array([float(x) for x in t[10:13]]); d['s2'] = np.array([float(x) for x in t[18:21]])
        elif line.startswith('PAIR CLOCK701 k'):
            d['k'].append({t[i]: float(t[i+1]) for i in range(2, len(t) - 1, 2)})
    return d


def means(d, cand):
    Jk, Lk = CAND[cand]; rows = d['k']
    Jr = [sum(v[a] for a in Jk) - sum(v[a] for a in Lk)*v['dphi'] for v in rows]
    L = [sum(v[a] for a in Lk) for v in rows]
    return dict(Jr=np.mean(Jr), L=np.mean(L), E=np.mean([v['E'] for v in rows]), T=np.mean([v['T'] for v in rows]),
                dphi=np.mean([v['dphi'] for v in rows]), n=np.mean([(jr/H + l/HBAR) - math.sqrt(d['binding']/abs(v['E'])) for jr, l, v in zip(Jr, L, rows)]))


def residual(R, tag, s, cand, de):
    m, p = means(parse(f'{R}/{tag}_{s}_e-{de}.err'), cand), means(parse(f'{R}/{tag}_{s}_e{de}.err'), cand)
    T = 0.5*(m['T'] + p['T']); dphi = 0.5*(m['dphi'] + p['dphi']); dE = p['E'] - m['E']
    return (p['Jr'] - m['Jr'] - T*dE + dphi*(p['L'] - m['L']))/abs(T*dE)


def strip(path):
    return ''.join(re.sub(r' Jq \S+ JTh \S+', '', re.sub(r' Lq \S+ LTh \S+', '', l)) for l in open(path))


def main():
    R = sys.argv[1]; out = []
    a = strip(f'{R}/ret/proton_0_e0.err') == open(f'{R}/ctl/ret_proton_0.err').read()
    b = re.sub(r' Lq \S+ LTh \S+', '', open(f'{R}/ins/proton_0_e0.err').read()) == open(f'{R}/ctl/ins_proton_0.err').read()
    out.append(f'H709a (files): eps = 0 with the switch, new fields removed = without: retarded {a}, instantaneous {b}')
    ok = {}
    for sumtag, D in (('instantaneous (weight 0)', 'ins'), ('retarded (weight 1)', 'ret')):
        for cand in CAND:
            spin = []
            for s in ('0', '1'):
                r = {}
                for tag in ('proton', 'antiproton'):
                    r1 = residual(f'{R}/{D}', tag, s, cand, '1e-3'); r2 = residual(f'{R}/{D}', tag, s, cand, '2e-3')
                    r[tag] = (r1, r2, (4*r1 - r2)/3)
                sp = r['proton'][2] - r['antiproton'][2]; spin.append(sp)
                out.append(f'{sumtag:24s} {cand:12s} seed {s}: R/|T dE| H {r["proton"][0]:+.3e} / {r["proton"][1]:+.3e} -> {r["proton"][2]:+.3e},'
                           f' Hbar {r["antiproton"][0]:+.3e} / {r["antiproton"][1]:+.3e} -> {r["antiproton"][2]:+.3e} | spin part (H - Hbar) {sp:+.3e}'
                           f' | mean {0.5*(r["proton"][2] + r["antiproton"][2]):+.3e}')
            ok[(D, cand)] = all(abs(x) < 1e-7 for x in spin)
    out.append('criterion |spin part| < 1e-7, both seeds: ' + ', '.join(f'{D} {c}: {v}' for (D, c), v in ok.items()))
    okB = ok[('ins', 'K3 J+Jq+JTh')] and not ok[('ins', 'K0 J')]
    winners = [c for c in CAND if ok[('ret', c)]]
    out.append(f'H709b (instantaneous: K3 passes, K0 fails): {okB}')
    out.append(f'H709c (retarded: exactly one candidate passes; predicted K2): {len(winners) == 1} -> {winners}')
    for D in ('ins', 'ret'):
        for cand in CAND:
            rs = []
            for s in ('0', '1'):
                h = parse(f'{R}/{D}/proton_{s}_e0.err'); a_ = parse(f'{R}/{D}/antiproton_{s}_e0.err')
                mh, ma = means(h, cand), means(a_, cand); cos = float(h['L'] @ h['s2']); n = math.sqrt(h['binding']/abs(mh['E']))
                rs.append(-(2*h['binding']/n**3)*(mh['n'] - ma['n'])/EV/(DE_H*cos))
            out.append(f'{D} {cand:12s} (H - Hbar)/(Delta E cos) at eps 0: ' + ' / '.join(f'{x:.4f}' for x in rs) + ' (BP x 8/9 = 0.8889)')
    okD = None
    if len(winners) == 1:
        okD = True
        for s in ('0', '1'):
            h = parse(f'{R}/ret/proton_{s}_e0.err'); a_ = parse(f'{R}/ret/antiproton_{s}_e0.err')
            mh, ma = means(h, winners[0]), means(a_, winners[0]); cos = float(h['L'] @ h['s2']); n = math.sqrt(h['binding']/abs(mh['E']))
            okD &= abs(-(2*h['binding']/n**3)*(mh['n'] - ma['n'])/EV/(DE_H*cos)/(8/9) - 1) <= 0.03
    out.append(f'H709d (the retarded winner with weight 1 = BP x 8/9 within 3 %): {okD}')
    text = '\n'.join(out); print(text)
    os.makedirs('out', exist_ok=True); open('out/ana709.txt', 'w').write(text + '\n')


if __name__ == '__main__':
    main()

# Audit 710: Hamilton's relation dJ_r = T dE - dphi dL on a family that moves L (velocities about the centre of mass x (1 + eps),
# CREM_PAIR_CLOCK_VSCALE), residual R = Delta J_r - <T> Delta E + <dphi> Delta L in units of |<dphi> Delta L|, spin part = hydrogen - antihydrogen;
# candidates and parsing as in audit 709.
# usage: python3 ana710.py <710 runs dir with ret/ ins/ ctl/> <709 runs dir>
import sys, os, math
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'section-709'))
from ana709 import parse, means, CAND, HBAR, EV, DE_H


def residual(R, tag, s, cand, de):
    m, p = means(parse(f'{R}/{tag}_{s}_e-{de}.err'), cand), means(parse(f'{R}/{tag}_{s}_e{de}.err'), cand)
    m0, p0 = means(parse(f'{R}/{tag}_{s}_e-{de}.err'), 'K0 J'), means(parse(f'{R}/{tag}_{s}_e{de}.err'), 'K0 J')
    T = 0.5*(m['T'] + p['T']); dphi = 0.5*(m['dphi'] + p['dphi'])
    return (p['Jr'] - m['Jr'] - T*(p['E'] - m['E']) + dphi*(p['L'] - m['L']))/abs(dphi*(p0['L'] - m0['L']))


def main():
    R, R9 = sys.argv[1:3]; out = []
    same = [open(f'{R}/ctl/{D}_proton_0.err').read() == open(f'{R9}/{D}/proton_0_e0.err').read() for D in ('ret', 'ins')]
    out.append(f'H710a (files): CREM_PAIR_CLOCK_VSCALE=0 = 709 eps = 0 byte for byte: retarded {same[0]}, instantaneous {same[1]}')
    spin = {}
    for sumtag, D in (('instantaneous (weight 0)', 'ins'), ('retarded (weight 1)', 'ret')):
        for cand in CAND:
            sp = []
            for s in ('0', '1'):
                r = {}
                for tag in ('proton', 'antiproton'):
                    r1 = residual(f'{R}/{D}', tag, s, cand, '1e-3'); r2 = residual(f'{R}/{D}', tag, s, cand, '2e-3')
                    r[tag] = (r1, r2, (4*r1 - r2)/3)
                sp.append(r['proton'][2] - r['antiproton'][2])
                out.append(f'{sumtag:24s} {cand:12s} seed {s}: R/|dphi dL| H {r["proton"][0]:+.3e} / {r["proton"][1]:+.3e} -> {r["proton"][2]:+.3e},'
                           f' Hbar {r["antiproton"][0]:+.3e} / {r["antiproton"][1]:+.3e} -> {r["antiproton"][2]:+.3e} | spin part (H - Hbar) {sp[-1]:+.3e}')
            spin[(D, cand)] = sp
    pas = {k: all(abs(x) < 1e-6 for x in v) for k, v in spin.items()}
    fail = {k: all(abs(x) > 1e-5 for x in v) for k, v in spin.items()}
    out.append('pass (< 1e-6 both seeds) / fail (> 1e-5 both seeds): ' + ', '.join(f'{D} {c}: {pas[(D, c)]}/{fail[(D, c)]}' for (D, c) in spin))
    okB = pas[('ins', 'K3 J+Jq+JTh')] and fail[('ins', 'K0 J')]
    winners = [c for c in CAND if pas[('ret', c)]]
    okC = len(winners) == 1 and all(fail[('ret', c)] for c in CAND if c not in winners)
    out.append(f'H710b (instantaneous: K3 < 1e-6, K0 > 1e-5): {okB}')
    out.append(f'H710c (retarded: exactly one < 1e-6, the others > 1e-5; predicted K2): {okC} -> {winners}')
    # precession: spin part of dphi along the family (the quantity the L derivative tests)
    for D in ('ins', 'ret'):
        for s in ('0', '1'):
            row = []
            for e in ('-2e-3', '-1e-3', '1e-3', '2e-3'):
                h = means(parse(f'{R}/{D}/proton_{s}_e{e}.err'), 'K0 J'); a = means(parse(f'{R}/{D}/antiproton_{s}_e{e}.err'), 'K0 J')
                row.append(f'{e} {h["dphi"] - a["dphi"]:+.4e}')
            out.append(f'{D} seed {s}: periapsis precession per period H - Hbar along the family: ' + ', '.join(row) + ' rad')
    okD = None
    if len(winners) == 1:
        okD = True
        for s in ('0', '1'):
            h = parse(f'{R9}/ret/proton_{s}_e0.err'); a_ = parse(f'{R9}/ret/antiproton_{s}_e0.err')
            mh, ma = means(h, winners[0]), means(a_, winners[0]); cos = float(h['L'] @ h['s2']); n = math.sqrt(h['binding']/abs(mh['E']))
            v = -(2*h['binding']/n**3)*(mh['n'] - ma['n'])/EV/(DE_H*cos); okD &= abs(v/(8/9) - 1) <= 0.03
            out.append(f'winner {winners[0]} seed {s} at eps 0 (709 runs): (H - Hbar)/(Delta E cos) {v:.4f} (BP x 8/9 = 0.8889)')
    out.append(f'H710d (the winner = BP x 8/9 within 3 %, only together with H710c): {okD}')
    text = '\n'.join(out); print(text)
    os.makedirs('out', exist_ok=True); open('out/ana710.txt', 'w').write(text + '\n')


if __name__ == '__main__':
    main()

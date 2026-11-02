# Audit 685: Delta m of the pair's 2p -> 1s Zeeman photon from Born's rule (CREM_PAIR_BORN_DELTA_M=1) against the default path (684 runs, same seeds).
# usage: python3 ana685.py <685 runs dir with born/ ps/> <684 runs dir with def/>
import sys, glob, math, os
R, R4 = sys.argv[1], sys.argv[2]
UEV = 1.602176634e-25
FIELDS = ('0.01', '0.1', '1.0329')


def vec(t, key, n=4):
    i = t.index(key); return [float(x) for x in t[i+1:i+1+n]]


def parse(path):
    d = {}
    for line in open(path):
        t = line.split()
        if line.startswith('PAIR BORNDM'):
            d['pBd'] = vec(t, 'pB', 3); d['cl'] = int(t[t.index('classical')+1]); d['born'] = int(t[t.index('born')+1])
            d['dJz'] = float(t[t.index('dJz')+1]); d['dE'] = float(t[t.index('dE')+1])/UEV
        elif line.startswith('PAIR QM1S684'):
            d['pB'] = vec(t, 'pB', 3); d['iL'] = int(t[t.index('iL')+1])
        elif line.startswith('PAIR QM1S path'):
            d['sc'], d['qm'] = vec(t, 'sc'), vec(t, 'qm'); d['m'] = [int(x) for x in vec(t, 'm')]
        elif line.startswith('PAIR eigen ledger'):
            e = {t[i]: float(t[i+1]) for i in range(3, len(t)-1, 2)}
            d['imp'] = (e['B'] + e['S2P'] - e['S1S'])/UEV
    return d


def gap(rows):
    n = len(rows)
    d = [[r['sc'][i] - r['qm'][i] for r in rows] for i in range(4)]
    mu = [sum(x)/n for x in d]
    sd = [2*math.sqrt(sum((x - m)**2 for x in xs)/(n - 1)/n) for xs, m in zip(d, mu)]
    i = max(range(4), key=lambda j: abs(mu[j])); return mu, abs(mu[i]), sd[i]


def rms(xs): return math.sqrt(sum(x*x for x in xs)/len(xs)) if xs else float('nan')


out = ['born = CREM_PAIR_BORN_DELTA_M=1, def = default path (684 runs, same seeds); 1S states in energy order']
okB1 = True; okC = True; okD = True; chi2 = 0.0; dof = 0; maxdp = 0.0; orphan = 0; qmsame = 0; nph = 0
for B in FIELDS:
    for P in ('1', '2'):
        born = [parse(f) for f in sorted(glob.glob(f'{R}/born/B{B}_A{P}_*.err'))]
        dflt = [parse(f.replace(R + '/born/', R4 + '/def/')) for f in sorted(glob.glob(f'{R}/born/B{B}_A{P}_*.err'))]
        rb = [r for r in born if 'sc' in r]; rd = [r for r in dflt if 'sc' in r]
        for b, d in zip(born, dflt):
            if 'pBd' in b and 'pB' not in b: orphan += 1
            if 'pBd' in b and 'pB' in b: maxdp = max(maxdp, max(abs(x - y) for x, y in zip(b['pBd'], b['pB'])))
            if 'qm' in b and 'qm' in d: nph += 1; qmsame += b['qm'] == d['qm']
        mb, gb, sb = gap(rb); md, gd, sd = gap(rd)
        out.append(f'B {B:>6} T PHEN {P} n born {len(rb)} def {len(rd)} m {rb[0]["m"]}')
        out.append('    QM   ' + ' '.join(f'{sum(r["qm"][i] for r in rb)/len(rb):.4f}' for i in range(4)))
        out.append('    born ' + ' '.join(f'{sum(r["sc"][i] for r in rb)/len(rb):.4f}' for i in range(4))
                   + '   diff ' + ' '.join(f'{x:+.4f}' for x in mb) + f'   max {gb:.4f} ({sb:.4f})')
        out.append('    def  ' + ' '.join(f'{sum(r["sc"][i] for r in rd)/len(rd):.4f}' for i in range(4))
                   + '   diff ' + ' '.join(f'{x:+.4f}' for x in md) + f'   max {gd:.4f} ({sd:.4f})')
        okC &= gb <= 0.02
        if B in ('0.01', '0.1') and P == '1': okC &= gb <= 0.5*gd
        ib = rms([r['imp'] for r in born if 'imp' in r]); idf = rms([r['imp'] for r in dflt if 'imp' in r])
        okD &= ib <= 1.2*idf
        dj = [r['dJz'] for r in born if 'dJz' in r]; de = [r['dE'] for r in born if 'dE' in r]
        changed = sum(1 for r in born if 'born' in r and r['born'] != r['cl'])
        out.append(f'    import rms born {ib:.2f} def {idf:.2f} ueV;  Delta m changed {changed}/{len(dj)}, dJz != 0 {sum(abs(x) > 1e-12 for x in dj)},'
                   f' mean |dJz| {sum(abs(x) for x in dj)/len(dj):.4f} hbar, rms dE {rms(de):.4f} ueV')
        bb = [r for r in born if 'pBd' in r]
        cnt = [sum(1 - r['born'] == i for r in bb) for i in range(3)]
        exp = [sum(r['pBd'][i] for r in bb) for i in range(3)]
        var = [sum(r['pBd'][i]*(1 - r['pBd'][i]) for r in bb) for i in range(3)]
        for i in range(3):
            if var[i] > 0: chi2 += (cnt[i] - exp[i])**2/var[i]; dof += 1
        out.append('    Born Delta m (+1/0/-1) count ' + ' '.join(map(str, cnt)) + '  expected ' + ' '.join(f'{x:.1f}' for x in exp))
p = math.exp(-chi2/2)*sum((chi2/2)**k/math.factorial(k) for k in range(dof//2))
okB1 = maxdp <= 1e-9 and orphan == 0
out.append(f'H685a: Ps with CREM_PAIR_BORN_DELTA_M=1 = without, byte for byte: '
           f'{sum(open(f, "rb").read() == open(f.replace("/on_", "/off_"), "rb").read() for f in sorted(glob.glob(R + "/ps/on_*")))}/{len(glob.glob(R + "/ps/on_*"))} files')
out.append(f'H685b: max |pB(draw) - pB(684 line)| {maxdp:.2e}, draws without a 661 photon {orphan}; Delta m frequency chi^2 {chi2:.1f} for {dof} cells, p {p:.3g}:'
           f' {okB1 and p > 0.05}')
out.append(f'  (QM reference identical to the default runs in {qmsame}/{nph} photons -- same trajectory up to the photon)')
out.append(f'H685c (born max gap <= 0.02 everywhere, <= half of def at 0.01 / 0.1 T PHEN 1): {okC}')
out.append(f'H685d (import rms born <= 1.2 x def everywhere): {okD}')
text = '\n'.join(out); print(text)
os.makedirs('out', exist_ok=True); open('out/ana685.txt', 'w').write(text + '\n')

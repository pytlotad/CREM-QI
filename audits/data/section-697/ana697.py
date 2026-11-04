# Audit 697: the pair projection at the photon from the state transported to the photon instant (CREM_PAIR_PHOTON_INSTANT=1)
# against the checkpoint-midpoint state of the default path (693 poi files, the same seeds), 1.0329 T, hydrogen and antihydrogen, PHEN 1 / 2.
# (A) = v3 - QM and the path gap sc - QM on the 1S states as in 693; Hotelling T^2 on 3 states ~ chi^2(3), a pair's 2 cells ~ chi^2(6).
# usage: python3 ana697.py <697 runs dir with on/ dry/ ctl/> <693 runs dir with poi/>
import sys, glob, math, os
import numpy as np
from scipy.stats import chi2
PAIRS = (('antiproton', 'antihydrogen'), ('proton', 'hydrogen'))


def vec(t, key, n=4):
    i = t.index(key); return [float(x) for x in t[i+1:i+1+n]]


def parse(path):
    d = {}
    for line in open(path):
        t = line.split()
        if line.startswith('PAIR QM1S684'): d['v3'] = vec(t, 'v3')
        elif line.startswith('PAIR QM1S path'): d['sc'], d['qm'] = vec(t, 'sc'), vec(t, 'qm')
        elif line.startswith('PAIR QM1S697'):
            d['delay'] = float(t[t.index('delay')+1]); d['done'] = int(t[t.index('done')+1])
            d['mid'] = np.array(vec(t, 'mid', 12)); d['inst'] = np.array(vec(t, 'inst', 12))
        elif line.startswith('END'): d['end'] = (t[t.index('stop')+1], t[t.index('photons')+1])
    return d


def hot(D):
    D = np.array(D); n = len(D); mu = D.mean(0)
    return mu, float(n*mu[:3] @ np.linalg.solve(np.cov(D[:, :3], rowvar=False), mu[:3]))



def main():
    R, R3 = sys.argv[1:3]
    out = []
    # H697a: controls
    ctl = sorted(glob.glob(f'{R}/ctl/*.err'))
    same = sum(open(f, 'rb').read() == open(f'{R3}/poi/{os.path.basename(f)}', 'rb').read() for f in ctl)
    dry = sorted(glob.glob(f'{R}/dry/*.err'))
    sdry = sum(''.join(l for l in open(f) if not l.startswith('PAIR QM1S697')) == open(f'{R}/ctl/{os.path.basename(f)}').read() for f in dry)
    out.append(f'H697a (files): no switch = 693 poi byte for byte {same}/{len(ctl)}; dry without QM1S697 = no switch {sdry}/{len(dry)}: '
               f'{same == len(ctl) == 16 and sdry == len(dry) == 16}')
    # H697b
    on = {os.path.basename(f): parse(f) for f in sorted(glob.glob(f'{R}/on/*.err'))}
    alldone = all(r.get('done') == 1 and r['delay'] >= 0 for r in on.values())
    dm = max(float(np.max(np.abs(on[os.path.basename(f)]['mid'] - parse(f)['mid']))) for f in dry)
    sm = max(abs(parse(f)['mid'].sum() - 1) for f in dry)
    ends = sorted(set(r.get('end') for r in on.values()))
    out.append(f'H697b: {len(on)} files, ends {ends}; done 1 and delay >= 0 in every photon {alldone}; max |mid(on) - mid(dry)| {dm:.1e};'
               f' max |sum mid - 1| (dry) {sm:.1e}: {alldone and dm <= 1e-9 and sm <= 1e-8}')
    dl = np.array([r['delay'] for r in on.values()])*1e12
    out.append('delay t_ph - T/2 [ps]: min {:.2f} quartiles {:.2f} {:.2f} {:.2f} max {:.2f} mean {:.2f}'.format(dl.min(), *np.percentile(dl, (25, 50, 75)), dl.max(), dl.mean()))
    okC = okD = okE = True; rmsd = {}
    for tag, name in PAIRS:
        T = 0.0
        for P in ('1', '2'):
            keys = sorted(k for k in on if k.startswith(f'{tag}_B1.0329_A{P}_'))
            a = [on[k] for k in keys]; b = [parse(f'{R3}/poi/{k}') for k in keys]
            qsame = max(max(abs(x - y) for x, y in zip(p['qm'], q['qm'])) for p, q in zip(a, b))
            Ai = np.array([[x - y for x, y in zip(p['v3'], p['qm'])] for p in a]); Am = np.array([[x - y for x, y in zip(q['v3'], q['qm'])] for q in b])
            dA = Ai - Am; n = len(keys)
            mi, mm = Ai.mean(0), Am.mean(0); s = int(np.argmax(np.abs(mm)))
            md, td = hot(dA); se = dA[:, s].std(ddof=1)/math.sqrt(n)
            r_d = math.sqrt((dA**2).sum(1).mean()); r_m = math.sqrt((Am**2).sum(1).mean()); rmsd[(tag, P)] = r_d
            G = np.array([[x - y for x, y in zip(p['sc'], p['qm'])] for p in a]); mg, tg = hot(G); T += tg
            G3 = np.array([[x - y for x, y in zip(q['sc'], q['qm'])] for q in b]); mg3, tg3 = hot(G3)
            i = int(np.argmax(np.abs(mg))); se2 = 2*G[:, i].std(ddof=1)/math.sqrt(n)
            p2 = np.array([float(np.sqrt(((p['inst'] - p['mid'])**2).sum())) for p in a]); dq = np.array([p['delay'] for p in a])
            qs = np.percentile(dq, (25, 50, 75)); qi = np.digitize(dq, qs)
            out.append(f'{name:12s} PHEN {P} n {n} (QM same as 693 to {qsame:.1e})')
            out.append(f'   (A) inst ' + ' '.join(f'{x:+.4f}' for x in mi) + f' | mid (693) ' + ' '.join(f'{x:+.4f}' for x in mm)
                       + f' | state {s}: inst {mi[s]:+.4f} mid {mm[s]:+.4f} diff {md[s]:+.4f} +- {se:.4f} (paired T2 {td:.2f} p {chi2.sf(td, 3):.2g})')
            out.append(f'   rms |(A)_inst - (A)_mid| {r_d:.4f} vs rms |(A)_mid| {r_m:.4f} (ratio {r_d/r_m:.2f}); max_s |mean (A)| inst {abs(mi).max():.4f} mid {abs(mm).max():.4f}'
                       f' (ratio {abs(mi).max()/abs(mm).max():.2f})')
            out.append(f'   gap sc - QM inst ' + ' '.join(f'{x:+.4f}' for x in mg) + f' max {abs(mg[i]):.4f} (2 sigma {se2:.4f}) T2 {tg:.2f} p {chi2.sf(tg, 3):.3g}'
                       f' | 693 max {abs(mg3).max():.4f} T2 {tg3:.2f} p {chi2.sf(tg3, 3):.3g}')
            out.append(f'   2P populations rms |inst - mid| {math.sqrt((p2**2).mean()):.4f}; by delay quartile ' + ' '.join(f'{math.sqrt((p2[qi == k]**2).mean()):.4f}' for k in range(4))
                       + ' (delay quartile edges ' + ' '.join(f'{x*1e12:.1f}' for x in qs) + ' ps)')
            if P == '2':
                okC &= abs(mi).max() >= 0.5*abs(mm).max()
                okD &= r_d >= 0.5*r_m and abs(md[s]) < 2*se
        out.append(f'{name:12s} pooled gap (inst) chi^2(6) {T:.1f} p {chi2.sf(T, 6):.3g}'); okE &= chi2.sf(T, 6) > 0.05
    out.append(f'H697c (PHEN 2: max |mean (A)| inst >= 0.5 x mid, both pairs): {okC}')
    out.append(f'H697d (PHEN 2: rms |inst - mid| >= 0.5 rms |mid| and |mean diff| < 2 SE in the largest state, both pairs): {okD}')
    out.append(f'H697e (gap inst pooled p > 0.05 both pairs): {okE}')
    out.append(f'H697f (rms |(A)_inst - (A)_mid| PHEN 1 < PHEN 2, both pairs): {all(rmsd[(t, "1")] < rmsd[(t, "2")] for t, _ in PAIRS)}')
    text = '\n'.join(out); print(text)
    os.makedirs('out', exist_ok=True); open('out/ana697.txt', 'w').write(text + '\n')


if __name__ == '__main__':
    main()

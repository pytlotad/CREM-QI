# Audit 708: the level n = sqrt(Ry/|E|) the mechanical path reads with the charge-dipole weight 1 against weight 0 (the conserved energy, audit 145):
# (b) range of both energies over each clock period on the retarded and the instantaneous sum, (c) photon energies at site A5 from both weights, and the
# per-site |n1 - n0| of the production statistical runs converted to a photon energy E(n) - E(n - 1).
# usage: python3 ana708.py <708 runs dir with ck/ ci/ a5/ ctl/ st/ stc/> <707 runs dir with ch/>
import sys, os, glob, re, filecmp
import numpy as np
EV = 1.602176634e-19; HEV = 4.135667696e-15; HC_EVNM = 1239.84198
RY_PS = 6.8028; RY_H = 13.5983   # eV, pair binding (reduced mass)
SITES = ('A1 ladder window', 'A3 E1 count factor', 'A4 rule rates', 'A5 photon', 'A6 window length')


def strip(path):
    out = []
    for line in open(path):
        if line.startswith('LEVELW708'): continue
        out.append(re.sub(r' rng1 \S+ rng0 \S+', '', line))
    return ''.join(out)


def periods(path):
    rows = []
    for line in open(path):
        if line.startswith('PAIR CLOCK701 k'):
            t = line.split(); v = {t[i]: float(t[i+1]) for i in range(2, len(t) - 1, 2)}; rows.append(v)
    return rows


def units(dE_eV, E_eV):
    lam = HC_EVNM/E_eV if E_eV > 0 else float('nan')
    return f'{dE_eV*1e6:+.4f} ueV = {dE_eV/HEV/1e9:+.5f} GHz (line {E_eV:.6f} eV, {lam:.4f} nm; shift {-lam*dE_eV/E_eV*1e6:+.4f} fm)'


def main():
    R, R7 = sys.argv[1:3]; out = []
    # H708a
    same = tot = 0
    for tag in ('1', '2', 'L2'):
        for f in ('out.txt', 'err.txt'):
            a, b = f'{R}/st/{tag}/{f}', f'{R}/stc/{tag}/{f}'
            ok = (open(a).read() == open(b).read()) if f == 'out.txt' else strip(a) == open(b).read()
            tot += 1; same += ok
            if f == 'out.txt' and not ok:   # POST HOC: wall-clock censoring (the traced run is slower) moves only the Kaplan-Meier lines
                wc = lambda p: [l for l in open(p) if 'Kaplan-Meier' not in l and 'is censored' not in l]
                out.append(f'   exp {tag} out.txt differs; without the Kaplan-Meier / censoring lines (wall-clock budget): {"identical" if wc(a) == wc(b) else "DIFFERENT"}')
    for f in sorted(glob.glob(f'{R}/ctl/*.err')):
        tot += 1; same += filecmp.cmp(f, f'{R7}/ch/{os.path.basename(f)}', shallow=False)
        tot += 1; same += strip(f'{R}/ck/{os.path.basename(f)}') == open(f).read()
    out.append(f'H708a (files): with the switch, LEVELW708 lines and rng fields removed = without; without = 707: {same}/{tot}: {same == tot}')
    # H708b
    okB = True
    for sumtag, D in (('retarded (default)', 'ck'), ('instantaneous', 'ci')):
        for f in sorted(glob.glob(f'{R}/{D}/*.err')):
            rows = periods(f)
            if not rows: out.append(f'{sumtag:18s} {os.path.basename(f)[:-4]:26s} no periods'); okB = False; continue
            r = np.array([v['rng0']/v['rng1'] for v in rows])
            okB &= bool((r < 0.1).all())
            out.append(f'{sumtag:18s} {os.path.basename(f)[:-4]:26s} periods {len(rows)} | range per period: weight 1 {np.mean([v["rng1"] for v in rows])/EV*1e6:.4f} ueV,'
                       f' weight 0 {np.mean([v["rng0"] for v in rows])/EV*1e6:.4f} ueV | ratio 0/1 mean {r.mean():.4f} max {r.max():.4f}')
    out.append(f'H708b (every clock period: range weight 0 < 0.1 x range weight 1, both sums, hydrogen and Ps): {okB}')
    # H708c
    okC = False; nph = 0
    for f in sorted(glob.glob(f'{R}/a5/*.err')):
        ph = []
        for line in open(f):
            if line.startswith('LEVELW708 photon'):
                t = line.split(); v = {t[i]: float(t[i+1]) for i in range(2, len(t) - 1, 2)}; ph.append(v)
        nph += len(ph)
        for v in ph:
            rel = abs(v['dph'])/v['ph1'] if v['ph1'] > 0 else 0.0; okC |= rel >= 1e-6
        if ph:
            d = np.array([v['dph'] for v in ph])/EV; e = np.array([v['ph1'] for v in ph])/EV
            i = int(np.argmax(np.abs(d)))
            out.append(f'A5 {os.path.basename(f)[:-4]:26s} photons {len(ph)} | n1 - n0 max {max(abs(v["n1"] - v["n0"]) for v in ph):.3e} | largest photon difference '
                       + units(d[i], e[i]) + f' | mean |d| {np.abs(d).mean()*1e6:.4f} ueV')
        else:
            out.append(f'A5 {os.path.basename(f)[:-4]:26s} photons 0')
    for tag, ry in (('1', RY_PS), ('2', RY_PS), ('L2', RY_PS)):
        cnt = np.zeros(5); mx = np.zeros(5); sm = np.zeros(5)
        for line in open(f'{R}/st/{tag}/err.txt'):
            if line.startswith('LEVELW708 run'):
                t = line.split()
                for k in range(5):
                    c, m, a = int(t[3 + 4*k]), float(t[4 + 4*k]), float(t[5 + 4*k]); cnt[k] += c; sm[k] += c*m; mx[k] = max(mx[k], a)
        for k in range(5):
            if cnt[k]:
                n = 2.0; dEdn = ry*(-2/(n - 1)**3 + 2/n**3)   # d/dn [Ry (1/(n-1)^2 - 1/n^2)] at n = 2
                out.append(f'statistical exp {tag:2s} {SITES[k]:20s} reads {int(cnt[k])} | mean n1 - n0 {sm[k]/cnt[k]:+.3e} max |n1 - n0| {mx[k]:.3e}'
                           f' -> equivalent on a 2 -> 1 photon (Ps; this site sets rates, not a photon): ' + units(mx[k]*dEdn, ry*(1 - 0.25)))
            else:
                out.append(f'statistical exp {tag:2s} {SITES[k]:20s} reads 0')
    out.append(f'H708c (A5: |photon(1) - photon(0)| >= 1e-6 x photon in at least one photon; {nph} photons): {okC if nph else "not decidable (no A5 photon)"}')
    text = '\n'.join(out); print(text)
    os.makedirs('out', exist_ok=True); open('out/ana708.txt', 'w').write(text + '\n')


if __name__ == '__main__':
    main()

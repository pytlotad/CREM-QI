# Audit 695 (print only): the 2P eigen-populations of the classical coherent state along the transport (PAIR QM1S695).
# D_k(t) = p_k(t) - p_k(formation); settling time, spectrum of sum |D|^2 against the quantum beat frequencies, the photon's spin instant.
# usage: python3 ana695.py <695 runs dir with tr/> <694 runs dir with poi/> <693 runs dir with poi/>
import sys, glob, os
import numpy as np
R, R4, R3 = sys.argv[1:4]
H = 6.62607015e-34
PAIRS = (('antiproton', 'antihydrogen'), ('proton', 'hydrogen'))


def parse(path):
    d = dict(t=[], p=[])
    for line in open(path):
        s = line.split()
        if line.startswith('PAIR QM1S695 form'):
            d['w'] = np.array([float(x) for x in s[4:16]]); d['pf'] = np.array([float(x) for x in s[17:29]])
        elif line.startswith('PAIR QM1S695 t'): d['t'].append(float(s[3])); d['p'].append([float(x) for x in s[5:17]])
        elif line.startswith('PAIR QM1S695 photon'): d['ts'] = float(s[4]); d['pp'] = np.array([float(x) for x in s[8:20]])
        elif line.startswith('CREM_LEDGER t'): d['tau'] = float(s[2])
    d['t'] = np.array(d['t']); d['p'] = np.array(d['p']); return d


def settle(t, D, win=50e-12):
    r = np.sqrt((D**2).sum(1)); full = np.sqrt((r**2).mean())
    for i in range(len(t)):
        m = (t >= t[i]) & (t < t[i] + win)
        if t[i] + win > t[-1]: break
        if np.sqrt((r[m]**2).mean()) >= 0.8*full: return t[i] + win, full
    return float('nan'), full


out = []
# H695a (file part): trace files without QM1S695 = 694 / 693 files
same = tot = 0
for f in sorted(glob.glob(f'{R}/tr/*.err')):
    b = os.path.basename(f); s = int(b.rsplit('_', 1)[1][:-4])
    ref = f'{R3}/poi/{b}' if '_B1.0329_' in b else f'{R4}/poi/{b}'
    if not os.path.exists(ref): continue
    tot += 1; same += ''.join(l for l in open(f) if not l.startswith('PAIR QM1S695')) == open(ref).read()
out.append(f'H695a (files): trace files without QM1S695 = 694 / 693 byte for byte {same}/{tot} (regression and validation: out/regression695.txt, out/validation.txt)')
okB1 = okB2 = okC = okD = okE = okF = True; sumerr = 0.0; firsterr = 0.0
for B in ('0.01', '0.1', '1.0329'):
    for tag, name in PAIRS:
        for P in ('1', '2'):
            for f in sorted(glob.glob(f'{R}/tr/{tag}_B{B}_A{P}_*.err')):
                d = parse(f); t, p = d['t'], d['p']; D = p - d['pf']
                se = np.abs(p.sum(1) - 1).max(); fe = np.abs(p[0] - d['pf']).max(); sumerr = max(sumerr, se); firsterr = max(firsterr, fe)
                okB1 &= se < 1e-9; okB2 &= fe < 1e-3
                mx = np.abs(D).max(); ts, full = settle(t, D)
                # spectrum of S = sum |D|^2 on a 1 ps grid
                g = np.arange(t[0], t[-1], 1e-12); S = np.interp(g, t, (D**2).sum(1)); S = S - S.mean()
                F = np.fft.rfftfreq(len(g), 1e-12); A = np.abs(np.fft.rfft(S)); A[0] = 0; fpk = F[np.argmax(A)] if len(A) > 1 else float('nan')
                beats = np.array(sorted({abs(a - b)/H for a in d['w'] for b in d['w'] if a != b}))
                near = beats[np.argmin(np.abs(beats - fpk))]; rel = abs(fpk - near)/near
                j = np.argmin(np.abs(t - d['ts'])); pm = np.abs(p[j] - d['pp']).max()
                okF &= pm < 1e-4 and t[j] <= d['tau']
                if P == '2' and B == '1.0329': okC &= mx > 0.05; okD &= ts < 150e-12; okE &= rel < 0.1
                out.append(f'{name:12s} B {B:>6} T PHEN {P} {os.path.basename(f)[-6:-4].strip("_")}: n {len(t)} t_end {t[-1]*1e9:.3f} ns | max|D| {mx:.4f} rms|D| {full:.4f}'
                           f' t_s {ts*1e12:.0f} ps | S peak {fpk/1e9:.2f} GHz nearest beat {near/1e9:.2f} GHz ({rel*100:.1f} %)'
                           f' | spins {d["ts"]*1e12:.1f} ps photon {d["tau"]*1e12:.1f} ps (lag {(d["tau"] - d["ts"])*1e12:.1f} ps) match {pm:.1e}')
out.append(f'H695b: max |sum p - 1| {sumerr:.1e} (< 1e-9: {okB1}); max |first trace point - p(form)| {firsterr:.1e} (< 1e-3: {okB2})')
out.append(f'H695c (PHEN 2, 1.0329 T: max |D| > 0.05 every run): {okC}')
out.append(f'H695d (PHEN 2, 1.0329 T: t_s < 150 ps every run): {okD}')
out.append(f'H695e (PHEN 2, 1.0329 T: S peak within 10 % of a beat (w_k - w_l)/h every run): {okE}')
out.append(f'H695f (photon state = nearest trace point to 1e-4 and before the photon, every run): {okF}')
text = '\n'.join(out); print(text)
os.makedirs('out', exist_ok=True); open('out/ana695.txt', 'w').write(text + '\n')

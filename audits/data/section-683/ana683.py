# Audit 683: 1S distribution of the default (photon-time) and formation paths against the QM distribution of the formation state.
# usage: python3 ana683.py <runs dir with def/ form/ ctl/ ps/>
import sys, glob, math, os, re
R = sys.argv[1]
UEV = 1.602176634e-25
FIELDS = ('0.01', '0.1', '1.0329', '3')


def parse(path):
    qm = sc = m = led = None; ovl = None; end = None
    for line in open(path):
        if line.startswith('PAIR QM1S'):
            t = line.split()
            i = t.index('qm'); qm = [float(x) for x in t[i+1:i+5]]
            i = t.index('sc'); sc = [float(x) for x in t[i+1:i+5]]
            i = t.index('m'); m = [int(x) for x in t[i+1:i+5]]
            ovl = float(t[t.index('ovl')+1])
        elif line.startswith('PAIR eigen ledger'):
            t = line.split(); d = {t[i]: float(t[i+1]) for i in range(3, len(t)-1, 2)}
            led = ((d['S2P']-d['Zsb'])/UEV, (d['S1S']-d['Zsa'])/UEV)
        elif line.startswith('END'):
            end = line
    return qm, sc, m, ovl, led, end


out = []
res = {}
for path in ('def', 'form'):
    for B in FIELDS:
        for P in ('1', '2'):
            rows = [parse(f) for f in sorted(glob.glob(f'{R}/{path}/B{B}_A{P}_*.err'))]
            rows = [r for r in rows if r[0] is not None]
            res[(path, B, P)] = rows
out.append('mean over photons; 1S states in energy order (m labels per field shown); diff = mean(sc - qm) with 2 sigma of the mean')
okC = True; okD = True; okB = True; worstD = {}
for B in FIELDS:
    for P in ('1', '2'):
        for path in ('def', 'form'):
            rows = res[(path, B, P)]; n = len(rows)
            if n == 0:
                out.append(f'B {B} T PHEN {P} {path}: no QM1S photons'); continue
            mq = [sum(r[0][i] for r in rows)/n for i in range(4)]
            ms = [sum(r[1][i] for r in rows)/n for i in range(4)]
            d = [[r[1][i]-r[0][i] for r in rows] for i in range(4)]
            md = [sum(x)/n for x in d]
            sd = [2*math.sqrt(sum((x-mu)**2 for x in xs)/max(1, n-1)/n) for xs, mu in zip(d, md)]
            labels = rows[0][2]
            s2 = math.sqrt(sum(r[4][0]**2 for r in rows if r[4])/max(1, sum(1 for r in rows if r[4])))
            s1 = math.sqrt(sum(r[4][1]**2 for r in rows if r[4])/max(1, sum(1 for r in rows if r[4])))
            out.append(f'B {B:>6} T PHEN {P} {path:4s} n {n:2d} m {labels}  QM ' + ' '.join(f'{x:.4f}' for x in mq)
                       + '  path ' + ' '.join(f'{x:.4f}' for x in ms) + '  diff ' + ' '.join(f'{a:+.4f}({b:.4f})' for a, b in zip(md, sd))
                       + f'  rms S2P* {s2:.2f} S1S* {s1:.2f} ueV')
            if path == 'form':
                okC &= all(abs(a) <= b for a, b in zip(md, sd))
                okB &= all(r[3] >= 0.999 for r in rows)
            else:
                w = max(abs(a) for a in md); worstD[(B, P)] = w; okD &= w <= 0.02
out.append(f'H683b (overlap >= 0.999 in every formation photon): {okB}')
out.append(f'H683c (formation path mean diff within 2 sigma, every state and field): {okC}')
out.append(f'H683d (default path max |mean diff| <= 0.02 at every field): {okD}  worst ' + ', '.join(f'{k[0]} T P{k[1]} {v:.4f}' for k, v in worstD.items()))
# H683a controls
ctl = sorted(glob.glob(f'{R}/ctl/*.err')); same = 0
for f in ctl:
    e1 = parse(f)[5]; e2 = parse(f.replace('/ctl/', '/def/'))[5]; same += e1 == e2 and e1 is not None
out.append(f'H683a: END lines with CREM_PAIR_QM1S=1 = without: {same}/{len(ctl)}')
ps = 0; tot = 0
for f in sorted(glob.glob(f'{R}/ps/on_*')):
    tot += 1; ps += open(f, 'rb').read() == open(f.replace('/on_', '/off_'), 'rb').read()
out.append(f'H683a: Ps with both variables = without, byte for byte: {ps}/{tot} files')
text = '\n'.join(out); print(text)
os.makedirs('out', exist_ok=True); open('out/ana683.txt', 'w').write(text + '\n')

# Audit 845: closure of the 840-844 chain from the committed reports (out/ana*.txt).  Written after prereg845.txt, before running.
import os, re
B = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')
rd = lambda n, f: open(os.path.join(B, f'section-{n}', 'out', f)).read()
pr = lambda *a: print(*a, flush=True)
num = r'([+-]\d+\.\d+)'
t840, t841, t843, t844 = rd(840, 'ana840.txt'), rd(841, 'ana841.txt'), rd(843, 'ana843.txt'), rd(844, 'ana844.txt')
m = re.search(r'Z   : F = 0 ' + num + ', F = 1 ' + num, t840); ZD840 = {'H0': float(m[1]), 'H1': float(m[2])}
SYS = ('ops', 'H0', 'H1')
dP841 = dict(zip(SYS, map(float, re.findall(r'Delta_P \(ZD - Z\) ' + num, t841))))
dM841 = dict(zip(SYS, map(float, re.findall(r'Delta_p814 ' + num, t841))))
dP843 = dict(zip(SYS, map(float, re.findall(r'Delta_C statycznie \(ZDC - ZD\) ' + num, t843))))
dM843 = dict(zip(SYS, map(float, re.findall(r'Delta_C p814 ' + num, t843))))
m = re.search(r'\(F = 0 ' + num + ', F = 1 ' + num + r' MHz\)', t843); S843 = {'H0': float(m[1]), 'H1': float(m[2])}
R844 = {}
for blk, sk in zip(re.split(r'\n(?=\S)', t844)[:4], ('ops', 'pps', 'H0', 'H1')):
    for q in ('Z', 'ZD', 'ZDC'):
        R844[sk, q] = float(re.search(r'\n    ' + q + r' *: nu_silnik .*?reszta alpha\^6 ' + num, blk)[1])
lab = {'ops': 'o-Ps (PHEN=2)', 'H0': 'wodor F = 0 (PHEN=1)', 'H1': 'wodor F = 1 (PHEN=2)'}
okc = okd = okb = True
for sk in ('H0', 'H1'):
    c = ZD840[sk] + dP841[sk] + dP843[sk]; d = ZD840[sk] + dM841[sk] + dM843[sk]
    okc &= abs(c - S843[sk]) <= 0.001; okd &= abs(d - S843[sk]) <= 0.001
    pr(f'{lab[sk]}: 840 Z - D {ZD840[sk]:+.4f} + 841 {dP841[sk]:+.4f} + 843 {dP843[sk]:+.4f} = {c:+.4f} (statycznie);  '
       f'p814: {ZD840[sk]:+.4f} {dM841[sk]:+.4f} {dM843[sk]:+.4f} = {d:+.4f};  843 S ZDC - D {S843[sk]:+.5f} MHz')
for sk in SYS:
    b = R844[sk, 'Z'] + dP841[sk] + dP843[sk]; okb &= abs(b - R844[sk, 'ZDC']) <= 0.001
    pr(f'{lab[sk]}: 844 (Z - BOC) {R844[sk, "Z"]:+.4f} + 841 {dP841[sk]:+.4f} + 843 {dP843[sk]:+.4f} = {b:+.4f};  844 (ZDC - BOC) {R844[sk, "ZDC"]:+.4f} MHz')
pr(f"H845c (domkniecie statyczne do 0,001 MHz): {'TAK' if okc else 'NIE'}")
pr(f"H845d (domkniecie dynamiczne p814 do 0,001 MHz): {'TAK' if okd else 'NIE'}")
pr(f"H845b (844 ZDC - BOC = Z - BOC + 841 + 843 do 0,001 MHz, z o-Ps): {'TAK' if okb else 'NIE'}")

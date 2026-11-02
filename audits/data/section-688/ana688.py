# Audit 688: Born Delta m on by default -- byte-for-byte checks of the new probe against the 686 / 685 files, and the regression of ba08e143.
# usage: python3 ana688.py <688 runs dir> <686 runs dir> <685 runs dir>
import sys, glob, os, hashlib
R, R6, R5 = sys.argv[1:4]


def same(a, b): return os.path.exists(b) and open(a, 'rb').read() == open(b, 'rb').read()


def count(pattern, ref):
    fs = sorted(glob.glob(pattern)); return sum(same(f, ref(f)) for f in fs), len(fs)


out = []
b = count(f'{R}/born/*.err', lambda f: f.replace(R + '/born/', R6 + '/born/'))
d = count(f'{R}/def/*.err', lambda f: f.replace(R + '/def/', R6 + '/def/'))
out.append(f'H688b: no variable = 686 born {b[0]}/{b[1]}; CREM_NO_PAIR_BORN_DELTA_M=1 = 686 def {d[0]}/{d[1]}: {b[0] == b[1] == d[0] == d[1] == 64}')
ps = [count(f'{R}/ps/off_A{P}_*.{x}', lambda f: f.replace(R + '/ps/', R5 + '/ps/')) for P in ('1', '2') for x in ('out', 'err')]
out.append('H688c: e+e- 1.0329 T, no variable = 685 off (PHEN 1 out, err; PHEN 2 out, err): ' + ', '.join(f'{a}/{n}' for a, n in ps)
           + f': {all(a == n == 4 for a, n in ps)}')
a = count(f'{R}/acc/*.err', lambda f: f.replace(R + '/acc/', R + '/born/'))
out.append(f'H688d: CREM_PAIR_BORN_DELTA_M=1 = no variable {a[0]}/{a[1]}: {a[0] == a[1] == 8}')
out.append('H688a regression: ba08e143 (old) vs working tree (new), no variables, sequential; cmp out.txt and err.txt')
ok = 0
for tag in ('1', '2', 'L2', '5Ps', '5H', '6'):
    s = [same(f'{R}/reg/old/{tag}/{x}', f'{R}/reg/new/{tag}/{x}') for x in ('out.txt', 'err.txt')]; ok += all(s)
    out.append(f'  {tag} out {"identical" if s[0] else "DIFFERENT"} err {"identical" if s[1] else "DIFFERENT"} md5 '
               + hashlib.md5(open(f'{R}/reg/new/{tag}/out.txt', 'rb').read()).hexdigest())
v = open(f'{R}/validation.txt').read()
out.append(f'  {ok}/6 identical; validation: ' + ' '.join(l.strip() for l in v.splitlines() if 'Enforced checks' in l or 'EXIT' in l))
text = '\n'.join(out); print(text)
os.makedirs('out', exist_ok=True); open('out/ana688.txt', 'w').write(text + '\n')

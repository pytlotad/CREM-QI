# Audit 804: apsis-detector hysteresis (engine print) + default per-orbit measure (per_orbit.py).  usage: python3 ana804.py
# 799 / 800 / 803 runs unpacked into $R799, $R800, $R803 (../section-79x|80x/runs/runs*.tar.xz); new runs in out/runs.
import os, sys, glob, math, shutil
import numpy as np
D = os.path.dirname(os.path.abspath(__file__)); O = os.path.join(D, 'out', 'runs'); sys.path.insert(0, D)
from per_orbit import om
SP = '/tmp/claude-1000/-home-teddy-Projekty-CREM-QI/4fa83d5c-8878-4cfd-898e-1a9b53a01ab2/scratchpad'
R799 = os.environ.get('R799', SP + '/r799x/runs'); R800 = os.environ.get('R800', SP + '/r800x/runs'); R803 = os.environ.get('R803', SP + '/r803x/runs')
pr = lambda *a: print(*a, flush=True)
def old_om(p, tag):   # measure of ana799-803
    L = []; S = []; per = None
    for ln in open(p):
        if not ln.startswith(f'PAIR TORQUE799 {tag} '): continue
        w = ln.split()
        if w[3] == 'apo':
            i = w.index('L'); L.append(np.array([float(x) for x in w[i + 1:i + 4]])); i = w.index('s2'); S.append(np.array([float(x) for x in w[i + 1:i + 4]]))
        elif w[3] == 'sum': per = float(w[w.index('period') + 1])
    d = [np.linalg.norm(S[j] - S[j - 1])/np.linalg.norm(np.cross(L[j - 1], S[j - 1])) for j in range(1, len(S))]
    return float(np.mean(d))/(2*math.pi*per)
R1, R1o = f'{O}/R1_F1_proton_n4L4_5.err', f'{R803}/M/F1_proton_n4L4_5.err'
a, b = open(R1o).read().splitlines(), open(R1).read().splitlines()
exp = [ln.replace(' apo 200', ' apo 199') if ln.startswith('PAIR TORQUE799 mechB0 sum') else ln for ln in a if not ln.startswith('PAIR TORQUE799 mechB0 apo 200 ')]
okK = b == exp and len(a) - len(b) == 1
pr(f"H804k (R1 = 803 bez linii mechB0 apo 200 i z 'apo 199' w linii sum, reszta bajt w bajt): {'TAK' if okK else 'NIE'} (linii 803 {len(a)}, R1 {len(b)})")
pairs = [(f'{O}/R2_F1_antiproton_n4L4_5.err', f'{R803}/M/F1_antiproton_n4L4_5.err'), (f'{O}/R3_proton_n3L1_5.err', f'{R799}/A/proton_n3L1_5.err'),
         (f'{O}/R4_proton_n2L2_7.err', f'{R800}/B/proton_n2L2_7.err')]
okC = True
for n, o in pairs:
    same = open(n, 'rb').read() == open(o, 'rb').read(); okC &= same; pr(f'  {os.path.basename(n)} = {o.split("/runs/")[1]}: {"bajt w bajt" if same else "ROZNY"}')
pr(f"H804c (R2, R3, R4 bajt w bajt): {'TAK' if okC else 'NIE'}")
files = sorted(glob.glob(f'{R799}/*/*.err') + glob.glob(f'{R800}/*/*.err') + glob.glob(f'{R803}/M/*.err'))
okM = True; nser = 0; rej = 0
for f in files:
    g = R1 if f == R1o else f
    for tag in ('mech', 'mechB0'):
        if not any(ln.startswith(f'PAIR TORQUE799 {tag} sum') for ln in open(g)): continue
        nser += 1; w, r, _, _ = om(g, tag); rej += r
        if f == R1o and tag == 'mechB0':
            T = [float(ln.split()[6]) for ln in open(f) if ln.startswith('PAIR TORQUE799 mechB0 apo')]
            ref = om(f, tag)[0]; q = abs(w/ref - 1); okM &= q <= 1e-12
            pr(f'  R1 mechB0: miara domyslna / POST HOC 803 - 1 = {q:.1e}; stara miara na pliku 803 / domyslna - 1 = {old_om(f, tag)/w - 1:+.3e}')
        else: okM &= w == old_om(g, tag)
okM &= rej == 0
pr(f"H804m (miara domyslna na {nser} seriach: odrzuconych odcinkow {rej}; = miara 799-803 bit w bit poza R1; R1 = POST HOC 803 do 1e-12): {'TAK' if okM else 'NIE'}")
# RAPORT: 803 scored with the default measure (R1 in place of the 803 file)
tmp = os.path.join(SP, 'r804eval'); shutil.rmtree(tmp, ignore_errors=True); os.makedirs(f'{tmp}/out/runs/M')
for f in glob.glob(f'{R803}/M/*.err'): os.symlink(R1 if f == R1o else f, f'{tmp}/out/runs/M/{os.path.basename(f)}')
src = open(os.path.join(D, '..', 'section-803', 'ana803.py')).read()
a0 = src.index('def om(p):'); a1 = src.index('def sec(p):')
src = src[:a0] + 'def om(p):\n    w, r, n, napo = _om(p); return w, napo\n' + src[a1:]
pr('=== RAPORT: ocena 803 z miara domyslna (wodor F = 1 (PHEN=2) i F = 0 (PHEN=1)) ===')
exec(compile(src, 'ana803-default', 'exec'), {'__file__': f'{tmp}/ana803.py', '_om': om, '__name__': '__main__'})

# Audit 866: which Thomas form agrees with hydrogen 2P -- indirectly through Dirac.  Reads the committed reports of 840 / 841 / 843;
# computes the alpha^6 Dirac term of 2P3/2 - 2P1/2.  Written after prereg866.txt.  Usage: LC_ALL=C python3 ana866.py -> out/ana866.txt.
import os, re
from decimal import Decimal as Dd, getcontext
getcontext().prec = 60
D = os.path.dirname(os.path.abspath(__file__)); A = os.path.join(D, '..')
pr = lambda *a: print(*a, flush=True)
HEV, C = 4.135667696e-15, 299792458.0
SIG = 0.015   # MHz, 2P measurement (Hagley i Pipkin 1994)
units = lambda f: f'{f:+.4f} MHz = {f*1e6*HEV*1e6:+.5f} ueV' + (f' (|lambda| {C/abs(f*1e6)*100:.3g} cm)' if f else '')
# --- reports
ph = open(f'{A}/section-840/out/posthoc840.txt').read().split('\n')
blk = None; ref = {'D': {}, 'H': {}}
for ln in ph:
    if ln.startswith('Dirac FW'): blk = 'D'
    elif ln.startswith('odniesienie z 1/2'): blk = 'H'
    m = re.match(r'\s+(Z|ZT|Z0|ZT0)\s+F0 ([+-][\d.]+) F1 ([+-][\d.]+)', ln)
    if m and blk: ref[blk][m.group(1)] = (float(m.group(2)), float(m.group(3)))
d41 = open(f'{A}/section-841/out/ana841.txt').read()
dP = [float(x) for x in re.findall(r'wodor F = [01] \(PHEN=[12]\): Delta_P \(ZD - Z\) ([+-][\d.]+)', d41)]
d43 = open(f'{A}/section-843/out/ana843.txt').read()
zd = re.search(r'nu_ZD - nu_D: F = 0 ([+-][\d.]+) MHz.*?F = 1 ([+-][\d.]+) MHz', d43)
zdc = re.search(r'F = 0 ([+-][\d.]+), F = 1 ([+-][\d.]+) MHz\)', d43)
dC = [float(x) for x in re.findall(r'wodor F = [01] \(PHEN=[12]\): Delta_C statycznie \(ZDC - ZD\) ([+-][\d.]+)', d43)]
ZD = (float(zd.group(1)), float(zd.group(2))); ZDC = (float(zdc.group(1)), float(zdc.group(2)))
pr('Wejscie (raporty 840 / 841 / 843):')
for b in ('D', 'H'):
    for k, v in ref[b].items(): pr(f'  {"Dirac FW" if b == "D" else "odn. 1/2 "} {k:4s} F0 {v[0]:+.4f} F1 {v[1]:+.4f}')
pr(f'  841 Delta_P (ZD - Z) F0 {dP[0]:+.4f} F1 {dP[1]:+.4f};  843 ZD - D {ZD[0]:+.4f} / {ZD[1]:+.4f}, Delta_C {dC[0]:+.4f} / {dC[1]:+.4f}, ZDC - D {ZDC[0]:+.5f} / {ZDC[1]:+.5f}')
# --- H866b
F = {'Z': ref['D']['Z'], 'ZT': ref['D']['ZT'], 'Z0': ref['D']['Z0'], 'ZT0': ref['D']['ZT0'],
     'Z0B': ref['H']['Z0'], 'ZT0B': ref['H']['ZT0'], 'ZDC': ZDC}
c1 = [abs(F['Z0B'][i] - (F['Z'][i] + dP[i])) for i in (0, 1)]
c2 = [abs(ZDC[i] - (ZD[i] + dC[i])) for i in (0, 1)]
pr(f'H866b: |Z0B - (Z - D + Delta841)| = {c1[0]:.4f} / {c1[1]:.4f};  |ZDC - (ZD + Delta843)| = {c2[0]:.4f} / {c2[1]:.4f} MHz -> '
   f'{"TAK" if max(c1 + c2) <= 0.001 else "NIE"}')
# --- H866e: Dirac interval, reduced mass
alpha = Dd(1) / Dd('137.035999177'); me, mp = Dd('9.1093837139e-31'), Dd('1.67262192595e-27'); mu = me*mp/(me + mp)
mec2_Hz = Dd('1.2355899865e20')   # m_e c^2 / h  [Hz] (CODATA 2022: 8.1871057880e-14 J / 6.62607015e-34)
mu_Hz = mec2_Hz*mu/me
def f(n, j):
    k = j + Dd('0.5'); dl = k - (k*k - alpha*alpha).sqrt(); x = alpha/(n - dl)
    return 1/(1 + x*x).sqrt() - 1
exact = (f(2, Dd('1.5')) - f(2, Dd('0.5')))*mu_Hz/Dd(1e6)
lead = alpha**4/32*mu_Hz/Dd(1e6)
a6 = exact - lead
th, me_ = 10969.0415, 10969.045
pr(f'  Dirac 2P3/2 - 2P1/2 (masa zredukowana, bez QED): dokladnie {float(exact):.4f} MHz, alpha^4 {float(lead):.4f}, alpha^6 i wyzej {float(a6):+.4f} MHz = {float(a6)/SIG:+.1f} sigma_2P')
pr(f'  teoria (Horbatsch i Hessels 2016) {th} wobec pomiaru {me_}(15): {th - me_:+.4f} MHz = {(th - me_)/SIG:+.2f} sigma')
okE = abs(th - me_) <= SIG and abs(float(a6)) >= 10*SIG
pr(f'H866e: {"TAK" if okE else "NIE"}')
# --- H866m, H866t, H866r
pr('Tabela nu_X - nu_D (wodor, czastka probna):')
for k, v in F.items():
    pr(f'  {k:4s} F = 0 {units(v[0])} [{v[0]/SIG:+6.1f} s];  F = 1 {units(v[1])} [{v[1]/SIG:+6.1f} s];  parzysta {(v[0]+v[1])/2:+.4f}, nieparzysta {(v[1]-v[0])/2:+.4f}')
agree = [k for k, v in F.items() if max(abs(v[0]), abs(v[1])) <= SIG]
pr(f'H866m (GLOWNA; zgodna tylko ZDC): zgodne = {agree} -> {"TAK" if agree == ["ZDC"] else "NIE"}')
mx = lambda k: max(abs(F[k][0]), abs(F[k][1]))
lv = (('1 dipol lab.', 'Z', 'ZT'), ('2 dipol wlasny', 'Z0', 'ZT0'), ('3 + czynnik BMT', 'Z0B', 'ZT0B'))
okT = True
for lab, a, b in lv:
    pr(f'  poziom {lab}: max|.| L_T {mx(a):.4f} ({mx(a)/SIG:.1f} s), Tulczyjew {mx(b):.4f} ({mx(b)/SIG:.1f} s) MHz')
    okT &= mx(a) < mx(b)
pr(f'H866t (L_T blizej na poziomach 1-3): {"TAK" if okT else "NIE"}')
odd = (F['ZT0B'][1] - F['ZT0B'][0])/2
pr(f'H866r (ZT0B nieparzysta {odd:+.4f} MHz = {abs(odd)/SIG:.1f} sigma >= 3): {"TAK" if abs(odd) >= 3*SIG else "NIE"}')

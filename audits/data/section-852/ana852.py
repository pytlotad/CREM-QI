# Audit 852: default stable photon increment.  Usage: LC_ALL=C python3 ana852.py -> out/ana852.txt.  Written during the runs, before any result.
# Reads out/runs/M (mechanical), out/regression852.txt, out/validation_{old,new}.txt and the reg files in $R852REG (scratchpad r852reg/reg).
import os, re, tarfile, math
D = os.path.dirname(os.path.abspath(__file__)); O = os.path.join(D, 'out'); pr = lambda *a: print(*a, flush=True)
SP = '/tmp/claude-1000/-home-teddy-Projekty-CREM-QI/4fa83d5c-8878-4cfd-898e-1a9b53a01ab2/scratchpad'
REG = os.environ.get('R852REG', SP + '/r852reg') + '/reg'
t850 = tarfile.open(os.path.join(D, '..', 'section-850', 'runs', 'runs850.tar.xz')); t851 = tarfile.open(os.path.join(D, '..', 'section-851', 'runs', 'runs851.tar.xz'))
okM = True
for s in ('H0', 'H1', 'M0', 'oPs'):
    for ext in ('err', 'out'):
        ref = t851.extractfile(f'runs/W/{s}_4f_0.{ext}').read().replace(b'PAIR PHOTONSTABLE851 on\n', b'')
        new = open(f'{O}/runs/M/{s}_4f_0.{ext}', 'rb').read(); same = ref == new; okM &= same
        pr(f'  M {s}_4f_0.{ext} = 851 (bez wiersza PHOTONSTABLE851): {same}')
for ext in ('err', 'out'):
    same = t850.extractfile(f'runs/W/H0_4f_0.{ext}').read() == open(f'{O}/runs/M/L_H0_4f_0.{ext}', 'rb').read(); okM &= same
    pr(f'  L H0_4f_0.{ext} (LEGACY) = 850: {same}')
pr(f"H852m (sciezka mechaniczna: domyslnie = 851, LEGACY = 850, bajt w bajt): {'TAK' if okM else 'NIE'}")
reg = open(f'{O}/regression852.txt').read(); pr(reg.strip())
TAGS = ('1', '2', 'L2', '5Ps', '5H', '6')
okO = all(f'old vs new_legacy {t} SAME' in reg for t in TAGS)
pr(f"H852o (new_legacy = old bajt w bajt na 6 zestawach): {'TAK' if okO else 'NIE'}")
def val(v):
    t = open(f'{O}/validation_{v}.txt').read(); m = re.search(r'Enforced checks:\s+(\d+)/(\d+)', t)
    return ('EXIT 0' in t, 'validation verdict: PASS' in t, (int(m.group(1)), int(m.group(2))) if m else None)
vo, vn = val('old'), val('new')
okV = vn[0] and vn[1] and vn[2] is not None and vo[2] is not None and vn[2][0] == vo[2][0] and vn[2][0] == vn[2][1]
pr(f"H852v (walidacja new EXIT 0, PASS, liczba kontroli jak old): {'TAK' if okV else 'NIE'}  (old {vo}, new {vn})")
PM = re.compile(r'^(.*?)([-+]?\d[\d.]*(?:e[-+]?\d+)?)\s*\+/-\s*([\d.]+(?:e[-+]?\d+)?)')
okZ = True
for t in TAGS:
    st = 'SAME' if f'old vs new {t} SAME' in reg else 'DIFF'
    o = open(f'{REG}/old/{t}/out.txt').read().split('\n'); n = open(f'{REG}/new/{t}/out.txt').read().split('\n')
    pr(f'  zestaw {t}: {st}; linii out.txt old {len(o)} new {len(n)}')
    if st == 'SAME': continue
    po = {}; 
    for ln in o:
        m = PM.match(ln.strip())
        if m: po.setdefault(m.group(1).strip(), []).append((float(m.group(2)), float(m.group(3))))
    pn = {}
    for ln in n:
        m = PM.match(ln.strip())
        if m: pn.setdefault(m.group(1).strip(), []).append((float(m.group(2)), float(m.group(3))))
    for lab in po:
        for (xo, so), (xn, sn) in zip(po[lab], pn.get(lab, [])):
            sg = math.hypot(so, sn); z = abs(xn - xo) / sg if sg > 0 else (0.0 if xn == xo else float('inf')); okZ &= z <= 3
            if xn != xo: pr(f'    {lab[:90]}: {xo:.6g} +/- {so:.3g} -> {xn:.6g} +/- {sn:.3g}  ({(xn / xo - 1) * 100 if xo else float("nan"):+.3f} %, {z:.2f} sigma)')
        if len(po[lab]) != len(pn.get(lab, [])): okZ = False; pr(f'    {lab[:90]}: liczba wierszy {len(po[lab])} -> {len(pn.get(lab, []))}')
    d = [(a, b) for a, b in zip(o, n) if a != b and not PM.match(a.strip())]
    pr(f'    inne rozne wiersze: {len(d)}' + ''.join(f'\n      old: {a.strip()[:150]}\n      new: {b.strip()[:150]}' for a, b in d[:8]))
pr(f"H852z (GLOWNA: kazde X +/- Y w zestawach DIFF w 3 sigma): {'TAK' if okZ else 'NIE'}")
pr(f"DECYZJA: {'sciezka stabilna domyslna' if okM and okO and okV and okZ else 'powrot do f2216d77 (przelacznik domyslnie wylaczony)'}")

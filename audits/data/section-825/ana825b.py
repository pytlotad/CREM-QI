# Audit 825, added AFTER reading ana825.txt (UZUPELNIENIE 1): the delay measure named in the prereg (MIARY: time from the first refused
# photon to the payout), not implemented in ana825.py.  Delay in units of the rule lifetime of the refusing state: tau = 1/(SCALE*downE1).
import os, re, glob
D = os.path.dirname(os.path.abspath(__file__)); SCALE = 8e5
kv = lambda s: dict(re.findall(r'([A-Za-z_][A-Za-z_0-9]*)=(\S+)', s))
for p in sorted(glob.glob(f'{D}/out/runs/*.err')):
    L = open(p).read().splitlines(); first = None; down = None; paid = None; nmin = None
    ph = [ln for ln in L if ln.startswith('PHASE787')]; lad = [ln for ln in L if ln.startswith('LADDER E_rel')]; rul = [ln for ln in L if ln.startswith('RULES')]
    for ln, lr, pl in zip(lad, rul, ph):
        d = kv(ln); t = float(pl.split()[3])
        if float(d['photon']) == 0.0 and first is None: first = t; down = float(kv(lr)['downE1'])
        if float(d['photon']) == 0.0: nmin = min(nmin or 9, float(kv(lr)['n']))
        if first is not None and float(d['photon']) > 0 and paid is None: paid = t
    if first is None: continue
    end = float([ln for ln in L if ln.startswith('END')][0].split()[6]); tau = 1/(SCALE*down)
    nref = sum(1 for ln in lad if float(kv(ln)['photon']) == 0.0)
    dl = f'{(paid - first)/tau:.1f} tau' if paid else f'> {(end - first)/tau:.1f} tau (cenzura OBS)'
    print(f'  {os.path.basename(p)[:-4]}: odrzucen {nref}, pierwsze t = {first:.3e} s przy n = {nmin:.6f} (najnizsze), tau(2p, skala) = {tau:.3e} s, opoznienie {dl}')

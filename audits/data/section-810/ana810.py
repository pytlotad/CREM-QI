# Audit 810: verdicts H810o, H810p, H810e, H810f, H810h, H810v, H810n, H810m -> out/ana810.txt
# Inputs (runs/all.sh): out/p810_{new,new_off,old}.txt, out/probes810.txt, out/probes_{old,new}/*.err, out/validation_{old,new}.txt, out/regression810.txt
import os, re, math
from mpmath import mp, mpf, sqrt as msqrt
mp.dps = 40
D = os.path.dirname(os.path.abspath(__file__)); O = os.path.join(D, 'out')
h = 6.62607015e-34; HEV = 4.135667696e-15; C = 299792458.0
NU = 10969.045; SIG = math.hypot(0.015, 0.020)
pr = lambda *a: print(*a, flush=True)
rd = lambda f: open(os.path.join(O, f)).read()
# --- model (mpmath, as prereg)
al = mpf(1)/mpf('137.035999177'); me = mpf('9.1093837139e-31'); mpp = mpf('1.67262192595e-27'); cc = mpf(299792458); hh = mpf('6.62607015e-34')
mu = me*mpp/(me + mpp); lam = mpf('2.00231930436092') - 1
def ES(n, L): return mu*cc**2*(1/msqrt(1 + al**2/(n - L + msqrt(L**2 - al**2))**2) - 1)
def E4(n, L): return -mu*cc**2*al**2/(2*n**2) - mu*cc**2*al**4/(2*n**4)*(n/L - mpf(3)/4)
def R(n, L): return float((ES(n, L) - E4(n, L))/hh/10**6)   # MHz
def lev(n, l, j2):   # shift of level (n, l, j = j2/2) in MHz: S via k = 1; l >= 1: R(l + 1/2) + lam [R(k) - R(l + 1/2)]
    if l == 0: return R(n, 1)
    k = (j2 + 1)//2; return R(n, l + 0.5) + float(lam)*(R(n, k) - R(n, l + 0.5))
# --- H810o / H810p
reg = dict((' '.join(l.split()[:4]), l.split()[4]) for l in rd('regression810.txt').splitlines())
prb = dict((' '.join(l.split()[:5]), l.split()[5]) for l in rd('probes810.txt').splitlines())
tags = ('1', '2', 'L2', '5Ps', '5H', '6'); probes = ('H2p', 'H3p', 'H4p', 'H3s', 'H3d', 'Ps2p')
o_reg = all(reg[f'old vs new_off {t}'] == 'SAME' for t in tags); o_prb = all(prb[f'probe {t} old vs new_off'] == 'SAME' for t in probes)
o_p810 = rd('p810_old.txt') == rd('p810_new_off.txt')
pr(f"H810o (new_off = old bajt w bajt): {'TAK' if o_reg and o_prb and o_p810 else 'NIE'} (regresja {sum(reg[f'old vs new_off {t}'] == 'SAME' for t in tags)}/6, "
   f"sondy p585 {sum(prb[f'probe {t} old vs new_off'] == 'SAME' for t in probes)}/6, p810 {'bajt w bajt' if o_p810 else 'ROZNE'})")
pr('  regresja old vs new: ' + ', '.join(f'{t} {reg[f"old vs new {t}"]}' for t in tags) + '; sondy old vs new: ' + ', '.join(f'{t} {prb[f"probe {t} old vs new"]}' for t in probes))
okP = all(reg[f'old vs new {t}'] == 'SAME' for t in ('1', '2', 'L2', '5Ps', '6')) and prb['probe Ps2p old vs new'] == 'SAME'
pr(f"H810p (Ps bez zmian: zestawy 1 p-Ps, 2 o-Ps, L2, 5Ps, 6 i sonda Ps2p bajt w bajt): {'TAK' if okP else 'NIE'}")
# --- H810e
def splits(f): return [float(m.group(1)) for m in (re.search(r'= (-?[\d.]+) MHz', x) for x in rd(f).splitlines() if x.startswith('split')) if m]
def grid(f): return {(int(w[2]), float(w[4])): float(w[6]) for w in (x.split() for x in rd(f).splitlines() if x.startswith('grid'))}
def levels(f): return {(w[1], int(w[3])): (float(w[7]), float(w[10])) for w in (x.split() for x in rd(f).splitlines() if x.startswith('level'))}
sN, sF = splits('p810_new.txt'), splits('p810_new_off.txt')
OLD763 = [10968.704747, 3249.986592, 1083.328864, 1371.088093, 457.029364, 228.514682]
NL = [(2, 1), (3, 1), (3, 2), (4, 1), (4, 2), (4, 3)]
pred = [OLD763[i] + float(lam)*(R(n, l + 1) - R(n, l)) for i, (n, l) in enumerate(NL)]
gN, gF = grid('p810_new.txt'), grid('p810_new_off.txt')
dg = 0.0
for (n, L), v in gN.items():
    eng = (v - gF[(n, L)])/h/1e6; mod = R(n, max(L, 1.0)) if L < 1 else R(n, L); dg = max(dg, abs(eng - mod))
ds = max(abs(a - b) for a, b in zip(sN, pred)); doff = max(abs(a - b) for a, b in zip(sF, OLD763))
lN, lF = levels('p810_new.txt'), levels('p810_new_off.txt')
dfine = {k: lN[k][1] - lF[k][1] for k in lN}
c12 = dfine[('2S1/2', 1)] - dfine[('1S1/2', 1)]; cL = dfine[('2S1/2', 1)] - dfine[('2P1/2', 1)]
p12 = R(2, 1) - R(1, 1); pL = R(2, 1) - lev(2, 1, 1)
okE = dg <= 1e-6 and ds <= 5e-4 and doff <= 1e-6 and abs(c12 - p12) <= 5e-4 and abs(cL - pL) <= 5e-4
pr(f"H810e: {'TAK' if okE else 'NIE'} -- R silnika vs mpmath max {dg:.2e} MHz ({len(gN)} punktow siatki); rozszczepienia new vs przewidywanie max {ds:.2e}; "
   f"new_off vs 763 max {doff:.2e}; 1S - 2S zmiana {c12:+.6f} (przew. {p12:+.6f}); Lamb 2S - 2P1/2 zmiana {-cL:+.6f}... (2S - 2P1/2 {cL:+.6f}, przew. {pL:+.6f}) MHz")
for i, (n, l) in enumerate(NL):
    pr(f'  n {n} l {l}: new {sN[i]:.6f}, przewidywanie {pred[i]:.6f}, new_off {sF[i]:.6f}, zmiana {sN[i] - sF[i]:+.6f} MHz')
# --- H810f
dF = max(abs((lN[(nm, 1)][0] - lF[(nm, 1)][0]) - (lN[(nm, 0)][0] - lF[(nm, 0)][0])) for nm in {k[0] for k in lN})
pr(f"H810f (zmiana poziomu F = 1 = F = 0 do 1e-9 MHz): {'TAK' if dF <= 1e-9 else 'NIE'} (max {dF:.2e})")
for nm in ('1S1/2', '2S1/2', '2P1/2', '2P3/2', '3S1/2', '3P1/2', '3P3/2', '3D3/2', '3D5/2'):
    pr(f'  {nm}: zmiana F = 1 {lN[(nm, 1)][0] - lF[(nm, 1)][0]:+.6f}, F = 0 {lN[(nm, 0)][0] - lF[(nm, 0)][0]:+.6f} MHz')
# --- H810h: photons of the hydrogen probes
def led(f): return [(float(m.group(1)), float(m.group(2))) for m in (re.search(r'CREM_LEDGER t (\S+) Eg (\S+)', l) for l in open(f)) if m]
LV = {(n, l, j2): lev(n, l, j2) for n in range(1, 5) for l in range(0, n) for j2 in ((1,) if l == 0 else (2*l - 1, 2*l + 1))}
okH = True
for t in ('H2p', 'H3p', 'H4p', 'H3s', 'H3d'):
    a = led(os.path.join(O, 'probes_old', t + '.err')); b = led(os.path.join(O, 'probes_new', t + '.err'))
    same_t = len(a) == len(b) and all(x[0] == y[0] for x, y in zip(a, b)); okH &= same_t
    for (ta, ea), (tb, eb) in zip(a, b):
        d = (eb - ea)/h/1e6
        cand = sorted({round(LV[u] - LV[v], 6) for u in LV for v in LV if u[0] > v[0] and abs(u[1] - v[1]) == 1 and abs(u[2] - v[2]) <= 2})
        best = min(cand, key=lambda x: abs(x - d)); okH &= abs(best - d) <= 0.003
        pr(f'  {t}: t {ta:.6e} (czasy {"identyczne" if same_t else "ROZNE"}), Eg old {ea*HEV/h:.9f} eV, zmiana {d:+.4f} MHz, najblizsze przejscie modelu {best:+.4f} MHz')
pr(f"H810h (czasy identyczne, zmiany = roznice R do 0,003 MHz): {'TAK' if okH else 'NIE'}")
# --- H810v
vN, vO = rd('validation_new.txt'), rd('validation_old.txt')
VOL = ('speedup', 'cost:', 'wall', 'rate:', 's/ps', 'worst |tau', 'tau/tau_1s')
ROWS = ('H 2P3/2 - 2P1/2', 'H 2p -> 1s')
fil = lambda s: [x for x in s.splitlines() if not any(k in x for k in VOL + ROWS)]
okV = 'EXIT 0' in vN and 'Enforced checks:     72/72' in vN and fil(vN) == fil(vO)
pr(f"H810v (walidacja new EXIT 0, 72/72, poza 3 wierszami wodoru i liniami zmiennymi = old): {'TAK' if okV else 'NIE'}")
for x in vN.splitlines():
    if any(k in x for k in ROWS): pr('  new ' + x.strip())
for x in vO.splitlines():
    if any(k in x for k in ROWS): pr('  old ' + x.strip())
# --- H810n, H810m
nu, nu0 = sN[0], sF[0]; z = (nu - NU)/SIG
units = lambda f: f'{f:.4f} MHz = {f*1e6*HEV*1e6:.4f} ueV = {C/(f*1e6)*100:.4f} cm = {f/NU:.6f} pomiaru'
okN = abs(z) <= 3 and abs(nu - NU) < abs(nu0 - NU)
for F in (1, 0):
    pr(f'=== F = {F} (PHEN={F + 1}) === 2P3/2 - 2P1/2 new = {units(nu)}; {nu - NU:+.4f} MHz = {z:+.2f} sigma; przed 810 {nu0:.4f} ({nu0 - NU:+.4f} MHz = {(nu0 - NU)/SIG:+.2f} sigma)')
    pr(f'  o-Ps przez czynnik z wodoru: {nu/NU*0.9998:.5f}-{nu/NU*0.9999:.5f} multipletu (Hagena i in., PRL 71, 2887 (1993)); p-Ps bez czesci spinowej (Ps bez zmian: H810p)')
pr(f"H810n (GLOWNA; F = 1 i F = 0 -- zmiana niezalezna od F, H810f): {'TAK' if okN and dF <= 1e-9 else 'NIE'}")
M809 = 10943.7009 + 25.3678
pr(f"H810m (|2P new - mechanika 809 {M809:.4f}| <= 0,010 MHz): {'TAK' if abs(nu - M809) <= 0.010 else 'NIE'} ({nu - M809:+.4f} MHz)")
pr(f'1S - 2S wodoru: zmiana {c12:+.4f} MHz; reszta wiersza [W] Model.md -3,94 -> {-3.94 + c12:+.2f} MHz (Parthey i in., PRL 107, 203001 (2011))')
pr('Pomiar 2P: 10 969,045(15) MHz (Hagley i Pipkin, PRL 72, 1172 (1994); Lundeen i Pipkin, PRL 46, 232 (1981)); import QED 0,000 +- 0,020 MHz '
   '(Blumer i in., arXiv:2412.19580, Tab. 1; Horbatsch i Hessels, PRA 93, 022513 (2016)).')
pr(f'DECYZJA: {"zmiana zostaje domyslna" if (o_reg and o_prb and o_p810 and okP and okE and okV and okN) else "zmiana NIE spelnia kryteriow -- do wycofania / poprawy"}')

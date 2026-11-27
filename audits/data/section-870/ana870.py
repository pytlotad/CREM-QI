# Audit 870: phase closure for a field without a mode basis (free / local field with a clock) -- verdicts H870v, d, p, m, f, w, e, r -> out/ana870.txt.
# Written after prereg870.txt, before any result.  hbar = mu = 1, e^2 = 1, beta = 0.  Closure: C = Delta phi_s + omega_f T + mu = 2 pi n.
import math, numpy as np
from decimal import Decimal as Dm, getcontext
getcontext().prec = 60
pr = lambda *a: print(*a, flush=True)
PI = math.pi
# ---------- (1) Kepler integrals in time (uniform mean anomaly, Kepler equation by Newton) ----------
def kepler_integrals(N, e, M=200000):
    a = N * N; Mg = 2 * PI * np.arange(M) / M; u = Mg.copy()
    for _ in range(60): u -= (u - e * np.sin(u) - Mg) / (1 - e * np.cos(u))
    r = a * (1 - e * np.cos(u)); dt = 2 * PI * N ** 3 / M
    Tk = 1 / r - 1 / (2 * a); V = -1 / r
    return Tk.sum() * dt, V.sum() * dt, (2 * Tk).sum() * dt
okv = True
pr("(1) calki po orbicie Keplera (numerycznie, czas jednostajny):")
for N in (1, 2):
    for e in (0.0, 0.5, 0.9):
        iT, iV, ipq = kepler_integrals(N, e)
        dv = max(abs(iT / (PI * N) - 1), abs(iV / (-2 * PI * N) - 1), abs(ipq / (2 * PI * N) - 1)); okv &= dv < 1e-9
        pr(f"  N {N} e {e:.1f}: int T dt / (pi N) {iT / (PI * N):.12f}, int V dt / (-2 pi N) {iV / (-2 * PI * N):.12f}, oint p dq / (2 pi N) {ipq / (2 * PI * N):.12f}  (max odch. {dv:.1e})")
pr(f"H870v (calki = analityczne do 1e-9): {'TAK' if okv else 'NIE'}")
# ---------- (2) spectra ----------
NMAX = 8.0; TOL = 1e-9
def classify(sols):
    sols = sorted(s for s in sols if 0 < s <= NMAX + 1e-12)
    ints = sorted({round(s) for s in sols if abs(s - round(s)) < TOL and round(s) >= 1})
    spur = [s for s in sols if abs(s - round(s)) >= TOL or round(s) < 1]
    full = ints == list(range(1, int(NMAX) + 1))
    cls = 'DOKLADNIE BOHR' if full and not spur else 'NADZBIOR' if full else 'PODZBIOR' if ints and not spur else 'OBOK'
    return cls, ints, spur, sols
def lin_solutions(c, mu):   # pi c N + mu = 2 pi k
    out = []
    for k in range(-200, 400):
        N = (2 * k - mu / PI) / c
        if 0 < N <= NMAX + 1e-12: out.append(N)
    return out
def cubic_solutions(cs, w0, mu):   # pi cs N + 2 pi w0 N^3 + mu = 2 pi k
    Ng = np.linspace(1e-9, NMAX, 20001); C = (cs * Ng + 2 * w0 * Ng ** 3 + mu / PI) / 2
    out = []
    for k in range(math.floor(C.min()) - 1, math.ceil(C.max()) + 2):
        for z in np.roots([2 * w0, 0, cs, mu / PI - 2 * k]):
            if abs(z.imag) < 1e-9 and 0 < z.real <= NMAX + 1e-12: out.append(float(z.real))
    return out
SRC = {'S1 lagranzjan': 3.0, 'S2 czas wlasny': 1.0, 'S3 akcja skrocona': 2.0}
MUS = {'0': 0.0, '-pi': -PI, '-2pi': -2 * PI}
res = {}
pr("\n(2) klasy widma (N w (0, 8]):")
pr("  zegar               pole  mu     c        klasa            calkowite                     falszywe (pierwsze)")
for sn, cs in SRC.items():
    for fn, cf in (('F0', 0.0), ('FE', -1.0), ('FD', cs * cs / 4)):
        for mn, mu in MUS.items():
            c = cs + cf
            if abs(c) < 1e-12:   # degenerate: C = mu for every orbit -- all N close (continuum) or none
                cls = 'CIAGLE (kazde N)' if abs(math.remainder(mu, 2 * PI)) < 1e-12 else 'BRAK'; res[(sn, fn, mn)] = (cls, [])
                pr(f"  {sn:19s} {fn:4s}  {mn:5s}  {c:7.4f}  {cls:15s}  (c = 0: faza niezalezna od N)"); continue
            cls, ints, spur, sols = classify(lin_solutions(c, mu)); res[(sn, fn, mn)] = (cls, sols)
            pr(f"  {sn:19s} {fn:4s}  {mn:5s}  {c:7.4f}  {cls:15s}  {str(ints):28s}  {', '.join(f'{s:.4f}' for s in spur[:5])}{' ...' if len(spur) > 5 else ''} ({len(spur)})")
for sn in ('S1 lagranzjan', 'S3 akcja skrocona'):
    for w0 in (-0.5, 0.5, 1.0):
        for mn, mu in MUS.items():
            cls, ints, spur, sols = classify(cubic_solutions(SRC[sn], w0, mu)); res[(sn, f'FW{w0:+.1f}', mn)] = (cls, sols)
            pr(f"  {sn:19s} FW w0 {w0:+.1f} {mn:5s}  {cls:15s}  {str(ints):28s}  falszywych {len(spur)} (np. {', '.join(f'{s:.4f}' for s in spur[:4])})")
# H870p
cls0 = res[('S1 lagranzjan', 'FE', '-2pi')][0]
ph = {}
for nm, mu in (('-5pi/3', -5 * PI / 3), ('-7pi/3', -7 * PI / 3)):
    cls, ints, spur, sols = classify(lin_solutions(2.0, mu)); ph[nm] = ints
    pr(f"  S1 FE z faza 869 mu = {nm}: {cls}, calkowite {ints}, rozwiazania {', '.join(f'{s:.4f}' for s in sols[:5])} ...")
okp = cls0 == 'DOKLADNIE BOHR' and all(not v for v in ph.values())
pr(f"H870p (S1 FE mu -2pi = Bohr; fazy 869 bez calkowitych N): {'TAK' if okp else 'NIE'}")
bohr = {k for k, v in res.items() if v[0] == 'DOKLADNIE BOHR'}
exp = {('S1 lagranzjan', 'FE', '0'), ('S1 lagranzjan', 'FE', '-2pi'), ('S3 akcja skrocona', 'F0', '0'), ('S3 akcja skrocona', 'F0', '-2pi')}
pr(f"  DOKLADNIE BOHR: {sorted(bohr)}")
pr(f"H870m (GLOWNA; Bohr tylko S1 FE i S3 F0 przy mu = 0 mod 2 pi): {'TAK' if bohr == exp else 'NIE'}")
csx = 2 * math.sqrt(3) - 2
okf = all(res[(sn, 'FD', mn)][0] != 'DOKLADNIE BOHR' for sn in SRC for mn in MUS)
pr(f"  FD: c = c_s + c_s^2 / 4 = 2 wymaga c_s* = 2 sqrt 3 - 2 = {csx:.4f} (sprawdzenie: {csx + csx * csx / 4:.12f})")
N, e = 1.0, 0.5; a = N * N; rr = np.array([a * (1 - e), a * (1 + e)]); wf = 1 / rr - 1 / (2 * a); Tor = 2 * PI * N ** 3
pr(f"  FD na orbicie e 0,5 (S3, k = p(x)): omega_f(x) T od {wf.min() * Tor:.3f} do {wf.max() * Tor:.3f} rad -- rozrzut {(wf.max() - wf.min()) * Tor:.3f} rad (> 2 pi: domkniecie punktowe niemozliwe)")
pr(f"H870f (FD nigdy Bohr dla c_s 1, 2, 3): {'TAK' if okf else 'NIE'}")
okw = all(v[0] != 'DOKLADNIE BOHR' for k, v in res.items() if k[1].startswith('FW'))
pr(f"H870w (FW nigdy dokladnie Bohr): {'TAK' if okw else 'NIE'}")
# ---------- (3) 1S-2S vs measurement ----------
RINF = 10973731.568157 * 299792458.0; MEMP = 5.446170214889e-4; ALPHA = 7.2973525643e-3; C0 = 299792458.0; H_EVS = 4.135667696e-15
meas = {'wodor': (2466061413.187e6, 1 / (1 + MEMP), 'Parthey 2011'), 'pozytonium': (1233607216.4e6, 0.5, 'Fee 1993')}
pr("\n(3) 1S-2S z dwoch najnizszych poziomow kazdej klasy (poziomy grube bez spinu: jednakowe dla F = 0 / F = 1 i p-Ps / o-Ps):")
for ch, (fm, red, cit) in meas.items():
    fb = RINF * red * 0.75
    pr(f"  {ch}: pomiar {fm / 1e6:.3f} MHz = {fm * H_EVS:.4f} eV, {C0 / fm * 1e9:.3f} nm ({cit}); Bohr {fb / 1e6:.3f} MHz, odch. wzgl. {fb / fm - 1:+.2e}")
oke = True
for k, (cls, sols) in sorted(res.items()):
    s = sorted(set(round(x, 9) for x in sols))
    if len(s) < 2: continue
    rel = {}
    for ch, (fm, red, _) in meas.items():
        f = RINF * red * (1 / s[0] ** 2 - 1 / s[1] ** 2); rel[ch] = f / fm - 1
    good = all(abs(v) < 1e-4 for v in rel.values()); far = all(abs(v) >= 0.10 for v in rel.values())
    oke &= good if cls == 'DOKLADNIE BOHR' else far
    if cls == 'DOKLADNIE BOHR' or k[2] == '0':
        pr(f"  {k[0]:19s} {k[1]:7s} mu {k[2]:5s} {cls:15s} N1 {s[0]:.4f} N2 {s[1]:.4f}: wodor {rel['wodor']:+.3e}, pozytonium {rel['pozytonium']:+.3e}")
pr(f"H870e (1S-2S: Bohr w 1e-4, pozostale >= 10 % od pomiaru, oba uklady): {'TAK' if oke else 'NIE'}")
# ---------- (5) relativistic n = 2 ----------
A = Dm('7.2973525643e-3'); R = Dm('10973731.568157') * Dm(299792458); mec2 = 2 * R / (A * A); muc2 = mec2 / (1 + Dm('5.446170214889e-4'))
def f(n, k): return (1 + A * A / (Dm(n) - k + (k * k - A * A).sqrt()) ** 2) ** Dm('-0.5')
som = muc2 * (f(2, Dm(2)) - f(2, Dm(1))) / Dm(10 ** 6); kg = muc2 * (f(2, Dm('1.5')) - f(2, Dm('0.5'))) / Dm(10 ** 6)
pr("\n(5) przedzial n = 2 z domkniecia S3 z pedem relatywistycznym (masa zredukowana, mu c^2 (f - 1)):")
for nm, v in (('k calkowite (Sommerfeld = Dirac, j = k - 1/2): 2P3/2 - 2P1/2', som), ('k polowkowe (Langer + Maslov = Klein-Gordon, bez spinu): 2S - 2P', kg)):
    fv = float(v)
    pr(f"  {nm}: {fv:.4f} MHz = {fv * 1e6 * H_EVS * 1e6:.3f} ueV, {C0 / (fv * 1e6) * 100:.4f} cm")
okd = abs(float(som) - 10943.6883) < 0.001
pr(f"H870d (Sommerfeld k = 1, 2 = Dirac 866 10 943,6883 MHz do 0,001): {'TAK' if okd else 'NIE'} ({float(som) - 10943.6883:+.5f} MHz)")
fm = 10969.045
pr(f"  pomiar 10 969,045(15) MHz (Hagley i Pipkin 1994) = {fm * 1e6 * H_EVS * 1e6:.3f} ueV, {C0 / (fm * 1e6) * 100:.4f} cm;  k calkowite: {float(som) - fm:+.4f} MHz = {(float(som) / fm - 1) * 100:+.3f} % "
   f"({(float(som) - fm) / 0.015:+.0f} sigma; reszta = QED, Horbatsch i Hessels 2016);  k polowkowe / pomiar = {float(kg) / fm:.4f}")
okr = abs(float(som) / fm - 1) < 0.005 and float(kg) / fm > 2
pr(f"H870r (k calkowite w 0,5 %; k polowkowe > 2 razy dalej): {'TAK' if okr else 'NIE'}")

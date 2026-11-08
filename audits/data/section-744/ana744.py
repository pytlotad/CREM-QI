# Audit 744: analysis of out/runs -> out/ana744.txt
import sys, os, glob, numpy as np
D = os.path.dirname(os.path.abspath(__file__)); O = os.path.join(D, 'out', 'runs')
sys.path.insert(0, os.path.join(D, '..', 'section-743'))
from modes743 import torus, modes2d, response
pr = lambda *a: print(*a, flush=True)
ld = lambda f: np.loadtxt(os.path.join(O, f), comments='#', ndmin=2)
cmean = lambda a: np.angle(np.mean(np.exp(1j * a)))
G = 0.1 / torus(0.5, 1.0, 0.2)['Tr']; M2 = modes2d(0.2, 12, 6)

pr("=== H744a: kappa = 0 (faza cechowania, pole przyczynowe) wobec stanu ustalonego 743 M ===")
w = 0
for J in ['0.30', '0.40', '0.50', '0.60', '0.70']:
    d = ld(f'k0_Jr{J}.txt'); P = d[-10:, 5].mean(); PM = response(torus(float(J), 1.0, 0.2), M2, G)[0]
    w = max(w, abs(P / PM - 1)); pr(f"  J_r {J}: P = {P:.6e}, 743 M = {PM:.6e}, P/M - 1 = {P/PM-1:+.2e}; max |slip| = {np.abs(d[:, 4]).max():.1e}")
pr(f"  max |P/M - 1| = {w:.2e} (kryterium 1e-2)")

def scan(kind, tau, sg):
    pat = f'scan{kind}_tau{tau}_s{sg}_' + ('Jr' if kind == 'J' else 'L')
    fs = sorted(glob.glob(os.path.join(O, pat + '*.txt')))
    xs = np.array([float(os.path.basename(f)[len(pat):-4]) for f in fs])
    P = np.array([np.loadtxt(f)[-10:, 5].mean() for f in fs]); S = np.array([cmean(np.loadtxt(f)[-20:, 4]) for f in fs])
    Sd = np.array([np.std(np.loadtxt(f)[-20:, 4]) for f in fs]); M = np.array([np.loadtxt(f)[-10:, 6].mean() for f in fs])
    return xs, P, S, Sd, M

for tau in ['0.5', '0.25', '1.0']:
    pr(f"=== tau = {tau} T_r0, kappa T_r0 = 10 ===")
    for sg, name in [('1', 'harmonia'), ('-1', 'antyharmonia')]:
        for kind, ebk in [('J', 0.5), ('L', 1.0)]:
            xs, P, S, Sd, M = scan(kind, tau, sg); i0 = np.argmin(abs(xs - ebk)); im = P.argmax()
            pr(f"  {name}, skan {'J_r (L = 1)' if kind == 'J' else 'L (J_r = 0,5)'}: max P przy {xs[im]:.2f} ({'w 0,02 od EBK' if abs(xs[im]-ebk) <= 0.02 + 1e-9 else 'NIE przy EBK'}); "
               f"P(EBK)/max = {P[i0]/P.max():.3f}; poslizg w EBK {S[i0]:+.4f} rad (|s| - pi = {abs(S[i0]) - np.pi:+.4f}, rozrzut {Sd[i0]:.3f})")
            pr("    x     P            slip     sd     <|psi_mem|>")
            for x, p, s, sd, m in zip(xs, P, S, Sd, M): pr(f"    {x:.2f}  {p:.4e}  {s:+.4f}  {sd:.3f}  {m:.4e}")

pr("=== H744e: sila zwrotna (kappa T 10, tau T/2, harmonia, eps 3e-5, 300 okresow) ===")
for f in sorted(glob.glob(os.path.join(O, 'br_*.txt'))):
    d = np.loadtxt(f, comments='#', ndmin=2); n = len(d) // 10 * 10
    J = d[:n, 2].reshape(-1, 10).mean(1); L = d[:n, 3].reshape(-1, 10).mean(1)
    rJ, rL = np.ptp(J[-20:]), np.ptp(L[-20:]); lock = rJ < 1e-3 and rL < 1e-3 and abs(J[-1] - .5) < 1e-3 and abs(L[-1] - 1) < 1e-3
    pr(f"  {os.path.basename(f)}: okresy {len(d)}; J_r " + ' '.join(f"{x:.4f}" for x in J[::5]) + "; L " + ' '.join(f"{x:.4f}" for x in L[::5])
       + f"; zakres ost. 200: {rJ:.2e} / {rL:.2e}; poslizg (ost. 20) {cmean(d[-20:, 4]):+.3f}; zatrzasniecie: {'TAK' if lock else 'NIE'}")

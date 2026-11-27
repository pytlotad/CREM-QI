# 871 POST HOC (po biegach, po ana871.py): kierunek i tempo dryfu J wobec najblizszego celu;
# porownanie z ksztaltem sily z wyprowadzenia preregu: dJ/dt = -(gamma_r) * S(C) / (C' * S2),
# S(C) = sum_k k x^k sin(kC), S2 = sum_k k^2 x^k, x = e^{-Gamma T}, gamma_r T = 0,2 / N^2 (dobor kappa).
import numpy as np, glob, os
O = 'out/runs'
def load(name):
    hdr = {}; rows = []
    for ln in open(f'{O}/{name}.txt'):
        if ln.startswith('# 871'):
            t = ln.split()[2:]; hdr = {t[i]: float(t[i + 1]) for i in range(0, len(t) - 1, 2)}
        elif not ln.startswith('#') and ln.strip(): rows.append([float(x) for x in ln.split()])
    return hdr, rows
GAM = 1e-3
def S(C, x):
    z = x * np.exp(1j * C)
    return (z / (1 - z) ** 2).imag, x * (1 + x) / (1 - x) ** 3
def target(J, clock):
    if clock == 3: return round(J)
    return round(J * 1.5) / 1.5
print('bieg          N0      J_pocz    J_kon     cel     dJ_obs       kierunek  rate_obs     rate_pred(ksztalt)  obs/pred')
stat = {}
for f in sorted(glob.glob('out/runs/*.txt')):
    name = os.path.basename(f)[:-4]
    grp = name.split('_')[0]
    if grp not in ('GH', 'GP', 'GL', 'GZH', 'GZP', 'RH', 'RN', 'RP', 'SH0', 'SH1', 'SP0', 'SP1'): continue
    hdr, rows = load(name)
    if hdr.get('KAPPA', 0) == 0 or len(rows) < 20: continue
    a = np.array(rows); t = a[:, 1]; J = a[:, 2]
    clock = int(hdr['CLOCK']); mu = hdr['M1'] * hdr['M2'] / (hdr['M1'] + hdr['M2'])
    J0 = J[:5].mean(); J1 = J[-max(5, len(J) // 20):].mean()
    tg = target(J0, clock); good = np.sign(J1 - J0) == np.sign(tg - J0)
    h = len(J) // 2; p = np.polyfit(t[h:], J[h:], 1); robs = p[0]
    # prognoza ksztaltu (tylko NR, orbita z J0)
    rpred = float('nan')
    if hdr['MODE'] == 0:
        N = J0; T = 2 * np.pi * N ** 3 / mu; x = np.exp(-GAM * T)
        Cp = 2 * np.pi if clock == 3 else 3 * np.pi
        C = Cp * N; s, s2 = S(C, x)
        rpred = -(0.2 / N ** 2 / T) * s / (Cp * s2)  # linearnie: -(gamma_r) (N - cel)
    stat.setdefault(grp, []).append(good)
    print(f'{name:13s} {hdr["N0"]:5.2f}  {J0:9.5f} {J1:9.5f} {tg:7.4f}  {J1 - J0:+.3e}  {"do celu" if good else "OD celu":8s}  {robs:+.3e}   {rpred:+.3e}          {robs / rpred if rpred == rpred else float("nan"):+.3f}')
print()
for g, v in stat.items(): print(f'{g}: dryf w strone najblizszego celu {sum(v)}/{len(v)}')

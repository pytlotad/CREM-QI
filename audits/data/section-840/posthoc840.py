# Audit 840 POST HOC (written after ana840): (1) secular rate of the spinless residual (recoil m/M of the x1e6 nucleus);
# (2) Dirac reference with the BMT factor 1/(gamma+1) frozen at 1/2 (P840_HALF=1).
import os, math, numpy as np
D = os.path.dirname(os.path.abspath(__file__)); P = os.path.join(D, 'out', 'posthoc')
def gauss(N):
    x, w = np.polynomial.legendre.leggauss(N); u = 0.75 + 0.25*x[::-1]; return 1/u, 0.25*w[::-1]
NODES = dict(zip([f'u5{c}' for c in 'abcde'], zip(*gauss(5))))
nu_of = lambda Om: sum(w*L**2*Om[k] for k, (L, w) in NODES.items())/(2*math.pi)
def rd(f, key):
    o = {}
    for ln in open(f):
        w = ln.split()
        if key == 'SPINLESS' and ln.startswith('SPINLESS'): o['s'] = -float(w[2])
        if key == 'OMX' and ln.startswith('OMX'): o[w[1]] = -float(w[3])
    return o
nu = {}
for sk in ('H0', 'H1'):
    sp = nu_of({nd: rd(f'{P}/{sk}_{nd}_spinless.txt', 'SPINLESS')['s'] for nd in NODES})
    print(f'{sk}: spinless residual (odrzut jadra x1e6) nu {sp/1e6:+.2e} MHz')
    for var in ('rebuilt', 'half'):
        O = {nd: rd(f'{P}/{sk}_{nd}_{var}.txt', 'OMX') for nd in NODES}
        for x in ('Z', 'ZT', 'Z0', 'ZT0'): nu[sk, var, x] = nu_of({nd: O[nd][x] for nd in NODES})/1e6
for var, lab in (('rebuilt', 'Dirac FW (1/(gamma+1))'), ('half', 'odniesienie z 1/2 zamiast 1/(gamma+1)')):
    print(lab)
    for x in ('Z', 'ZT', 'Z0', 'ZT0'):
        a, b = nu['H0', var, x], nu['H1', var, x]
        print(f'  {x:4s} F0 {a:+.4f} F1 {b:+.4f} parzysta {(a+b)/2:+.4f} nieparzysta {(b-a)/2:+.4f} MHz')
print('udzial czynnika 1/(gamma+1) (Dirac - odniesienie 1/2), nieparzysta: ' +
      f"{((nu['H1','rebuilt','Z']-nu['H0','rebuilt','Z']) - (nu['H1','half','Z']-nu['H0','half','Z']))/2:+.4f} MHz, parzysta "
      f"{((nu['H1','rebuilt','Z']+nu['H0','rebuilt','Z']) - (nu['H1','half','Z']+nu['H0','half','Z']))/2:+.4f} MHz")

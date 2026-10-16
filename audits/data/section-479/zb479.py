# Audit 479: near-zone Zitterbewegung echo -- phase selection rule, implied pair states, and time available in one pass.
import math, random
hbar=1.054571817e-34; c=299792458.0; me=9.1093837015e-31; e=1.602176634e-19; eps0=8.8541878128e-12; mu0=1.25663706212e-6; muB=9.2740100783e-24
k=2*me*c/hbar; wZB=2*me*c*c/hbar; kC=e*e/(4*math.pi*eps0); rstar=1/k; mperp2=2*muB**2
# (1) near-zone (static) coupling of two circularly rotating transverse moments, separation at angle chi to the cone axis:
# U(dphi) = -(mu0/4pi)(m_perp^2/2)(1/r^3) * (1/2)[(3 sin^2 chi) - 2] cos dphi  -> sign flips at sin^2 chi = 2/3 (chi = 54.74 deg)
def sel(chi):
    f=3*math.sin(chi)**2-2     # in-plane (chi = 90): +1 -> min at dphi = 0;  along axis (chi = 0): -2 -> min at pi
    return 0.0 if f>0 else math.pi
print(f'(1) selection: chi < {math.degrees(math.asin(math.sqrt(2/3))):.2f} deg -> dphi = pi; otherwise 0;  fraction pi for isotropic axis = 1 - 1/sqrt3 = {1-1/math.sqrt(3):.4f}')
random.seed(479); N=400000; cnt={}
for _ in range(N):
    cz=2*random.random()-1; chi=math.acos(cz); d=sel(chi); p1=random.choice((0.5,-0.5)); p2=random.choice((0.5,-0.5))
    s12=p1*p2+0.5*math.cos(d); key=(('opposite' if p1!=p2 else 'same'),round(d,3),round(s12,3)); cnt[key]=cnt.get(key,0)+1
print('(2) implied states (projections random +-1/2 on the common axis):')
for kk,v in sorted(cnt.items()): print(f'    projections {kk[0]:8s} dphi {kk[1]:.3f}  S1.S2 = {kk[2]:+.3f}  fraction {v/N:.4f}')
print('    QM: singlet -3/4 (0.25), T0 +1/4 (0.25), T+-1 +1/4 (0.50)')
# (3) time in one pass: relativistic hyperbola approximated by energy conservation along a Kepler hyperbola (nonrel geometry, speed from energy)
print('(3) one pass, E = 10 keV relative:')
mu=me/2; E=10e3*e
for rmin_fm in (190,100,50,20,10):
    rmin=rmin_fm*1e-15
    # nonrel hyperbola: a = kC/(2E), r_min = a(ecc - 1)
    a=kC/(2*E); ecc=1+rmin/a; L=math.sqrt(mu*kC*a*(ecc*ecc-1))
    # integrate over true anomaly while r < r*: dt = mu r^2 dtheta / L ; speed from energy (relativistic gamma for each particle)
    th_max=math.acos((a*(ecc*ecc-1)/rstar-1)/ecc) if rstar>rmin else 0.0
    n=20000; phase=0.0; t=0.0
    for i in range(n):
        th=-th_max+2*th_max*(i+0.5)/n; r=a*(ecc*ecc-1)/(1+ecc*math.cos(th)); dt=mu*r*r*(2*th_max/n)/L
        U=kC/(2*r)        # near-zone amplitude ~ half Coulomb (474 identity), upper bound of the coupling
        phase+=U/hbar*dt; t+=dt
    print(f'    r_min {rmin_fm:4d} fm: time inside r* {t:.3e} s = {t*wZB/(2*math.pi):.2f} ZB periods; accumulated coupling phase <= {phase:.3f} rad')

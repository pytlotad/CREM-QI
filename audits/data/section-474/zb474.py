# Audit 474: Zitterbewegung-frequency echo between two transverse spin moments.
import cmath, math
hbar=1.054571817e-34; c=299792458.0; me=9.1093837015e-31; mu0=1.25663706212e-6; muB=9.2740100783e-24; alpha=7.2973525693e-3
h=2*math.pi*hbar; eV=1.602176634e-19
wZB=2*me*c*c/hbar; k=wZB/c; lbarC=hbar/(me*c)
a0=lbarC/alpha; aPs=2*a0
print(f'omega_ZB {wZB:.4e} rad/s, k {k:.4e} 1/m (= 2/lambdabar_C), k a0 {k*a0:.2f}, k a_Ps {k*aPs:.2f}')
mperp=muB*math.sqrt(2)   # transverse moment for length sqrt(3)/2 hbar, projection 1/2: t = sqrt(1/2) hbar -> |m_perp| = 2 muB t/hbar
def G(r):   # magnetic dipole Green dyadic coefficients: B = (mu0/4pi) e^{ikr} [ (3nn - I)(1/r^3 - ik/r^2) + (I - nn) k^2/r ] . m
    e=cmath.exp(1j*k*r); return (mu0/(4*math.pi))*e*(1/r**3-1j*k/r**2), (mu0/(4*math.pi))*e*k*k/r
# transverse moments rotating in the plane normal to the common axis z; separation n in that plane (worst case for coupling)
# rotating vectors m_i = m_perp Re[(x + i y) e^{-i(w t + phi_i)}]; time-averaged energy U = -1/2 Re[m1* . B2] with B2 = G . m2
def U(r,dphi):
    A,Bf=G(r); # n along x: (3nn - I) = diag(2,-1,-1); (I - nn) = diag(0,1,1)
    # complex amplitudes: m1 = m_perp (1, i, 0)/sqrt2 ... use circular: m = m_perp/ sqrt(2) (x + i y)
    v1=[1,1j,0]; v2=[cmath.exp(-1j*dphi),1j*cmath.exp(-1j*dphi),0]
    T=[[2*A,0,0],[0,-A+Bf,0],[0,0,-A+Bf]]
    Bv=[sum(T[i][j]*v2[j] for j in range(3)) for i in range(3)]
    return -0.5*(mperp**2/2)*sum((v1[i].conjugate()*Bv[i]).real for i in range(3))
def scan(r):
    us=[(U(r,2*math.pi*j/3600),2*math.pi*j/3600) for j in range(3600)]
    umin=min(us); umax=max(us); return umin[1], (umax[0]-umin[0])/2
HFS_Ps=203.3942e9*h; HFS_H=1420.405751768e6*h
print('\ncircular separations (r = a): selected relative phase and implied cos th12 (opposite projections: S1.S2 = -1/4 + 1/2 cos dphi)')
for name,r in (('Ps 1s r = a_Ps',aPs),('Ps 2p r = 4 a_Ps',4*aPs),('H 1s r = a0',a0)):
    ph,amp=scan(r); s12=-0.25+0.5*math.cos(ph); c12=s12/0.75
    print(f'  {name:18s} k r = {k*r:9.3f} (mod 2pi {math.fmod(k*r,2*math.pi):.3f}); dphi* = {ph:.4f} rad; cos th12 = {c12:+.4f} (P2: singlet -1, T0 +1/3);'
          f' half-amplitude of U(dphi) = {amp/eV:.3e} eV = {amp/HFS_Ps:.2e} x HFS(Ps), {amp/HFS_H:.2e} x HFS(H)')
# sensitivity: how much r must change to move dphi* by 2 pi
print(f'\nr shift for a full turn of dphi*: 2 pi / k = {2*math.pi/k:.4e} m = {2*math.pi/k/a0:.4e} a0')
# orbit averaging of the far-field phase factor over a Kepler orbit, weight 1/r (far field) and 1/r^3
def orbit_avg(a,e,N=200000):
    s1=0; s2=0; w1=0; s3=0; w3=0
    for i in range(N):
        E=2*math.pi*(i+0.5)/N; r=a*(1-e*math.cos(E)); wt=(1-e*math.cos(E))
        s1+=wt*cmath.exp(1j*k*r)/r; w1+=wt/r; s3+=wt*cmath.exp(1j*k*r)/r**3; w3+=wt/r**3
    return abs(s1)/w1, abs(s3)/w3
print('\norbit averages |<e^{ikr} r^-p>| / <r^-p> (time average over a Kepler orbit):')
for name,a,e in (('Ps 1s (L = hbar/2, e = 0.968)',aPs,math.sqrt(1-0.25)), ('Ps 2p (L = 1.5 hbar, e = 0.661)',4*aPs,math.sqrt(1-(1.5/2)**2)), ('H 1s (L = hbar/2, e = 0.968)',a0,math.sqrt(1-0.25)), ('Ps circular n = 2 (e = 0)',4*aPs,0.0)):
    p1,p3=orbit_avg(a,e); print(f'  {name:32s} p = 1: {p1:.3e}   p = 3: {p3:.3e}')
# radiation of the transverse moment at omega_ZB
P=mu0*wZB**4*mperp**2/(6*math.pi*c**3)
print(f'\nradiated power of m_perp = sqrt(2) mu_B rotating at omega_ZB: {P:.3e} W; m c^2 / P = {me*c*c/P:.3e} s; per period of omega_ZB the energy lost / (hbar omega_ZB) = {P*2*math.pi/wZB/(hbar*wZB):.3e}')

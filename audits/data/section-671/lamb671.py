# Audit 671: can the model's classical 1S orbit get a Lamb shift from radiation reaction or from isotropic field fluctuations?  usage: python3 lamb671.py
import numpy as np, math
erfc=np.vectorize(math.erfc)
hb=1.054571817e-34; me=9.1093837015e-31; mp=1.67262192369e-27; e0=1.602176634e-19; eps=8.8541878128e-12; c=299792458.0; h=2*math.pi*hb
al=7.2973525693e-3; k=e0**2/(4*math.pi*eps); mu=me*mp/(me+mp)
MHz=lambda E: E/h/1e6; ueV=lambda E: E/e0*1e6
LAMB=8172.9
E=-2.178685816657088e-18; L=0.5*hb                       # engine 1S (probe 662)
a=k/(2*-E); ecc=math.sqrt(1+2*E*L**2/(mu*k**2)); rp=a*(1-ecc); ra=a*(1+ecc)
vp=math.sqrt(k/mu*(2/rp-1/a))
print(f'1S orbit: a {a*1e12:.4f} pm, e {ecc:.6f}, r_p {rp*1e12:.3f} pm, r_a {ra*1e12:.3f} pm, v_p {vp:.4e} m/s')
# (a) radiation reaction
tau0=2*k/(3*me*c**3); wp=vp/rp
dEa=(tau0*wp)**2*abs(E)
print(f'(a) tau0 {tau0:.4e} s, omega_p {wp:.4e} 1/s, tau0 omega_p {tau0*wp:.3e}; bound (tau0 omega_p)^2 |E| = {MHz(dEa):.4f} MHz ({ueV(dEa):.3e} ueV)'
      f' -> H671a {MHz(dEa)<=1.0}  (Lamb / bound = {LAMB/MHz(dEa):.2e})')
# (b) Welton smearing on the classical orbit, time average over one Kepler period
lc=hb/(me*c); Ry=13.605693122994*e0
dr2=(2*al/math.pi)*lc**2*math.log(me*c**2/Ry); sig=math.sqrt(dr2/3)
M=np.linspace(0,2*math.pi,200001)[:-1]; u=M.copy()
for _ in range(100): u-=(u-ecc*np.sin(u)-M)/(1-ecc*np.cos(u))
r=a*(1-ecc*np.cos(u))
x=r/(math.sqrt(2)*sig); dV=k*erfc(x)/r
dEb=float(np.mean(dV))
print(f'(b) Welton <dr^2>^1/2 {math.sqrt(dr2)*1e15:.2f} fm, sigma {sig*1e15:.2f} fm; closest approach r_p / (sqrt2 sigma) = {rp/(math.sqrt(2)*sig):.1f};'
      f' orbit-averaged shift {MHz(dEb):.3e} MHz -> H671b {MHz(dEb)<=1.0}')
for f in (10,100,1000):
    s=sig*f; dEf=float(np.mean(k*erfc(r/(math.sqrt(2)*s))/r))
    print(f'    with sigma x {f}: {s*1e15:.0f} fm, shift {MHz(dEf):.3e} MHz')
# (c) the same mechanism with the QM density at the nucleus (Bethe)
redm=mu/me
bethe=(4/(3*math.pi))*al*al**4*me*c**2*(math.log(al**-2)-2.984128+19/30)*redm**3
print(f'(c) Bethe 1S with QM |psi(0)|^2 (reduced-mass factor (mu/m)^3 {redm**3:.6f}): {MHz(bethe):.1f} MHz ({ueV(bethe):.3f} ueV); measured {LAMB} MHz:'
      f' {100*(MHz(bethe)/LAMB-1):+.2f} % -> H671c {abs(MHz(bethe)/LAMB-1)<=0.02}')
print(f'    Welton with QM density: (2/3) <dr^2> k 4 pi |psi(0)|^2 ... = {MHz((2*math.pi/3)*dr2*k*(1/(math.pi*(a)**3))):.1f} MHz (same order; log cut-off conventional)')

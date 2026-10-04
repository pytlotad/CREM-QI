# Audit 324: ground-state hyperfine splitting, model vs measurement.
import numpy as np
mu0=1.25663706212e-6; mu=9.2847647043e-24; e=1.602176634e-19; c=299792458.0
h=6.62607015e-34; a0=5.29177210903e-11; alpha=7.2973525693e-3; mec2=8.1871057769e-14
rstar=mu/(e*c); eps=0.96682*rstar; a=2*a0
MEAS=203.3942e9; LO=(7/12)*alpha**4*mec2/h
print(f"r* = {rstar:.6e} m, eps = {eps:.6e} m, a_pair = {a:.6e} m, a/eps = {a/eps:.1f}")
print(f"measured 203.3942 GHz; LO QED 7/12 alpha^4 mc^2 = {LO/1e9:.4f} GHz (theory)")
def n(r): return 3*eps**2/(4*np.pi*(r*r+eps*eps)**2.5)
def B_field(rvec,m):   # code: plummerDipoleField + plummerMagnetizationField
    rho2=(rvec**2).sum(-1)+eps**2; ir=rho2**-0.5
    pole=(rvec*(3*(rvec@m)*ir**5)[:,None]-m*(ir**3)[:,None])*mu0/(4*np.pi)
    mag=m*(mu0*3*eps**2/(4*np.pi*rho2**2.5))[:,None]
    return pole+mag
# R324a: spherical average identity
rng=np.random.default_rng(1); dirs=rng.normal(size=(400000,3)); dirs/=np.linalg.norm(dirs,axis=1)[:,None]
m=np.array([0,0,1.0]); worst=0
for r in (0.1*eps,eps,3*eps,30*eps):
    Bav=B_field(dirs*r,m).mean(0)[2]; want=(2*mu0/3)*n(r); worst=max(worst,abs(Bav/want-1))
print(f"R324a: spherical average / (2mu0/3) m n(r): worst rel. dev {worst:.1e} (Monte Carlo, 4e5 directions)")
def dE(nbar): return 2*(2*mu0/3)*mu**2*nbar     # E_ortho - E_para, J
def report(name,nbar):
    d=dE(nbar)/h; print(f"  {name:38s} <n> = {nbar:.4e} m^-3  dE = {d/1e9:12.5f} GHz = {d/MEAS:10.4e} x pomiar")
    return d
# QM 1s density (import)
r=np.geomspace(1e-6*eps,40*a,2_000_001); psi2=np.exp(-2*r/a)/(np.pi*a**3)
nq=np.trapezoid(psi2*n(r)*4*np.pi*r*r,r)
print(f"|psi(0)|^2 = {1/(np.pi*a**3):.4e}; smeared <n>_1s / |psi(0)|^2 = {nq*np.pi*a**3:.6f}")
res={}
res['qm']=report("gestosc kwantowa 1s (import)",nq)
def orbit_nbar(J0):
    ecc=np.sqrt(max(0.0,1-J0*J0))
    # eccentric anomaly grid clustered at E=0 (periapsis)
    u=np.linspace(0,1,2_000_001); E=np.pi*u**4
    w=1-ecc*np.cos(E); f=n(a*w)*w
    return np.trapezoid(f,E)/np.pi       # (1/2pi) * 2 * int_0^pi
for J0,name in ((1.0,"kolo L = hbar (Bohr, domyslny start)"),(0.5,"L = hbar/2 (Langer, l = 0)"),
                (0.07,"L = 0,07 hbar (eksperyment 6)"),(0.0,"radialna L = 0")):
    res[J0]=report(name,orbit_nbar(J0))
# microcanonical: e^2 uniform on [0,1]  <=>  J0^2 = 1 - e^2 uniform
J=np.sqrt(np.linspace(0,1,4001)); vals=np.array([orbit_nbar(j) for j in J[::40]])
nmc=np.trapezoid(vals,J[::40]**2)
res['mc']=report("mikrokanoniczny (e^2 jednorodne)",nmc)
print(f"  mikrokanoniczny <n>/|psi(0)|^2 = {nmc*np.pi*a**3:.2f}  (a/eps)^0.5 = {(a/eps)**0.5:.2f}")

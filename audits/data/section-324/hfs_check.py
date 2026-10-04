exec(open(__file__.replace('hfs_check.py','hfs.py')).read().split('# R324a')[0])
# R324a exact: azimuthal average kills x,y; B_z depends on cos(theta) only -> Gauss-Legendre
x,w=np.polynomial.legendre.leggauss(64); worst=0
for r in (0.01*eps,0.1*eps,eps,3*eps,30*eps,300*eps):
    rho2=r*r+eps*eps
    Bz=mu0/(4*np.pi)*(3*r*r*x*x/rho2**2.5-1/rho2**1.5)+mu0*3*eps**2/(4*np.pi*rho2**2.5)
    avg=(w*Bz).sum()/2; want=(2*mu0/3)*n(r); worst=max(worst,abs(avg/want-1))
print(f"R324a (Gauss-Legendre, 6 radii 0.01-300 eps): worst rel. dev {worst:.1e}")
def dE(nbar): return 2*(2*mu0/3)*mu**2*nbar
def orbit_nbar(J0,N=2_000_001):
    ecc=np.sqrt(max(0.0,1-J0*J0)); u=np.linspace(0,1,N); E=np.pi*u**4
    wgt=1-ecc*np.cos(E); return np.trapezoid(n(a*wgt)*wgt,E)/np.pi
print("radial convergence:", orbit_nbar(0.0,500_001), orbit_nbar(0.0,2_000_001), orbit_nbar(0.0,8_000_001))
for M in (200,400,800):
    q=np.concatenate(([0.0],np.geomspace(1e-7,1.0,M)))      # q = J0^2
    vals=np.array([orbit_nbar(np.sqrt(v),400_001) for v in q])
    nmc=np.trapezoid(vals,q)
    print(f"microcanonical, {M} points in J0^2 (geometric from 1e-7): <n>/|psi0|^2 = {nmc*np.pi*a**3:.4f}, dE = {dE(nmc)/h/1e9:.3f} GHz = {dE(nmc)/h/MEAS:.4f} x pomiar")

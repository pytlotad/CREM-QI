exec(open(__file__.replace('hfs_micro_direct.py','hfs.py')).read().split('# R324a')[0])
r=np.geomspace(1e-9*eps,2*a*(1-1e-15),4_000_001)
rho=np.sqrt(np.clip(1/r-1/(2*a),0,None)); norm=np.trapezoid(rho*4*np.pi*r*r,r)
nbar=np.trapezoid(rho*n(r)*4*np.pi*r*r,r)/norm
d=2*(2*mu0/3)*mu**2*nbar/h
print(f"direct microcanonical density: <n>/|psi0|^2 = {nbar*np.pi*a**3:.4f}, dE = {d/1e9:.3f} GHz = {d/MEAS:.4f} x pomiar")

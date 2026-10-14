# Audit 441: two clocks on a FROZEN Kepler orbit (hbar = mu = k = 1).  theta_r = psi_r + int p_r rdot dt = psi_r + J(E + e sinE) - L f,
# theta_phi = psi_phi + L f (f true anomaly).  Couplings U = -eps g(r) cos theta_r - eps h(f) cos theta_phi, g = exp(-r), h = exp(cos f).
# Time averages over one orbit at resonance (J_r = k integer, L = m integer):
#   radial drive  dJ_r/dt|g   = <cos theta_r dg/dt> / omega  ->  -(sin psi_r) Ir      (radial force changes E only; dJ_r = dE/omega)
#   angular drive dL/dt|h     = <h'(f) cos theta_phi>       ->  -(sin psi_phi) Iphi   (torque)
#   torque also changes E: dE/dt|h = <h'(f) fdot cos theta_phi>  ->  dJ_r/dt|h = dE/omega - dL   (cross term, reported)
import math
def simpson(fun,a,b,n=8000):
    h=(b-a)/n; s=fun(a)+fun(b)
    for i in range(1,n): s+=(4 if i%2 else 2)*fun(a+i*h)
    return s*h/3
for (k,m) in [(1,1),(1,2),(2,1),(1,3),(2,2)]:
    J=k+m; L=float(m); e=math.sqrt(1-(L/J)**2); a=J*J; w=1/J**3; T=2*math.pi/w
    def fofE(E): return 2*math.atan2(math.sqrt(1+e)*math.sin(E/2), math.sqrt(1-e)*math.cos(E/2))
    dtdE=lambda E:(1-e*math.cos(E))*J**3
    r=lambda E: a*(1-e*math.cos(E))
    th_r=lambda E: J*(E+e*math.sin(E))-L*fofE(E)
    gdot_dE=lambda E: math.exp(-r(E))*(-1)*a*e*math.sin(E)          # dg/dE
    Ir=simpson(lambda E: gdot_dE(E)*math.sin(th_r(E)),-math.pi,math.pi)/T/w     # radial drive coefficient (dJ_r/dt per eps)
    hp=lambda f: -math.sin(f)*math.exp(math.cos(f))
    Iphi=simpson(lambda E: hp(fofE(E))*math.sin(m*fofE(E))*dtdE(E),-math.pi,math.pi)/T   # <h' sin(m f)>_t
    fdot=lambda E: L/r(E)**2
    Ecross=simpson(lambda E: hp(fofE(E))*fdot(E)*math.sin(m*fofE(E))*dtdE(E),-math.pi,math.pi)/T   # <h' fdot sin mf>_t
    print(f"J_r={k} L={m} (J={J}, e={e:.4f}): radial drive coeff {Ir:+.4e}; angular drive coeff <h' sin mf> {Iphi:+.4e}; torque power coeff {Ecross:+.4e} -> dJ_r/dt cross {Ecross/w - Iphi:+.4e}")

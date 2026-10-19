# Audit 537 (side question): is B ~ 1 T near a 2:1 resonance of the classical two-spin model (relative Larmor 2 g mu_B B/h vs magnetic contact 4K_c/h)?
import numpy as np, sys
sys.path.insert(0,'.')
h=6.62607015e-34; hbar=h/(2*np.pi); muB=9.2740100783e-24; g=2.00231930436092
dE=203.38910e9*h; Kc=204.6577e9*h/1.748264/4      # engine's magnetic contact (4K_c = HFS/f, audit 537)
def run(K,B,T=1e-9,dt=2.5e-14):
    Z=(g/2)*muB*B; n1=np.array([1.0,0,0]); n2=n1.copy(); zz=np.array([0,0,Z]); wsum=0.0; N=int(T/dt)
    def f(n1,n2): return np.cross((2/hbar)*(K*n2+zz),n1),np.cross((2/hbar)*(K*n1-zz),n2)
    for i in range(N):
        a1,a2=f(n1,n2); b1,b2=f(n1+0.5*dt*a1,n2+0.5*dt*a2); c1,c2=f(n1+0.5*dt*b1,n2+0.5*dt*b2); d1,d2=f(n1+dt*c1,n2+dt*c2)
        n1=n1+dt/6*(a1+2*b1+2*c1+d1); n2=n2+dt/6*(a2+2*b2+2*c2+d2); n1/=np.linalg.norm(n1); n2/=np.linalg.norm(n2); wsum+=(1-n1@n2)/2
    return wsum/N
print("4K_c/h = %.3f GHz; relative Larmor 2 g mu_B B/h per tesla = %.3f GHz"%(4*Kc/h*1e-9,2*g*muB/h*1e-9))
for B in (0.6,0.8,0.9,0.95,1.0,1.05,1.1,1.2,1.5):
    x=2*g*muB*B/dE; y=x/(1+np.sqrt(1+x*x)); PS=y*y/(1+y*y); w=run(Kc,B); wr=run(dE/4,B)
    print(f"B {B:.2f} T: ratio (2 g mu_B B)/(4 K_c) = {2*g*muB*B/(4*Kc):.3f}; magnetic K_c: <w> {w:.4e} ({w/PS:.2f} P_S); K = dE/4: <w> {wr:.4e} ({wr/PS:.2f} P_S)")

# Audit 558: QM {S, T0} with decay from the Hermitian S' -- singlet weight of the surviving (normalized) state over 0-2 ns, against the engine's w range (p558.txt).
# Model numbers as p558: Gamma_p 8.03250e9/s, Gamma_o 7.21117e6/s, HFS 204.6577 GHz, Z = (g/2) mu_B B.
import numpy as np
h=6.62607015e-34; hbar=h/(2*np.pi); muB=9.2740100783e-24; g=2.00231930436; Gp=8.03250e9; Go=7.21117e6; K=204.6577e9*h/4
for B in (0.1,1.0,3.0):
    Z=(g/2)*muB*B
    H=np.array([[-3*K-0.5j*hbar*Gp,2*Z],[2*Z,K-0.5j*hbar*Go]])
    th=0.5*np.arctan2(4*Z,4*K); psi0=np.array([np.cos(th),-np.sin(th)],dtype=complex)
    w,V=np.linalg.eig(H); c=np.linalg.solve(V,psi0)
    ws=[]
    for t in np.linspace(0,2e-9,4001):
        psi=V@(c*np.exp(-1j*w*t/hbar)); ws.append(abs(psi[0])**2/np.vdot(psi,psi).real)
    ws=np.array(ws)
    print(f"B {B} T: QM conditioned singlet weight from S': start {ws[0]:.6f}, range {ws.min():.6f}-{ws.max():.6f}, at 2 ns {ws[-1]:.6f}")

# Czy kontakt stanow S (import 17) to Sommerfeld k = j + 1/2 = 1 takze dla Ps?
# Klasyczna srednia hamiltonianu Darwina (dwa ciala) po orbicie Keplera przy (J = n hbar, L):
#   H2 = -(p^4/8c^2)(1/m1^3 + 1/m2^3) - (k/(2 m1 m2 c^2 r)) (p^2 + p_r^2)
# Jednostki: hbar = c = m_e = 1, k = alpha. Szukamy L_eff: E_rel(L_eff) = E_rel(L=1/2) + C_n (import 17).
import numpy as np
alpha=7.2973525693e-3; k=alpha
def Erel(L,n,m1,m2):
    mu=m1*m2/(m1+m2); J=n
    E=-mu*k*k/(2*J*J); a=J*J/(mu*k); b=a*L/J
    r1=1/a; r2=1/(a*b); r3=1/b**3
    p4=4*mu*mu*(E*E+2*E*k*r1+k*k*r2)
    p2r=2*mu*(E*r1+k*r2); pr2r=p2r-L*L*r3
    return -(p4/8)*(1/m1**3+1/m2**3)-(k/(2*m1*m2))*(p2r+pr2r)
def C(n,m1,m2):  # import 17: pi k hbar^2 n_QR (1/2m1^2+1/2m2^2+1/(m1 m2)), n_QR = |psi_ns(0)|^2 = 1/(pi a^3 n^3)
    mu=m1*m2/(m1+m2); a=1/(mu*k)
    return np.pi*k/(np.pi*a**3*n**3)*(1/(2*m1**2)+1/(2*m2**2)+(1/(m1*m2) if m2<1e6 else 0))
def solve(n,m1,m2):
    tgt=Erel(0.5,n,m1,m2)+C(n,m1,m2); lo,hi=0.3,5.0
    for _ in range(200):
        mid=(lo+hi)/2
        if (Erel(mid,n,m1,m2)-tgt)*(Erel(lo,n,m1,m2)-tgt)<=0: hi=mid
        else: lo=mid
    return (lo+hi)/2
for lab,m2 in [("wodor (m2 -> inf)",1e12),("muonium",206.768),("Ps",1.0)]:
    print(lab, ["n=%d: L_eff=%.6f"%(n,solve(n,1.0,m2)) for n in (1,2,3,5)])

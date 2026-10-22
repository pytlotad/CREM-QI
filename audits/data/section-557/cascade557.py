# Audit 557: exact hydrogen cascade (Case A, no field): Einstein A_{nl -> n'l'} from bound-bound radial integrals (atomic units, reduced mass ~ m_e), capture shares (n, l) exact at
# K/B = 0.05 (gaunt553 functions), cascade probabilities of Lyman-alpha, 2s, direct Lyman-beta+, Halpha, Hbeta per capture, n <= 8.
import mpmath as mp, numpy as np
src=open('../section-553/gaunt553.py').read().split('print("check')[0]; ns={}; exec(src,ns)
mp.mp.dps=20
alpha=float(ns['alpha'])
mp.mp.dps=60
def coeffs(n,l):   # R_nl = N sum_j a_j r^(l+j) e^(-r/n)
    N=mp.sqrt((mp.mpf(2)/n)**3*mp.factorial(n-l-1)/(2*n*mp.factorial(n+l)))
    return [N*(-1)**j*mp.binomial(n+l,n-l-1-j)/mp.factorial(j)*(mp.mpf(2)/n)**(l+j) for j in range(n-l)]
NMAX=8
rad={}
def R(n,l,m,lp):   # exact: int R_nl r R_m,lp r^2 dr
    key=(n,l,m,lp)
    if key not in rad:
        beta=mp.mpf(1)/n+mp.mpf(1)/m; s=0
        for j,a in enumerate(coeffs(n,l)):
            for jj,b in enumerate(coeffs(m,lp)):
                p=l+j+lp+jj+3; s+=a*b*mp.factorial(p)/beta**(p+1)
        rad[key]=float(s)
    return rad[key]
# A (atomic units, 1/t_au): A = (4/3) alpha^3 omega^3 (l_max/(2l+1)) R^2
A={}
for n in range(2,NMAX+1):
    for l in range(n):
        for m in range(1,n):
            for lp in (l-1,l+1):
                if lp<0 or lp>m-1: continue
                w=0.5*(1/m**2-1/n**2)
                A[(n,l,m,lp)]=(4/3)*alpha**3*w**3*max(l,lp)/(2*l+1)*R(n,l,m,lp)**2
tot={(n,l):sum(v for k,v in A.items() if k[:2]==(n,l)) for n in range(2,NMAX+1) for l in range(n)}
tau_au=2.4188843265857e-17
print("lifetimes [ns]: "+' '.join(f"{n}{'spdfghi'[l]} {1/tot[(n,l)]*tau_au*1e9:.3f}" for (n,l) in [(2,1),(3,0),(3,1),(3,2),(4,1)]))
# capture shares at K/B = 0.05, n <= NMAX
k=mp.sqrt(mp.mpf('0.05')); cap={}
for n in range(1,NMAX+1):
    omega=k*k/2+mp.mpf(1)/(2*n*n)
    for l in range(n):
        cap[(n,l)]=float(2*(2*l+1)*(omega*ns['alpha']/k)**2*ns['sigma_PI'](n,l,k))
Z=sum(cap.values())
# propagate: probability flux from top
flux={key:cap[key]/Z for key in cap}; lines={}
for n in range(NMAX,1,-1):
    for l in range(n):
        f=flux[(n,l)]
        if f==0 or (n,l)==(2,0): continue
        for (a,b,m,lp),v in A.items():
            if (a,b)!=(n,l): continue
            p=f*v/tot[(n,l)]; flux[(m,lp)]=flux.get((m,lp),0)+p; lines[(n,m)]=lines.get((n,m),0)+p
end1s=flux[(1,0)]; end2s=flux[(2,0)]
lya=lines.get((2,1),0); direct=sum(v for (a,m),v in lines.items() if m==1 and a>2)
print(f"capture shares (n <= {NMAX}) normalized; shares n = 1..4: "+' '.join(f"{sum(v for (n,l),v in cap.items() if n==m)/Z:.3f}" for m in range(1,5)))
print(f"QM cascade per capture: captured directly into 1s {cap[(1,0)]/Z:.3f}; Lyman-alpha {lya:.3f}; ends in 2s {end2s:.3f} (incl. direct capture {cap[(2,0)]/Z:.3f}); Lyman-beta+ (n >= 3 -> 1s) {direct:.3f}")
print(f"  Balmer per capture: Halpha (3->2) {lines.get((3,2),0):.4f}, Hbeta (4->2) {lines.get((4,2),0):.4f}, Halpha/Hbeta {lines.get((3,2),0)/lines.get((4,2),1e-300):.3f}; Paschen-alpha (4->3) {lines.get((4,3),0):.4f}")

# Audit 464: De Luca (PRE 73, 026221, 2006) resonance -- PP modes of eqs (42), (90), condition eps1 - eps2 + p theta = 0.
import cmath, math, sys
def fxy(l,th,r):   # eq (42): planar, r = mu/M
    return (1+2/l+7/l**2+10/l**3-5/l**4)*r*th**4*cmath.exp(-2*l)-(1-2/3*th**2*l+th**4*l**2/9)
def fz(l,th,r):    # eq (90): perpendicular
    return 1-2/3*th**2*l+4*r/9*th**4*l**2-r*th**4*(1+cmath.exp(-2*l))*(1+2/l-1/l**2-1/l**3+1/l**4)
def root(f,th,r,q):
    sig=0.5*math.log(1/(r*th**4)); l=complex(-sig,math.pi*q)
    for _ in range(200):
        h=1e-7*abs(l); d=(f(l+h,th,r)-f(l-h,th,r))/(2*h); s=f(l,th,r)/d; l-=s
        if abs(s)<1e-14*abs(l): return l
    return None
def g(th,r,q,p):
    a=root(fxy,th,r,q); b=root(fz,th,r,q)
    if a is None or b is None: return None
    return (a.imag-math.pi*q)-(b.imag-math.pi*q)-p*th   # printed (56) has +; the sign that reproduces Table 1
def solve(r,q,p):
    # scan l_z = 1/theta from 20 to 6000 (log grid), bisect sign changes
    out=[]; N=3000; prev=None
    for i in range(N+1):
        lz=20*(300**(i/N)); th=1/lz; v=g(th,r,q,p)
        if v is not None and prev is not None and prev[1] is not None and (v>0)!=(prev[1]>0):
            a,b=prev[0],th; fa=prev[1]
            for _ in range(80):
                m=0.5*(a+b); fm=g(m,r,q,p)
                if fm is None: break
                if (fm>0)==(fa>0): a,fa=m,fm
                else: b=m
            out.append(1/(0.5*(a+b)))
        prev=(th,v)
    return out
if __name__=='__main__':
    m2=float(sys.argv[1]); p=float(sys.argv[2]); qmax=int(sys.argv[3]) if len(sys.argv)>3 else 21
    mu=m2/(1+m2); M=1+m2; r=mu/M; H=137.036
    print(f'# m1 = 1, m2 = {m2}, mu/M = {r:.6e}, p = {p}')
    print('# q  l_z roots (e^2/c)  ->  per root: l_z/137.036, frequency mu(137.036 theta)^3 [a.u.]')
    for q in range(1,qmax+1):
        rs=solve(r,q,p)
        print(q,'  '.join(f'{lz:.2f} (k~{lz/H:.3f}, w {mu*(H/lz)**3:.4e})' for lz in rs) if rs else 'none')

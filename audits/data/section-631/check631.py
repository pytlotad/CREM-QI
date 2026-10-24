"""Audit 631, H631a: p631 output against exact 6j (Fraction, Racah) -- weights, and chi^2 of the drawn counts.  usage: python3 check631.py p631.txt"""
import math, sys
from fractions import Fraction as F
def fac(n): return math.factorial(n)
def tri(a,b,c): return abs(a-b)<=c<=a+b
def d2(a,b,c): return F(fac(a+b-c)*fac(a-b+c)*fac(-a+b+c),fac(a+b+c+1))
def sixj2(a,b,c,d,e,f):
    if not(tri(a,b,c) and tri(a,e,f) and tri(d,b,f) and tri(d,e,c)): return F(0)
    pre=d2(a,b,c)*d2(a,e,f)*d2(d,b,f)*d2(d,e,c); s=0
    for t in range(max(a+b+c,a+e+f,d+b+f,d+e+c),min(a+b+d+e,b+c+e+f,c+a+f+d)+1):
        s+=(-1)**t*F(fac(t+1),fac(t-a-b-c)*fac(t-a-e-f)*fac(t-d-b-f)*fac(t-d-e-c)*fac(a+b+d+e-t)*fac(b+c+e+f-t)*fac(c+a+f+d-t))
    return pre*s*s
def qm(l,J,lp):
    w=[(2*Jp+1)*sixj2(l,J,1,Jp,lp,1) for Jp in range(abs(lp-1),lp+2)]; z=sum(w); return [float(x/z) for x in w]
def gammaq(a,x):
    if x<=0: return 1.0
    if x<a+1:
        s=t=1/a; ap=a
        for _ in range(2000):
            ap+=1; t*=x/ap; s+=t
            if abs(t)<abs(s)*1e-15: break
        return 1-s*math.exp(-x+a*math.log(x)-math.lgamma(a))
    b=x+1-a; c=1e300; d=1/b; h=d
    for i in range(1,2000):
        an=-i*(i-a); b+=2; d=an*d+b; d=d if abs(d)>1e-300 else 1e-300; c=b+an/c; c=c if abs(c)>1e-300 else 1e-300; d=1/d; h*=d*c
        if abs(d*c-1)<1e-15: break
    return math.exp(-x+a*math.log(x)-math.lgamma(a))*h
worst=0; chi=0; dof=0
for line in open(sys.argv[1]):
    w=line.split()
    if w[0]=='W':
        l,J,lp=int(w[2]),int(w[4]),int(w[6]); eng=[float(x) for x in w[7:]]; ref=qm(l,J,lp)
        worst=max(worst,max(abs(a-b) for a,b in zip(eng,ref)))
    elif w[0]=='C':
        l,up,J,n=int(w[2]),int(w[4]),int(w[6]),int(w[8]); obs=[int(x) for x in w[9:12]]; lp=l+(1 if up else -1)
        ref=qm(l,J,lp)
        for o,p in zip(obs,ref):
            if p*n>=5: chi+=(o-p*n)**2/(p*n); dof+=1
            elif p==0 and o>0: print("zdarzenie w galezi o wadze 0:",line.strip())
        dof-=1
    elif w[0]=='TURN': print(line.strip())
print("H631a: max |w_silnik - w_QM| = %.3e; chi^2 losowan = %.1f przy %d stopniach swobody, p = %.3f" % (worst,chi,dof,gammaq(dof/2,chi/2)))

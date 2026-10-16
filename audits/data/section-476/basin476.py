# Audit 476: rule Q on the precession loop of a spin pair about S: basin weights and per-pair channel vs the Born rule.
import math, random
random.seed(476)
def iso():
    z=2*random.random()-1; ph=2*math.pi*random.random(); r=math.sqrt(1-z*z); return (r*math.cos(ph),r*math.sin(ph),z)
def vmf(P,axis):  # von Mises-Fisher direction with mean cosine P about axis (+z or -z)
    if P<=1e-9: return iso()
    # solve coth(k) - 1/k = P
    lo,hi=1e-9,1e6
    for _ in range(200):
        k=0.5*(lo+hi); m=(1/math.tanh(k)-1/k) if k<700 else 1-1/k
        if m<P: lo=k
        else: hi=k
    u=random.random(); w=1+math.log(u+(1-u)*math.exp(-2*k))/k
    ph=2*math.pi*random.random(); r=math.sqrt(max(0,1-w*w)); v=(r*math.cos(ph),r*math.sin(ph),w)
    return v if axis>0 else (v[0],v[1],-v[2])
def dot(a,b): return sum(x*y for x,y in zip(a,b))
N=1000000
print('(1) isotropic ensemble, basin of rule Q on |S| (S = 0 if |S| < boundary):')
for name,l,bnd in (('hbar/2, Langer rungs 1/2, 3/2',0.5,1.0),('sqrt(3)/2, Langer rungs',math.sqrt(0.75),1.0),('sqrt(3)/2, sqrt rungs 0, sqrt2',math.sqrt(0.75),math.sqrt(2)/2),('hbar (rule Q), Langer rungs',1.0,1.0)):
    s=0
    for _ in range(N):
        c=dot(iso(),iso()); S=math.sqrt(max(0,2*l*l*(1+c)))
        s+=S<bnd
    ps=s/N; print(f'  length {name:32s}: P(S=0) = {ps:.5f} (+-{math.sqrt(ps*(1-ps)/N):.5f}); singlet : triplet = 1 : {(1-ps)/ps if ps>0 else float("inf"):.3f}')
print('\n(2) single pair, length hbar: P_S basin = step(cos th12 < -1/2); Born = (1 - cos)/4')
for c in (-1,-0.75,-0.5001,-0.4999,-0.25,0,0.5,1): print(f'  cos {c:+.4f}: basin {1.0 if c<-0.5 else 0.0:.2f}  Born {(1-c)/4:.4f}')
print('\n(3) partially polarized ensembles (vMF, mean cos P), basin (length hbar) vs Born <(1 - s1.s2)/4>:')
M=200000
for Pp,Pe in ((0.3,0.03),(0.3,0.1),(0.3,0.3),(0.6,0.6),(0.9,0.9),(1.0,1.0)):
    for ax in (+1,-1):
        b=0; br=0
        for _ in range(M):
            s1=vmf(Pp,1); s2=vmf(Pe,ax); c=dot(s1,s2); b+=c<-0.5; br+=(1-c)/4
        print(f'  P+ {Pp:.2f}, P- {Pe:.2f}, axes {"parallel" if ax>0 else "antiparallel"}: P_S basin {b/M:.4f}  Born {br/M:.4f}  diff {b/M-br/M:+.4f} (MC +-{math.sqrt(0.25*0.75/M):.4f})')
# small-polarization slopes vs product P+P-
print('\n(4) slopes dP_S/d(P+P-) at small polarization (P+ = P- = p, antiparallel):')
for p in (0.05,0.1,0.2):
    b=0;br=0
    for _ in range(M):
        s1=vmf(p,1); s2=vmf(p,-1); c=dot(s1,s2); b+=c<-0.5; br+=(1-c)/4
    print(f'  p {p:.2f}: (P_S - 1/4)/p^2  basin {(b/M-0.25)/p**2:+.3f}  Born {(br/M-0.25)/p**2:+.3f}')

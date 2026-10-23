# Audit 595: Heisenberg-Kramers sum over intermediate l' for <(S.r^)^2>, classical Fourier amplitudes between l and l'.
import math
def qm(l,S,J):
    if S==0.5: return 0.25
    return {l+1:(l+2)/(2*l+3), l:1.0, l-1:(l-1)/(2*l-1)}[J]
def sperp2(J2,S2,L2):
    if L2<=0: return None
    SL=(J2-L2-S2)/2; v=S2-SL*SL/L2
    return v if v>=0 else 0.0
def W(v,l,S,J):
    tot=0.0
    for lp in (l-1,l+1):
        if lp<0: continue
        Lm=max(l,lp)
        if v=='W1': x=sperp2((J+.5)**2,(S+.5)**2,Lm**2); a=x/4 if x is not None else 0
        elif v=='W2': x=sperp2(J*(J+1),S*(S+1),Lm**2-0.25); a=x/4 if x is not None else 0
        elif v=='W3':
            x1=sperp2((J+.5)**2,(S+.5)**2,(l+.5)**2); x2=sperp2((J+.5)**2,(S+.5)**2,(lp+.5)**2)
            a=math.sqrt(x1*x2)/4 if (x1 is not None and x2 is not None) else 0
        elif v=='W4':
            x1=sperp2(J*(J+1),S*(S+1),l*(l+1)); x2=sperp2(J*(J+1),S*(S+1),lp*(lp+1))
            a=math.sqrt(x1*x2)/4 if (x1 is not None and x2 is not None) else 0
        tot+=a
    return tot
cases=[(l,S,J) for l in range(1,7) for S in (0.5,1) for J in ([l-0.5,l+0.5] if S==0.5 else [l-1,l,l+1]) if J>=0]
for v in ('W1','W2','W3','W4'):
    d=[abs(W(v,*c)-qm(*c)) for c in cases]
    print(f"{v}: max |W - QM| = {max(d):.3e};  S=1/2 l=1: "+", ".join(f"{W(v,1,.5,J):.4f}" for J in (.5,1.5))+
          ";  S=1 l=1 J=0,1,2: "+", ".join(f"{W(v,1,1,J):.4f}/{qm(1,1,J):.4f}" for J in (0,1,2))+
          ";  S=1 l=3 J=2,3,4: "+", ".join(f"{W(v,3,1,J):.4f}/{qm(3,1,J):.4f}" for J in (2,3,4)))
print("\nW1 convergence (raport, nie regula): max_J |W1 - QM| per l, and times l^2")
for l in (1,2,3,4,6,10,20,40):
    cs=[(l,S,J) for S in (0.5,1) for J in ([l-0.5,l+0.5] if S==0.5 else [l-1,l,l+1])]
    d=max(abs(W('W1',*c)-qm(*c)) for c in cs); print(f"  l={l:3d}  {d:.3e}   x l^2 = {d*l*l:.3f}")

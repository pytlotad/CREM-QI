# Audit 594: can the ordering term (S.L)^2 -> (S.L)(S.L + 1/2) come from the midpoint principle (two dipole channels of r^)?
from fractions import Fraction as Fr
def qm(l,S,J):
    if S==0.5: return 0.25
    return {l+1:(l+2)/(2*l+3), l:1.0, l-1:(l-1)/(2*l-1)}[J]
def C(S2,SL,L2): return S2/3-(SL*SL-S2*L2/3)/(2*L2)
def variant(v,l,S,J):
    S2=S*(S+1); SL=(J*(J+1)-l*(l+1)-S2)/2; Ll=l+0.5; tot=0.0
    for lp in (l-1,l+1):
        if lp<0: continue
        lg=max(l,lp); w=2*lg/(2*l+1); Lm=float(lg)
        if v=='V1': L2=Lm*Lm; sl=SL
        elif v=='V2': L2=(Lm-0.5)*(Lm+0.5); sl=SL
        elif v=='V3': L2=Lm*Ll; sl=SL
        elif v=='V4': L2=Lm*Lm; sl=SL*Lm/Ll
        if L2<=0: continue
        tot+=w/2*C(S2,sl,L2)
    return tot
def rule(l,S,J):   # the engine's rule (521)
    S2=S*(S+1); SL=(J*(J+1)-l*(l+1)-S2)/2; L2=l*(l+1.0); Ll=l+0.5
    return S2/3-(SL*(SL+0.5)-S2*L2/3)/(2*(Ll-1)*(Ll+1))
cases=[(l,S,J) for l in range(1,7) for S in (0.5,1) for J in ([l-0.5,l+0.5] if S==0.5 else [l-1,l,l+1]) if J>=0]
print("engine rule (521) vs QM: max dev", max(abs(rule(*c)-qm(*c)) for c in cases))
for v in ('V1','V2','V3','V4'):
    devs=[abs(variant(v,*c)-qm(*c)) for c in cases]
    print(f"{v}: max |variant - QM| = {max(devs):.3e};  l=1 S=1 J=0,1,2: "+", ".join(f"{variant(v,1,1,J):.4f}/{qm(1,1,J):.4f}" for J in (0,1,2)))

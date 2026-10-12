# Audit 427: B_1/2 of magnetic quenching (audit 422 formula, variant b: Born projection on field eigenstates)
# with the model's HFS after 426 and 427, against the same QM curve with the measured HFS.
import math
h=6.62607015e-34; muB=9.2740100783e-24; g=2.00231930436
Gp=1/125.14e-12; Go=1/142.04e-9
def y_of(B,hfs): x=2*g*muB*B/(h*hfs); return x/(1+math.sqrt(1+x*x))
def G_T0(B,hfs): y=y_of(B,hfs); return (Go+y*y*Gp)/(1+y*y)
def Bhalf(hfs):
    target=Go+0.5*((Gp+Go)/2-Go); lo,hi=1e-4,1e3
    for _ in range(200):
        mid=math.sqrt(lo*hi)
        if G_T0(mid,hfs)<target: lo=mid
        else: hi=mid
    return math.sqrt(lo*hi)
for name,hfs in [("measured 203.3942",203.3942e9),("model 422 (58.53)",58.53e9),("model 426 (117.063)",117.063395e9),("model 427 (204.658)",204.657666e9)]:
    print(f"{name:22s} B_1/2 = {Bhalf(hfs):.4f} T")

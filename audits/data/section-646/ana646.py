import re,sys
hb=1.054571817e-34; eV=1.602176634e-19
def info(f):
    t=open(f).read().splitlines()
    wl=[l for l in t if l.startswith("WAVE_LOCK t")]; wlL=[l for l in t if l.startswith("WAVE_LOCK_L")]
    led=[l for l in t if l.startswith("CREM_LEDGER t")]
    born=[l.replace("BORN ","") for l in t if l.startswith("BORN formation")]
    n0=re.search(r"n_before=([\d.]+)",wl[0]).group(1) if wl else "-"
    Lafter=[float(re.search(r"L_after=([\d.]+)",l).group(1)) for l in wlL]
    lines=[]; Lb=[]; nb=[]
    for l in led:
        w=l.split(); lines.append(float(w[4])/eV); B=l.split("|")[1].split(); Lb.append(float(B[B.index("L")+1])/hb)
    return n0,len(wlL),Lafter[-1] if Lafter else None,Lb,lines,born
for sysn in ["oPs","pPs","H"]:
    print("=====",sysn)
    for st in ["S1","S2","S3","RL","R0"]:
        for i in range(4):
            n0,nl,La,Lb,lines,born=info(f"runs/{sysn}_{st}_{i}.txt")
            print("%s %d n0 %s locks %d L before 1st photon %s lines %s %s"%(st,i,n0,nl,"%.9f"%Lb[0] if Lb else "-"," ".join("%.10f"%x for x in lines)," ".join(born)))

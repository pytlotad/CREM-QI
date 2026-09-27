import sys,glob,math,statistics as st
S=sys.argv[1]
def load(pat):
    rows=[]
    for f in sorted(glob.glob(S+"/"+pat)):
        for ln in open(f):
            t=ln.split()
            if len(t)==8 and t[0]=="ROW":
                rows.append(dict(seed=int(t[1]),life=float(t[2]),
                    ph=int(t[3]),ceil=int(t[4]),kin=int(t[5]),
                    rec=int(t[6]),out=int(t[7])))
    return rows
ax=load("p?.txt"); na=load("na?.txt")
dl={}
for ln in open(S+"/delta100.txt"):
    t=ln.split()
    if len(t)==9 and t[0].isdigit(): dl[int(t[0])]=float(t[8])
def summarise(rows,label):
    if not rows: print("%s: brak danych"%label); return None
    R=sum(r['ceil']+r['kin']+r['rec'] for r in rows)
    N=sum(r['ph'] for r in rows)
    cens=sum(1 for r in rows if r['out']!=0)
    p=R/(R+N) if R+N else float('nan')
    se=math.sqrt(p*(1-p)/(R+N)) if R+N else float('nan')
    print("%-12s n=%3d  ocenzurowane %3d  fotony %4d  odmowy %4d"
          %(label,len(rows),cens,N,R))
    print("%-12s   sufit %d, kinematyka %d, odrzut %d"
          %("",sum(r['ceil'] for r in rows),sum(r['kin'] for r in rows),
            sum(r['rec'] for r in rows)))
    print("%-12s   p = R/(R+N) = %.4f +- %.4f"%("",p,se))
    return p
pa=summarise(ax,"osiowe")
print()
pn=summarise(na,"bez osiowego")
if pa and pn:
    print()
    print("roznica p: %+.4f"%(pn-pa))
    print("(bez osiowego p jest oszacowaniem OD GORY: ostatnia sekwencja")
    print(" odmow jest ucieta cenzura bez akceptacji)")
# orientation test on the non-axial set
sel=[r for r in na if r['seed'] in dl]
if len(sel)>=20:
    print()
    print("%14s %5s %7s %7s %9s"%("przedzial delta","n","fotony","odmowy","p"))
    edges=[0,60,90,120,180]
    for i in range(len(edges)-1):
        b=[r for r in sel if edges[i]<=dl[r['seed']]<edges[i+1]]
        if not b: continue
        R=sum(r['ceil']+r['kin']+r['rec'] for r in b); N=sum(r['ph'] for r in b)
        print("%6d-%-7d %5d %7d %7d %9.4f"
              %(edges[i],edges[i+1],len(b),N,R,R/(R+N) if R+N else float('nan')))

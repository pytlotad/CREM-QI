import sys,glob,math,collections
S=sys.argv[1]
rows=[]
for f in sorted(glob.glob(S+"/p?.txt")):
    for ln in open(f):
        t=ln.split()
        if len(t)==8 and t[0]=="ROW":
            rows.append((int(t[1]),float(t[2]),int(t[3]),int(t[4]),
                         int(t[5]),int(t[6]),int(t[7])))
ok=[r for r in rows if r[6]==0]
print("trajektorie: %d, nieocenzurowane: %d"%(len(rows),len(ok)))
ceil=sum(r[3] for r in ok); kin=sum(r[4] for r in ok); rec=sum(r[5] for r in ok)
print("odmowy wg powodu: sufit %d, kinematyka %d, odrzut %d"%(ceil,kin,rec))
ks=[r[3]+r[4]+r[5] for r in ok]
h=collections.Counter(ks)
n=len(ks); kbar=sum(ks)/n
p=kbar/(1.0+kbar)
# Fisher information for geometric P(k)=p^k(1-p): I(p)=1/(p(1-p)^2) per obs
se=math.sqrt(p*(1-p)**2/n)
print()
print("kbar = %.4f   MLE p = %.4f +- %.4f"%(kbar,p,se))
print()
print("%4s %10s %12s"%("k","obserw.","geom. oczek."))
chi=0.0; used=0
for k in range(0,max(h)+1):
    obs=h.get(k,0); exp=n*(p**k)*(1-p)
    print("%4d %10d %12.2f"%(k,obs,exp))
    if exp>=5: chi+=(obs-exp)**2/exp; used+=1
print()
print("chi^2 na komórkach z oczek.>=5: %.3f na %d st. swobody"%(chi,max(0,used-1)))
# lifetime check
base=186.72
bad=[(r[0],r[1],r[3]+r[4]+r[5]) for r in ok
     if abs(r[1]-(1+r[3]+r[4]+r[5])*base)>0.5]
print()
print("czas = (1+k)*186.72 ps: %d z %d zgodnych do 0.5 ps"%(n-len(bad),n))
if bad: print("  odstające:",bad[:6])

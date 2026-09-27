import sys,glob,math,statistics as st
S=sys.argv[1]
def load(pats):
    d={}
    for pat in pats:
        for f in sorted(glob.glob(S+"/"+pat)):
            for ln in open(f):
                t=ln.split()
                if len(t)==4 and t[0].isdigit():
                    d[int(t[0])]=(float(t[1]),float(t[2]),float(t[3]))
    return d
a=load(["seeds.txt"]); b=load(["s15_?.txt"]); c=load(["s75_*.txt"])
seeds=sorted(set(a)&set(b)&set(c))
print("ziaren z trzema krokami: %d"%len(seeds))
def aitken(v1,v2,v3):
    d1=v2-v1; d2=v3-v2
    if d2==d1 or d2==0: return float('nan'),float('nan')
    p=math.log(abs(d1/d2))/math.log(2.0) if d1*d2>0 else float('nan')
    return v3-d2*d2/(d2-d1), p
print()
print("%5s %11s %11s %11s %12s %8s"%("ziarno","s=0.30","s=0.15","s=0.075",
      "granica","rzad"))
lims=[];ords=[]
for s in seeds:
    v=[a[s][2],b[s][2],c[s][2]]
    L,p=aitken(*v)
    print("%5d %11.4e %11.4e %11.4e %12.4e %8.3f"%(s,v[0],v[1],v[2],L,p))
    if math.isfinite(L): lims.append(L)
    if math.isfinite(p) and 0.1<p<6: ords.append(p)
print()
m=[st.mean([a[s][2] for s in seeds]),st.mean([b[s][2] for s in seeds]),
   st.mean([c[s][2] for s in seeds])]
L,p=aitken(*m)
print("SREDNIA ZESPOLU")
print("  s=0.30 %.4e   s=0.15 %.4e   s=0.075 %.4e"%(m[0],m[1],m[2]))
print("  rzad zespolu p = %.4f   granica = %.4e"%(p,L))
print("  (214b pozyczylo p = 0.746 z jednego ziarna i dalo 5.73e-04)")
if lims:
    print("  srednia granic po ziarnach = %.4e  sd %.4e  sem %.4e"
          %(st.mean(lims),st.stdev(lims),st.stdev(lims)/math.sqrt(len(lims))))
if ords:
    print("  rzedy po ziarnach: mediana %.3f, n=%d z %d sensownych"
          %(st.median(ords),len(ords),len(seeds)))
# absolute para lifetime, same treatment
mp=[st.mean([a[s][0] for s in seeds]),st.mean([b[s][0] for s in seeds]),
    st.mean([c[s][0] for s in seeds])]
Lp,pp=aitken(*mp)
print()
print("CZAS PARA (bezwzgledny)")
print("  %.6e  %.6e  %.6e"%(mp[0],mp[1],mp[2]))
print("  rzad %.4f  granica %.6e s = %.4f ps"%(pp,Lp,Lp*1e12))
print("  postac zamknieta 223b: 31.123549 ps")

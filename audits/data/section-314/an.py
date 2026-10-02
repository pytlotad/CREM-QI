import glob,math,statistics as st
D='/tmp/claude-1000/-home-teddy-Projekty-Positronium/5981d9ee-dfbe-485a-8f93-5c6d87fd7653/scratchpad/c314/'
def load(prefix,J):
    d={}
    for f in glob.glob(D+'%s_%s_*.txt'%(prefix,J)):
        t=open(f).read().split()
        if len(t)<7: continue
        s,ph=int(t[0]),int(t[1]); d.setdefault(s,{})[ph]=(float(t[2]),int(t[3]),float(t[4]))
    return d
def binom(k,n): return min(1.0,sum(math.comb(n,i) for i in range(n+1) if abs(i-n/2)>=abs(k-n/2))/2**n)
for J in ('0.06','0.08','0.10'):
    d=load('r',J); pairs=[(s,v[1],v[2]) for s,v in sorted(d.items()) if 1 in v and 2 in v]
    if not pairs: print(J,'brak'); continue
    closer=sum(1 for s,p,o in pairs if p[0]<o[0]); rel=[(o[0]-p[0])/p[0] for s,p,o in pairs]
    cutp=sum(1 for s,p,o in pairs if p[1]==0); cuto=sum(1 for s,p,o in pairs if o[1]==0)
    print('J0=%s: par %d, para blizej w %d (%.0f%%), p=%.2g; (orto-para)/para srednio %+.3f, mediana %+.3f; odciecie para %d orto %d; min para %.3f..%.3f, orto %.3f..%.3f r*'%(
        J,len(pairs),closer,100*closer/len(pairs),binom(closer,len(pairs)),st.mean(rel),st.median(rel),cutp,cuto,
        min(p[0] for s,p,o in pairs),max(p[0] for s,p,o in pairs),min(o[0] for s,p,o in pairs),max(o[0] for s,p,o in pairs)))
    print('        cos12 para srednio %+.3f, orto %+.3f'%(st.mean(p[2] for s,p,o in pairs),st.mean(o[2] for s,p,o in pairs)))
k=load('k','0.08'); pk=[(s,v[1],v[2]) for s,v in sorted(k.items()) if 1 in v and 2 in v]
if pk: print('KONTROLA (sila wyl., J0=0,08): %d par, identycznych %d, max |roznica| %.2e r*'%(len(pk),sum(1 for s,p,o in pk if p[0]==o[0]),max(abs(p[0]-o[0]) for s,p,o in pk)))

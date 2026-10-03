import glob,math
D='/tmp/claude-1000/-home-teddy-Projekty-Positronium/5981d9ee-dfbe-485a-8f93-5c6d87fd7653/scratchpad/a316/'
W=40.0
def load(prefix,J):
    d={}
    for f in glob.glob(D+'%s_%s_*.txt'%(prefix,J)):
        t=open(f).read().split()
        if len(t)<7: continue
        s,ph,out,orb=int(t[0]),int(t[1]),int(t[2]),float(t[4])
        d.setdefault(s,{})[ph]=(out==0,min(orb,W),float(t[6]))
    return d
def binom(k,n): return min(1.0,sum(math.comb(n,i) for i in range(n+1) if abs(i-n/2)>=abs(k-n/2))/2**n) if n else 1.0
for J in ('0.07','0.08'):
    d=load('r',J); pairs=[(v[1],v[2]) for s,v in sorted(d.items()) if 1 in v and 2 in v]
    if not pairs: continue
    for name,i in (('para',0),('orto',1)):
        ev=[p[i] for p in pairs]
        n_ann=sum(1 for e in ev if e[0]); orbs=sum(e[1] for e in ev)
        rate=n_ann/orbs; rmst=sum(e[1] for e in ev)/len(ev)
        print('J0=%s %s: anihiluje %2d/%d, RMST(40 orb) = %6.2f orbit, tempo = %.4f /orbite, cos12 sr %+.2f'%(J,name,n_ann,len(ev),rmst,rate,sum(e[2] for e in ev)/len(ev)))
    both=[(p,o) for p,o in pairs if p[0] and o[0]]
    pe=sum(1 for p,o in both if p[1]<o[1])
    # para wczesniej takze gdy para anihiluje, a orto nie (rozstrzygalne z cenzura)
    dec=[(p,o) for p,o in pairs if p[0] or o[0]]
    pfirst=sum(1 for p,o in dec if (p[0] and (not o[0] or p[1]<o[1])))
    ofirst=sum(1 for p,o in dec if (o[0] and (not p[0] or o[1]<p[1])))
    rp=sum(1 for p,o in pairs if p[0])/sum(p[1] for p,o in pairs); ro=sum(1 for p,o in pairs if o[0])/sum(o[1] for p,o in pairs)
    print('   pary rozstrzygalne %d: para pierwsza %d, orto pierwsze %d, p = %.2g; oba anihiluja %d (para wczesniej %d)'%(len(dec),pfirst,ofirst,binom(pfirst,pfirst+ofirst),len(both),pe))
    print('   stosunek temp para/orto = %s  (wykladniczy odpowiednik tau_orto/tau_para)'%('%.2f'%(rp/ro) if ro>0 else 'nieskonczony (orto 0 wejsc)'))
k=load('k','0.08'); pk=[(v[1],v[2]) for s,v in sorted(k.items()) if 1 in v and 2 in v]
if pk: print('KONTROLA sila wyl.: %d par, identycznych %d'%(len(pk),sum(1 for p,o in pk if p[0]==o[0] and p[1]==o[1])))

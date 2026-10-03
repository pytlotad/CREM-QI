import glob,math,statistics as st
D='/tmp/claude-1000/-home-teddy-Projekty-Positronium/5981d9ee-dfbe-485a-8f93-5c6d87fd7653/scratchpad/w317/'
W=40.0; eps=4*(math.pi**2-9)*7.2973525693e-3/(9*math.pi)
print('eps (Ore-Powell) = %.6e = 1/%.1f'%(eps,1/eps))
for q in ('0','1'):
    rows={}
    for f in glob.glob(D+'q%s_*.txt'%q):
        t=open(f).read().split()
        if len(t)<10: continue
        rows[(int(t[0]),int(t[1]))]=dict(ann=int(t[2])==0,orb=min(float(t[4]),W),w0=float(t[7]),w=float(t[8]),c=float(t[9]))
    print('== %s'%('Q0 bez kwantyzacji' if q=='0' else 'Q1 z kwantyzacja'))
    res={}
    for ph,name in ((1,'para'),(2,'orto')):
        r=[v for (s,p),v in rows.items() if p==ph]
        if not r: continue
        ann=[v for v in r if v['ann']]
        nu=len(ann)/sum(v['orb'] for v in r)
        wann=[v['w'] for v in ann]
        pann=st.mean([w+(1-w)*eps for w in wann]) if wann else float('nan')
        res[ph]=(nu,pann)
        print('  %s: n=%d, wejsc %d, tempo wejsc %.4f/orb, <w0> start %.3f, <w> przy wejsciu %s, <w+(1-w)eps> %.4e'%(
            name,len(r),len(ann),nu,st.mean(v['w0'] for v in r),'%.4f (zakres %.4f..%.4f)'%(st.mean(wann),min(wann),max(wann)) if wann else '-',pann))
    if 1 in res and 2 in res:
        (nup,pp),(nuo,po)=res[1],res[2]
        geo=nup/nuo if nuo>0 else float('inf')
        tot=(nup*pp)/(nuo*po) if nuo>0 and po>0 else float('inf')
        print('  stosunek tau_o/tau_p: sama geometria %.2f, geometria x regula wyboru %.1f'%(geo,tot))

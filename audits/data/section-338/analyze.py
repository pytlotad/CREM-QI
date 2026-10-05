import sys
for base in sys.argv[1:]:
    c1w=[];m=[];Sp=[]
    for l in open(base+'.err'):
        if not l.startswith('PROJ'): continue
        d={}
        for tok in l.split()[1:]:
            if '=' in tok:
                k,v=tok.split('=',1)
                try: d[k]=float(v)
                except: pass
        c1w.append(d['cos1w']); Sp.append(d['Spair_hbar'])
        m.append(0.5*(-d['cos1L']+d['cos2L']))   # S.Lhat/hbar: electron gamma<0, positron gamma>0
    n=len(c1w)
    band=sum(abs(abs(x)-1)<0.05 for x in c1w)/n
    near=sum(min(abs(x-k) for k in (-1,0,1))<0.05 for x in m)/n
    lo,hi=min(m),max(m); span=hi-lo
    ints=[k for k in (-1,0,1) if lo-0.05<=k<=hi+0.05]
    uniform=min(1.0,0.1*len(ints)/span) if span>0 else float('nan')
    print(f"{base.split('/')[-1]}: substeps {n}")
    print(f"   R2 cos(mu1,omega1): start {c1w[0]:+.3f} end {c1w[-1]:+.3f} min {min(c1w):+.3f} max {max(c1w):+.3f}; within 0.05 of +-1: {100*band:.1f}%")
    print(f"   R3 S.Lhat/hbar: start {m[0]:+.4f} end {m[-1]:+.4f} range [{lo:+.4f}, {hi:+.4f}]; within 0.05 of an integer: {100*near:.1f}% (uniform on range: {100*uniform:.1f}%)")
    print(f"      |S|/hbar range [{min(Sp):.4f}, {max(Sp):.4f}]")

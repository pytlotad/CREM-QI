# Audit 650 post hoc: classify first-photon lines by multiplet component (model zero-field 2P level -> 1S level, + model m=0 field offset).
import re,sys,math,random
d=sys.argv[1]; eV=1.602176634e-19; N=32
P={'3P0':5.1017085689,'3P1':5.1017312538,'3P2':5.1017493913,'1P1':5.1017388096}; HFS=5.1025799222-5.1017388096
M={f'{k}->1{"3S" if s==0 else "1S"}':v+s*HFS for k,v in P.items() for s in (0,1)}
shift={'0.1':5.1017494372-5.1017493913,'1':5.1017537412-5.1017493913}
def parse(f):
    t=open(f).read(); E=float(re.search(r'Eg (\S+)',t).group(1))/eV
    k=re.search(r'\(2P\) k=(\d+)',t); img=re.findall(r'\(1S\) q=\S+ m=(-?\d+) P_S=\S+ -> (\S+)',t)
    return E,(k.group(1) if k else '-'),(img[-1][1]+' m='+img[-1][0] if img else '-')
def comp(E,B):
    for n,v in M.items():
        for s,tag in ((0,''),(shift[B],'+m0')):
            if abs(E-(v+s))<1e-9: return n+tag  # constants rounded to 1e-10 eV; smallest spacing 4.6e-8 eV
    return 'off'
def fisher(a,b,c,d_):
    n=a+b+c+d_; r1=a+b; c1=a+c; p=lambda x: math.comb(r1,x)*math.comb(n-r1,c1-x)/math.comb(n,c1); p0=p(a)
    return min(1.0,sum(p(x) for x in range(max(0,c1-(n-r1)),min(r1,c1)+1) if p(x)<=p0*(1+1e-9)))
def chi(a,b):
    cats=sorted(set(a)|set(b)); s=0.0
    for c in cats:
        o1,o2=a.count(c),b.count(c); e=(o1+o2)/2
        s+=((o1-e)**2+(o2-e)**2)/e
    return s
random.seed(650)
for B in ['0.1','1']:
    R={(c,s):[parse(f'{d}/B{B}_{c}_{s}_{i}.err') for i in range(N)] for c,s in [('O1','L2'),('O1','c2'),('O1','n23'),('P1','L2')]}
    print(f'=== B = {B} T (model m=0 offset {shift[B]*1e6:+.4f} ueV)')
    C={k:[comp(E,B) for E,_,_ in v] for k,v in R.items()}
    for k,v in R.items():
        offs=[E for (E,_,_),c in zip(v,C[k]) if c=='off']
        dev=[min(abs(E-x) for x in M.values())*1e6 for E in offs]
        hist={}
        for c in C[k]: hist[c]=hist.get(c,0)+1
        print(f'  {k[0]} {k[1]:3s}: off-multiplet {len(offs)}/{N}'+(f' (nearest-component distance {min(dev):.1f}-{max(dev):.1f} ueV)' if dev else '')+'  components '+
              ', '.join(f'{c} {n}' for c,n in sorted(hist.items())))
    L=C[('O1','L2')]
    for s in ['c2','n23']:
        Pc=C[('O1',s)]; a=Pc.count('off'); b=L.count('off')
        obs=chi(L,Pc); pool=L+Pc; cnt=0
        for _ in range(20000):
            random.shuffle(pool); cnt+=chi(pool[:N],pool[N:])>=obs-1e-12
        print(f'  {s} vs L2: (a) off {a}/{N} vs {b}/{N}, Fisher p = {fisher(a,N-a,b,N-b):.3f}; (b) component chi2 {obs:.2f}, permutation p = {cnt/20000:.3f}')
    km={}
    for (E,k,img),c in zip(R[('O1','L2')],L): km.setdefault((k,img),set()).add(c)
    multi=[(k,v) for k,v in km.items() if len(v)>1]
    print(f'  (c) L2 O1: {len(km)} distinct (k, 1S image) labels, {len(multi)} with >1 component:',' ; '.join(f'k={k[0]} {k[1]}: {sorted(v)}' for k,v in multi))
    seen={}
    for (E,k,img),c in zip(R[('O1','L2')],L): seen.setdefault(c,set()).add(k)
    print('      component -> k values (L2 O1):',' ; '.join(f'{c}: {sorted(v,key=lambda x:int(x) if x!="-" else -1)}' for c,v in sorted(seen.items())))

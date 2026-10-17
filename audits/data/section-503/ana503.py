# Audit 503: Zeeman structure of the engine's 2p -> 1s photons in B = 10 T (random direction per run) vs B = 0 (same seeds).
# python3 ana503.py B10.err B0.err
import math, sys
hbar=1.054571817e-34; eV=1.602176634e-19; muB=9.2740100783e-24
def V(a): return [float(x) for x in a[:3]]
def dot(a,b): return sum(x*y for x,y in zip(a,b))
def nrm(a): return math.sqrt(dot(a,a))
def unit(a): n=nrm(a); return [x/n for x in a]
def load(f):
    out={}; idx=None; B=None
    for l in open(f):
        if l.startswith('START'): idx=int(l.split()[2])
        elif l.startswith('BFIELD'): t=l.split(); B=[float(x) for x in t[3:6]]
        elif l.startswith('CREM_LEDGER ') and idx is not None and idx not in out:
            t=l.split(); b=[k for k,x in enumerate(t) if x=='|']; sb=t[b[0]+2:]; sa=t[b[1]+2:]
            g=lambda s,k: s[s.index(k)+1:]
            out[idx]=dict(B=B,Eg=float(t[t.index('Eg')+1]),dir=V(t[t.index('dir')+1:]),hel=float(t[t.index('hel')+1]),
                          sph=V(t[t.index('sphoton')+1:]),LB=float(g(sb,'L')[0]),axB=V(g(sb,'ax')),LA=float(g(sa,'L')[0]),axA=V(g(sa,'ax')),
                          UFB=float(g(sb,'UF')[0]),UFA=float(g(sa,'UF')[0]))
    return out
F=load(sys.argv[1]); Z=load(sys.argv[2]); keys=sorted(set(F)&set(Z))
print(f'paired photons (2p -> 1s, first photon): {len(keys)}')
dE=[]; mg=[]; mL=[]; cth=[]; rows=[]
for k in keys:
    f,z=F[k],Z[k]; B=nrm(f['B']); bh=unit(f['B'])
    dE.append((f['Eg']-z['Eg'])/(muB*B))
    mg.append(dot(f['sph'],bh)/hbar)
    LBv=[x*f['LB'] for x in unit(f['axB'])]; LAv=[x*f['LA'] for x in unit(f['axA'])]
    mL.append((dot(LBv,bh)-dot(LAv,bh))/hbar)
    cth.append(dot(unit(f['dir']),bh)); rows.append((k,dE[-1],mg[-1],mL[-1],cth[-1],f['hel'],(f['UFB']-f['UFA']-(z['UFB']-z['UFA']))/(muB*B)))
print('P1 energy shift (E(10 T) - E(0)) / (mu_B B): min %+.3e max %+.3e mean %+.3e;  of which level-energy (UF) part: min %+.3e max %+.3e'%(min(dE),max(dE),sum(dE)/len(dE),min(r[6] for r in rows),max(r[6] for r in rows)))
h=[0]*10
for m in mg: h[min(9,int((m+1)/0.2))]+=1
print('P2 photon angular momentum along B, m_gamma = s_photon.B^/hbar: histogram over [-1,1] (10 bins):',h,' near -1/0/+1 (|m - k| < 0.1):',
      sum(1 for m in mg if abs(m+1)<0.1),sum(1 for m in mg if abs(m)<0.1),sum(1 for m in mg if abs(m-1)<0.1),' (QM: all at -1/0/+1, 1/3 each)')
print('    orbital Delta m_L = Delta(L.B^)/hbar equals m_gamma: max |Delta m_L - m_gamma| %.2e'%max(abs(a-b) for a,b in zip(mL,mg)))
hc=[0]*5
for c in cth: hc[min(4,int((c+1)/0.4))]+=1
print('P3 photon direction cos(k, B): histogram (5 bins over [-1,1]):',hc)
print('    per photon: idx, dE/muB B, m_gamma, cos(k,B), helicity')
for r in rows: print('     %2d %+.2e %+.3f %+.3f %+.0f'%(r[0],r[1],r[2],r[4],r[5]))

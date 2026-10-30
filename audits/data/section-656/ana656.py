# Audit 656: which energy the spin-orbit transport conserves (TRANSPORT_E half-step prints), formation to first photon.
# usage: python3 ana656.py <656 runs> <654 runs>
import re,sys,math,os
d,base=sys.argv[1],sys.argv[2]; u=1e6/1.602176634e-19
def halfsteps(f):
    L=open(f).read().splitlines(); rows=[]
    for l in L:
        if not l.startswith('TRANSPORT_E'): continue
        t=l.split(); kv={t[i]:t[i+1] for i in range(1,len(t)-1,2)}
        if float(kv['L'])<0.99: break          # after the photon (1S)
        rows.append({k:float(v) for k,v in kv.items()})
    return rows,L
cands={'Ps':[('U_eng',['dUeng']),('U_eng + Z_s',['dUeng','dZs']),('U_cl + Z_s + Z_orb',['dUcl','dZs','dZorb']),('F_tr',['dFtr'])],
       'H':[('U_eng + Z_s + Z_orb',['dUeng','dZs','dZorb']),('U_cl + Z_s + Z_orb',['dUcl','dZs','dZorb']),('F_sc + Z_s + Z_orb',['dFsc','dZs','dZorb'])]}
for B in ('0','0.1','1'):
    for c,name,kind in (('O1','o-Ps','Ps'),('P1','p-Ps','Ps'),('H1','H PHEN 1','H'),('H2','H PHEN 2','H')):
        R=[halfsteps(f'{d}/B{B}_{c}_L2_{i}.err') for i in range(8)]
        n=[len(r[0]) for r in R]
        print(f'=== B {B} T, {name}: half-steps before the photon per seed {min(n)}-{max(n)}')
        for lab,keys in cands[kind]:
            ab=[]; net=[]; bad=0
            for rows,_ in R:
                v=[sum(r[k] for k in keys) for r in rows]
                if any(math.isnan(x) for x in v): bad+=1; continue
                ab.append(sum(abs(x) for x in v)*u); net.append(abs(sum(v))*u)
            if not ab: print(f'  {lab:22s}: rule not acting ({bad}/8 seeds NaN)'); continue
            print(f'  {lab:22s}: sum|dX| mean {sum(ab)/len(ab):10.4e} max {max(ab):10.4e}   |net| mean {sum(net)/len(net):10.4e} max {max(net):10.4e} ueV'
                  +(f'  ({bad} seeds NaN)' if bad else ''))
        if kind=='Ps':
            rel=[abs(r['dFtr']/r['Ftr']) for rows,_ in R for r in rows if r['Ftr']==r['Ftr'] and r['Ftr']!=0]
            dj=[abs(r['dJ']) for rows,_ in R for r in rows]
            print(f'        max |dF_tr/F_tr| per half-step {max(rel,default=float("nan")):.2e};  max |d|J_e|| per half-step {max(dj):.3f} hbar')
        if B!='0' and kind=='Ps':
            same=0; m=0
            for i,(rows,L) in enumerate(R):
                g=f'{base}/B{B}_{c}_L2_{i}.err'
                if not os.path.exists(g): continue
                a=[re.sub(r' Uk .*','',l) for l in L if not l.startswith('TRANSPORT_E')]
                m+=1; same+=a==open(g).read().splitlines()
            print(f'  H656d other lines = 654 (Uk/Zk fields of 655 stripped): {same}/{m}')

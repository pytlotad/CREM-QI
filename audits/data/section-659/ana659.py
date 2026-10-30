# Audit 659: one transport energy for fermion pairs (CREM_PAIR_TRANSPORT=1).
# usage: python3 ana659.py <659 runs> <657 runs> <658 runs B>
import re,sys,math
d,r657,r658=sys.argv[1:4]; eV=1.602176634e-19; u=1e6/eV; h=6.62607015e-34
def lines(f): return open(f).read().splitlines()
def kvline(L,tag):
    for l in L:
        if l.startswith(tag):
            t=l.split(); return {t[i]:float(t[i+1]) for i in range(2,len(t)-1,2) if re.match(r'^-?[0-9.]',t[i+1])}
    return None
def led(L):
    for l in L:
        m=re.match(r'CREM_LEDGER t (\S+) Eg (\S+)',l)
        if m: return float(m.group(1)),float(m.group(2))
    return None
mx=lambda x: max((abs(v) for v in x),default=float('nan'))
mean=lambda x: sum(abs(v) for v in x)/len(x) if x else float('nan')
print('=== H659a: Ps with CREM_PAIR_TRANSPORT=1 against 657 (CREM_PS_TRANSPORT_LEVEL=1), the added print fields stripped')
for B in ('0','0.1','1'):
    for c in ('O1','P1'):
        n=0; same=0
        for i in range(8):
            a=[re.sub(r' CZ .*','',l) for l in lines(f'{d}/B{B}_{c}_L2_{i}.err')]
            b=lines(f'{r657}/B{B}_{c}_L2_{i}.err'); n+=1; same+=a==b
        print(f'  B {B} T {c}: {same}/{n} byte-identical')
for pair,tag in (('H','hydrogen (proton,electron)'),('M','muonium (antimuon,electron)')):
    print(f'=== {tag}')
    for B in ('0','0.1','1'):
        for ph in ('1','2'):
            R=[lines(f'{d}/B{B}_{pair}{ph}_L2_{i}.err') for i in range(8)]
            T=[kvline(L,'CREM_TRANSPORT_LEVEL photon') for L in R]
            ok=[t for t in T if t]
            orb=[(t['Eorb']-t['rung'])*u for t in ok]; W=[t['W']*u for t in ok]
            dF=[sum(float(l.split()[l.split().index('dFH')+1]) for l in L if l.startswith('TRANSPORT_H'))*u for L in R]
            iv=[float(l.split()[l.split().index('interval')+1])/h/1e6 for L in R for l in L if l.startswith('TRANSPORT_H')]
            print(f'  B {B} T PHEN {ph}: photons {len(ok)}/8; dF net mean {mean(dF):.2e} max {mx(dF):.2e} ueV; |E_orb - E_2| max {mx(orb):.2e}, |W| max {mx(W):.2e} ueV;'
                  f' C rms {math.sqrt(sum((t["C"]*u)**2 for t in ok)/len(ok)):.2f} max {mx([t["C"]*u for t in ok]):.2f} ueV; dZ used max {mx([t["used"]*u for t in ok]):.2f} ueV;'
                  f' interval {min(iv):.3f} .. {max(iv):.3f} MHz')
            if pair=='H':
                P=[(led(L),led(lines(f'{r658}/B{B}_H{ph}_L2_{i}.err'))) for i,L in enumerate(R)]
                dE=[(a[1]-b[1])/h/1e6 for a,b in P if a and b]
                print(f'      line vs 658 B (MHz): '+' '.join(f'{v:+.1f}' for v in dE)+f';  rel. <= 1e-9: {sum(abs(v)*1e6*h/b[1]<=1e-9 for v,(a,b) in zip(dE,P))}/{len(dE)}')
# post hoc (after the first output): the switch's own effect on the line is -(C + dZ used)/h at the photon state; the rest of the
# change against 658 B is a different trajectory (time to photon, Zeeman channel, image).
print('=== post hoc: H line change against 658 B = -(C + dZ used)/h + rest (MHz)')
for B in ('0.1','1'):
    for ph in ('1','2'):
        out=[]
        for i in range(8):
            L=lines(f'{d}/B{B}_H{ph}_L2_{i}.err'); t=kvline(L,'CREM_TRANSPORT_LEVEL photon'); a=led(L); b=led(lines(f'{r658}/B{B}_H{ph}_L2_{i}.err'))
            direct=-(t['C']+t['used'])/h/1e6; tot=(a[1]-b[1])/h/1e6
            out.append(f'{direct:+.0f}/{tot-direct:+.0f}(t {a[0]/b[0]-1:+.0e})')
        print(f'  B {B} T PHEN {ph}: direct/rest: '+' '.join(out))

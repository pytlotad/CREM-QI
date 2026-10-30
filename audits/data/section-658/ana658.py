# Audit 658 stage A: is the hydrogen transport the gradient of F_H = Phi(|J_e|) + a (J_e . I_e) + U_Z + omega_L . L ?
# usage: python3 ana658.py <runs>
import sys,math
d=sys.argv[1]; eV=1.602176634e-19; u=1e6/eV; h=6.62607015e-34
def parse(f):
    E=[];H=[]
    for l in open(f):
        if l.startswith('TRANSPORT_E'): E.append(l.split())
        elif l.startswith('TRANSPORT_H'): H.append(l.split())
    g=lambda L,k:[float(x[x.index(k)+1]) for x in L]
    r={k:g(H,k) for k in ('dFH','dFfs','dFhfs','dUZ','dZorb','perp','Ffs')}
    r['dZs']=g(E,'dZs'); r['dUeng']=g(E,'dUeng'); r['J']=g(E,'J')
    r['iv']=g(H,'interval') if H and 'interval' in H[0] else []
    r['kj']=sum(1 for x in H if x[x.index('kb')+1]!=x[x.index('ka')+1]); r['n']=len(H)
    return r
net=lambda x: sum(x)*u
ab=lambda x: sum(abs(v) for v in x)*u
for B in ('0','0.1','1'):
    for ph,name in (('1','PHEN 1 (singlet start)'),('2','PHEN 2 (triplet start)')):
        R=[parse(f'{d}/B{B}_H{ph}_L2_{i}.err') for i in range(8)]
        print(f'=== B {B} T, H {name}: half-steps {sum(r["n"] for r in R)}, k jumps {sum(r["kj"] for r in R)}')
        for k in ('dFH','dFfs','dFhfs','dUZ','dZorb','dZs','dUeng'):
            N=[net(r[k]) for r in R]; A=[ab(r[k]) for r in R]
            print(f'  {k:6s} net mean|.| {sum(abs(v) for v in N)/8:.3e} max|.| {max(abs(v) for v in N):.3e}   abs-sum mean {sum(A)/8:.3e} ueV')
        print(f'  max |w_SO perp|/|w_SO| {max(max(r["perp"]) for r in R):.2e}')
        if R[0]['iv']:
            I=[v/h/1e6 for r in R for v in r['iv']]
            print(f'  fine-structure interval of the generator Phi\'*hbar, all half-steps: {min(I):.3f} .. {max(I):.3f} MHz'
                  f' (measured 10969.045 MHz: {min(I)-10969.045:+.3f} .. {max(I)-10969.045:+.3f})')
# Stage B: with CREM_H_HAMILTONIAN_TRANSPORT=1 (runs B) against stage A (same seeds, without the switch): lines and time to photon.
# usage: python3 ana658.py <runs B> <runs A>
import re
if len(sys.argv)>2:
    A=sys.argv[2]
    def led(f):
        m=re.search(r'CREM_LEDGER t (\S+) Eg (\S+)',open(f).read()); return (float(m.group(1)),float(m.group(2))) if m else None
    print('=== stage B vs stage A (same seeds): photon line and time to photon')
    for B in ('0','0.1','1'):
        for ph in ('1','2'):
            P=[(led(f'{d}/B{B}_H{ph}_L2_{i}.err'),led(f'{A}/B{B}_H{ph}_L2_{i}.err')) for i in range(8)]
            P=[(b,a) for b,a in P if b and a]
            dE=[b[1]/a[1]-1 for b,a in P]; dt=[b[0]/a[0]-1 for b,a in P]
            dMHz=[(b[1]-a[1])/6.62607015e-34/1e6 for b,a in P]
            print(f'  B {B} T PHEN {ph}: photons {len(P)}/8; line rel. change max {max(map(abs,dE)):.2e} (max |dE| {max(map(abs,dMHz)):.3f} MHz); time to photon rel. change max {max(map(abs,dt)):.2e}')

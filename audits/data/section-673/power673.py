# Audit 673: classical dipole power of the Langer orbit against the QM emitted power sum_f hbar omega A, state by state (offline, hydrogen).
# usage: python3 power673.py   (run from this directory)
import numpy as np
src=open('../section-560/tempo_e1_punkt_srodkowy.py').read().split('if __name__ == "__main__":')[0]
g={}; exec(src,g)
rate_qm,channels,omega,C=g['rate_qm'],g['channels'],g['omega'],g['C']
def spectra(n,L,N=1<<20):
    a=n*n; e=np.sqrt(max(0.0,1-(L/n)**2)); b=a*np.sqrt(1-e*e)
    M=(np.arange(N))*2*np.pi/N; u=M.copy() if e<1e-12 else np.full_like(M,np.pi)
    for _ in range(200): u-= (u-e*np.sin(u)-M)/(1-e*np.cos(u))
    z=a*(np.cos(u)-e)+1j*b*np.sin(u); c=np.fft.fft(z)/N
    k=np.arange(1,N//2); Om=1.0/n**3; cp=np.abs(c[k])**2; cm=np.abs(c[N-k])**2
    f=2/3/C**3
    return dict(Pp=np.sum(f*cp*(k*Om)**4),Pm=np.sum(f*cm*(k*Om)**4),Np=np.sum(f*cp*(k*Om)**3))
names="spdfghiklm"; rows=[]
print('state   R_P=P+cl/P(-1)QM   R_N=C0/A(-1)QM   E-factor   P-cl/P(+1)QM')
for n in range(2,11):
    for l in range(1,n):
        s=spectra(n,l+0.5)
        dn=[(n2,l2) for n2,l2 in channels(n,l) if l2==l-1]; up=[(n2,l2) for n2,l2 in channels(n,l) if l2==l+1]
        Pd=sum(omega(n,n2)*rate_qm(n,l,n2,l2) for n2,l2 in dn); Ad=sum(rate_qm(n,l,n2,l2) for n2,l2 in dn)
        Pu=sum(omega(n,n2)*rate_qm(n,l,n2,l2) for n2,l2 in up)
        RP=s['Pp']/Pd; RN=s['Np']/Ad; rows.append((n,l,RP,RN))
        print(f'{n}{names[l]:<4}   {RP:10.4f}        {RN:10.4f}     {RN/RP:8.4f}   '+(f'{s["Pm"]/Pu:10.4f}' if Pu>0 else '       -'),flush=True)
p=[r for r in rows if r[1]==1]
da=max(abs(r[2]-r[3]) for r in p)
print(f'\nH673a (p states, |R_P - R_N| <= 0.01): max {da:.4f} -> {da<=0.01}')
print(f'H673b (p states R_P in [0.92, 0.96]): {min(r[2] for r in p):.4f} .. {max(r[2] for r in p):.4f} -> {all(0.92<=r[2]<=0.96 for r in p)}')
ok=True
for n in range(6,11):
    v=[abs(1-r[2]) for r in rows if r[0]==n]
    mono=all(v[i]>=v[i+1]-1e-12 for i in range(len(v)-1)); ok&=mono
    print(f'  n={n}: |1-R_P| by l = '+' '.join(f'{x:.4f}' for x in v)+f'  decreasing {mono}')
print(f'H673c: {ok}')

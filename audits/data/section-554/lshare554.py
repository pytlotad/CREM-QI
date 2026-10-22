# Audit 554: exact nonrelativistic dipole recombination into (n, l) (gaunt553.py functions) -- l shares within each n, at K/B of the experiment-5 samples.
import sys, mpmath as mp
sys.path.insert(0,'../section-553')
import importlib.util
spec=importlib.util.spec_from_file_location('g553','../section-553/gaunt553.py')
src=open('../section-553/gaunt553.py').read().split('print("check')[0]
ns={}; exec(src,ns)
for KoverB in ('0.05','0.18'):
    k=mp.sqrt(mp.mpf(KoverB)); tot=0; rows=[]
    for n in range(1,5):
        omega=k*k/2+mp.mpf(1)/(2*n*n)
        per=[2*(2*l+1)*(omega*ns['alpha']/k)**2*ns['sigma_PI'](n,l,k) for l in range(n)]
        s=sum(per); rows.append((n,per,s)); tot+=s
    print(f"K/B = {KoverB}:")
    for n,per,s in rows:
        print(f"  n = {n}: share of n {mp.nstr(s/tot,4)}; l shares "+' '.join(f"l={l}: {mp.nstr(x/s,4)}" for l,x in enumerate(per)))

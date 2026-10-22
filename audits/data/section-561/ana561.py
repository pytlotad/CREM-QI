# Audit 561: first-photon times of an ensemble (CREM_LEDGER after START idx), MLE tau, KS against Exp(tau), and against the measured lifetime.  python3 ana561.py label tau_meas_ns sigma_ns files...
import sys, math
label, tm, sm = sys.argv[1], float(sys.argv[2]), float(sys.argv[3])
times=[]; cens=0
for f in sys.argv[4:]:
    cur=None; got=False
    for l in open(f):
        if l.startswith('START idx'):
            if cur is not None and not got: cens+=1
            cur=int(l.split()[2]); got=False
        elif l.startswith('CREM_LEDGER t') and cur is not None and not got:
            times.append(float(l.split()[2])*1e9); got=True
    if cur is not None and not got: cens+=1
n=len(times); tau=sum(times)/n; se=tau/math.sqrt(n)
xs=sorted(times); d=max(max((i+1)/n-(1-math.exp(-x/tau)), (1-math.exp(-x/tau))-i/n) for i,x in enumerate(xs))
lam=(math.sqrt(n)+0.12+0.11/math.sqrt(n))*d; p=max(0.0,min(1.0,2*sum((-1)**(k-1)*math.exp(-2*k*k*lam*lam) for k in range(1,100))))
cv=math.sqrt(sum((x-tau)**2 for x in times)/(n-1))/tau
print(f"{label}: N {n} (no photon: {cens}); MLE tau {tau:.4f} +- {se:.4f} ns; sigma/mean {cv:.3f} (Exp: 1); KS vs Exp(tau) D {d:.4f} p {p:.3f} (fitted mean -- p conservative); "
      f"measured {tm} +- {sm} ns: ratio {tau/tm:.4f} +- {se/tm:.4f}")

# Audit 575: correction variants (prereg575.txt); measured: p-Ps 125,142 ps (Al-Ramadhan, Gidley 1994),
# o-Ps 142,037 ns (Vallery, Zitzewitz, Gidley 2003, as in the program), HFS 203,3942 GHz (Ishida et al. 2014).
import math
a=7.2973525693e-3; api=a/math.pi; L=math.log(1/a); mc2_hbar=510998.95/6.582119569e-16
Gp0=a**5*mc2_hbar/2; Go0=2*(math.pi**2-9)/(9*math.pi)*a**6*mc2_hbar
var={'one loop only':(1-(5-math.pi**2/4)*api,1-10.286606*api),
 'one loop + log, literature signs':(1-(5-math.pi**2/4)*api+2*a*a*L,1-10.286606*api-a*a*L/3),
 'author script (o-Ps log sign reversed)':(1-(5-math.pi**2/4)*api+2*a*a*L,1-10.286606*api+a*a*L/3),
 'order alpha^2 (B_p 5.1243, B_o 45.06) -- adopted':(1-(5-math.pi**2/4)*api+2*a*a*L+5.1243*api**2,1-10.286606*api-a*a*L/3+45.06*api**2)}
for k,(fp,fo) in var.items():
    tp=1e12/(Gp0*fp); to=1e9/(Go0*fo)
    print(f"{k:48s} tau_p {tp:8.3f} ps ({100*(tp/125.142-1):+.3f}%)  tau_o {to:8.3f} ns ({100*(to/142.037-1):+.3f}%)  ratio {to*1e3/tp:7.1f}")
E4=a**4*mc2_hbar/(2*math.pi)/1e9; anom=E4/3*((1+0.00115965218)**2-1); one=-api*(8/9+0.5*math.log(2))*E4
print("HFS: LO %.4f, anomaly %.4f, one loop %.4f, 204.66 + one loop - anomaly = %.4f GHz"%(7/12*E4,anom,one,204.66+one-anom))

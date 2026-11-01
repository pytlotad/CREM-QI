# Audit 674: Ps n = 2 fine structure nu_J = E(2^3S1) - E(2^3P_J) of the current engine against the measurements.  usage: python3 nu674.py out/p640_now.txt
import re,sys,math
t=open(sys.argv[1]).read(); E={m.group(1):float(m.group(2)) for m in re.finditer(r'(2\^\dP\d|2\^3S1) (\S+)',t)}
h=6.62607015e-34; e0=1.602176634e-19; ueV=lambda f: h*f*1e6/e0*1e6
c0=299792458.0
meas={0:(18501.02,0.61,'Gurung i in., PRL 125, 073002 (2020)',18498.25,0.08),
      1:(13012.42,math.hypot(0.65,1.54),'Hagena i in., PRL 71, 2887 (1993)',13012.41,0.08),
      2:(8627.94,0.95,'Sheldon i in., PRL 131, 043001 (2023)',8626.71,0.08)}
for J in (0,1,2):
    nu=E['2^3S1']-E[f'2^3P{J}']; m,s,src,th,ths=meas[J]
    res=nu-m; pred=-266.04+(th-m); unc=2*math.hypot(s,ths); ok=abs(res-pred)<=unc
    print(f'nu{J}: engine {nu:.2f} MHz ({c0/(nu*1e6)*1e3:.4f} mm; {ueV(nu):.3f} ueV); measured {m:.2f}({s:.2f}) MHz [{src}]; theory {th:.2f}({ths:.2f})')
    print(f'      engine - measured {res:+.2f} MHz ({ueV(res):+.4f} ueV); predicted from 640 block -266.04 + (theory - measured) = {pred:+.2f} +- {unc:.2f}: H674b {ok};'
          f' theory - measured {th-m:+.2f} MHz ({(th-m)/s:+.1f} sigma)')
print(f"2^3S1 - 2^1P1: engine {E['2^3S1']-E['2^1P1']:.2f} MHz (Ley i in. 1994: 11 180 MHz; teoria 11 185,37)")

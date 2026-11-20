import sys; import os; sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import mpmath as mp; from cdarw816 import Pair, HBAR
m1=mp.mpf('9.1093837139e-31'); kq=mp.mpf('-2.3070775508e-28'); c=mp.mpf(299792458)
PS=Pair(m1,1e12*m1,kq,c,dar=1); k=-kq
Es=lambda Jr,L: m1*c**2/mp.sqrt(1+(k/(c*(Jr+mp.sqrt(L**2-k**2/c**2))))**2)
for L,Jr in ((2*HBAR,0),(HBAR,HBAR)):
    E=PS.etorus(2,L); print('code E-Mc2',mp.nstr(E-PS.M*c**2,15),' Sommerfeld',mp.nstr(Es(Jr,L)-m1*c**2,15))
def ecirc2(P,L):
    r0=L**2/(P.mu*(-P.kq)); f=lambda r: P.e_rs(r,0,L)[0]
    d=lambda r: (f(r*(1+mp.mpf('1e-12')))-f(r*(1-mp.mpf('1e-12'))))/(2*r*mp.mpf('1e-12'))
    rc=mp.findroot(d,(r0*mp.mpf('0.9'),r0*mp.mpf('1.1')),solver='anderson'); return f(rc)
print('fixed one-body circ', mp.nstr(ecirc2(PS,2*HBAR)-PS.M*c**2,15))
P=Pair(m1,m1,kq,c,dar=1); print('Ps circ old',mp.nstr(P.ecirc(2*HBAR)-P.M*c**2,20),' new',mp.nstr(ecirc2(P,2*HBAR)-P.M*c**2,20))

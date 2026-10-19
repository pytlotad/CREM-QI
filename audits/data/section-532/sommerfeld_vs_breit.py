# Import 21 (Sommerfeld k = j+1/2, postac F(k) - F(L_l)) a spin-orbita Breita-Pauliego xi<L.S>
# z <r^-3> regula Borna (519) i <L.S> z trojkata Borna. Wodor (m_p -> inf), g = 2, jednostki m c^2.
# F(L) = A + B/L, B = -alpha^4/(2 n^3) (klasyczna srednia Darwina, 517). Wynik: rownosc dokladna.
from fractions import Fraction as Fr
for n in range(2,9):
    B=Fr(-1,2*n**3)                         # w jednostkach alpha^4 m c^2
    for l in range(1,n):
        F=lambda k: B/k
        mean=(F(l+1)*(2*l+2)+F(l)*2*l)/(4*l+2)   # srednia (2j+1)-wazona = F(L_l) importu 21
        r3=Fr(1,n**3*l*(2*l+1)*(l+1))*2           # <r^-3> Borna: 1/(n^3 l(l+1/2)(l+1)) w a0^-3
        for j2 in (2*l+1,2*l-1):
            j=Fr(j2,2); k=j+Fr(1,2)
            LS=(j*(j+1)-l*(l+1)-Fr(3,4))/2
            som=F(k)-mean
            bp=Fr(1,2)*r3*LS                      # xi = alpha^4 m c^2 /2 * <r^-3> (z Thomasem)
            assert som==bp,(n,l,j,som,bp)
print("Sommerfeld F(k)-F(L_l) == Breit-Pauli xi<L.S> (regula Borna) dla n=2..8, wszystkich l, j: OK")
# Energia splatania singletu (import 8) a trojkat Borna: <sigma1.sigma2> = 2[S(S+1) - 3/2]
for S,cl in ((0,-1),(1,1)):
    print("S=%d: trojkat Borna %+d, klasyczne wektory %+d"%(S,2*(S*(S+1))-3,cl))

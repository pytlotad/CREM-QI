// Znak sily dipol-dipol w modelu (pairDipoleForce): ujemna skladowa radialna
// sily na cel = PRZYCIAGANIE.  Separacja wzdluz x; zrodlo w zerze.
#include "modules/crem_collapse.hpp"
#include <cstdio>
int main(){
    const double rs=comptonBarrierRadius, mu=firstMagneticMoment;
    const double coul=pairCoulombStrength;   // k e^2 [J m]
    std::printf("eps (Plummer) = %.4f r*\n",magneticDipoleRadius()/rs);
    std::printf("%8s | %14s %14s | %14s %14s | %14s %14s\n","r/r*",
        "para wzdl. L","orto wzdl. L","para wzdl. r","orto wzdl. r","para 45deg","orto 45deg");
    const Vec3 z{0,0,1}, x{1,0,0}, d45=Vec3{1,0,1}/std::sqrt(2.0);
    for(double q:{0.25,0.5,1.0,2.0,5.0,20.0,100.0,547.5}){
        const Vec3 sep=x*(q*rs);
        auto F=[&](Vec3 a,Vec3 b){ // a: cel, b: zrodlo; skladowa radialna / Coulomb
            const Vec3 f=pairDipoleForce(sep,a*mu,b*mu);
            return dot(f,x)/(coul/(q*rs*q*rs)); };
        std::printf("%8.2f | %+14.4e %+14.4e | %+14.4e %+14.4e | %+14.4e %+14.4e\n",q,
            F(z,z),F(z,z*-1.0),F(x,x),F(x,x*-1.0),F(d45,d45),F(d45,d45*-1.0));
    }
    std::printf("(wartosci = radialna sila dipolowa / sila Coulomba; < 0 przyciaganie)\n");
}

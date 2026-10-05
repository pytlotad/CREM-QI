// E_SO = sum_i S_i . omega_i(orbital fields only), orbit-averaged, from CREM_SPINVEC lines; and the
// energy change a circular orbit would need for its |L| (from CREM_SPINENERGY/SKIP not needed: use L from vector state).
#include "modules/secular_spin_orbit.hpp"
#include <cstdio>
#include <iostream>
#include <string>
static Vec3 v3(const char* s){ Vec3 v; std::sscanf(s,"%lf,%lf,%lf",&v.x,&v.y,&v.z); return v; }
int main(){
    const double M=firstMass+secondMass, mu=firstMass*secondMass/M, eV=1.602176634e-19;
    std::string line; int run=0; double prevT=-1;
    while(std::getline(std::cin,line)){
        if(line.rfind("idx",0)==0){ std::printf("%s\n",line.c_str()); continue; }
        if(line.rfind("CREM_SPINVEC",0)!=0) continue;
        double t,a,e; char m1s[200],m2s[200],Ls[200],Ps[200];
        if(std::sscanf(line.c_str(),"CREM_SPINVEC t=%lf a=%lf e=%lf m1=%199s m2=%199s L=%199s P=%199s",&t,&a,&e,m1s,m2s,Ls,Ps)!=7) continue;
        Vec3 m1=v3(m1s), m2=v3(m2s), L=v3(Ls), P=v3(Ps); L=L*(1/L.norm()); P=P-L*dot(P,L); P=P*(1/P.norm()); const Vec3 Q=cross(L,P);
        const double ec=std::sqrt(std::max(0.0,1-e*e)), n=std::sqrt(pairCoulombStrength/(mu*a*a*a));
        const int N=4096; double E=0,W=0;
        const Vec3 S1=m1*(1.0/firstGyromagneticRatioOf()), S2=m2*(1.0/secondGyromagneticRatioOf());
        for(int k=0;k<N;++k){ const double Ea=2*pi*(k+0.5)/N, cs=std::cos(Ea), sn=std::sin(Ea), tw=1-e*cs;
            const Vec3 r=P*(a*(cs-e))+Q*(a*ec*sn), v=(P*(-a*n*sn)+Q*(a*n*ec*cs))*(1/tw);
            State s{}; s.firstPosition=r*(secondMass/M); s.secondPosition=r*(-firstMass/M);
            s.firstVelocity=v*(secondMass/M); s.secondVelocity=v*(-firstMass/M);
            s.firstProperDipole=m1*1e-30; s.secondProperDipole=m2*1e-30; synchronizeCovariantDipoles(s);  // partner dipole ~0: orbital fields only
            const StateHistory h{State{s}}; const LocalElectromagneticFields f=localRelativisticFields(s,h);
            const Vec3 w1=thomasBmtEffectiveField(s.firstVelocity,f.atFirst,firstGFactor)*(-firstCharge/firstMass);
            const Vec3 w2=thomasBmtEffectiveField(s.secondVelocity,f.atSecond,secondGFactor)*(-secondCharge/secondMass);
            E+=tw*(dot(S1,w1)+dot(S2,w2)); W+=tw; }
        E/=W;
        std::printf("t=%9.4f ps  E_SO=%+.6e eV  |S|/hbar=%.5f\n",t*1e12,E/eV,(S1+S2).norm()/hbar);
    }
}
